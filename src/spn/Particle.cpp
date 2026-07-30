
#include "Particle.h"

#include "Spring.h"
#include "SpringNetwork.h"
#include "forcefield/constants.hpp"
#include "logging.h"
#include "measure.hpp"

#include <sstream>

namespace logging = biospring::logging;

namespace biospring
{
namespace spn
{

void Particle::updateFromForceField(const biospring::forcefield::ForceField & ff)
{
    if (getCGName() == "")
    {
        logging::die("particle %d (%s:%d:%s) has no CG name", getExtid(), getResName().c_str(), getResId(),
                     getName().c_str());
    }
    if (!ff.hasProperty(getCGName()))
    {
        logging::warning("No Forcefield found for '%s'!", getCGName().c_str());
    }
    else
    {
        ParticleProperty pp = ff.getPropertiesFromName(getCGName());
        setCharge(pp.getCharge());
        setRadius(pp.getRadius());
        setMass(pp.getMass());
        setEpsilon(pp.getEpsilon());
        setHydrophobicity(pp.getHydrophobicity());
        setTransferEnergyByAccessibleSurface(pp.getTransferEnergyByAccessibleSurface());
    }
}

void Particle::addToSpringNeighbors(unsigned index, Spring * spring)
{
    _springneighbors.insert(unordered_map<unsigned, Spring *>::value_type(index, spring));
}

string Particle::tostr() const
{
    ostringstream oss;
    oss << _chainname << "::" << _resname << "::" << _resid << "::" << _name;
    return oss.str();
}

float Particle::distance(const Particle & p1, const Particle & p2) { return biospring::measure::distance(p1, p2); }

float Particle::distance(const Particle & p) const { return distance(*this, p); }

// ======================================================================================
//
// Setters
//
// ======================================================================================

void Particle::setPosition(const Vector3f & p)
{
    _previousposition = _position;
    _position = p;
}

void Particle::resetForce()
{
    _force.reset();
    _electrostaticenergy = 0.0;
    _stericenergy = 0.0;
    _kineticenergy = 0.0;
    _impenergy = 0.0;
    _hydrophobicityenergy = 0.0;
    _hydrogenbondcorerepulsionenergy = 0.0;
}

void Particle::applyViscosity(float viscosity)
{
    Vector3f visc = _velocity * viscosity;
    _force = _force - visc;
}

// ======================================================================================
//
// Integration methods
//
// ======================================================================================

void Particle::IntegrateVelocityVerlet(float timestep)
{
    // A configured mass of 0 (e.g. a haptic probe meant to be positioned
    // externally rather than driven by F=ma) would otherwise divide by
    // zero and send velocity/position to NaN; treat it as immovable by
    // force instead.
    if (getMass() > 0.0f)
    {
        Vector3f var = (_force / getMass()) * 0.5f * timestep;
        _velocity = _velocity + var;
        _position = _position + _velocity * timestep;
        _velocity = _velocity + var;
    }
    float vitesse = _velocity.norm();
    _kineticenergy = 0.5f * getMass() * (vitesse * vitesse) * biospring::forcefield::GLOBAL_KINETIC_ENERGY_CONVERT;
}

void Particle::IntegrateEuler(float timestep)
{
    _integrateForce(timestep);
    _integrateVelocity(timestep);
}

void Particle::_integrateForce(float timestep)
{
    // See IntegrateVelocityVerlet: guard against a configured mass of 0.
    if (getMass() > 0.0f)
        _velocity = _velocity + (_force / getMass()) * timestep;
    float vitesse = _velocity.norm();
    _kineticenergy = 0.5f * getMass() * (vitesse * vitesse) * biospring::forcefield::GLOBAL_KINETIC_ENERGY_CONVERT;
}

void Particle::_integrateVelocity(float timestep)
{
    _previousposition = _position;
    _position = _position + _velocity * timestep;
}

// ======================================================================================
//
// Add force methods
//
// ======================================================================================

void Particle::addDensityFieldForce()
{
    // Own scale (densitygrid.scale), independent from the steric force's
    // steric.gridscale -- see SpringNetwork::getDensityGridScale.
    float gridscale = _springnetwork->getDensityGridScale();
    const biospring::grid::PotentialGrid & potentialgrid = _springnetwork->getDensityGrid();

    // Off-grid guard: a steered or free-moving particle can leave the density
    // grid. DenseGrid::get -> at() throws std::out_of_range for an out-of-bounds
    // cell, which aborts the whole run. The JAX port (potential/density_field.py)
    // instead contributes ZERO force for out-of-grid particles (it clamps the
    // index only to keep the gather safe, then masks the value to zero). Mirror
    // that here: skip the lookup and add nothing.
    if (potentialgrid.is_out_of_grid(biospring::grid::real_coordinates(getX(), getY(), getZ())))
    {
        BIOSPRING_WARN_ONCE("particle %u left the density grid: contributing zero density force", getId());
        return;
    }

    Vector3f force = potentialgrid.get(getX(), getY(), getZ()).vector;
    force = force * getBurying() * gridscale;

    addForce(force);
}

void Particle::addElectrostaticFieldForce()
{
    const biospring::forcefield::ForceField * ff = _springnetwork->getForceField();
    // Same scale (potentialgrid.scale, via setForceFieldScale) as the energy
    // computed below (computeElectrostaticFieldEnergy): force and energy must
    // agree on what they are scaling. Previously read steric.gridscale here,
    // an unrelated setting shared with the steric force.
    float gridscale = ff->getForceFieldScale();
    const biospring::grid::PotentialGrid & potentialgrid = _springnetwork->getPotentialGrid();

    // Off-grid guard (see addDensityFieldForce): mirror the JAX port's zero-force
    // out-of-bounds behaviour instead of throwing std::out_of_range and crashing.
    // Note this also skips the energy accumulation below, so an off-grid particle
    // stops contributing to the reported electrostatic energy.
    if (potentialgrid.is_out_of_grid(biospring::grid::real_coordinates(getX(), getY(), getZ())))
    {
        BIOSPRING_WARN_ONCE("particle %u left the electrostatic potential grid: contributing zero field force "
                            "and zero field energy",
                            getId());
        return;
    }

    const auto & cell = potentialgrid.get(getX(), getY(), getZ());

    Vector3f force = cell.vector * getCharge() * gridscale;
    addForce(force);

    _electrostaticenergy += ff->computeElectrostaticFieldEnergy(cell.scalar, getCharge());
}

void Particle::addElectrostaticForce(std::vector<DeferredNonbondedContribution> & deferred)
{
    Vector3f f = Vector3f();
    float distance = 0.0;
    float cutoff = _springnetwork->getElectrostaticCutoff();
    bool apply = true;

    const biospring::forcefield::ForceField * ff = _springnetwork->getForceField();
    _springnetwork->getNeighborSearch().electrostatic->for_each_neighbor(*this, [&](size_t neighbor_index) {
        if (_springnetwork->isProbeParticle(neighbor_index))
            return;

        const Particle & p = _springnetwork->getParticle(neighbor_index);

        // A dynamic neighbor will independently visit this same pair from its
        // own side; process it only once, from the lower-id side, instead of
        // recomputing it twice. A static neighbor never visits any pair on
        // its own, so it must always be processed here.
        if (p.isDynamic() && neighbor_index < static_cast<size_t>(getId()))
            return;

        apply = true;
        if (_springnetwork->isSpringEnabled())
            // getId() stays signed (used as an "unassigned" sentinel elsewhere,
            // e.g. SpringNetwork::isProbeParticle), so cast explicitly here
            // where it is consumed as a spring-neighbor index (always >= 0
            // once assigned by SpringNetwork).
            apply = !isInSpringNeighbors(static_cast<unsigned>(p.getId()));

        if (apply)
        {
            f = p.getPosition() - getPosition();
            distance = f.norm();
            if (distance < cutoff && distance != 0.0)
            {
                const float pair_energy = ff->computeElectrostaticEnergy(p.getCharge(), getCharge(), distance);
                f.normalize();
                f = f * ff->computeElectrostaticForceModule(p.getCharge(), getCharge(), distance);
                addForce(f);

                if (p.isDynamic())
                {
                    // Combining rules are symmetric, so a dynamic neighbor
                    // would have computed exactly `-f` and `0.5 *
                    // pair_energy` had it visited this pair on its own
                    // (Newton's third law, applied explicitly instead of
                    // through a redundant computation).
                    _electrostaticenergy += 0.5f * pair_energy;
                    deferred.push_back({static_cast<unsigned>(p.getId()), -f, 0.5f * pair_energy});
                }
                else
                {
                    // A static neighbor never visits this pair on its own, so
                    // it never contributes its own half: crediting only half
                    // here would silently drop the other half of the pair's
                    // energy from the system total. It never receives a
                    // force either way (it is never integrated or reset).
                    _electrostaticenergy += pair_energy;
                }
            }
        }
    });
}

/// @brief Add IMPALA force to the particle. 
/// See @ref biospring::forcefield::ForceField::computeIMPEnergy and 
/// @ref biospring::forcefield::ForceField::computeIMPForceVector for more informations.
/// @callergraph @callgraph
void Particle::addIMPForce()
{
    const biospring::forcefield::ForceField * ff = _springnetwork->getForceField();
    _impenergy = ff->computeIMPEnergy(getX(), getY(), getZ(), getSolventAccessibilitySurface(),
                                      getTransferEnergyByAccessibleSurface());
    Vector3f f = ff->computeIMPForceVector(getX(), getY(), getZ(), getSolventAccessibilitySurface(),
                                             getTransferEnergyByAccessibleSurface());
    addForce(f);
}

void Particle::addHydrophobicityForce(std::vector<DeferredNonbondedContribution> & deferred)
{
    Vector3f f = Vector3f();
    float distance = 0.0;
    float cutoff = _springnetwork->getHydrophobicCutoff();
    bool apply = true;

    if (_springnetwork->isHydrophobicityEnabled())
    {
        const biospring::forcefield::ForceField * ff = _springnetwork->getForceField();
        _springnetwork->getNeighborSearch().hydrophobic->for_each_neighbor(*this, [&](size_t neighbor_index) {
            if (_springnetwork->isProbeParticle(neighbor_index))
                return;

            const Particle & p = _springnetwork->getParticle(neighbor_index);

            // See addElectrostaticForce: process each pair once, from the
            // lower-id side when the neighbor is dynamic; a static neighbor
            // never visits any pair on its own, so it is always processed.
            if (p.isDynamic() && neighbor_index < static_cast<size_t>(getId()))
                return;

            apply = true;
            if (_springnetwork->isSpringEnabled())
                // getId() stays signed (used as an "unassigned" sentinel elsewhere,
            // e.g. SpringNetwork::isProbeParticle), so cast explicitly here
            // where it is consumed as a spring-neighbor index (always >= 0
            // once assigned by SpringNetwork).
            apply = !isInSpringNeighbors(static_cast<unsigned>(p.getId()));

            if (apply)
            {
                f = p.getPosition() - getPosition();
                distance = f.norm();
                if (distance < cutoff && distance != 0.0)
                {
                    const float pair_energy =
                        ff->computeHydrophobicityEnergy(p.getHydrophobicity(), getHydrophobicity(), distance);
                    f.normalize();
                    f = f * ff->computeHydrophobicityForceModule(p.getHydrophobicity(), getHydrophobicity(), distance);
                    addForce(f);

                    if (p.isDynamic())
                    {
                        // See addElectrostaticForce for why `-f` and half the
                        // energy are exactly what a dynamic neighbor would
                        // have computed on its own.
                        _hydrophobicityenergy += pair_energy / 2.0f;
                        deferred.push_back({static_cast<unsigned>(p.getId()), -f, pair_energy / 2.0f});
                    }
                    else
                    {
                        // A static neighbor never contributes its own half;
                        // see addElectrostaticForce for why the full energy
                        // must be credited here instead.
                        _hydrophobicityenergy += pair_energy;
                    }
                }
            }
        });
    }
}

void Particle::addHydrogenBondCoreRepulsion(std::vector<DeferredNonbondedContribution> & deferred)
{
    Vector3f f = Vector3f();
    float distance = 0.0;
    bool apply = true;

    if (_springnetwork->isHydrogenBondEnabled() && _springnetwork->getNeighborSearch().hbond)
    {
        const biospring::forcefield::ForceField * ff = _springnetwork->getForceField();
        const float equilibrium = ff->getHydrogenBondEquilibrium();
        const int my_partner = _springnetwork->getHydrogenBondPartner(static_cast<size_t>(getId()));
        const bool self_is_donor = isDonor();
        const bool self_is_acceptor = isAcceptor();

        _springnetwork->getNeighborSearch().hbond->for_each_neighbor(*this, [&](size_t neighbor_index) {
            if (_springnetwork->isProbeParticle(neighbor_index))
                return;

            // Already fully handled (attraction and repulsion both) by
            // computeHydrogenBondForces for this particle's own current
            // exclusive partner -- applying this term too would double-count.
            if (static_cast<int>(neighbor_index) == my_partner)
                return;

            const Particle & p = _springnetwork->getParticle(neighbor_index);

            // No explicit H, so only the opposite role is meaningful (same
            // rule as the exclusive mechanism -- see _assignHydrogenBondPairs).
            const bool roles_match = (self_is_donor && p.isAcceptor()) || (self_is_acceptor && p.isDonor());
            if (!roles_match)
                return;

            // See addElectrostaticForce: process each pair once, from the
            // lower-id side when the neighbor is dynamic; a static neighbor
            // never visits any pair on its own, so it is always processed.
            if (p.isDynamic() && neighbor_index < static_cast<size_t>(getId()))
                return;

            apply = true;
            if (_springnetwork->isSpringEnabled())
                apply = !isInSpringNeighbors(static_cast<unsigned>(p.getId()));

            if (apply)
            {
                f = p.getPosition() - getPosition();
                distance = f.norm();
                // Repulsive-only regime: the true Morse force is already
                // zero exactly at distance=equilibrium, so cutting off here
                // introduces no discontinuity, and leaves the attractive
                // range (distance > equilibrium) exclusively to whichever
                // pair is actually engaged via _assignHydrogenBondPairs.
                if (distance < equilibrium && distance != 0.0)
                {
                    const float pair_energy = ff->computeHydrogenBondEnergy(distance);
                    f.normalize();
                    f = f * ff->computeHydrogenBondForceModule(distance);
                    addForce(f);

                    if (p.isDynamic())
                    {
                        _hydrogenbondcorerepulsionenergy += pair_energy / 2.0f;
                        deferred.push_back({static_cast<unsigned>(p.getId()), -f, pair_energy / 2.0f});
                    }
                    else
                    {
                        _hydrogenbondcorerepulsionenergy += pair_energy;
                    }
                }
            }
        });
    }
}

void Particle::addElectrostaticForceNoGrid(float cutoff)
{
    Vector3f f = Vector3f();
    float distance = 0.0;
    bool apply = true;

    const biospring::forcefield::ForceField * ff = _springnetwork->getForceField();
    vector<unsigned> chargedparticules = _springnetwork->getChargedParticles();
    for (unsigned i = 0; i < chargedparticules.size(); i++)
    {
        const Particle & p = _springnetwork->getParticle(chargedparticules[i]);
        apply = true;
        if (_springnetwork->isSpringEnabled())
        {
            // getId() stays signed (used as an "unassigned" sentinel elsewhere,
            // e.g. SpringNetwork::isProbeParticle), so cast explicitly here
            // where it is consumed as a spring-neighbor index (always >= 0
            // once assigned by SpringNetwork).
            apply = !isInSpringNeighbors(static_cast<unsigned>(p.getId()));
        }

        if (apply && p.getId() != getId())
        {
            f = p.getPosition() - getPosition();
            distance = f.norm();
            if (distance < cutoff && distance != 0.0f)
            {
                _electrostaticenergy += ff->computeElectrostaticEnergy(p.getCharge(), getCharge(), distance);
                f.normalize();
                f = f * ff->computeElectrostaticForceModule(p.getCharge(), getCharge(), distance);
                addForce(f);
            }
        }
    }
}

void Particle::addStericForce(std::vector<DeferredNonbondedContribution> & deferred)
{
    Vector3f f = Vector3f();
    float distance = 0.0;

    float cutoff = _springnetwork->getStericCutoff();

    bool apply = true;

    const biospring::forcefield::ForceField * ff = _springnetwork->getForceField();
    _springnetwork->getNeighborSearch().steric->for_each_neighbor(*this, [&](size_t neighbor_index) {
        if (_springnetwork->isProbeParticle(neighbor_index))
            return;

        const Particle & p = _springnetwork->getParticle(neighbor_index);

        // See addElectrostaticForce: process each pair once, from the
        // lower-id side when the neighbor is dynamic; a static neighbor never
        // visits any pair on its own, so it is always processed.
        if (p.isDynamic() && neighbor_index < static_cast<size_t>(getId()))
            return;

        apply = true;
        if (_springnetwork->isSpringEnabled())
            // getId() stays signed (used as an "unassigned" sentinel elsewhere,
            // e.g. SpringNetwork::isProbeParticle), so cast explicitly here
            // where it is consumed as a spring-neighbor index (always >= 0
            // once assigned by SpringNetwork).
            apply = !isInSpringNeighbors(static_cast<unsigned>(p.getId()));

        if (apply && p.getId() != getId())
        {
            f = p.getPosition() - getPosition();
            distance = f.norm();
            if (distance < cutoff)
            {
                const float pair_energy =
                    ff->computeStericEnergy(p.getRadius(), getRadius(), p.getEpsilon(), getEpsilon(), distance);
                f.normalize();
                f = f *
                    ff->computeStericForceModule(p.getRadius(), getRadius(), p.getEpsilon(), getEpsilon(), distance);
                addForce(f);

                if (p.isDynamic())
                {
                    // See addElectrostaticForce for why `-f` and half the
                    // energy are exactly what a dynamic neighbor would have
                    // computed on its own.
                    _stericenergy += pair_energy / 2.0f;
                    deferred.push_back({static_cast<unsigned>(p.getId()), -f, pair_energy / 2.0f});
                }
                else
                {
                    // A static neighbor never contributes its own half; see
                    // addElectrostaticForce for why the full energy must be
                    // credited here instead.
                    _stericenergy += pair_energy;
                }
            }
        }
    });
}

// ======================================================================================
// Add forces to probe

float Particle::addStericProbeForce(Particle & probe)
{
    const biospring::forcefield::ForceField * ff = _springnetwork->getForceField();
    Vector3f f = probe.getPosition() - getPosition();
    float distance = f.norm();

    const float pair_energy =
        ff->computeStericEnergy(probe.getRadius(), getRadius(), probe.getEpsilon(), getEpsilon(), distance);
    _stericenergy += pair_energy * 0.5f;
    probe.setStericEnergy(probe.getStericEnergy() + pair_energy * 0.5f);
    f.normalize();
    f = f * ff->computeStericForceModule(probe.getRadius(), getRadius(), probe.getEpsilon(), getEpsilon(), distance);
    addForce(f);
    probe.addForce(-f);
    return pair_energy;
}

float Particle::addElectrostaticProbeForce(Particle & probe)
{
    const biospring::forcefield::ForceField * ff = _springnetwork->getForceField();
    Vector3f f = probe.getPosition() - getPosition();
    float distance = f.norm();

    const float pair_energy = ff->computeElectrostaticEnergy(probe.getCharge(), getCharge(), distance);
    _electrostaticenergy += pair_energy * 0.5f;
    probe.setElectrostaticEnergy(probe.getElectrostaticEnergy() + pair_energy * 0.5f);
    f.normalize();
    f = f * ff->computeElectrostaticForceModule(probe.getCharge(), getCharge(), distance);
    addForce(f);
    probe.addForce(-f);
    return pair_energy;
}

} // namespace spn
} // namespace biospring
