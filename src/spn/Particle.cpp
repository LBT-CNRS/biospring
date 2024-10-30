
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
    Vector3f var = (_force / getMass()) * 0.5f * timestep;
    _velocity = _velocity + var;
    _position = _position + _velocity * timestep;
    _velocity = _velocity + var;
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
    float gridscale = _springnetwork->getGridScale();
    const biospring::grid::PotentialGrid & potentialgrid = _springnetwork->getDensityGrid();

    Vector3f force = potentialgrid.get(getX(), getY(), getZ()).vector;
    force = force * getBurying() * gridscale;

    addForce(force);
}

void Particle::addElectrostaticFieldForce()
{
    const biospring::forcefield::ForceField * ff = _springnetwork->getForceField();
    float gridscale = _springnetwork->getGridScale();
    const biospring::grid::PotentialGrid & potentialgrid = _springnetwork->getPotentialGrid();

    const auto & cell = potentialgrid.get(getX(), getY(), getZ());

    Vector3f force = cell.vector * getCharge() * gridscale;
    addForce(force);

    _electrostaticenergy += ff->computeElectrostaticFieldEnergy(cell.scalar, getCharge());
}

void Particle::addElectrostaticForce()
{
    Vector3f f = Vector3f();
    float distance = 0.0;
    float cutoff = _springnetwork->getElectrostaticCutoff();
    bool apply = true;

    const auto & coulombneighbors = _springnetwork->getNeighborSearch().electrostatic->get_neighbors(*this);
    const biospring::forcefield::ForceField * ff = _springnetwork->getForceField();
    for (auto neighbor_index : coulombneighbors)
    {
        const Particle & p = _springnetwork->getParticle(neighbor_index);
        apply = true;
        if (_springnetwork->isSpringEnabled())
            apply = !isInSpringNeighbors(p.getId());

        if (apply)
        {
            f = p.getPosition() - getPosition();
            distance = f.norm();
            if (distance < cutoff && distance != 0.0)
            {
                _electrostaticenergy += 0.5 * ff->computeElectrostaticEnergy(p.getCharge(), getCharge(), distance);
                f.normalize();
                f = f * ff->computeElectrostaticForceModule(p.getCharge(), getCharge(), distance);
                addForce(f);
            }
        }
    }
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

void Particle::addHydrophobicityForce()
{
    Vector3f f = Vector3f();
    float distance = 0.0;
    float cutoff = _springnetwork->getHydrophobicCutoff();
    bool apply = true;

    if (_springnetwork->isHydrophobicityEnabled())
    {
        const auto & hydrophobicneighbors = _springnetwork->getNeighborSearch().hydrophobic->get_neighbors(*this);
        const biospring::forcefield::ForceField * ff = _springnetwork->getForceField();
        for (auto neighbor_index : hydrophobicneighbors)
        {
            const Particle & p = _springnetwork->getParticle(neighbor_index);
            apply = true;
            if (_springnetwork->isSpringEnabled())
                apply = !isInSpringNeighbors(p.getId());

            if (apply)
            {
                f = p.getPosition() - getPosition();
                distance = f.norm();
                if (distance < cutoff && distance != 0.0)
                {
                    _hydrophobicityenergy +=
                        ff->computeHydrophobicityEnergy(p.getHydrophobicity(), getHydrophobicity(), distance) / 2.0;
                    f.normalize();
                    f = f * ff->computeHydrophobicityForceModule(p.getHydrophobicity(), getHydrophobicity(), distance);
                    addForce(f);
                }
            }
        }
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
            apply = !isInSpringNeighbors(p.getId());
        }

        if (apply && p.getId() != getId())
        {
            f = p.getPosition() - getPosition();
            distance = f.norm();
            _electrostaticenergy += ff->computeElectrostaticEnergy(p.getCharge(), getCharge(), distance);
            f.normalize();
            f = f * ff->computeElectrostaticForceModule(p.getCharge(), getCharge(), distance);
            addForce(f);
        }
    }
}

void Particle::addStericForce()
{
    Vector3f f = Vector3f();
    float distance = 0.0;

    float cutoff = _springnetwork->getStericCutoff();

    bool apply = true;

    const auto & vanderwaalsneighbors = _springnetwork->getNeighborSearch().steric->get_neighbors(*this);
    const biospring::forcefield::ForceField * ff = _springnetwork->getForceField();
    for (auto neighbor_index : vanderwaalsneighbors)
    {
        const Particle & p = _springnetwork->getParticle(neighbor_index);
        apply = true;
        if (_springnetwork->isSpringEnabled())
            apply = !isInSpringNeighbors(p.getId());

        if (apply && p.getId() != getId())
        {
            f = p.getPosition() - getPosition();
            distance = f.norm();
            if (distance < cutoff)
            {
                _stericenergy +=
                    ff->computeStericEnergy(p.getRadius(), getRadius(), p.getEpsilon(), getEpsilon(), distance) / 2.0;
                f.normalize();
                f = f *
                    ff->computeStericForceModule(p.getRadius(), getRadius(), p.getEpsilon(), getEpsilon(), distance);
                addForce(f);
            }
        }
    }
}

// ======================================================================================
// Add forces to probe

void Particle::addStericProbeForce(Particle & probe)
{
    const biospring::forcefield::ForceField * ff = _springnetwork->getForceField();
    Vector3f f = probe.getPosition() - getPosition();
    float distance = f.norm();

    _stericenergy +=
        ff->computeStericEnergy(probe.getRadius(), getRadius(), probe.getEpsilon(), getEpsilon(), distance) / 2.0;
    f.normalize();
    f = f * ff->computeStericForceModule(probe.getRadius(), getRadius(), probe.getEpsilon(), getEpsilon(), distance);
    addForce(f);
    probe.addForce(-f);
}

void Particle::addElectrostaticProbeForce(Particle & probe)
{
    const biospring::forcefield::ForceField * ff = _springnetwork->getForceField();
    Vector3f f = probe.getPosition() - getPosition();
    float distance = f.norm();

    _electrostaticenergy = ff->computeElectrostaticEnergy(probe.getCharge(), getCharge(), distance);
    probe.setElectrostaticEnergy(probe.getElectrostaticEnergy() + _electrostaticenergy);
    f.normalize();
    f = f * ff->computeElectrostaticForceModule(probe.getCharge(), getCharge(), distance);
    addForce(f);
    probe.addForce(-f);
}

} // namespace spn
} // namespace biospring
