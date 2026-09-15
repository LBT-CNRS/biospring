#include "SpringNetwork.h"

#include <fstream>
#include "logging.h"
#include "measure.hpp"

#include "forcefield/constants.hpp"
#include "forcefield/ForceField.h"
#include "forcefield/ForceFieldElectrostaticCoulombAndStericLennardJones_12_6Amber.h"
#include "forcefield/ForceFieldElectrostaticCoulombAndStericLennardJones_8_6Lewitt.h"
#include "forcefield/ForceFieldElectrostaticCoulombAndStericLennardJones_8_6Zacharias.h"
#include "forcefield/ForceFieldElectrostaticCoulombAndStericLinear.h"

#include "IO/CSVSampleWriter.h"
#include "IO/OpenDXReader.h"
#include "IO/PDBTrajectoryWriter.h"
#include "IO/XTCTrajWriter.h"
#include "interactor/Interactor.h"
#ifdef MDDRIVER_SUPPORT
	#include "interactor/mddriver/InteractorMDDriver.h"
#endif
#ifdef FREESASA_SUPPORT
    #include "interactor/freesasa/InteractorFreeSASA.h"
#endif

#include <algorithm>
#include <sstream>
#include <map>
#include <set>
#include <utility>
#include <iostream>
#include <limits>
#include <math.h>
#include <memory>
#include <stdlib.h>
#include <string.h>
#include <utility>
#include <chrono>
#include <thread>

#ifdef OPENMP_SUPPORT
#include <omp.h>
#endif

#include "rigidbody/RigidBodiesManager.h"

namespace biospring
{
namespace spn
{

// The .bi.ff spelling of each family, in DihedralFamilyIndex order. The same
// list BondedForceFieldReader keeps; duplicated rather than shared because that
// reader is a build-time component and this is the runtime.
static constexpr const char * DIHEDRAL_FAMILY_KEYWORDS_RT[SpringNetwork::DIHEDRAL_FAMILY_COUNT] = {
    "PHI", "PSI", "OMEGA", "SIDECHAIN", "PLANARITY", "NUCLEIC_BACKBONE", "NUCLEIC_CHI", "NUCLEIC_SUGAR"};

unsigned SpringNetwork::_currentstructid = 0;

SpringNetwork::~SpringNetwork() {}

void SpringNetwork::_resetEnergies()
{
    _energies.reset();
}

void SpringNetwork::_updateInsertionVector()
{
    {
        _insertionVector->computeVector();
        _insertionVector->computeAngle();
        _insertionVector->computeRollAngle();
    }
}

// Calculates spring forces and applies them to the particles.
// Updates global `_energies.spring` variable.
void SpringNetwork::computeSpringForces()
{
    float springenergy = 0.0f;
    _springForceScratch.resize(_dynamicsprings.size());

#ifdef OPENMP_SUPPORT
#pragma omp parallel for schedule(static)
#endif
    for (size_t i = 0; i < _dynamicsprings.size(); ++i)
    {
        Spring & spring = getSpring(_dynamicsprings[i]);
        _springForceScratch[i] = spring.computeForce(*_ff);
    }

    // Particle forces and energies are accumulated in a deterministic serial
    // pass. The expensive spring evaluation remains parallel, while concurrent
    // writes and thread-count-dependent floating-point reductions are avoided.
    for (size_t i = 0; i < _dynamicsprings.size(); ++i)
    {
        Spring & spring = getSpring(_dynamicsprings[i]);
        const Vector3f & force = _springForceScratch[i];
        spring.getParticle1().addForce(force);
        spring.getParticle2().addForce(-force);
        springenergy += spring.getEnergy();
    }

    _energies.spring = springenergy;
}

// Shared parallel-compute/serial-accumulate loop for the DIHEDRAL
// spring collections -- the exact pattern computeSpringForces
// established: the spring evaluation (the expensive part) runs in parallel
// into _springForceScratch, then a deterministic serial pass applies the
// forces (particles are shared between springs, so they must not be
// written concurrently) and accumulates the energy in a fixed order,
// keeping the result independent of the thread count. Returns the summed
// energy. The scratch buffer is shared with computeSpringForces: the four
// loops run sequentially within one step, and resize never shrinks
// capacity, so no per-step allocation happens either way.
float SpringNetwork::_computeSpringCollectionForces(std::vector<Spring> & springs, bool ignoreDynamicState,
                                                    bool subtractDcOffset, const std::vector<unsigned> * axes)
{
    // The torsional evaluation is not a force on the pair distance at all, so
    // it does not go through computeForce: it produces its own two forces from
    // the energy's exact gradient in the azimuths. It runs serially -- the same
    // particles are shared between springs -- which is what the parallel pass
    // below avoids, and is affordable because it replaces that pass rather than
    // adding to it.
    if (axes != nullptr && _config.dihedral.torsionalonly)
    {
        float torsional = 0.0f;
        for (size_t i = 0; i < springs.size(); ++i)
        {
            const unsigned axisIndex = (*axes)[i];
            if (axisIndex == NO_AXIS || !springs[i].hasTorsionalFrame())
            {
                const Vector3f force = springs[i].computeForce(*_ff, ignoreDynamicState);
                springs[i].getParticle1().addForce(force);
                springs[i].getParticle2().addForce(-force);
                torsional += springs[i].getEnergy();
                continue;
            }
            const float e = computeTorsionalDihedral(springs[i], axisIndex);
            torsional += subtractDcOffset ? e - springs[i].getDcOffset() : e;
        }
        return torsional;
    }

    _springForceScratch.resize(springs.size());

#ifdef OPENMP_SUPPORT
#pragma omp parallel for schedule(static)
#endif
    for (size_t i = 0; i < springs.size(); ++i)
        _springForceScratch[i] = springs[i].computeForce(*_ff, ignoreDynamicState);

    float energy = 0.0f;
    for (size_t i = 0; i < springs.size(); ++i)
    {
        Spring & spring = springs[i];
        const Vector3f & force = _springForceScratch[i];
        const unsigned axisIndex = axes != nullptr ? (*axes)[i] : NO_AXIS;
        if (axes != nullptr && axisIndex != NO_AXIS && _config.dihedral.tangentialonly)
        {
            // Keep only what turns the axis, at the point the torsion term
            // produces it. The axis belongs to the SPRING, not to either
            // endpoint: a ghost hangs off exactly one axis, but a real atom
            // used as an endpoint can serve several, so asking the particle
            // gives an ambiguous answer while asking the spring cannot.
            applyProjectedDihedralForce(spring.getParticle1(), force, axisIndex);
            applyProjectedDihedralForce(spring.getParticle2(), -force, axisIndex);
        }
        else
        {
            spring.getParticle1().addForce(force);
            spring.getParticle2().addForce(-force);
        }
        energy += subtractDcOffset ? spring.getEnergy() - spring.getDcOffset() : spring.getEnergy();
    }
    return energy;
}

// Component of `f` that turns `p` about the ghost axis it belongs to.
//
// A torsion should push a substituent ONLY around its axis. The ring does
// not: measured on ubiquitin, 78 % of what it applies (86 % median) is
// radial or axial, which carries no torque about the axis at all and simply
// deforms bonds and angles. Radial and axial components contribute exactly
// zero to the axial torque, so dropping them keeps the torque exactly.
//
// This used to live in redistributeGhostForces, applied to the force AFTER
// it had been rotated onto the real atom the ghost images. The two are
// algebraically the same -- rotating about the axis maps a tangential
// direction to a tangential direction, so projecting then rotating and
// rotating then projecting agree -- but doing it here makes the filter a
// property of the TORSION TERM rather than of the ghost mechanism, so it
// applies just as well to an endpoint that is a real atom rather than a
// ghost. A quarter of all ring ghosts sit at azimuth 0, where the placement
// is the identity and the ghost is an exact copy of its reference atom; they
// can only be replaced by that atom once the filter no longer hangs off the
// ghost.
//
// A particle with no axis, or sitting on the axis (no lever arm, hence no
// torsional role), gets nothing rather than an arbitrary direction.
Vector3f SpringNetwork::tangentialAboutAxis(const Particle & p, const Vector3f & f, unsigned axisIndex) const
{
    const GhostAxis & axis = _ghostaxes[axisIndex];
    const Vector3f & B = getParticle(axis.anchorBIndex).getPosition();
    Vector3f ahat = getParticle(axis.anchorCIndex).getPosition() - B;
    ahat.normalize();

    const Vector3f rel = p.getPosition() - B;
    const Vector3f radial = rel - ahat * rel.dot(ahat);
    const float rn = radial.norm();
    if (rn <= 1e-6f)
        return Vector3f();

    const Vector3f that = ahat ^ (radial / rn);
    return that * f.dot(that);
}

// Applies one dihedral spring's share of force to one of its endpoints, kept
// tangential, and books what a real atom receives against its axis.
//
// A ghost needs no booking here: its force accumulates on the ghost and
// redistributeGhostForces transfers the total to the real atom it images,
// booking it there. An endpoint that IS a real atom has no such pass -- the
// force lands on it directly -- so without this it would never appear in
// axis.sumAtomForces and the axis reaction would balance against a total that
// is missing it. Momentum would leak, silently, in proportion to how many
// ghosts were replaced by the atom they sit on.
void SpringNetwork::applyProjectedDihedralForce(Particle & p, const Vector3f & f, unsigned axisIndex)
{
    const Vector3f projected = tangentialAboutAxis(p, f, axisIndex);
    p.addForce(projected);

    // A ghost normally needs no booking here: its force accumulates on it and
    // redistributeGhostForces books the total where it lands. Spring-held, no
    // such pass exists -- the ghost keeps what it is given, exactly like a real
    // endpoint -- so it has to be booked here or the axis reaction would
    // balance against a total that is missing it.
    const unsigned particleIndex = static_cast<unsigned>(p.getId());
    if (!_ghostsAreSpringHeld && particleIndex < _isGhost.size() && _isGhost[particleIndex])
        return;

    GhostAxis & axis = _ghostaxes[axisIndex];
    const Vector3f & B = getParticle(axis.anchorBIndex).getPosition();
    axis.sumAtomForces += projected;
    axis.sumAtomTorquesAboutB += (p.getPosition() - B) ^ projected;
}

// One dihedral spring, evaluated on IDEALISED positions: each endpoint keeps
// the azimuth it really has about the axis, but is put back at the radius and
// axial offset it had when the model was built. Their distance is then
//
//     d^2 = rho1^2 + rho2^2 + (z1 - z2)^2 - 2*rho1*rho2*cos(dpsi)
//
// the same closed form the ring construction is derived from, and a function of
// the azimuth difference ALONE.
//
// That is what makes this conservative where the tangential filter is not. The
// filter throws away the radial part of a force while the energy keeps its
// radial dependence, so an atom moving radially changes the energy with no work
// done against it -- measured on ubiquitin, 124152 kJ/mol of kinetic energy
// conjured out of a standing start. Here the energy cannot see radial motion at
// all, so there is nothing to throw away and nothing to leak.
//
// The force is the exact gradient of that energy in the azimuths: grad(psi) is
// t/rho with rho the atom's REAL radius, which is the only place the real
// radius still enters -- converting a torque back into a force, not deciding
// the energy. The reaction on the two axis atoms is left to the same
// closed-form pass that already serves the rings.
float SpringNetwork::computeTorsionalDihedral(Spring & spring, unsigned axisIndex)
{
    const GhostAxis & axis = _ghostaxes[axisIndex];
    const Vector3f & B = getParticle(axis.anchorBIndex).getPosition();
    Vector3f ahat = getParticle(axis.anchorCIndex).getPosition() - B;
    ahat.normalize();

    Particle & p1 = spring.getParticle1();
    Particle & p2 = spring.getParticle2();

    Vector3f r1 = p1.getPosition() - B;
    r1 = r1 - ahat * r1.dot(ahat);
    Vector3f r2 = p2.getPosition() - B;
    r2 = r2 - ahat * r2.dot(ahat);
    const float n1 = r1.norm();
    const float n2 = r2.norm();
    if (n1 <= 1e-6f || n2 <= 1e-6f)
        return 0.0f; // on the axis: no azimuth, hence no torsional role
    r1 = r1 / n1;
    r2 = r2 / n2;

    const Vector3f t1 = ahat ^ r1;
    const Vector3f t2 = ahat ^ r2;
    const float cos_dpsi = r1.dot(r2);
    const float sin_dpsi = t1.dot(r2); // = ahat . (r1 x r2)

    const float rho1 = spring.getRho1();
    const float rho2 = spring.getRho2();
    const float dz = spring.getZ1() - spring.getZ2();
    const float d = std::sqrt(std::max(rho1 * rho1 + rho2 * rho2 + dz * dz - 2.0f * rho1 * rho2 * cos_dpsi, 1e-12f));

    const float k = spring.getStiffness();
    const float ext = d - spring.getEquilibrium();
    // GLOBAL_SPRING_FORCE_CONVERT, and springscale, are what Spring::computeForce
    // applies on the way out: forces are integrated in Da.A.fs-2 while k is a
    // kJ.mol-1.A-2 number. Bypassing computeForce means applying them here --
    // without it the force is 1e4 too large, which is exactly what two springs
    // needed to reach 1.9e8 kJ/mol of kinetic energy.
    const float dEdpsi = _ff->getSpringScale() * forcefield::GLOBAL_SPRING_FORCE_CONVERT * k * ext *
                         (rho1 * rho2 * sin_dpsi) / d;

    // dpsi is measured from endpoint 1 to endpoint 2, so it grows with psi2 and
    // shrinks with psi1.
    //
    // The lever arm is the REFERENCE radius, not the current one. grad(psi) is
    // t/rho_real exactly, but that is singular: let an atom drift towards the
    // axis and the force diverges, which feeds back and blew the network up at
    // 4e9 kJ/mol. The reference radius is the arm the atom is meant to turn on
    // -- it is the arm a ghost would have had -- and it keeps the force bounded
    // by the same token that keeps the energy free of radial dependence.
    const float arm1 = std::max(rho1, 0.1f);
    const float arm2 = std::max(rho2, 0.1f);
    applyAxisBookedForce(p1, t1 * (dEdpsi / arm1), axisIndex);
    applyAxisBookedForce(p2, t2 * (-dEdpsi / arm2), axisIndex);
    return 0.5f * k * ext * ext;
}

// Adds a force to a dihedral endpoint and books it against its axis, so the
// closed-form reaction has the whole total to balance against.
void SpringNetwork::applyAxisBookedForce(Particle & p, const Vector3f & f, unsigned axisIndex)
{
    p.addForce(f);
    GhostAxis & axis = _ghostaxes[axisIndex];
    const Vector3f & B = getParticle(axis.anchorBIndex).getPosition();
    axis.sumAtomForces += f;
    axis.sumAtomTorquesAboutB += (p.getPosition() - B) ^ f;
}

void SpringNetwork::resetGhostAxisSums()
{
    for (GhostAxis & axis : _ghostaxes)
    {
        axis.sumGhostForces = Vector3f();
        axis.sumGhostTorquesAboutB = Vector3f();
        axis.sumAtomForces = Vector3f();
        axis.sumAtomTorquesAboutB = Vector3f();
    }
}

// Gives every dihedral spring endpoint its axis, including the ones that are
// real atoms rather than ghosts.
//
// Both endpoints of a dihedral spring hang off the same axis -- that is what
// the ring construction means -- so an endpoint with no axis of its own takes
// the one its partner knows. No file format carries this: addGhostParticle
// knows the axis because it created the ghost, and this pass propagates it
// across each spring afterwards.
//
// Run once, lazily, because it has to happen after every ghost is registered
// AND every dihedral spring is built, and no single construction path
// guarantees an ordering of those two.
// Linear search is fine: this runs at build time, and the number of distinct
// axes is small next to the number of springs hanging off them (example 072:
// 6484 ghosts, a few hundred axes).
unsigned SpringNetwork::findOrCreateGhostAxis(unsigned anchorBIndex, unsigned anchorCIndex)
{
    for (unsigned i = 0; i < _ghostaxes.size(); ++i)
        if (_ghostaxes[i].anchorBIndex == anchorBIndex && _ghostaxes[i].anchorCIndex == anchorCIndex)
            return i;
    _ghostaxes.push_back(GhostAxis{anchorBIndex, anchorCIndex, Vector3f(), Vector3f(), Vector3f(), Vector3f()});
    return static_cast<unsigned>(_ghostaxes.size() - 1);
}

void SpringNetwork::bindDihedralEndpointsToAxes()
{
    if (_dihedralAxesBound)
        return;
    _dihedralAxesBound = true;

    // The axis is a property of the SPRING. Both its endpoints hang off the
    // same one -- that is what a ring is -- so whichever endpoint is a ghost
    // names it. Asking the particle instead would be ambiguous: a ghost has
    // exactly one axis, but a real atom used as an endpoint (every azimuth-0
    // placement is its own reference atom) can serve several.
    //
    // Run once, lazily, because it must happen after every ghost is registered
    // AND every dihedral spring built, and no construction path guarantees an
    // order between those two.
    auto ghost_axis = [this](const Particle & p) {
        const unsigned i = static_cast<unsigned>(p.getId());
        if (i >= _isGhost.size() || !_isGhost[i])
            return NO_AXIS;
        return i < _axisOfParticle.size() ? _axisOfParticle[i] : NO_AXIS;
    };

    unsigned unresolved = 0;
    for (unsigned family = 0; family < DIHEDRAL_FAMILY_COUNT; ++family)
    {
        std::vector<Spring> & springs = _dihedralsprings[family];
        _dihedralAxis[family].assign(springs.size(), NO_AXIS);
        for (size_t i = 0; i < springs.size(); ++i)
        {
            // The spring's own axis first, when it carries one: a spring
            // between two real substituents has no ghost to ask, and asking
            // its endpoints would answer for some other torsion they happen
            // to serve. Ring springs carry none and keep answering through
            // their ghost, exactly as before.
            unsigned a = NO_AXIS;
            if (springs[i].hasAxis())
                a = findOrCreateGhostAxis(springs[i].getAxisB(), springs[i].getAxisC());
            if (a == NO_AXIS)
            {
                const unsigned a1 = ghost_axis(springs[i].getParticle1());
                const unsigned a2 = ghost_axis(springs[i].getParticle2());
                a = a1 != NO_AXIS ? a1 : a2;
            }
            _dihedralAxis[family][i] = a;
            if (a == NO_AXIS)
                ++unresolved;
        }
    }

    if (unresolved > 0)
        logging::warning("SpringNetwork: %u dihedral spring(s) have a real atom at BOTH ends, so no ghost "
                         "names their axis; they are applied unfiltered. Keep one ghost per ring spring.",
                         unresolved);
}

void SpringNetwork::computeDihedralForces()
{
    bindDihedralEndpointsToAxes();

    float dihedralenergy = 0.0f;

    auto accumulate = [&](std::vector<Spring> & springs, const std::vector<unsigned> * axes) {
        dihedralenergy += _computeSpringCollectionForces(springs, /*ignoreDynamicState=*/true,
                                                         /*subtractDcOffset=*/true, axes);
    };

    const bool enabled[DIHEDRAL_FAMILY_COUNT] = {
        isDihedralPhiEnabled(),             isDihedralPsiEnabled(),        isDihedralOmegaEnabled(),
        isDihedralChiEnabled(),             isDihedralPlanarityEnabled(),  isDihedralNucleicBackboneEnabled(),
        isDihedralNucleicChiEnabled(),      isDihedralNucleicSugarEnabled()};
    for (unsigned family = 0; family < DIHEDRAL_FAMILY_COUNT; ++family)
        if (enabled[family])
            // Both treatments need the axis: one to project onto it, the other
            // to measure azimuths about it.
            accumulate(_dihedralsprings[family],
                       (_config.dihedral.tangentialonly || _config.dihedral.torsionalonly)
                           ? &_dihedralAxis[family]
                           : nullptr);

    _energies.dihedral = dihedralenergy + computeTorsionForces();
}

// Calculate forces that apply on dynamic particles.
void SpringNetwork::computeParticleForces()
{
    float electrostatic_energy = 0.0f;
    float steric_energy = 0.0f;
    float imp_energy = 0.0f;
    float hydrophobic_energy = 0.0f;

    _resizeNonbondedPairScratch();

#ifdef OPENMP_SUPPORT
#pragma omp parallel for schedule(static)
#endif
    for (size_t i = 0; i < _dynamicparticules.size(); ++i)
    {
        Particle & p = getParticle(_dynamicparticules[i]);

        if (isElectrostaticEnabled())
        {
            if (isElectrostaticCoulombEnabled() && p.isCharged() && _nsearch.electrostatic)
                p.addElectrostaticForce(_electrostaticPairScratch[i]);
            if (isElectrostaticFieldEnabled())
                p.addElectrostaticFieldForce();
        }

        if (isDensityGridEnabled())
            p.addDensityFieldForce();

        if (isStericEnabled())
            p.addStericForce(_stericPairScratch[i]);

        if (isViscosityEnabled())
            p.applyViscosity(getViscosity());

        if (isIMPEnabled())
            p.addIMPForce();

        if (isHydrophobicityEnabled() && p.isHydrophobic() && _nsearch.hydrophobic)
            p.addHydrophobicityForce(_hydrophobicPairScratch[i]);
    }

    // Applies the deferred "other side" of each unique nonbonded pair
    // (Newton's third law) serially, since two threads may have deferred a
    // contribution to the same target particle. Must run before the force
    // read by rigid-body torque aggregation and setPreviousForce() below, and
    // before summing per-particle energies, since it feeds both.
    _applyNonbondedPairScratch(_stericPairScratch, steric_energy);
    _applyNonbondedPairScratch(_electrostaticPairScratch, electrostatic_energy);
    _applyNonbondedPairScratch(_hydrophobicPairScratch, hydrophobic_energy);

    // Sum per-particle energies in particle order to keep results reproducible
    // across OpenMP thread counts.
    for (const unsigned particle_id : _dynamicparticules)
    {
        const Particle & p = getParticle(particle_id);
        electrostatic_energy += p.getElectrostaticEnergy();
        steric_energy += p.getStericEnergy();
        imp_energy += p.getIMPEnergy();
        hydrophobic_energy += p.getHydrophobicityEnergy();
    }

    // The probe is shared by every particle, therefore probe interactions must
    // not mutate it from an OpenMP loop. It is integrated exactly once per step.
    if (isProbeEnabled())
    {
        _probeparticule.resetForce();

        for (const unsigned particle_id : _dynamicparticules)
        {
            Particle & p = getParticle(particle_id);
            if (isProbeStericEnabled())
                steric_energy += p.addStericProbeForce(_probeparticule);
            if (isProbeElectrostaticEnabled())
                electrostatic_energy += p.addElectrostaticProbeForce(_probeparticule);
        }

        if (isProbeElectrostaticFieldEnabled())
        {
            const float previous_energy = _probeparticule.getElectrostaticEnergy();
            _probeparticule.addElectrostaticFieldForce();
            electrostatic_energy += _probeparticule.getElectrostaticEnergy() - previous_energy;
        }

        if (isViscosityEnabled())
            _probeparticule.applyViscosity(getViscosity());

        _probeparticule.IntegrateEuler(getTimeStep());
        _energies.kinetic += _probeparticule.getKineticEnergy();
        _syncProbeParticle();
    }

    // Rigid-body accumulators are shared between their particles. Aggregate
    // them serially after all per-particle forces are complete.
    for (const unsigned particle_id : _dynamicparticules)
    {
        Particle & p = getParticle(particle_id);
        if (isRigidBodyEnabled() && p.isRigid() && !isImpalaSamplingEnabled() && !isMonteCarloEnabled())
            rigidbody::RigidBody::computeParticleForceAndTorque(p);
        p.setPreviousForce();
    }

    _energies.electrostatic = electrostatic_energy;
    _energies.steric = steric_energy;
    _energies.imp = imp_energy;
    _energies.hydrophobic = hydrophobic_energy;
}

// Update the positions of the particles.
// Update global `_energies.kinetic` variable.
void SpringNetwork::updateParticlePositions()
{
    float kinetic_energy_particle = 0.0;
#ifdef OPENMP_SUPPORT
#pragma omp parallel default(shared)
#endif
    {
#ifdef OPENMP_SUPPORT
#pragma omp for reduction(+ : kinetic_energy_particle) schedule(static)
#endif
        // i stays a signed int: MSVC only supports OpenMP 2.0, which requires
        // a signed loop counter for #pragma omp parallel for.
        for (int i = 0; i < (int)_dynamicparticules.size(); i++)
        {
            Particle & p = getParticle(_dynamicparticules[static_cast<size_t>(i)]);
            if (p.isRigid())
                rigidbody::RigidBody::integrateParticleVelocity(p, i, getTimeStep());
            else
                p.IntegrateEuler(getTimeStep());

            // Check if position exploses, if one of float is NaN
            float x = p.getPosition().getX();
            float y = p.getPosition().getY();
            float z = p.getPosition().getZ();
            if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z))
            {
                logging::die("Found non-finite position for particle %d.", p.getId());
            }

            kinetic_energy_particle += p.getKineticEnergy();
            p.resetForce();
        } // omp for loop
    }     // omp parallel
    _energies.kinetic += kinetic_energy_particle;

    // The spatial grids must follow particle motion. Rebuild them once here,
    // after all particle positions have been integrated for the current step.
    _markNeighborSearchesDirty();
    _updateNeighborSearches();
}

/// @brief Compute particles force and, if activated, springs forces.
void SpringNetwork::computeForces()
{
    if (isSpringEnabled())
    {
        computeSpringForces();
        computeDihedralForces();
    }
    computeParticleForces();
}

void SpringNetwork::computeStep()
{
    idleRun();
    _meanConstraintsDistances = 0.0;

    // Before any force is produced, not inside redistributeGhostForces: with
    // the tangential filter on, a dihedral endpoint that is a real atom
    // contributes to its axis's totals while the spring loop runs, which is
    // earlier than the redistribution pass that used to clear them.
    resetGhostAxisSums();

    computeForces();
    redistributeGhostForces();
    applyGhostDamping();

    if (isConstraintEnabled())
        applyConstraints();

    if (isRigidBodyEnabled())
        rigidbody::RigidBodiesManager::SolveRigidBodiesDynamic();

    updateParticlePositions();
    updateGhostPositions();

    if (isInsertionVectorEnabled())
        _updateInsertionVector();
}

void SpringNetwork::run()
{
    initRun();

#ifdef MDDRIVER_SUPPORT
    // getInteractorInstance returns null when no interactor of that type was
    // registered, which is every caller that drives a SpringNetwork directly
    // rather than through biospring-cli -- the OpenCL parity tests, for one.
    // Dereferencing it unconditionally segfaulted all of them as soon as
    // MDDriver support was compiled in, and only then, which is why it went
    // unnoticed: the builds those tests were written against had it off.
    // The wait belongs with the announcement: it exists to let the server
    // finish opening its port, and there is no port without an interactor.
    if (interactor::InteractorMDDriver * interactormddriver =
            getInteractorInstance<interactor::InteractorMDDriver>())
    {
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        logging::info("    MDDriver parameters:");
        logging::info("      port: %d is open for connection.", interactormddriver->getPort());
    }
#endif

    if (isRigidBodyEnabled())
    {
        // Set all structure to rigid.
        rigidbody::RigidBodiesManager::InitRigidBodies(this, getDynamicParticles());

        // If write csv and automatic sampling enabled, do not write relative 
        // to the option csvsampling.frequency (set to a big number)
        // The information will be written at the end of every protein rotation
        // around the insertion vector axis
        if (isImpalaSamplingEnabled() && _config.csvsample.enable)
            _config.csvsample.frequency = 1000000;
    }

    while (!isEnd())
    {
        computeStep();
        if (_isTimeToLogData())
        {
            _updateFrameRate();
            _displayFrameData();
        }
    }

    // Stop measuring time and calculate the elapsed time.
    _profiler["main"].stop();
    float elapsed = _profiler["main"].elapsed_seconds();
    logging::info("Total time measured: %5.2f seconds.", elapsed);
    logging::info("Average framerate: %5.2f sec-1.", getMaxIteration() / elapsed);

    endRun();
}

void SpringNetwork::initRun()
{
    // Starts measuring time.
    _profiler["main"].start();
    _profiler["samplerate"].start();

    _nbiter = 0;

    for (Interactor* interactor : getInteractors()) 
    {
        if (interactor != nullptr)
        {
            interactor->startInteractionThread();
        }
    }
}

void SpringNetwork::endRun()
{
    for (Interactor* interactor : getInteractors()) 
    {
        if (interactor != nullptr)
            interactor->stopInteractionThread();
    }
    for (Interactor* interactor : getInteractors())
    {
        if (interactor != nullptr)
            interactor->waitForInteractionThread();
    }
}

void SpringNetwork::idleRun()
{
    while (_pause)
        std::this_thread::sleep_for(std::chrono::milliseconds(10));

    for (Interactor* interactor : getInteractors()) 
    {
        if (interactor != nullptr)
        {
            interactor->syncSystemStateData();
        }
    }

    _writeNextStep();

    _resetEnergies();

    _nbiter++;

    if (_hasReachedEndOfRun())
        setEnd(true);
}

void SpringNetwork::_updateFrameRate()
{
    _framerate = getSampleRate() / _profiler["samplerate"].elapsed_seconds();
    _profiler["samplerate"].reset();
}

void SpringNetwork::_displayFrameData()
{
    logging::info("Step: %5d", _nbiter);
    logging::info("Framerate: %5.2f", _framerate);
    logging::info("Kinetic energy: %5.2f kJ.mol-1", _energies.kinetic);
    if (isSpringEnabled())
    {
        logging::info("Spring energy: %5.2f kJ.mol-1", _energies.spring);
        logging::info("Dihedral energy: %5.2f kJ.mol-1", _energies.dihedral);
    }
    if (isElectrostaticEnabled())
        logging::info("Electrostatic energy: %5.2f kJ.mol-1", _energies.electrostatic);
    if (isStericEnabled())
        logging::info("Steric energy: %5.2f kJ.mol-1", _energies.steric);
    if (isIMPEnabled())
        logging::info("IMP energy: %5.2f kJ.mol-1", _energies.imp);
    if (isHydrophobicityEnabled())
        logging::info("Hydrophobic energy: %5.2f kJ.mol-1", _energies.hydrophobic);
    if (isInsertionVectorEnabled())
    {
        logging::info("Insertion angle: %5.2lf °", _insertionVector->getAngle());
        logging::info("Roll angle: %5.2lf °", _insertionVector->getRollAngle());
        logging::info("Insertion depth: %5.2lf.", _insertionVector->getInsertionDepth());
    }
    if (isConstraintEnabled())
        logging::info("Constraints mean distance: %5.2f A", _meanConstraintsDistances);
#ifdef FREESASA_SUPPORT
    // Only when a FreeSASA interactor was actually registered. Printed
    // unconditionally it reported 0 on every run that has none, and on the
    // runs that do have one but are not dynamic it reports a figure computed
    // once at setup -- so say which, rather than letting a frozen number look
    // like a live measurement.
    if (getInteractorInstance<interactor::InteractorFreeSASA>())
        logging::info("Total SASA (%s): %5.2lf A2",
                      _freesasaState.isDynamic ? "updated" : "computed once at setup",
                      _freesasaState.sasaTotal);
#endif

}

// _nbiter stays a signed int: it is compared against getMaxIteration()
// (config's simulation.nbsteps, whose -1 means "infinite run"), and logged
// with %d in several places. Cast explicitly here (always >= 0 in practice,
// only incremented from 0).
void SpringNetwork::_writeNextStep() { _trajectories.write_step(static_cast<size_t>(_nbiter)); }

void SpringNetwork::writeNextStepNow() { _trajectories.write_step(); }

std::vector<Particle>::const_reference SpringNetwork::getParticleFromId(unsigned id) const
{
    for (unsigned i = 0; i < _particles.size(); i++)
    {
        if (_particles[i].getId() == static_cast<int>(id))
            return _particles[i];
    }
    throw std::out_of_range("SpringNetwork::getParticleFromId: Particle id not found.");
}

std::vector<Particle>::reference SpringNetwork::getParticleFromId(unsigned id)
{
    for (unsigned i = 0; i < _particles.size(); i++)
    {
        if (_particles[i].getId() == static_cast<int>(id))
            return _particles[i];
    }
    throw std::out_of_range("SpringNetwork::getParticleFromId: Particle id not found.");
}

void SpringNetwork::getParticlePosition(unsigned i, float position[3]) const
{
    const Particle & p = getParticle(i);
    position[0] = p.getX();
    position[1] = p.getY();
    position[2] = p.getZ();
}

void SpringNetwork::setForce(unsigned i, float force[3])
{
    Vector3f f(force[0], force[1], force[2]);
    getParticle(i).addForce(f);
}

// =====================================================================================
//
// Modification Methods.
//
// Should be used only to convert `Topology` to `SpringNetwork`.
//
// =====================================================================================

void SpringNetwork::addSpring(unsigned id1, unsigned id2, float equilibrium, float stiffness)
{
    if (id1 >= _particles.size() || id2 >= _particles.size())
        throw std::out_of_range("SpringNetwork::addSpring: particle index out of range");
    if (id1 == id2)
        throw std::invalid_argument("SpringNetwork::addSpring: a spring requires two distinct particles");

    // If the spring does not already exist.
    if (!_particles[id1].isInSpringNeighbors(id2))
    {
        // setup() may add one probe particle after springs have been created.
        // Reserve that slot before references to particles are stored in Spring.
        if (_springs.empty() && _particles.capacity() == _particles.size())
            _particles.reserve(_particles.size() + 1);

        Spring * old_storage = _springs.data();
        Particle & p1 = _particles[id1];
        Particle & p2 = _particles[id2];
        _springs.emplace_back(p1, p2, equilibrium, stiffness);
        _springs.back().setId(static_cast<unsigned>(_springs.size() - 1));

        if (old_storage != _springs.data())
            _rebuildSpringNeighbors();
        else
        {
            p1.addToSpringNeighbors(id2, &_springs.back());
            p2.addToSpringNeighbors(id1, &_springs.back());
        }

        if (p1.isStatic() && p2.isStatic())
            addStaticSpring(_springs.back().getId());
        else
            addDynamicSpring(_springs.back().getId());
    }
}

// Shared implementation behind addDihedralSpring:
// unlike addSpring, a ghost spring is always a new addition (no
// spring-neighbour check), never registered in the particles'
// spring-neighbour map (see the header comment on addDihedralSpring for
// why), and not split into static/dynamic subsets (always fully iterated in
// computeDihedralForces).
static void addDihedralSpringTo(std::vector<Spring> & collection, std::vector<Particle> & particles, unsigned id1,
                                unsigned id2, float equilibrium, float stiffness, float dcOffset)
{
    if (id1 >= particles.size() || id2 >= particles.size())
        throw std::out_of_range("SpringNetwork::addDihedral*Spring: particle index out of range");
    if (id1 == id2)
        throw std::invalid_argument("SpringNetwork::addDihedral*Spring: a spring requires two distinct particles");

    collection.emplace_back(particles[id1], particles[id2], equilibrium, stiffness);
    collection.back().setId(static_cast<unsigned>(collection.size() - 1));
    collection.back().setDcOffset(dcOffset);
}

void SpringNetwork::addDihedralSpring(unsigned family, unsigned id1, unsigned id2, float equilibrium, float stiffness,
                                      float dcOffset, unsigned axisB, unsigned axisC)
{
    if (family >= DIHEDRAL_FAMILY_COUNT)
        throw std::out_of_range("SpringNetwork::addDihedralSpring: dihedral family index out of range");
    addDihedralSpringTo(_dihedralsprings[family], _particles, id1, id2, equilibrium, stiffness, dcOffset);
    _dihedralsprings[family].back().setAxis(axisB, axisC);
}

void SpringNetwork::updateSpringState(unsigned id, bool isStatic) {
    if (isStatic) {
        removeDynamicSpring(id);
        addStaticSpring(id);
    } else {
        removeStaticSpring(id);
        addDynamicSpring(id);
    }
}

void SpringNetwork::addParticle(const Particle & source)
{
    if (!_springs.empty())
        throw std::logic_error("SpringNetwork::addParticle: particles must be added before springs");

    Particle p = Particle(source);
    p.setSpringNetwork(this);
    // Particle::getId()/setId() stay signed (used as an "unassigned" sentinel
    // for the probe particle, see isProbeParticle()), so cast explicitly at
    // this array-index assignment and below where it's consumed as an
    // unsigned index (always >= 0 once assigned here).
    p.setId(static_cast<int>(_particles.size()));

    p.setInternalStructId(_structid);

    if (p.isStatic())
        addStaticParticle(static_cast<unsigned>(p.getId()));
    else
        addDynamicParticle(static_cast<unsigned>(p.getId()));

    if (p.isCharged())
        _chargedparticules.push_back(static_cast<unsigned>(p.getId()));

    if (p.isHydrophobic())
        _hydrophobicparticules.push_back(static_cast<unsigned>(p.getId()));

    _particles.push_back(p);
    _initparticles.push_back(p);
    _markNeighborSearchesDirty();
}

unsigned SpringNetwork::addGhostParticle(unsigned placementValue, unsigned anchorBIndex, unsigned anchorCIndex,
                                         unsigned anchorRefIndex, float r, float theta_deg, float delta_deg)
{
    const GhostPlacement placement = static_cast<GhostPlacement>(placementValue);

    const float delta_rad_init = delta_deg * static_cast<float>(M_PI) / 180.0f;
    const Vector3f position =
        GhostParticle::computePositionByRotation(
                  getParticle(anchorBIndex).getPosition(), getParticle(anchorCIndex).getPosition(),
                  getParticle(anchorRefIndex).getPosition(), std::cos(delta_rad_init), std::sin(delta_rad_init));

    Particle p;
    p.setPosition(position);
    p.setStatic(true);
    p.setMass(0.0f);
    addParticle(p);

    const unsigned ownIndex = static_cast<unsigned>(_particles.size() - 1);

    // Every ghost of one ring hangs off the same axis; find or create that
    // axis's accumulator (see GhostAxis). Linear search is fine: this runs
    // once at build time, and the number of distinct axes is small next to
    // the number of ghosts (example 072: 6484 ghosts, a few hundred axes).
    const unsigned axisIndex = findOrCreateGhostAxis(anchorBIndex, anchorCIndex);

    // The tangential filter needs to reach an endpoint's axis from the
    // endpoint alone (see tangentialAboutAxis), so record it per particle
    // rather than per ghost binding: a dihedral spring endpoint need not be
    // a ghost for the filter to apply to it.
    if (_axisOfParticle.size() <= ownIndex)
        _axisOfParticle.resize(ownIndex + 1, NO_AXIS);
    _axisOfParticle[ownIndex] = axisIndex;
    if (_isGhost.size() <= ownIndex)
        _isGhost.resize(ownIndex + 1, false);
    _isGhost[ownIndex] = true;

    const float delta_rad = delta_deg * static_cast<float>(M_PI) / 180.0f;
    _ghostparticles.push_back(GhostParticleBinding{ownIndex, anchorBIndex, anchorCIndex, anchorRefIndex, r, theta_deg,
                                                   delta_deg, std::cos(delta_rad), std::sin(delta_rad), axisIndex,
                                                   placement});
    return ownIndex;
}

// Anchors are heavily shared -- on example 072, 6484 ghosts hang off well
// under 2000 distinct anchors, dozens on a single one -- so the anchors
// cannot be written concurrently. Same split as computeSpringForces: the
// expensive part (the placement Jacobian, one per ghost) runs in parallel
// into a scratch buffer, then a deterministic serial pass accumulates,
// which also keeps the result independent of the thread count.
void SpringNetwork::redistributeGhostForces()
{
    _ghostForceScratch.resize(_ghostparticles.size());

    // Pass 1, parallel: rotate each ghost's force back onto the real atom
    // it images. That single Rodrigues application IS the exact Jacobian
    // transpose for this placement (see GhostParticle.h), so nothing is
    // approximated and no frame is built.
#ifdef OPENMP_SUPPORT
#pragma omp parallel for schedule(static)
#endif
    for (size_t i = 0; i < _ghostparticles.size(); ++i)
    {
        const GhostParticleBinding & binding = _ghostparticles[i];
        const Vector3f & B = getParticle(binding.anchorBIndex).getPosition();
        const Vector3f & C = getParticle(binding.anchorCIndex).getPosition();
        GhostForceContribution & out = _ghostForceScratch[i];
        out.F_Ref = GhostParticle::rotateForceToAtom(B, C, getParticle(binding.ownIndex).getForce(),
                                                     binding.cos_delta, binding.sin_delta);
    }

    // Pass 2, serial: apply to the real atoms (heavily shared, so not
    // concurrently writable) while accumulating each axis's force and
    // torque totals.
    // Totals are cleared once a step in computeStep, before any force exists:
    // a real dihedral endpoint books into them while the spring loop runs,
    // which is earlier than this pass.
    for (size_t i = 0; i < _ghostparticles.size(); ++i)
    {
        const GhostParticleBinding & binding = _ghostparticles[i];
        const GhostForceContribution & in = _ghostForceScratch[i];
        Particle & ghost = getParticle(binding.ownIndex);

        Particle & ref = getParticle(binding.anchorRefIndex);
        GhostAxis & axis = _ghostaxes[binding.axisIndex];
        const Vector3f & B = getParticle(axis.anchorBIndex).getPosition();

        // What acted on the ghost has to be transferred, not cancelled, so
        // both sides of the balance are accumulated.
        axis.sumGhostForces += ghost.getForce();
        axis.sumGhostTorquesAboutB += (ghost.getPosition() - B) ^ ghost.getForce();

        // A torsion should push a substituent ONLY around its axis. The ring
        // does not: measured on ubiquitin, 78 % of what it applies (86 %
        // median) is radial or axial, which carries no torque about the axis
        // at all and simply deforms bonds and angles -- the reason the mesh
        // has to be held at a stiffness of 8000. Projecting onto the
        // tangential direction keeps the torque exactly, since radial and
        // axial components contribute exactly zero to it, and drops the rest.
        // No projection here any more: the tangential filter is applied where
        // the dihedral spring produces its force (see tangentialAboutAxis),
        // so what arrives is already tangential and rotating it onto the real
        // atom keeps it tangential. Same numbers, but the filter no longer
        // depends on the endpoint being a ghost.
        const Vector3f & F_ref = in.F_Ref;

        ref.addForce(F_ref);
        axis.sumAtomForces += F_ref;
        axis.sumAtomTorquesAboutB += (ref.getPosition() - B) ^ F_ref;

        ghost.setForce(Vector3f(0.0f, 0.0f, 0.0f));
    }

    // Pass 3: one closed-form reaction per axis, restoring global force and
    // torque balance without ever differentiating the placement. The torsional
    // evaluation needs it just as much: its two forces turn on different lever
    // arms, so they do not cancel, and without the reaction the same network
    // ends 3 ps at 20031 kJ/mol of kinetic energy instead of 29.6.
    for (const GhostAxis & axis : _ghostaxes)
    {
        Vector3f F_B, F_C;
        // Filtering changes what the reaction has to balance against, and
        // getting this wrong would defeat the whole point. Unfiltered, the
        // reaction TRANSFERS the ghosts' force and torque onto the four real
        // atoms. Filtered, the discarded radial and axial parts carry no
        // torque about the axis but do carry torque about the point B, so
        // balancing against the ghost totals would hand exactly what was
        // removed back to B and C -- the leak relocated, not gone. Balancing
        // the filtered set against zero instead makes it self-cancelling in
        // force and torque, which is precisely AMBER's structure; and since
        // F_B and F_C act on the axis itself they cannot disturb the axial
        // torque the substituents carry.
        const Vector3f zero;
        GhostParticle::redistributeAxisReaction(
            getParticle(axis.anchorBIndex).getPosition(), getParticle(axis.anchorCIndex).getPosition(),
            _config.dihedral.tangentialonly ? zero : axis.sumGhostForces,
            _config.dihedral.tangentialonly ? zero : axis.sumGhostTorquesAboutB,
            axis.sumAtomForces, axis.sumAtomTorquesAboutB, F_B, F_C);
        getParticle(axis.anchorBIndex).addForce(F_B);
        getParticle(axis.anchorCIndex).addForce(F_C);
    }
}

// Embarrassingly parallel, unlike redistributeGhostForces: each iteration
// writes only its own ghost, and a ghost is never itself an anchor (checked
// on example 072: no anchor slot points at a ghost), so no
// iteration can depend on another's result and no write is shared.
void SpringNetwork::updateGhostPositions()
{
#ifdef OPENMP_SUPPORT
#pragma omp parallel for schedule(static)
#endif
    for (size_t i = 0; i < _ghostparticles.size(); ++i)
    {
        const GhostParticleBinding & binding = _ghostparticles[i];
        const Vector3f & B = getParticle(binding.anchorBIndex).getPosition();
        const Vector3f & C = getParticle(binding.anchorCIndex).getPosition();
        getParticle(binding.ownIndex)
            .setPosition(GhostParticle::computePositionByRotation(
                B, C, getParticle(binding.anchorRefIndex).getPosition(), binding.cos_delta, binding.sin_delta));
    }
}

// Turns every ring ghost into an ordinary dynamical particle held by springs,
// instead of an algebraic virtual site re-derived from its anchors at every
// step (dihedral.ghostsprings).
//
// The ring geometry is untouched. The ghost keeps the position
// addGhostParticle already gave it -- the same rotation of its reference atom
// about the B->C axis, at that atom's own radius and axial offset -- and every
// d0/k the generator calibrated is used exactly as generated. What changes is
// only what holds it there afterwards: springs to its rigid body instead of a
// formula re-applied every step.
//
// Three anchor springs per ghost, to B, C and the reference atom, each at the
// length it has at load. Those are precisely the three points the algebraic
// placement read, so the spring network is given the same information; what it
// is not given is the guarantee, and that is the whole of the trade.
//
// Their one blind spot is a ghost lying IN the B/C/reference plane: three
// distances to three coplanar points leave the mirror image just as valid. An
// n=2 ring is exactly that case -- its two b-side azimuths are 0 and 180, the
// first is the reference atom itself and the second is in the plane. The
// intra-ring chords close it, being purely tangential, which is the direction
// the three anchors constrain least.
// Records, for every dihedral spring that knows its axis, the radius and axial
// offset each endpoint has in the loaded structure. That is the geometry the
// torsional evaluation puts them back to, and it is taken from the structure
// for the same reason the mesh takes its own equilibrium lengths there: it is
// the conformation the model was built to hold.
void SpringNetwork::_setupTorsionalFrames()
{
    if (!_config.dihedral.torsionalonly)
        return;

    bindDihedralEndpointsToAxes();

    unsigned framed = 0;
    for (unsigned family = 0; family < DIHEDRAL_FAMILY_COUNT; ++family)
    {
        std::vector<Spring> & springs = _dihedralsprings[family];
        for (size_t i = 0; i < springs.size(); ++i)
        {
            const unsigned axisIndex = _dihedralAxis[family][i];
            if (axisIndex == NO_AXIS)
                continue;
            const GhostAxis & axis = _ghostaxes[axisIndex];
            const Vector3f & B = getParticle(axis.anchorBIndex).getPosition();
            Vector3f ahat = getParticle(axis.anchorCIndex).getPosition() - B;
            ahat.normalize();

            auto frame = [&](const Particle & p) {
                const Vector3f rel = p.getPosition() - B;
                const float z = rel.dot(ahat);
                return std::make_pair((rel - ahat * z).norm(), z);
            };
            const auto f1 = frame(springs[i].getParticle1());
            const auto f2 = frame(springs[i].getParticle2());
            springs[i].setTorsionalFrame(f1.first, f1.second, f2.first, f2.second);
            ++framed;
        }
    }
    logging::info("SpringNetwork: %u dihedral spring(s) evaluated on their torsion angle alone; "
                  "the tangential filter has nothing left to discard.",
                  framed);
}

// Reads a TORSION force field and resolves its four atoms against the loaded
// network, by residue and name, "+"/"-" naming the next/previous residue just
// as everywhere else. Resolution is by residue index within a chain, so a
// torsion missing an atom at a terminus is skipped rather than mis-bound.
void SpringNetwork::_setupTorsions()
{
    const std::string & path = _config.dihedral.torsionfile;
    if (path.empty())
        return;

    // residue id -> (atom name -> particle index). Built once.
    // The reducer renames every atom to a grain name -- the residue's
    // one-letter code followed by the atom name, MET's N becoming "MN" and
    // GLN's "QN" -- so both spellings are indexed and either resolves. The
    // build-time reader never meets this: it works on the topology before the
    // reduction.
    std::map<unsigned, std::map<std::string, unsigned>> byResidue;
    std::map<unsigned, std::string> resName;
    for (unsigned i = 0; i < getNumberOfParticles(); ++i)
    {
        const Particle & p = getParticle(i);
        const std::string & n = p.getName();
        byResidue[p.getResId()][n] = i;
        if (n.size() > 1)
            byResidue[p.getResId()].emplace(n.substr(1), i);
        resName[p.getResId()] = p.getResName();
    }

    std::ifstream in(path);
    if (!in)
        logging::die("SpringNetwork: cannot open torsion file '%s'", path.c_str());

    auto resolve = [&](unsigned resid, const std::string & name, unsigned & out) {
        int offset = 0;
        std::string bare = name;
        if (!name.empty() && (name[0] == '+' || name[0] == '-'))
        {
            offset = name[0] == '+' ? 1 : -1;
            bare = name.substr(1);
        }
        const auto r = byResidue.find(static_cast<unsigned>(static_cast<int>(resid) + offset));
        if (r == byResidue.end())
            return false;
        const auto a = r->second.find(bare);
        if (a == r->second.end())
            return false;
        out = a->second;
        return true;
    };

    unsigned applied = 0, skipped = 0, lineno = 0;
    std::string line;
    while (std::getline(in, line))
    {
        ++lineno;
        if (line.empty() || line[0] == '#')
            continue;
        std::istringstream ss(line);
        std::string kind;
        if (!(ss >> kind))
            continue;

        if (kind == "TORSIONTABLE")
        {
            unsigned id = 0, bins = 0;
            if (!(ss >> id >> bins) || bins == 0)
                logging::die("SpringNetwork: torsion file line %u: malformed TORSIONTABLE", lineno);
            if (_torsiontables.size() <= id)
                _torsiontables.resize(id + 1);
            TorsionTable & tab = _torsiontables[id];
            tab.bins = bins;
            tab.energy.resize(bins + 1);
            tab.torque.resize(bins + 1);
            for (unsigned b = 0; b <= bins; ++b)
                if (!(ss >> tab.energy[b] >> tab.torque[b]))
                    logging::die("SpringNetwork: torsion file line %u: TORSIONTABLE %u is short at sample %u",
                                 lineno, id, b);
            continue;
        }
        if (kind != "TORSION")
            continue;

        std::string resname, family, names[4];
        unsigned table = 0;
        if (!(ss >> resname >> family >> names[0] >> names[1] >> names[2] >> names[3] >> table))
            logging::die("SpringNetwork: torsion file line %u is malformed", lineno);
        if (table >= _torsiontables.size() || _torsiontables[table].bins == 0)
            logging::die("SpringNetwork: torsion file line %u refers to table %u, which was never given", lineno,
                         table);

        unsigned fam = DIHEDRAL_FAMILY_COUNT;
        for (unsigned f = 0; f < DIHEDRAL_FAMILY_COUNT; ++f)
            if (family == DIHEDRAL_FAMILY_KEYWORDS_RT[f])
                fam = f;
        if (fam == DIHEDRAL_FAMILY_COUNT)
            logging::die("SpringNetwork: torsion file line %u names an unknown family '%s'", lineno, family.c_str());

        for (const auto & r : resName)
        {
            if (r.second != resname)
                continue;
            Torsion t;
            t.family = fam;
            t.table = table;
            bool ok = true;
            for (unsigned k = 0; k < 4 && ok; ++k)
                ok = resolve(r.first, names[k], t.atoms[k]);
            if (!ok)
            {
                ++skipped;
                continue;
            }
            // A torsion spanning two residues is written under each of them,
            // so the same four atoms can arrive twice; applying it twice would
            // double its torque. Keyed on the quadruplet itself, either way
            // round, since a torsion and its reverse are one torsion.
            std::array<unsigned, 4> q{t.atoms[0], t.atoms[1], t.atoms[2], t.atoms[3]};
            const std::array<unsigned, 4> rev{t.atoms[3], t.atoms[2], t.atoms[1], t.atoms[0]};
            if (rev < q)
                q = rev;
            if (!_seenTorsions.insert(q).second)
                continue;

            t.axisIndex = findOrCreateGhostAxis(t.atoms[1], t.atoms[2]);
            _torsions.push_back(t);
            ++applied;
        }
    }
    logging::info("SpringNetwork: %u torsion(s) applied as a couple about their own axis, from %zu "
                  "tabulated parameter set(s) in '%s' (%u skipped for a missing atom); no ring, no "
                  "ghost, no spring.",
                  applied, _torsiontables.size(), path.c_str(), skipped);
}

// AMBER's own V(phi), applied as a couple.
//
// V(phi) = sum_n V_n * (1 + cos(n*phi - gamma_n)), so the torque about the axis
// is -dV/dphi = sum_n n*V_n*sin(n*phi - gamma_n). A rigid side answers to that
// axial torque and to nothing else, which is why handing it to one atom per
// side is exact here rather than a distribution choice: the mesh has already
// removed every other freedom. The reaction on the two axis atoms is left to
// the same closed-form pass that serves the rings.
float SpringNetwork::computeTorsionForces()
{
    if (_torsions.empty())
        return 0.0f;

    const bool enabled[DIHEDRAL_FAMILY_COUNT] = {
        isDihedralPhiEnabled(),        isDihedralPsiEnabled(),           isDihedralOmegaEnabled(),
        isDihedralChiEnabled(),        isDihedralPlanarityEnabled(),     isDihedralNucleicBackboneEnabled(),
        isDihedralNucleicChiEnabled(), isDihedralNucleicSugarEnabled()};

    const float PI = static_cast<float>(M_PI);
    const float unit = _ff->getSpringScale() * static_cast<float>(forcefield::GLOBAL_SPRING_FORCE_CONVERT);
    float energy = 0.0f;

    for (const Torsion & t : _torsions)
    {
        if (!enabled[t.family])
            continue;

        Particle & p1 = getParticle(t.atoms[0]);
        Particle & p2 = getParticle(t.atoms[1]);
        Particle & p3 = getParticle(t.atoms[2]);
        Particle & p4 = getParticle(t.atoms[3]);

        const Vector3f b1 = p2.getPosition() - p1.getPosition();
        const Vector3f b2 = p3.getPosition() - p2.getPosition();
        const Vector3f b3 = p4.getPosition() - p3.getPosition();

        const Vector3f n1 = b1 ^ b2;
        const Vector3f n2 = b2 ^ b3;
        const float n1sq = n1.dot(n1);
        const float n2sq = n2.dot(n2);
        const float b2len = b2.norm();
        if (n1sq < 1e-12f || n2sq < 1e-12f || b2len < 1e-6f)
            continue; // three atoms in line: the dihedral is not defined

        const float phi = std::atan2(b2len * b1.dot(n2), n1.dot(n2));

        // One atan2, one index, two lerps: the harmonics were summed once,
        // offline, and the table holds both what the energy is worth and how
        // hard it pulls.
        const TorsionTable & tab = _torsiontables[t.table];
        const float x = (phi + PI) / (2.0f * PI) * static_cast<float>(tab.bins);
        const unsigned b = std::min(static_cast<unsigned>(std::max(x, 0.0f)), tab.bins - 1);
        const float f = x - static_cast<float>(b);
        energy += tab.energy[b] + f * (tab.energy[b + 1] - tab.energy[b]);
        const float torque = tab.torque[b] + f * (tab.torque[b + 1] - tab.torque[b]);
        if (torque == 0.0f)
            continue;

        // The exact gradient of phi, on all four atoms. A couple on the two
        // outer ones carries the same AXIAL torque and would be exact if the
        // two sides were rigid -- the mesh at 500 kJ.mol-1.A-2 is not, and the
        // difference measurably displaced the backbone by 3 degrees. This
        // costs a handful of cross products, needs no reaction on the axis and
        // conserves momentum by construction: the four forces sum to zero
        // identically, since F2 and F3 are built from F1 and F4.
        const float scale = unit * torque;
        const Vector3f F1 = n1 * (-scale * b2len / n1sq);
        const Vector3f F4 = n2 * (scale * b2len / n2sq);
        const float inv = 1.0f / (b2len * b2len);
        const float c1 = b1.dot(b2) * inv;
        const float c3 = b3.dot(b2) * inv;
        // Signs matter here and are not guessable: the decomposition is written
        // in terms of r_ij = r_i - r_j and r_kl = r_k - r_l, which are the
        // NEGATIVES of b1 and b3 as spelled above. Checked against a finite
        // difference of phi itself -- 2.7e-10 -- after the other sign
        // convention passed the sum-to-zero test while being wrong by 5 rad/A.
        const Vector3f F2 = F1 * (-(c1 + 1.0f)) + F4 * c3;
        const Vector3f F3 = F1 * c1 - F4 * (c3 + 1.0f);

        p1.addForce(F1);
        p2.addForce(F2);
        p3.addForce(F3);
        p4.addForce(F4);
    }
    return energy;
}

void SpringNetwork::_setupGhostSprings()
{
    if (!_config.dihedral.ghostsprings || _ghostparticles.empty())
        return;

    const float k = static_cast<float>(_config.dihedral.ghostspringstiffness);
    const float mass = static_cast<float>(_config.dihedral.ghostmass);
    const size_t nghosts = _ghostparticles.size();

    // Flip every ghost first, in one pass, then rebuild the two membership
    // lists once. updateParticleState erases from a vector by value, which is
    // linear, so calling it per ghost would be quadratic in the ring count.
    for (const GhostParticleBinding & binding : _ghostparticles)
    {
        Particle & ghost = getParticle(binding.ownIndex);
        ghost.setMass(mass);
        ghost.setStatic(false);
    }
    _staticparticules.clear();
    _dynamicparticules.clear();
    for (const Particle & p : _particles)
    {
        const unsigned id = static_cast<unsigned>(p.getId());
        if (p.isStatic())
            addStaticParticle(id);
        else
            addDynamicParticle(id);
    }

    // addSpring rebuilds every particle's spring-neighbour cache whenever the
    // spring vector moves; make it move once rather than once per doubling.
    _springs.reserve(_springs.size() + 4 * nghosts);
    _rebuildSpringNeighbors();

    auto tie = [this, k](unsigned a, unsigned b) -> unsigned {
        if (a == b || getParticle(a).isInSpringNeighbors(b))
            return 0;
        const float d = (getParticle(a).getPosition() - getParticle(b).getPosition()).norm();
        if (d < 1e-3f) // coincident: no direction to hold, so no spring to build
            return 0;
        addSpring(a, b, d, k);
        return 1;
    };

    unsigned nanchor = 0;
    for (const GhostParticleBinding & binding : _ghostparticles)
    {
        nanchor += tie(binding.ownIndex, binding.anchorRefIndex);
        nanchor += tie(binding.ownIndex, binding.anchorBIndex);
        nanchor += tie(binding.ownIndex, binding.anchorCIndex);
    }

    unsigned nchord = 0;
    if (_config.dihedral.ghostringchords)
    {
        // One ring side is the set of ghosts sharing an axis AND a reference
        // atom: they are that one atom replicated at M evenly spaced azimuths.
        std::map<std::pair<unsigned, unsigned>, std::vector<size_t>> sides;
        for (size_t i = 0; i < nghosts; ++i)
            sides[{_ghostparticles[i].axisIndex, _ghostparticles[i].anchorRefIndex}].push_back(i);

        for (auto & entry : sides)
        {
            std::vector<size_t> & side = entry.second;
            if (side.size() < 2)
                continue;
            std::sort(side.begin(), side.end(), [this](size_t a, size_t b)
                      { return _ghostparticles[a].delta_deg < _ghostparticles[b].delta_deg; });
            for (size_t i = 0; i + 1 < side.size(); ++i)
                nchord += tie(_ghostparticles[side[i]].ownIndex, _ghostparticles[side[i + 1]].ownIndex);
            // Close the ring, but only where that is a third distinct chord.
            if (side.size() > 2)
                nchord += tie(_ghostparticles[side.back()].ownIndex, _ghostparticles[side.front()].ownIndex);
        }
    }

    // The bindings have done their work. Dropping them is what actually stops
    // updateGhostPositions and the Rodrigues redistribution, since both iterate
    // this list -- while _ghostaxes, which the tangential filter and the axis
    // reaction both still need, is a separate vector and stays.
    _springHeldGhosts.reserve(nghosts);
    for (const GhostParticleBinding & binding : _ghostparticles)
        _springHeldGhosts.push_back(binding.ownIndex);

    _ghostparticles.clear();
    _ghostForceScratch.clear();
    _ghostsAreSpringHeld = true;
    _ghostDamping = static_cast<float>(_config.dihedral.ghostdamping);

    logging::info("SpringNetwork: %zu ghost(s) now spring-held at %.1f kJ.mol-1.A-2 and %.2f Da "
                  "(%u anchor spring(s), %u ring chord(s)); per-step placement and force "
                  "redistribution are off.",
                  nghosts, static_cast<double>(k), static_cast<double>(mass), nanchor, nchord);
}

// Friction on the ghosts alone. Applied after redistributeGhostForces, so it
// acts on everything a ghost has been given -- its ring spring's filtered
// share included, which is what the filter pumps through.
//
// Friction on a fictitious particle is not friction on the protein: no real
// atom is touched here, and the mesh keeps whatever global viscosity.value it
// was given. What leaves is the work the projection put in.
void SpringNetwork::applyGhostDamping()
{
    if (_ghostDamping <= 0.0f)
        return;
    for (const unsigned id : _springHeldGhosts)
        getParticle(id).applyViscosity(_ghostDamping);
}

void SpringNetwork::updateParticleState(unsigned id, bool isStatic) {
    if (isStatic) {
        removeDynamicParticle(id);
        addStaticParticle(id);
    } else {
        removeStaticParticle(id);
        addDynamicParticle(id);
    }
    _particles[id].setStatic(isStatic);
}

void SpringNetwork::clear()
{
    _initparticles.clear();
    _particles.clear();
    _staticparticules.clear();
    _dynamicparticules.clear();
    _chargedparticules.clear();
    _hydrophobicparticules.clear();
    _springs.clear();
    _staticsprings.clear();
    _dynamicsprings.clear();
    for (auto & family : _dihedralsprings)
        family.clear();
    _ghostparticles.clear();
    _springForceScratch.clear();
    _nsearch.steric.reset();
    _nsearch.electrostatic.reset();
    _nsearch.hydrophobic.reset();
    _neighborSearchesDirty = false;
    _insertionVector.reset();
    _probeparticule = Particle();
}

// TODO: implement this function in `Topology`
// int SpringNetwork::removeSpringsBetweenSelections(std::vector<int> sel1, std::vector<int> sel2)
// {
//     int nbRemove = 0;
//     for (vector<int>::iterator it1 = sel1.begin(); it1 != sel1.end(); it1++)
//     {
//         for (vector<int>::iterator it2 = sel2.begin(); it2 != sel2.end(); it2++)
//         {
//             Spring * s = getSpringFromParticlesIds(*it1, *it2);
//             if (s)
//             {
//                 if (s->getParticle1().getId() == *it1)
//                     s->getParticle1().removeSpringNeighbor(*it2);
//                 else
//                     s->getParticle2().removeSpringNeighbor(*it2);

//                 if (s->getParticle1().getId() == *it2)
//                     s->getParticle1().removeSpringNeighbor(*it1);
//                 else
//                     s->getParticle2().removeSpringNeighbor(*it1);
//                 _springs.erase(std::find(_springs.begin(), _springs.end(), s));
//                 nbRemove++;
//             }
//         }
//     }

//     return nbRemove;
// }

// =====================================================================================

void SpringNetwork::setInsertionVector(unsigned aa1, unsigned aa2)
{
    try
    {
        getParticleFromId(aa1);
        getParticleFromId(aa2);
    }
    catch (const std::out_of_range & e)
    {
        logging::die("Insertion vector could not be initialized. One or both of the pair of particles "
                     "defining the insertion vector have an invalid identifier.");
    }

    //   Check that the two elements of the vector are identical
    if (aa1 == aa2)
    {
        logging::warning("The two values defining the insertion vector are identical: %zu and %zu",
                         aa1, aa2);
    }

    _insertionVector =
        std::make_unique<InsertionVector>(*this, getParticleFromId(aa1), getParticleFromId(aa2));
    
    const Particle p1 = _insertionVector->getParticle(0);
    const Particle p2 = _insertionVector->getParticle(1);
    logging::info("Insertion vector set between %s%d:%s and %s%d:%s",
        p1.getResName().c_str(), p1.getResId(), p1.getName().c_str(),
        p2.getResName().c_str(), p2.getResId(), p2.getName().c_str());
}

void SpringNetwork::applyConstraints()
{
    if (_constraints.empty())
    {
        _meanConstraintsDistances = 0.0f;
        return;
    }

    float sumDistances = 0.0;
    for (unsigned i = 0; i < _constraints.size(); i++)
    {
        _constraints[i]->apply();
        sumDistances += _constraints[i]->getDistance();
    }
    _meanConstraintsDistances = sumDistances / _constraints.size();
}

void SpringNetwork::clearParticles(void) { clear(); }

// =====================================================================================
//
// Setup methods
//
// =====================================================================================

void SpringNetwork::setup(const configuration::Configuration & conf)
{
    _config = conf;

    _setupForceField();
    _setupProbe();
    _setupSteric();
    _setupElectrostatic();
    _setupHydrophobic();
    _setupDensityGrid();
    _setupInsertionVector();
    _setupTrajectories();
    // After everything else: it converts particles and adds springs, so it
    // needs the network whole and the configuration already stored.
    _setupTorsions();
    _setupTorsionalFrames();
    _setupGhostSprings();
    _neighborSearchesDirty = false;
    // _setupConstraints();
    // _setupSelections();
}

void SpringNetwork::_setupSteric()
{
    if (isStericEnabled())
    {
        if (getStericCutoff() < 1e-6)
            throw std::runtime_error("Steric cutoff must be > 0");
        _nsearch.steric = make_nsearch(_particles, getStericCutoff(), getNeighborSkin());
        _excludeProbeFromNeighborSearch(*_nsearch.steric);
    }
}

void SpringNetwork::_setupHydrophobic()
{
    if (isHydrophobicityEnabled())
    {
        if (getHydrophobicCutoff() < 1e-6)
            throw std::runtime_error("Hydrophobic cutoff must be > 0");
        const std::vector<size_t> hydrophobic_particles = _hydrophobicParticleIndexes();
        if (!hydrophobic_particles.empty())
        {
            _nsearch.hydrophobic =
                make_nsearch(_particles, getHydrophobicCutoff(), hydrophobic_particles, getNeighborSkin());
            _excludeProbeFromNeighborSearch(*_nsearch.hydrophobic);
        }
    }
}

void SpringNetwork::_setupForceField()
{
    const std::string steric = _config.steric.mode;

    if (steric == "lennard-jones-8-6Lewitt")
        _ff = std::make_unique<forcefield::ForceFieldElectrostaticCoulombAndStericLennardJones_8_6Lewitt>();
    else if (steric == "lennard-jones-8-6Zacharias")
        _ff = std::make_unique<forcefield::ForceFieldElectrostaticCoulombAndStericLennardJones_8_6Zacharias>();
    else if (steric == "lennard-jones-12-6Amber")
        _ff = std::make_unique<forcefield::ForceFieldElectrostaticCoulombAndStericLennardJones_12_6Amber>();
    else
        _ff = std::make_unique<forcefield::ForceFieldElectrostaticCoulombAndStericLinear>();

    _ff->setStericScale(_config.steric.gridscale);
    _ff->setCoulombScale(_config.electrostatic.scale);
    _ff->setDielectric(_config.electrostatic.dielectric);
    _ff->setForceFieldScale(_config.potentialgrid.scale);
    _ff->setSpringScale(_config.spring.scale);
    _ff->setIMPScale(_config.imp.scale);
    _ff->setHydrophobicityScale(_config.hydrophobicity.scale);
}

void SpringNetwork::_setupElectrostatic()
{
    if (!isElectrostaticEnabled())
        return;

    if (isElectrostaticFieldEnabled())
    {
        const std::string dxpath = _config.potentialgrid.path;
        logging::info("Reading electrostatic map from DX file '%s'", dxpath.c_str());
        _grids.potential = opendx::readGrid(dxpath);
    }

    if (isElectrostaticCoulombEnabled())
    {
        if (getElectrostaticCutoff() < 1e-6)
            throw std::runtime_error("Electrostatic cutoff must be > 0");

        const std::vector<size_t> charged_particles = _chargedParticleIndexes();
        if (!charged_particles.empty())
        {
            _nsearch.electrostatic =
                make_nsearch(_particles, getElectrostaticCutoff(), charged_particles, getNeighborSkin());
            _excludeProbeFromNeighborSearch(*_nsearch.electrostatic);
        }
    }
}

void SpringNetwork::_setupDensityGrid()
{
    if (isDensityGridEnabled())
    {
        const std::string dxpath = _config.densitygrid.path;
        logging::info("Reading density grid from DX file '%s'", dxpath.c_str());
        _grids.density = opendx::readGrid(dxpath);
    }
}

void SpringNetwork::_setupProbe()
{
    if (isProbeEnabled())
    {
        _probeparticule.setSpringNetwork(this);
        _probeparticule.setName("PRB");
        _probeparticule.setResName("PRB");
        _probeparticule.setCharge(_config.probe.charge);
        _probeparticule.setRadius(_config.probe.radius);
        _probeparticule.setEpsilon(_config.probe.epsilon);
        _probeparticule.setMass(_config.probe.mass);
        _probeparticule.setX(_config.probe.x);
        _probeparticule.setY(_config.probe.y);
        _probeparticule.setZ(_config.probe.z);
        _probeparticule.setId(static_cast<int>(_particles.size()));
        _particles.push_back(_probeparticule);
    }
}

void SpringNetwork::_setupTrajectories()
{
    if (_config.pdbtraj.enable)
        _trajectories.add_writer(
            std::make_unique<io::modern::PDBTrajectoryWriter>(_config.pdbtraj.path, *this, _config.pdbtraj.frequency));
    if (_config.xtctraj.enable)
        _trajectories.add_writer(
            std::make_unique<io::modern::XTCTrajectoryWriter>(_config.xtctraj.path, *this, _config.xtctraj.frequency));
    if (_config.csvsample.enable)
        _trajectories.add_writer(std::make_unique<io::modern::CSVTrajectoryWriter>(_config.csvsample.path, *this,
                                                                                   _config.csvsample.frequency));
}

void SpringNetwork::_setupInsertionVector()
{
    if (_config.ivector.enable)
        setInsertionVector(_config.ivector.vector[0], _config.ivector.vector[1]);
}


std::vector<size_t> SpringNetwork::_chargedParticleIndexes() const
{
    std::vector<size_t> indexes;
    indexes.reserve(_particles.size());

    for (size_t i = 0; i < _particles.size(); ++i)
    {
        if (_particles[i].isCharged())
            indexes.push_back(i);
    }

    return indexes;
}

std::vector<size_t> SpringNetwork::_hydrophobicParticleIndexes() const
{
    std::vector<size_t> indexes;
    indexes.reserve(_particles.size());

    for (size_t i = 0; i < _particles.size(); ++i)
    {
        if (_particles[i].isHydrophobic())
            indexes.push_back(i);
    }

    return indexes;
}

void SpringNetwork::_excludeProbeFromNeighborSearch(NeighborSearch::Searcher & searcher)
{
    if (!isProbeEnabled())
        return;

    const int probe_id = _probeparticule.getId();
    if (probe_id >= 0)
        searcher.exclude_index(static_cast<size_t>(probe_id));
}

void SpringNetwork::_markNeighborSearchesDirty()
{
    if (_nsearch.steric || _nsearch.electrostatic || _nsearch.hydrophobic)
        _neighborSearchesDirty = true;
}

void SpringNetwork::_updateNeighborSearches()
{
    if (!_neighborSearchesDirty)
        return;

    if (_nsearch.steric)
        _nsearch.steric->update();
    if (_nsearch.electrostatic)
        _nsearch.electrostatic->update();
    if (_nsearch.hydrophobic)
        _nsearch.hydrophobic->update();

    _neighborSearchesDirty = false;
}

void SpringNetwork::_resizeNonbondedPairScratch()
{
    const size_t n = _dynamicparticules.size();

    _stericPairScratch.resize(n);
    _electrostaticPairScratch.resize(n);
    _hydrophobicPairScratch.resize(n);

    // Clears logical contents but keeps each bucket's capacity, so the
    // simulation loop does not reallocate every step.
    for (auto & bucket : _stericPairScratch)
        bucket.clear();
    for (auto & bucket : _electrostaticPairScratch)
        bucket.clear();
    for (auto & bucket : _hydrophobicPairScratch)
        bucket.clear();
}

void SpringNetwork::_applyNonbondedPairScratch(
    const std::vector<std::vector<spn::DeferredNonbondedContribution>> & scratch, float & energy)
{
    for (const auto & bucket : scratch)
    {
        for (const auto & contribution : bucket)
        {
            getParticle(contribution.target).addForce(contribution.force);
            energy += contribution.energy;
        }
    }
}

void SpringNetwork::_syncProbeParticle()
{
    if (!isProbeEnabled())
        return;

    const int probe_id = _probeparticule.getId();
    if (probe_id < 0 || static_cast<size_t>(probe_id) >= _particles.size())
        return;

    _particles[static_cast<size_t>(probe_id)] = _probeparticule;
    _particles[static_cast<size_t>(probe_id)].setSpringNetwork(this);
}

void SpringNetwork::_rebuildSpringNeighbors()
{
    for (Particle & particle : _particles)
        particle.clearSpringNeighbors();

    for (Spring & spring : _springs)
    {
        Particle & p1 = spring.getParticle1();
        Particle & p2 = spring.getParticle2();
        p1.addToSpringNeighbors(static_cast<unsigned>(p2.getId()), &spring);
        p2.addToSpringNeighbors(static_cast<unsigned>(p1.getId()), &spring);
    }
}

void SpringNetwork::_setupSelections() { throw "Not Implemented Error"; }
void SpringNetwork::_setupConstraints() { throw "Not Implemented Error"; }

} // namespace spn
} // namespace biospring
