#include "SpringNetwork.h"

#include <fstream>
#include "logging.h"
#include "measure.hpp"

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

// Calculates STRETCH spring forces and applies them to the particles.
// Updates global `_energies.stretch`. A plain real-atom spring (no ghosts,
// unlike BEND/DIHEDRAL: a 1-2 bond length has no 3-body conflation to
// avoid), so no ignoreDynamicState override -- a stretch spring between two
// particles the user froze with --static correctly contributes nothing,
// same as it did before this spring lived in its own collection. Not
// parallelised, like computeDihedralForces: kept simple until profiling
// ever shows a need otherwise.
void SpringNetwork::computeStretchForces()
{
    _energies.stretch = _computeSpringCollectionForces(_stretchsprings, /*ignoreDynamicState=*/false,
                                                       /*subtractDcOffset=*/false);
}

// Shared parallel-compute/serial-accumulate loop for the STRETCH, BEND and
// DIHEDRAL spring collections -- the exact pattern computeSpringForces
// established: the spring evaluation (the expensive part) runs in parallel
// into _springForceScratch, then a deterministic serial pass applies the
// forces (particles are shared between springs, so they must not be
// written concurrently) and accumulates the energy in a fixed order,
// keeping the result independent of the thread count. Returns the summed
// energy. The scratch buffer is shared with computeSpringForces: the four
// loops run sequentially within one step, and resize never shrinks
// capacity, so no per-step allocation happens either way.
float SpringNetwork::_computeSpringCollectionForces(std::vector<Spring> & springs, bool ignoreDynamicState,
                                                    bool subtractDcOffset)
{
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
        spring.getParticle1().addForce(force);
        spring.getParticle2().addForce(-force);
        energy += subtractDcOffset ? spring.getEnergy() - spring.getDcOffset() : spring.getEnergy();
    }
    return energy;
}

// Calculates BEND ghost-ghost spring forces and applies them to the
// particles. Updates global `_energies.bend`. ignoreDynamicState=true for
// the same reason as computeDihedralForces: both ghosts are always static.
// Gated by isBendingEnabled() (a runtime debug on/off for a family that
// -bending already decided, at build time, to actually create springs
// for -- see Configuration.hpp's own comment) on top of the isSpringEnabled()
// master switch the caller (computeForces) already applies: when disabled,
// BOTH force and energy are skipped (a pure reporting-only toggle would
// leave the dynamics unchanged, defeating the point of isolating this
// family's contribution).
void SpringNetwork::computeBendForces()
{
    _energies.bend = isBendingEnabled() ? _computeSpringCollectionForces(_bendsprings, /*ignoreDynamicState=*/true,
                                                                         /*subtractDcOffset=*/false)
                                        : 0.0f;
}

// Calculates dihedral ghost-spring forces and applies them to the
// particles. Updates global `_energies.dihedral`. Which families were
// actually BUILT is a build-time decision (see -dihedral/--dihedral in
// pdb2spn-cli.cpp); which of the built ones are actually APPLIED here is
// independently gated per-family (isDihedralPhi/Psi/Omega/ChiEnabled(),
// see Configuration.hpp's own comment) on top of the isSpringEnabled()
// master switch the caller (computeForces) already applies -- same
// both-force-and-energy-skipped reasoning as computeBendForces above (a
// disabled family must not silently keep influencing the dynamics).
//
// ignoreDynamicState=true: both endpoints of a dihedral ghost spring are
// always static virtual sites (see Spring::computeForce's own doc) --
// without this, the force/energy here would silently stay zero for every
// single one of these springs. subtractDcOffset=true: each spring's own
// share of its axis's exact dihedral-energy correction (see
// Spring::getDcOffset's own comment) is subtracted from the reported
// total; it never affects the forces.
void SpringNetwork::computeDihedralForces()
{
    float dihedralenergy = 0.0f;

    auto accumulate = [&](std::vector<Spring> & springs) {
        dihedralenergy += _computeSpringCollectionForces(springs, /*ignoreDynamicState=*/true,
                                                         /*subtractDcOffset=*/true);
    };

    const bool enabled[DIHEDRAL_FAMILY_COUNT] = {
        isDihedralPhiEnabled(),             isDihedralPsiEnabled(),        isDihedralOmegaEnabled(),
        isDihedralChiEnabled(),             isDihedralPlanarityEnabled(),  isDihedralNucleicBackboneEnabled(),
        isDihedralNucleicChiEnabled(),      isDihedralNucleicSugarEnabled()};
    for (unsigned family = 0; family < DIHEDRAL_FAMILY_COUNT; ++family)
        if (enabled[family])
            accumulate(_dihedralsprings[family]);

    _energies.dihedral = dihedralenergy;
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
        computeStretchForces();
        computeBendForces();
        computeDihedralForces();
    }
    computeParticleForces();
}

void SpringNetwork::computeStep()
{
    idleRun();
    _meanConstraintsDistances = 0.0;

    computeForces();
    redistributeGhostForces();

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
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
    interactor::InteractorMDDriver* interactormddriver = getInteractorInstance<interactor::InteractorMDDriver>();
    logging::info("    MDDriver parameters:");
    logging::info("      port: %d is open for connection.", interactormddriver->getPort());
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
        logging::info("Stretch energy: %5.2f kJ.mol-1", _energies.stretch);
        logging::info("Bend energy: %5.2f kJ.mol-1", _energies.bend);
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
    logging::info("Total SASA (dynamic=%s): %5.2lf A2", _freesasaState.isDynamic? "true" : "false", _freesasaState.sasaTotal);
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

// Shared implementation for addDihedralSpring and the STRETCH/BEND adds:
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
                                      float dcOffset)
{
    if (family >= DIHEDRAL_FAMILY_COUNT)
        throw std::out_of_range("SpringNetwork::addDihedralSpring: dihedral family index out of range");
    addDihedralSpringTo(_dihedralsprings[family], _particles, id1, id2, equilibrium, stiffness, dcOffset);
}

void SpringNetwork::addStretchSpring(unsigned id1, unsigned id2, float equilibrium, float stiffness)
{
    addDihedralSpringTo(_stretchsprings, _particles, id1, id2, equilibrium, stiffness, /*dcOffset=*/0.0f);
}

void SpringNetwork::addBendSpring(unsigned id1, unsigned id2, float equilibrium, float stiffness)
{
    addDihedralSpringTo(_bendsprings, _particles, id1, id2, equilibrium, stiffness, /*dcOffset=*/0.0f);
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
        (placement == GhostPlacement::AxialOffset)
            ? GhostParticle::computePositionAxial(getParticle(anchorBIndex).getPosition(),
                                                  getParticle(anchorCIndex).getPosition(), r)
            : GhostParticle::computePositionByRotation(
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
    unsigned axisIndex = 0;
    for (; axisIndex < _ghostaxes.size(); ++axisIndex)
        if (_ghostaxes[axisIndex].anchorBIndex == anchorBIndex && _ghostaxes[axisIndex].anchorCIndex == anchorCIndex)
            break;
    if (axisIndex == _ghostaxes.size())
        _ghostaxes.push_back(GhostAxis{anchorBIndex, anchorCIndex, Vector3f(), Vector3f(), Vector3f(), Vector3f(),
                                       {}, {}, {}, {}, 0.0f, 0.0f});

    const float delta_rad = delta_deg * static_cast<float>(M_PI) / 180.0f;
    _ghostparticles.push_back(GhostParticleBinding{ownIndex, anchorBIndex, anchorCIndex, anchorRefIndex, r, theta_deg,
                                                   delta_deg, std::cos(delta_rad), std::sin(delta_rad), axisIndex,
                                                   placement});
    return ownIndex;
}

void SpringNetwork::setAxisSubstituents(unsigned anchorBIndex, unsigned anchorCIndex,
                                        const std::vector<unsigned> & substituentsB,
                                        const std::vector<unsigned> & substituentsC,
                                        const std::vector<double> & weightsB,
                                        const std::vector<double> & weightsC)
{
    // Stored as float, and only when the force field supplied one per atom --
    // a partial list would silently mix weighted and equal shares on the same
    // side, which is worse than either.
    std::vector<float> wB, wC;
    if (weightsB.size() == substituentsB.size())
        for (double w : weightsB)
            wB.push_back(static_cast<float>(w));
    if (weightsC.size() == substituentsC.size())
        for (double w : weightsC)
            wC.push_back(static_cast<float>(w));

    // BOTH ORIENTATIONS, and missing this leaves exactly half the model
    // uncovered. One torsion builds two GhostAxis entries, (B, C) and
    // (C, B): a ring's ghosts rotate about the axis with their OWN anchor as
    // B, so the substituents of each side hang off a different entry. Setting
    // only (B, C) covered the B side and left every C-side ring pushing its
    // reference atom alone (measured on the RNA hairpin: 82 of 160 axes, and
    // the 78 missing were precisely the reversed ones).
    //
    // The reversed entry takes the lists swapped, which is what it means: for
    // axis (C, B), the atoms bonded to its own B -- that is, to C -- are the
    // original C side.
    //
    // The force field also writes records for axes no ghost was ever built on
    // (a family switched off at build time, an unresolved anchor at a chain
    // terminus). Those match nothing here and are dropped: no torsion turns
    // about them, so there is nothing to share out.
    // A (B, C) pair can also exist as a BEND axis -- BEND anchors its two
    // axial ghosts on the same two atoms a torsion may turn about -- and such
    // an axis carries no torsion to share out. Its ghosts are AxialOffset and
    // are self-balancing, so a list there would never be read; it would only
    // inflate the coverage count into claiming more axes than exist.
    for (size_t i = 0; i < _ghostaxes.size(); ++i)
    {
        GhostAxis & axis = _ghostaxes[i];
        const bool matches = axis.anchorBIndex == anchorBIndex && axis.anchorCIndex == anchorCIndex;
        const bool reversed = axis.anchorBIndex == anchorCIndex && axis.anchorCIndex == anchorBIndex;
        if (!matches && !reversed)
            continue;

        bool carriesRing = false;
        for (const GhostParticleBinding & binding : _ghostparticles)
            if (binding.axisIndex == i && binding.placement == GhostPlacement::AxisRotation)
            {
                carriesRing = true;
                break;
            }
        if (!carriesRing)
            continue;

        axis.substituentsB = matches ? substituentsB : substituentsC;
        axis.substituentsC = matches ? substituentsC : substituentsB;
        axis.weightsB = matches ? wB : wC;
        axis.weightsC = matches ? wC : wB;
    }
}

size_t SpringNetwork::getNumberOfAxesWithSubstituents() const
{
    size_t n = 0;
    for (const GhostAxis & axis : _ghostaxes)
        if (!axis.substituentsB.empty() || !axis.substituentsC.empty())
            n++;
    return n;
}

void SpringNetwork::getAxisSubstituents(std::vector<std::array<unsigned, 2>> & anchors,
                                        std::vector<std::array<unsigned, 2>> & counts, std::vector<unsigned> & atoms,
                                        std::vector<float> * weights) const
{
    anchors.clear();
    counts.clear();
    atoms.clear();
    if (weights != nullptr)
        weights->clear();
    for (const GhostAxis & axis : _ghostaxes)
    {
        if (axis.substituentsB.empty() && axis.substituentsC.empty())
            continue;
        anchors.push_back({axis.anchorBIndex, axis.anchorCIndex});
        counts.push_back({static_cast<unsigned>(axis.substituentsB.size()),
                          static_cast<unsigned>(axis.substituentsC.size())});
        atoms.insert(atoms.end(), axis.substituentsB.begin(), axis.substituentsB.end());
        atoms.insert(atoms.end(), axis.substituentsC.begin(), axis.substituentsC.end());
        if (weights == nullptr)
            continue;
        // A side with no weights is written as -1 per atom, the same "no
        // weight given" marker the .bi.ff's unweighted form produces, so the
        // array stays parallel to `atoms` and one length check covers both.
        for (int side = 0; side < 2; ++side)
        {
            const std::vector<unsigned> & subs = side == 0 ? axis.substituentsB : axis.substituentsC;
            const std::vector<float> & w = side == 0 ? axis.weightsB : axis.weightsC;
            for (size_t k = 0; k < subs.size(); ++k)
                weights->push_back(w.size() == subs.size() ? w[k] : -1.0f);
        }
    }
}

// Shares `torque` (signed, about the axis) equally over `substituents`, as
// AMBER does: it splits a torsion into one term per (substituent,
// substituent) pair with an equal barrier each, so every substituent on a
// side ends up carrying 1/N of that side's torque. Converting a share back
// to a force is one division, |dphi/dr_i| = 1/(|r_i| sin theta) being the
// inverse perpendicular radius -- no Jacobian, only positions relative to
// the axis, which the projection already computes.
bool SpringNetwork::_distributeAxialTorque(GhostAxis & axis, const std::vector<unsigned> & substituents,
                                           const std::vector<float> & weights, float torque)
{
    const Vector3f & B = getParticle(axis.anchorBIndex).getPosition();
    const Vector3f & C = getParticle(axis.anchorCIndex).getPosition();
    Vector3f ahat = C - B;
    const float axisLength = ahat.norm();
    if (axisLength < 1e-6f)
        return false;
    ahat = ahat / axisLength;

    // Two passes over the same small list. An atom sitting ON the axis has no
    // lever arm and cannot carry any share, so it has to be excluded BEFORE
    // the shares are settled -- including it and then skipping it would
    // quietly lose its share and shrink the total torque. Same reason the
    // weights are renormalised here rather than trusted as written: the
    // force field lists every spelling of an atom and every substituent
    // including cross-residue ones, and what survives resolution is a subset.
    const bool weighted = weights.size() == substituents.size();
    std::vector<Vector3f> tangents;
    std::vector<float> radii, shares;
    tangents.reserve(substituents.size());
    radii.reserve(substituents.size());
    shares.reserve(substituents.size());
    float total = 0.0f;
    for (size_t k = 0; k < substituents.size(); ++k)
    {
        const Vector3f rel = getParticle(substituents[k]).getPosition() - B;
        const Vector3f radial = rel - ahat * rel.dot(ahat);
        const float rn = radial.norm();
        if (rn <= 1e-6f)
            continue;
        // A negative weight is the force field saying "no weight given".
        const float share = (weighted && weights[k] >= 0.0f) ? weights[k] : 1.0f;
        tangents.push_back(ahat ^ (radial / rn));
        radii.push_back(rn);
        shares.push_back(share);
        total += share;
    }
    // Every carrier weighted zero: AMBER really does assign this side no
    // barrier through these atoms (phi's amide H is the standard case), so
    // there is nothing to hand out rather than something to spread evenly.
    if (tangents.empty() || total <= 0.0f)
        return false;

    size_t slot = 0;
    for (unsigned index : substituents)
    {
        const Vector3f rel = getParticle(index).getPosition() - B;
        const Vector3f radial = rel - ahat * rel.dot(ahat);
        if (radial.norm() <= 1e-6f)
            continue;
        const Vector3f F = tangents[slot] * (torque * shares[slot] / total / radii[slot]);
        slot++;
        if (F.norm() <= 0.0f)
            continue;

        Particle & p = getParticle(index);
        p.addForce(F);
        axis.sumAtomForces += F;
        axis.sumAtomTorquesAboutB += (p.getPosition() - B) ^ F;
    }
    return true;
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
        if (binding.placement == GhostPlacement::AxialOffset)
        {
            // Self-balancing, no reference atom involved, no accumulation.
            GhostParticle::redistributeForceAxial(B, C, binding.r, getParticle(binding.ownIndex).getForce(), out.F_B,
                                                  out.F_C);
            out.F_Ref = Vector3f();
        }
        else
        {
            out.F_Ref = GhostParticle::rotateForceToAtom(B, C, getParticle(binding.ownIndex).getForce(),
                                                         binding.cos_delta, binding.sin_delta);
        }
    }

    // Pass 2, serial: apply to the real atoms (heavily shared, so not
    // concurrently writable) while accumulating each axis's force and
    // torque totals.
    for (GhostAxis & axis : _ghostaxes)
    {
        axis.sumGhostForces = Vector3f();
        axis.sumGhostTorquesAboutB = Vector3f();
        axis.sumAtomForces = Vector3f();
        axis.sumAtomTorquesAboutB = Vector3f();
        axis.axialTorqueB = 0.0f;
        axis.axialTorqueC = 0.0f;
    }
    for (size_t i = 0; i < _ghostparticles.size(); ++i)
    {
        const GhostParticleBinding & binding = _ghostparticles[i];
        const GhostForceContribution & in = _ghostForceScratch[i];
        Particle & ghost = getParticle(binding.ownIndex);

        if (binding.placement == GhostPlacement::AxialOffset)
        {
            getParticle(binding.anchorBIndex).addForce(in.F_B);
            getParticle(binding.anchorCIndex).addForce(in.F_C);
            ghost.setForce(Vector3f(0.0f, 0.0f, 0.0f));
            continue;
        }

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
        Vector3f F_ref = in.F_Ref;
        if (_config.dihedral.tangentialonly)
        {
            const Vector3f & C = getParticle(axis.anchorCIndex).getPosition();
            Vector3f ahat = C - B;
            ahat.normalize();
            const Vector3f rel = ref.getPosition() - B;
            Vector3f radial = rel - ahat * rel.dot(ahat);
            const float rn = radial.norm();
            // An atom sitting on the axis has no lever arm and no torsional
            // role: it gets nothing rather than an arbitrary direction.
            F_ref = Vector3f();
            if (rn > 1e-6f)
            {
                const Vector3f that = ahat ^ (radial / rn);
                F_ref = that * in.F_Ref.dot(that);
            }
        }

        // With distribution on, this force is not applied here: it is banked
        // as the torque it carries and shared out per axis in pass 2b, so
        // every bonded substituent gets its AMBER share instead of the ring's
        // one reference atom taking the lot. Only the axial component is
        // banked, which is all a filtered force has -- and it is independent
        // of where along the axis it is measured from, so B serves for both
        // sides (a shift d along the axis changes the torque by (d a) ^ F,
        // whose axial part is a . (a ^ F) = 0).
        bool banked = false;
        if (_distributetorque)
        {
            const Vector3f & C = getParticle(axis.anchorCIndex).getPosition();
            Vector3f ahat = C - B;
            if (ahat.norm() > 1e-6f)
            {
                ahat.normalize();
                const float tau = ahat.dot((ref.getPosition() - B) ^ F_ref);
                // Which side the reference atom is on decides which share it
                // joins. An atom the force field did not list -- an axis with
                // no records, a substituent that failed to resolve -- is not
                // banked, and falls through to the direct path below.
                for (unsigned index : axis.substituentsB)
                    if (index == binding.anchorRefIndex)
                    {
                        axis.axialTorqueB += tau;
                        banked = true;
                        break;
                    }
                if (!banked)
                    for (unsigned index : axis.substituentsC)
                        if (index == binding.anchorRefIndex)
                        {
                            axis.axialTorqueC += tau;
                            banked = true;
                            break;
                        }
            }
        }

        if (!banked)
        {
            ref.addForce(F_ref);
            axis.sumAtomForces += F_ref;
            axis.sumAtomTorquesAboutB += (ref.getPosition() - B) ^ F_ref;
        }

        ghost.setForce(Vector3f(0.0f, 0.0f, 0.0f));
    }

    // Pass 2b: share each side's banked torque over its bonded substituents.
    // Runs between the transfer and the axis reaction because pass 3 balances
    // against what actually reached the atoms, and that is only settled here.
    if (_distributetorque)
    {
        // The return value is deliberately ignored: a side only holds banked
        // torque because one of ITS OWN listed atoms carried a tangential
        // force in pass 2, which required a non-zero lever arm on the same
        // positions this pass reads. So a banked side always has at least one
        // carrier and the call cannot fail here. It can only return false on
        // a degenerate geometry, where nothing is left to conserve anyway.
        for (GhostAxis & axis : _ghostaxes)
        {
            if (axis.axialTorqueB != 0.0f)
                (void)_distributeAxialTorque(axis, axis.substituentsB, axis.weightsB, axis.axialTorqueB);
            if (axis.axialTorqueC != 0.0f)
                (void)_distributeAxialTorque(axis, axis.substituentsC, axis.weightsC, axis.axialTorqueC);
        }
    }

    // Pass 3: one closed-form reaction per axis, restoring global force and
    // torque balance without ever differentiating the placement.
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
            .setPosition(binding.placement == GhostPlacement::AxialOffset
                             ? GhostParticle::computePositionAxial(B, C, binding.r)
                             : GhostParticle::computePositionByRotation(
                                   B, C, getParticle(binding.anchorRefIndex).getPosition(), binding.cos_delta,
                                   binding.sin_delta));
    }
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
    _stretchsprings.clear();
    _bendsprings.clear();
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
    _setupDihedralTorqueDistribution();
    _neighborSearchesDirty = false;
    // _setupConstraints();
    // _setupSelections();
}

// Decides once whether the torque distribution can run at all, so the
// per-step path is a single bool and any complaint is made once instead of
// on every one of 20000 steps.
//
// Two things have to hold, and neither is the user's mistake to discover at
// step 4000. The model must carry DIHEDRALAXIS records -- a .nc built before
// they existed, or from a force field that writes none, has nothing to
// distribute over. And the projection must be on: converting a force to a
// torque and back keeps only what turns about the axis, so with the raw
// force it would silently throw away the radial and axial parts rather than
// distribute them, which is a different model than the one asked for.
void SpringNetwork::_setupDihedralTorqueDistribution()
{
    _distributetorque = false;
    if (!_config.dihedral.distributetorque)
        return;

    if (!_config.dihedral.tangentialonly)
    {
        logging::warning("dihedral.distributetorque needs dihedral.tangentialonly: sharing out a torque keeps only "
                         "what turns about the axis, so the unprojected radial and axial force would be dropped "
                         "rather than distributed. Distribution stays off.");
        return;
    }

    const size_t withLists = getNumberOfAxesWithSubstituents();
    if (withLists == 0)
    {
        logging::warning("dihedral.distributetorque is on but no torsion axis carries substituents: this model was "
                         "built without DIHEDRALAXIS records (an older .nc, or a force field that writes none). "
                         "Distribution stays off; the rings keep pushing their own reference atoms.");
        return;
    }

    // Counted against the axes a PROPER torsion actually turns about, which
    // is neither of the two obvious denominators. _ghostaxes counts every
    // axis, including one per valence angle from BEND's axial ghosts. And
    // filtering those out still leaves the PLANARITY impropers, which are
    // ghost rings on an axis too but restrain a hub's flatness rather than
    // share a barrier over substituent pairs -- they have no substituent
    // lists and want none. Reporting against either would read as poor
    // coverage when coverage is in fact near complete.
    std::vector<bool> isProperTorsionAxis(_ghostaxes.size(), false);
    std::unordered_map<unsigned, unsigned> axisOfGhost;
    for (const GhostParticleBinding & binding : _ghostparticles)
        if (binding.placement == GhostPlacement::AxisRotation)
            axisOfGhost[binding.ownIndex] = binding.axisIndex;
    for (unsigned family = 0; family < DIHEDRAL_FAMILY_COUNT; ++family)
    {
        if (family == DIHEDRAL_PLANARITY)
            continue;
        for (const Spring & spring : _dihedralsprings[family])
            for (unsigned index : {spring.getParticle1().getId(), spring.getParticle2().getId()})
            {
                const auto found = axisOfGhost.find(index);
                if (found != axisOfGhost.end())
                    isProperTorsionAxis[found->second] = true;
            }
    }
    // Both numbers come from this one set, which is the point: computed two
    // different ways they disagreed, and a coverage figure whose numerator
    // can exceed its denominator tells the reader nothing.
    size_t torsionAxes = 0, covered = 0;
    for (size_t i = 0; i < _ghostaxes.size(); ++i)
    {
        if (!isProperTorsionAxis[i])
            continue;
        torsionAxes++;
        if (!_ghostaxes[i].substituentsB.empty() || !_ghostaxes[i].substituentsC.empty())
            covered++;
    }

    _distributetorque = true;
    logging::info("Dihedral torque distributed over bonded substituents on %zu of %zu torsion axes.", covered,
                  torsionAxes);
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
