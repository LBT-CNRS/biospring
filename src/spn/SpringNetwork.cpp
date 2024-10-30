#include "SpringNetwork.h"
#include "logging.h"
#include "measure.hpp"

#include "forcefield/ForceField.h"
#include "forcefield/ForceFieldElectrostaticCoulombAndStericLennardJones_8_6Amber.h"
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
#include <math.h>
#include <memory>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>

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
    float springenergy = 0.0;
#ifdef OPENMP_SUPPORT
#pragma omp parallel default(shared)
#endif
    {
#ifdef OPENMP_SUPPORT
#pragma omp for reduction(+ : springenergy)
#endif
        for (size_t i = 0; i < _dynamicsprings.size(); i++)
        {
            Spring & s = getSpring(_dynamicsprings[i]);
            s.applyForceToParticle(*_ff);
            springenergy += s.getEnergy();
        }
    }
    _energies.spring = springenergy;
}

// Calculate forces that applies on dynamic particles.

/// @brief Compute Particles Forces
void SpringNetwork::computeParticleForces()
{
    float electrostatic_energy = 0.0;
    float steric_energy = 0.0;
    float kinetic_energy_probe = 0.0;
    float imp_energy = 0.0;
    // float hydrophobicity_energy = 0.0;

    // =========================================================================
    // Calculate forces on dynamic particles.
    // =========================================================================

#ifdef OPENMP_SUPPORT
#pragma omp parallel default(shared)
#endif
    {
#ifdef OPENMP_SUPPORT
#pragma omp for reduction(+ : electrostatic_energy, steric_energy, kinetic_energy_probe, imp_energy) schedule(static)
#endif
        for (size_t i = 0; i < _dynamicparticules.size(); i++)
        {
            Particle & p = getParticle(_dynamicparticules[i]);

            // ===================
            // Compute electrostatic forces.
            if (isElectrostaticEnabled())
            {
                if (isElectrostaticCoulombEnabled())
                {
                    // p->addElectrostaticForceNoGrid(getElectrostaticCutoff());
                    if (p.isCharged())
                    {
                        p.addElectrostaticForce();
                        electrostatic_energy += p.getElectrostaticEnergy();
                    }
                }
                if (isElectrostaticFieldEnabled())
                {
                    p.addElectrostaticFieldForce();
                    electrostatic_energy += p.getElectrostaticEnergy();
                }
            }

            // ===================
            // Compute density field forces.
            if (isDensityGridEnabled())
            {
                p.addDensityFieldForce();
            }

            // ===================
            // Compute steric forces.
            if (isStericEnabled())
            {
                p.addStericForce();
                steric_energy += p.getStericEnergy();
            }

            // ===================
            // Applies viscosity forces.
            if (isViscosityEnabled())
            {
                p.applyViscosity(getViscosity());
            }

            // ===================
            // Compute interactions with the probe.
            if (isProbeEnabled())
            {
                if (isProbeStericEnabled())
                {
                    p.addStericProbeForce(_probeparticule);
                    steric_energy += _probeparticule.getStericEnergy();
                    steric_energy += p.getStericEnergy();
                }
                if (isProbeElectrostaticEnabled())
                {
                    p.addElectrostaticProbeForce(_probeparticule);
                    electrostatic_energy += _probeparticule.getElectrostaticEnergy();
                    electrostatic_energy += p.getElectrostaticEnergy();
                }
                if (isProbeElectrostaticFieldEnabled())
                {
                    _probeparticule.addElectrostaticFieldForce();
                    electrostatic_energy += _probeparticule.getElectrostaticEnergy();
                }
                if (isViscosityEnabled())
                    _probeparticule.applyViscosity(getViscosity());

                _probeparticule.IntegrateEuler(getTimeStep());

                kinetic_energy_probe += _probeparticule.getKineticEnergy();
                p.resetForce();
            }

            // ===================
            // Compute Impala forces.
            if (isIMPEnabled())
            {
                p.addIMPForce();
                imp_energy += p.getIMPEnergy();
            }

            if (isRigidBodyEnabled() && p.isRigid() && !isImpalaSamplingEnabled() && !isMonteCarloEnabled())
                rigidbody::RigidBody::computeParticleForceAndTorque(p);

            // ===================
            // Compute hydrophobicity forces.
            if (isHydrophobicityEnabled())
            {
                p.addHydrophobicityForce();
            }

            p.setPreviousForce();
        } // omp for loop
    }     // omp parallel

    

    _energies.electrostatic = electrostatic_energy;
    _energies.steric = steric_energy;
    _energies.kinetic += kinetic_energy_probe;
    _energies.imp = imp_energy;
    // TODO: update hydrophobicity energy ? If yes, remember to update the line #pragma omp for reduction.
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
        for (int i = 0; i < (int)_dynamicparticules.size(); i++)
        {
            Particle & p = getParticle(_dynamicparticules[i]);
            if (p.isRigid())
                rigidbody::RigidBody::integrateParticleVelocity(p, i, getTimeStep());
            else
                p.IntegrateEuler(getTimeStep());

            // Check if position exploses, if one of float is NaN
            float x = p.getPosition().getX();
            float y = p.getPosition().getY();
            float z = p.getPosition().getZ();
            if (isnan(x) || isnan(y) || isnan(z))
            {
                logging::die("Found undefined position.");
            }

            kinetic_energy_particle += p.getKineticEnergy();
            p.resetForce();
        } // omp for loop
    }     // omp parallel
    _energies.kinetic += kinetic_energy_particle;
}

/// @brief Compute particles force and, if activated, springs forces.
void SpringNetwork::computeForces()
{
    if (isSpringEnabled())
        computeSpringForces();
    computeParticleForces();
}

void SpringNetwork::computeStep()
{
    idleRun();
    _meanConstraintsDistances = 0.0;

    computeForces();

    if (isConstraintEnabled())
        applyConstraints();

    if (isRigidBodyEnabled())
        rigidbody::RigidBodiesManager::SolveRigidBodiesDynamic(); 

    updateParticlePositions();

    if (isInsertionVectorEnabled())
        _updateInsertionVector();
}

void SpringNetwork::run()
{
    initRun();

#ifdef MDDRIVER_SUPPORT
    usleep(100000);
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
        {
            interactor->stopInteractionThread();
        }
    }
}

void SpringNetwork::idleRun()
{
    while (_pause)
        usleep(10000);

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
        logging::info("Spring energy: %5.2f kJ.mol-1", _energies.spring);
    if (isElectrostaticEnabled())
        logging::info("Electrostatic energy: %5.2f kJ.mol-1", _energies.electrostatic);
    if (isStericEnabled())
        logging::info("Steric energy: %5.2f kJ.mol-1", _energies.steric);
    if (isIMPEnabled())
        logging::info("IMP energy: %5.2f kJ.mol-1", _energies.imp);
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

void SpringNetwork::_writeNextStep() { _trajectories.write_step(_nbiter); }

void SpringNetwork::writeNextStepNow() { _trajectories.write_step(); }

std::vector<Particle>::const_reference SpringNetwork::getParticleFromId(unsigned id) const
{
    for (unsigned i = 0; i < _particles.size(); i++)
    {
        if (_particles[i].getId() == id)
            return _particles[i];
    }
    throw std::out_of_range("SpringNetwork::getParticleFromId: Particle id not found.");
}

std::vector<Particle>::reference SpringNetwork::getParticleFromId(unsigned id)
{
    for (unsigned i = 0; i < _particles.size(); i++)
    {
        if (_particles[i].getId() == id)
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
    // If the spring does not already exist.
    if (!_particles[id1].isInSpringNeighbors(id2))
    {
        Particle & p1 = _particles[id1];
        Particle & p2 = _particles[id2];
        Spring s = Spring(p1, p2, equilibrium, stiffness);
        s.setId(_springs.size());
        _springs.push_back(s);
        p1.addToSpringNeighbors(id2, &_springs.back());
        p2.addToSpringNeighbors(id1, &_springs.back());

        if (p1.isStatic() && p2.isStatic())
            addStaticSpring(s.getId());
        else
            addDynamicSpring(s.getId());
        
        
    }
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
    Particle p = Particle(source);
    p.setSpringNetwork(this);
    p.setId(_particles.size());

    p.setInternalStructId(_structid);

    if (p.isStatic())
        addStaticParticle(p.getId());
    else
        addDynamicParticle(p.getId());

    if (p.isCharged())
        _chargedparticules.push_back(p.getId());

    if (p.isHydrophobic())
        _hydrophobicparticules.push_back(p.getId());

    _particles.push_back(p);
    _initparticles.push_back(p);
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
    _particles.clear();
    _springs.clear();
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
        getParticleFromId(_config.ivector.vector[0]);
        getParticleFromId(_config.ivector.vector[1]);
    }
    catch (const std::out_of_range & e)
    {
        logging::die("Insertion vector could not be initialized. One or both of the pair of particles "
                     "defining the insertion vector have an invalid identifier.");
    }

    //   Check that the two elements of the vector are identical
    if (_config.ivector.vector[0] == _config.ivector.vector[1])
    {
        logging::warning("The two values defining the insertion vector are identical: %zu and %zu",
                         _config.ivector.vector[0], _config.ivector.vector[1]);
    }

    _insertionVector = std::make_unique<InsertionVector>(*this, getParticleFromId(_config.ivector.vector[0]),
                                                         getParticleFromId(_config.ivector.vector[1]));
    
    const Particle p1 = _insertionVector->getParticle(0);
    const Particle p2 = _insertionVector->getParticle(1);
    logging::info("Insertion vector set between %s%d:%s and %s%d:%s",
        p1.getResName().c_str(), p1.getResId(), p1.getName().c_str(),
        p2.getResName().c_str(), p2.getResId(), p2.getName().c_str());
}

void SpringNetwork::applyConstraints()
{
    float sumDistances = 0.0;
    for (unsigned i = 0; i < _constraints.size(); i++)
    {
        _constraints[i]->apply();
        sumDistances += _constraints[i]->getDistance();
    }
    _meanConstraintsDistances = sumDistances / _constraints.size();
}

void SpringNetwork::clearParticles(void) { _particles.clear(); }

// =====================================================================================
//
// Setup methods
//
// =====================================================================================

void SpringNetwork::setup(const configuration::Configuration & conf)
{
    _config = conf;

    _setupForceField();
    _setupSteric();
    _setupElectrostatic();
    _setupDensityGrid();
    _setupProbe();
    _setupInsertionVector();
    _setupTrajectories();
    // _setupConstraints();
    // _setupSelections();
}

void SpringNetwork::_setupSteric()
{
    if (isStericEnabled())
    {
        if (getStericCutoff() < 1e-6)
            throw std::runtime_error("Steric cutoff must be > 0");
        _nsearch.steric = make_nsearch(_particles, getStericCutoff());
    }
}

void SpringNetwork::_setupHydrophobic()
{
    if (isHydrophobicityEnabled())
    {
        if (getHydrophobicCutoff() < 1e-6)
            throw std::runtime_error("Hydrophobic cutoff must be > 0");
        _nsearch.hydrophobic = make_nsearch(_particles, getHydrophobicCutoff());
    }
}

void SpringNetwork::_setupForceField()
{
    const std::string steric = _config.steric.mode;

    if (steric == "lennard-jones-8-6Lewitt")
        _ff = std::make_unique<forcefield::ForceFieldElectrostaticCoulombAndStericLennardJones_8_6Lewitt>();
    else if (steric == "lennard-jones-8-6Zacharias")
        _ff = std::make_unique<forcefield::ForceFieldElectrostaticCoulombAndStericLennardJones_8_6Zacharias>();
    else if (steric == "lennard-jones-8-6Amber")
        _ff = std::make_unique<forcefield::ForceFieldElectrostaticCoulombAndStericLennardJones_8_6Amber>();
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
    if (isElectrostaticEnabled())
    {
        if (isElectrostaticFieldEnabled())
        {
            const std::string dxpath = _config.potentialgrid.path;
            logging::info("Reading electrostatic map from DX file '%s'", dxpath.c_str());
            _grids.potential = opendx::readGrid(dxpath);
        }
        else
        {
            if (getElectrostaticCutoff() < 1e-6)
                throw std::runtime_error("Electrostatic cutoff must be > 0");

            // TODO:: can we initialize it with all particles?
            _nsearch.electrostatic = make_nsearch(_particles, getElectrostaticCutoff());
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
        _probeparticule.setId(_particles.size());
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

void SpringNetwork::_setupSelections() { throw "Not Implemented Error"; }
void SpringNetwork::_setupConstraints() { throw "Not Implemented Error"; }

} // namespace spn
} // namespace biospring
