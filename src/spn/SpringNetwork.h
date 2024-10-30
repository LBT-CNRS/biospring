#ifndef _SPRINGNETWORK_H_
#define _SPRINGNETWORK_H_

#include "configuration/Configuration.hpp"

#include "forcefield/ForceField.h"
#include "reduce/Reduce.h"

#include "IO/modern.hpp"

#include "grid/PotentialGrid.hpp"

#include "measure.hpp"
#include "nsearch.hpp"
#include "timeit.hpp"

#include "Constraint.h"
#include "InsertionVector.h"
#include "interactor/Interactor.h"
#include "Particle.h"
#include "Selection.h"
#include "Spring.h"
#include "Vector3f.h"
#ifdef OPENGL_SUPPORT
#include "viewer/SpringNetworkViewer.h"
#endif
#include <cstdlib>
#include <cstring>
#include <stdio.h>
#include <vector>

#include <iostream>
#include <memory>

class Interactor;
namespace biospring
{
namespace spn
{

class Spring;
class SpringNetworkViewer;

class SpringNetwork
{
  private:
    struct Grids
    {
        grid::PotentialGrid potential;
        grid::PotentialGrid density;
    };

    struct Energies
    {
        float spring;
        float electrostatic;
        float steric;
        float kinetic;
        float imp;

        void reset()
        {
            spring = 0.0;
            electrostatic = 0.0;
            steric = 0.0;
            kinetic = 0.0;
            imp = 0.0;
        }
    };

    // Stores FreeSASA parameters.
    struct FreeSASAState
    {
        bool isDynamic = false;
        unsigned step = 1000; // for dynamic FreeSASA
        double sasaTotal = 0.0;
    };

    struct NeighborSearch
    {
        using Container = std::vector<Particle>;
        using Searcher = nsearch::NeighborSearchDynamic<Container>;
        using SearcherPtr = std::unique_ptr<Searcher>;

        SearcherPtr steric;
        SearcherPtr electrostatic;
        SearcherPtr hydrophobic;
    };

    constexpr static auto make_nsearch = [](const NeighborSearch::Container & particles, float cutoff) {
        return std::make_unique<NeighborSearch::Searcher>(particles, cutoff);
    };

  public:
    SpringNetwork()
        : _interactors(), _initparticles(), _particles(), _staticparticules(), _dynamicparticules(),
          _chargedparticules(), _hydrophobicparticules(), _probeparticule(), _springs(), _nbiter(0), _end(false),
          _pause(false), _constraintenabled(false), _framerate(0.0), _ff(nullptr), _insertionVector(nullptr),
          _constraints(), _meanConstraintsDistances(0.0), _structid(_currentstructid++), _config(), _profiler()
    {
        _profiler.create_timer("main");
        _profiler.create_timer("samplerate");
    }

    virtual ~SpringNetwork();

    // ================================================================================

    const NeighborSearch & getNeighborSearch() { return _nsearch; }

    // ================================================================================

    // Gets/Sets interator.
    void addInteractor(Interactor* interactor) { _interactors.push_back(interactor);}
    Interactor* getInteractor(size_t index) {if (index < _interactors.size()) { return _interactors[index]; } return nullptr;}
    void removeInteractor(size_t index) { if (index < _interactors.size()) { _interactors.erase(_interactors.begin() + index); } }
    const std::vector<Interactor*>& getInteractors() const { return _interactors; }
    template <typename T>
    T* getInteractorInstance() const { return Interactor::getInteractorInstance<T>(_interactors); }

    // Returns current number of iterations so far (read-only).
    int getNbIterations() const { return _nbiter; }

    // Gets/Sets force field.
    // Not declared with const because need to modify the instance of forcefield
    // in InteractorMDDriver::processIMDInteractions to call
    // ForceField::setImpDoubleMembraneOffset 
    biospring::forcefield::ForceField * getForceField() const { return _ff.get(); };

    // Gets/Sets insertion vector.
    void setInsertionVector(unsigned aa1, unsigned aa2);
    const InsertionVector & getInsertionVector() const { return *_insertionVector; }
    InsertionVector & getInsertionVector() { return *_insertionVector; } // Non const usefull to modify their particles!

    // Gets/Sets total SASA value
    void setSASATotal(double sasa_total) {_freesasaState.sasaTotal = sasa_total;}
    double getSASATotal() const { return _freesasaState.sasaTotal; }
    void isFreeSASADynamic(bool isdyn) {_freesasaState.isDynamic = isdyn;}

    // ================================================================================
    // Returns frame rate (read-only).
    float getFrameRate() const { return _framerate; }

    // ================================================================================
    // Returns energies (read-only).
    float getKineticEnergy() const { return _energies.kinetic; }
    float getSpringEnergy() const { return _energies.spring; }
    float getStericEnergy() const { return _energies.steric; }
    float getElectrostaticEnergy() const { return _energies.electrostatic; }
    float getIMPEnergy() const { return _energies.imp; }

    // ================================================================================
    // Used to define the minimum IMP energy of all possible conformations at a 
    // given insertion angle (rotation around the axis of the insertion vector).
    // Used in the specific case of automatic sampling in rigid body mode. 
    // This energy must be defined just before calling the write_step function
    //  of a trajectory writer like the CSV one
    void setIMPEnergy(float ener) { _energies.imp = ener; }

    void writeNextStepNow();

    // ================================================================================
    // Grid getters & setters.
    // ================================================================================

    // Gets/Sets potential grid.
    biospring::grid::PotentialGrid & getPotentialGrid() { return _grids.potential; }
    const biospring::grid::PotentialGrid & getPotentialGrid() const { return _grids.potential; }

    // Gets/Sets density grid.
    biospring::grid::PotentialGrid & getDensityGrid() { return _grids.density; }
    const biospring::grid::PotentialGrid & getDensityGrid() const { return _grids.density; }

    // ================================================================================
    // Modification Methods.
    // Should be used only to convert `Topology` to `SpringNetwork`.

    // Adds a spring to the network.
    void addSpring(unsigned id1, unsigned id2, float equilibrium, float stiffness);

    void updateSpringState(unsigned id, bool isStatic);
    void addStaticSpring(unsigned id) { _staticsprings.push_back(id); }
    void addDynamicSpring(unsigned id) { _dynamicsprings.push_back(id); }
    void removeStaticSpring(unsigned id) { _staticsprings.erase(std::remove(_staticsprings.begin(), _staticsprings.end(), id), _staticsprings.end()); }
    void removeDynamicSpring(unsigned id) {_dynamicsprings.erase(std::remove(_dynamicsprings.begin(), _dynamicsprings.end(), id), _dynamicsprings.end()); }

    // Adds a particle to the network.
    void addParticle(const Particle & p);

    void updateParticleState(unsigned id, bool isStatic);
    void addStaticParticle(unsigned id) { _staticparticules.push_back(id); }
    void addDynamicParticle(unsigned id) { _dynamicparticules.push_back(id); }
    void removeStaticParticle(unsigned id) { _staticparticules.erase(std::remove(_staticparticules.begin(), _staticparticules.end(), id), _staticparticules.end()); }
    void removeDynamicParticle(unsigned id) {_dynamicparticules.erase(std::remove(_dynamicparticules.begin(), _dynamicparticules.end(), id), _dynamicparticules.end()); }

    // Removes all particles.
    void clearParticles(void);

    // Empties the network.
    void clear();

    // ================================================================================
    // Gets springs/particles.

    // Returns the list of springs.
    const std::vector<Spring> & getSprings() const { return _springs; }

    // Returns ith spring in spring list.
    std::vector<Spring>::const_reference getSpring(unsigned index) const { return _springs[index]; }
    std::vector<Spring>::reference getSpring(unsigned index) { return _springs[index]; }

    // Returns ith particle in particle list.
    std::vector<Particle>::const_reference getParticle(unsigned index) const { return _particles[index]; }
    std::vector<Particle>::reference getParticle(unsigned index) { return _particles[index]; }

    // Returns a particle using its id (see Particle::getId)
    std::vector<Particle>::const_reference getParticleFromId(unsigned extid) const;
    std::vector<Particle>::reference getParticleFromId(unsigned extid);

    // ================================================================================
    // Shortcuts to configuration values.

    int getMaxIteration() const { return _config.sim.nbsteps; }

    float getSampleRate() const { return _config.sim.samplerate; }

    float getTimeStep() const { return _config.sim.timestep; }

    float getGridScale() const { return _config.steric.gridscale; }

    float getViscosity() const { return _config.viscosity.value; }
    void setViscosity(float visc) { _config.viscosity.value = visc; } // Viscosity able to be updated during simulation

    float getStericCutoff() const { return _config.steric.cutoff; }
    float getElectrostaticCutoff() const { return _config.electrostatic.cutoff; }
    float getHydrophobicCutoff() const { return _config.hydrophobicity.cutoff; }

    bool isSpringEnabled() const { return _config.spring.enable; }
    bool isViscosityEnabled() const { return _config.viscosity.enable; }
    bool isStericEnabled() const { return _config.steric.enable; }
    bool isElectrostaticEnabled() const { return _config.electrostatic.enable; }
    bool isElectrostaticCoulombEnabled() const { return _config.electrostatic.enable; }
    bool isElectrostaticFieldEnabled() const { return _config.potentialgrid.enable; }
    bool isIMPEnabled() const { return _config.imp.enable; }
    bool isDensityGridEnabled() const { return _config.densitygrid.enable; }
    bool isConstraintEnabled() const { return _constraintenabled; }
    bool isHydrophobicityEnabled() const { return _config.hydrophobicity.enable; }
    bool isProbeEnabled() const { return _config.probe.enable; }
    bool isProbeElectrostaticEnabled() const { return _config.probe.enableelectrostatic; }
    bool isProbeStericEnabled() const { return _config.probe.enablesteric; }
    bool isProbeElectrostaticFieldEnabled() const { return false; }

    bool isPDBTrajectoryWriterEnabled() const { return _config.pdbtraj.enable; }
    int getPDBTrajectoryWriterFreq() const { return _config.pdbtraj.frequency; }
    bool isXTCTrajectoryWriterEnabled() const { return _config.xtctraj.enable; }
    int getXTCTrajectoryWriterFreq() const { return _config.xtctraj.frequency; }
    bool isCSVTrajectoryWriterEnabled() const { return _config.csvsample.enable; }
    int getCSVTrajectoryWriterFreq() const { return _config.csvsample.frequency; }

    // ================================================================================
    bool isInsertionVectorEnabled() const { return _insertionVector != nullptr; }

    bool isRigidBodyEnabled() const { return _config.rigidbody.enable; }
    bool isImpalaSamplingEnabled() const { return _config.rigidbody.enablesampling; }
    bool isMonteCarloEnabled() const { return _config.rigidbody.enablemontecarlo; }
    double getMonteCarloTemperature() const { return _config.rigidbody.montecarlo_temperature; }
    void setMonteCarloTemperature(float temp) { _config.rigidbody.montecarlo_temperature = temp; }
    double getMonteCarloTranslationNorm() const { return _config.rigidbody.montecarlo_translation_norm; }
    void setMonteCarloTranslationNorm(float norm) { _config.rigidbody.montecarlo_translation_norm = norm; }
    double getMonteCarloRotationNorm() const { return _config.rigidbody.montecarlo_rotation_norm; }
    void setMonteCarloRotationNorm(float norm) { _config.rigidbody.montecarlo_rotation_norm = norm; }

    // ================================================================================

    bool writePDB(const char * pdbin, const char * pdbout);

    virtual void run();
    virtual void computeStep();
    virtual void computeForces();
    virtual void computeSpringForces();
    virtual void computeParticleForces();
    virtual void updateParticlePositions();

    unsigned getNumberOfSprings() const { return _springs.size(); }
    unsigned getNumberOfParticles() const { return _particles.size(); }

    // ================================================================================
    // Returns subsets of particles.

    const vector<Particle> & getParticles() const { return _particles; }
    const vector<Particle> & getInitParticles() const { return _initparticles; }

    // Returns subsets of particle ids.
    vector<unsigned> getDynamicParticles() const { return _dynamicparticules; }
    vector<unsigned> getChargedParticles() const { return _chargedparticules; }
    vector<unsigned> getStaticParticles() const { return _staticparticules; }
    vector<unsigned> getHydrophobicParticles() const { return _hydrophobicparticules; }

    // Returns the particle's centroid.
    auto getCentroid() const { return biospring::measure::centroid(_particles); }

    // ================================================================================
    // Run-related methods.
    virtual void initRun();
    virtual void endRun();
    virtual void idleRun();

    void setPause(bool pause) { _pause = pause; }
    bool getPause() const { return _pause; }

    bool isEnd() const { return _end; }
    void setEnd(bool end) { _end = end; }

    // ================================================================================

    // ================================================================================
    //
    // Setup methods.
    //
    // ================================================================================

    // Sets up the spring network.
    // Copies the configuration so it has ownership of it.
    void setup(const configuration::Configuration & conf);

  protected:
    void _setupSteric();
    void _setupHydrophobic();
    void _setupForceField();
    void _setupElectrostatic();
    void _setupDensityGrid();
    void _setupProbe();
    void _setupTrajectories();
    void _setupInsertionVector();
    void _setupSelections();
    void _setupConstraints();

    // ================================================================================
    //
    // Constraint methods.
    //
    // ================================================================================

  public:
    std::vector<Constraint *> getConstraints(void) { return _constraints; }
    void addConstraint(Constraint * constraint) { _constraints.push_back(constraint); }
    void applyConstraints();

    // ================================================================================

    virtual void getParticlePosition(unsigned i, float position[3]) const;

    virtual void setForce(unsigned i, float force[3]);

#ifdef OPENGL_SUPPORT
    SpringNetworkViewer * _viewer;
    void setSpringNetworkViewer(SpringNetworkViewer * viewer) { _viewer = viewer; }
    SpringNetworkViewer * getSpringNetworkViewer() const { return _viewer; }
#endif

  protected:
    std::vector<Interactor*> _interactors;
    std::vector<Particle> _initparticles;
    std::vector<Particle> _particles;

    std::vector<unsigned> _staticparticules;
    std::vector<unsigned> _dynamicparticules;
    std::vector<unsigned> _chargedparticules;
    std::vector<unsigned> _hydrophobicparticules;

    Particle _probeparticule;

    std::vector<Spring> _springs;
    std::vector<unsigned> _staticsprings;
    std::vector<unsigned> _dynamicsprings;

    Energies _energies;
    NeighborSearch _nsearch;

    int _nbiter;
    bool _end;
    bool _pause;

    Grids _grids;

    bool _constraintenabled;

    double _framerate;

    FreeSASAState _freesasaState;

    std::unique_ptr<forcefield::ForceField> _ff;

    io::modern::TrajectoryManager _trajectories;

    std::unique_ptr<InsertionVector> _insertionVector;

    vector<Constraint *> _constraints;
    float _meanConstraintsDistances;

    static unsigned _currentstructid;
    unsigned _structid;

    configuration::Configuration _config;

    timeit::Profiler _profiler;

    // ================================================================================

    // Computes insertion vector's angle and updates barycentre.
    void _updateInsertionVector();

    // Resets energy values to 0.
    void _resetEnergies();

    // Writes current positions using trajectory writers.
    void _writeNextStep();

    // Returns true if it's time to display current frame data.
    bool _isTimeToLogData() const { return _nbiter % (static_cast<unsigned>(getSampleRate())) == 0; }

    // Displays current frame informations on logging channel.
    void _displayFrameData();

    // Calculates frame rate.
    void _updateFrameRate();

    // Returns true if the run is endless (i.e. no maximum number of iterations).
    bool _isInfiniteRun() { return getMaxIteration() < 0; }

    // Returns true if reached the maximum number of iterations.
    bool _hasReachedEndOfRun() { return !_isInfiniteRun() && _nbiter == getMaxIteration(); }
};

} // namespace spn
} // namespace biospring

#endif
