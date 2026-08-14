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
#include "GhostParticle.h"
#include "InsertionVector.h"
#include "interactor/Interactor.h"
#include "Particle.h"
#include "Selection.h"
#include "Spring.h"
#include "Vector3f.h"
#include <cstdlib>
#include <cstring>
#include <stdio.h>
#include <vector>

#include <iostream>
#include <memory>
#include <utility>

class Interactor;
class SpringNetworkViewer;
namespace biospring
{
namespace spn
{

class Spring;

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
        // Only ever the untouched --rigidbody/ENM mesh's own energy now:
        // stretch/bend/dihedral each have their own dedicated channel below
        // instead of being folded in here, since "spring energy" is only a
        // meaningful, self-contained quantity for a plain ENM/rigid-body
        // network -- inflating it with the real force-field corrections
        // made it impossible to read any of the three back out again.
        float spring = 0.0f;
        float stretch = 0.0f;
        float bend = 0.0f;
        float dihedral = 0.0f;
        float electrostatic = 0.0f;
        float steric = 0.0f;
        float kinetic = 0.0f;
        float imp = 0.0f;
        float hydrophobic = 0.0f;
        float hbond = 0.0f;

        void reset()
        {
            spring = 0.0;
            stretch = 0.0;
            bend = 0.0;
            dihedral = 0.0;
            electrostatic = 0.0;
            steric = 0.0;
            kinetic = 0.0;
            imp = 0.0;
            hydrophobic = 0.0;
            hbond = 0.0;
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
        using Searcher = nsearch::NeighborSearch<Container>;
        using SearcherPtr = std::unique_ptr<Searcher>;

        SearcherPtr steric;
        SearcherPtr electrostatic;
        SearcherPtr hydrophobic;

        // Shared grid holding every donor AND acceptor particle together
        // (not split in two): a donor queries it and filters the returned
        // candidates by isAcceptor() itself, mirroring how steric's single
        // "all particles" grid is filtered per-pair rather than pre-split.
        SearcherPtr hbond;
    };

    static NeighborSearch::SearcherPtr make_nsearch(const NeighborSearch::Container & particles, float cutoff,
                                                     float skin)
    {
        return std::make_unique<NeighborSearch::Searcher>(particles, cutoff, skin);
    }

    static NeighborSearch::SearcherPtr make_nsearch(const NeighborSearch::Container & particles, float cutoff,
                                                     std::vector<size_t> included_indices, float skin)
    {
        return std::make_unique<NeighborSearch::Searcher>(particles, cutoff, std::move(included_indices), skin);
    }

  public:
    SpringNetwork()
        : _viewer(nullptr), _interactors(), _initparticles(), _particles(), _staticparticules(), _dynamicparticules(),
          _chargedparticules(), _hydrophobicparticules(), _probeparticule(), _springs(), _staticsprings(),
          _dynamicsprings(), _dihedralphisprings(), _dihedralpsisprings(), _dihedralomegasprings(),
          _dihedralsidechainsprings(), _dihedralplanaritysprings(), _stretchsprings(), _bendsprings(),
          _ghostparticles(),
          _springForceScratch(), _stericPairScratch(), _electrostaticPairScratch(),
          _hydrophobicPairScratch(), _hydrogenBondCoreRepulsionPairScratch(), _hydrogenBondPartner(),
          _hydrogenBondForceScratch(), _energies(), _nsearch(),
          _neighborSearchesDirty(false),
          _nbiter(0), _end(false), _pause(false), _grids(), _constraintenabled(false), _framerate(0.0),
          _freesasaState(), _ff(nullptr), _trajectories(), _insertionVector(nullptr), _constraints(),
          _meanConstraintsDistances(0.0), _structid(_currentstructid++), _config(), _profiler()
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
    void removeInteractor(size_t index) { if (index < _interactors.size()) { _interactors.erase(_interactors.begin() + static_cast<std::ptrdiff_t>(index)); } }
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
    // The untouched --rigidbody/ENM mesh only -- see Energies' own comment.
    float getSpringEnergy() const { return _energies.spring; }
    float getStretchEnergy() const { return _energies.stretch; }
    float getBendEnergy() const { return _energies.bend; }
    float getDihedralEnergy() const { return _energies.dihedral; }
    float getStericEnergy() const { return _energies.steric; }
    float getElectrostaticEnergy() const { return _energies.electrostatic; }
    float getIMPEnergy() const { return _energies.imp; }
    float getHydrophobicEnergy() const { return _energies.hydrophobic; }
    float getHydrogenBondEnergy() const { return _energies.hbond; }

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

    // Adds a dihedral ghost spring to the network -- always a new spring
    // (unlike addSpring, never checked against existing spring-neighbours),
    // and deliberately not registered as a spring-neighbour: a ghost spring
    // connects two real substituent atoms that are a real 1-4 pair (never a
    // real 1-2/1-3 bond), and BioSpring's nonbonded exclusion is driven
    // entirely by spring-neighbour membership (see Particle::isInSpringNeighbors) --
    // registering ghost springs there would silently exclude 1-4 pairs from
    // nonbonded forces as a side effect, which is a separate physical
    // modelling decision (matching AMBER's scaled 1-4 nonbonded convention)
    // that this feature does not make. See doc/BondedForceFieldSprings.md.
    void addDihedralPhiSpring(unsigned id1, unsigned id2, float equilibrium, float stiffness, float dcOffset = 0.0f);
    void addDihedralPsiSpring(unsigned id1, unsigned id2, float equilibrium, float stiffness, float dcOffset = 0.0f);
    void addDihedralOmegaSpring(unsigned id1, unsigned id2, float equilibrium, float stiffness,
                                float dcOffset = 0.0f);
    void addDihedralSidechainSpring(unsigned id1, unsigned id2, float equilibrium, float stiffness,
                                    float dcOffset = 0.0f);
    void addDihedralPlanaritySpring(unsigned id1, unsigned id2, float equilibrium, float stiffness,
                                    float dcOffset = 0.0f);

    // Adds a STRETCH spring (see _stretchsprings' own comment) -- a new
    // spring between the same 2 real atoms as the existing (now-zeroed)
    // --rigidbody spring, so also not registered as a spring-neighbour
    // itself: that role stays with the original --rigidbody spring, which
    // is kept (just zeroed) specifically so nonbonded exclusion is
    // unaffected by this split.
    void addStretchSpring(unsigned id1, unsigned id2, float equilibrium, float stiffness);

    // Adds a BEND ghost-ghost spring (see _bendsprings' own comment) --
    // always new, like the dihedral springs, never a spring-neighbour
    // either (the ghosts are virtual sites, not real 1-3 atoms; nonbonded
    // exclusion for the real 1-3 pair is still handled by the existing,
    // now-zeroed --rigidbody spring between the real atoms themselves).
    void addBendSpring(unsigned id1, unsigned id2, float equilibrium, float stiffness);

    void updateSpringState(unsigned id, bool isStatic);
    void addStaticSpring(unsigned id) { _staticsprings.push_back(id); }
    void addDynamicSpring(unsigned id) { _dynamicsprings.push_back(id); }
    void removeStaticSpring(unsigned id) { _staticsprings.erase(std::remove(_staticsprings.begin(), _staticsprings.end(), id), _staticsprings.end()); }
    void removeDynamicSpring(unsigned id) {_dynamicsprings.erase(std::remove(_dynamicsprings.begin(), _dynamicsprings.end(), id), _dynamicsprings.end()); }

    // Adds a particle to the network.
    void addParticle(const Particle & p);

    // Adds a massless virtual-site ("ghost") particle to the network,
    // bound to 3 already-added real anchor particles (by index -- see
    // GhostParticle.h for why an index, not a pointer/reference, is used).
    // Must be called during the same "add every particle" phase as
    // addParticle (before any spring is added -- enforced the same way,
    // via addParticle's own check), and anchorBIndex/anchorCIndex/
    // anchorRefIndex must already refer to previously-added particles.
    // Creates the ghost's own Particle entry (isStatic()=true, mass=0,
    // initial position placed from the anchors' current positions) and
    // returns its index. See redistributeGhostForces/updateGhostPositions
    // for the two per-step operations this binding drives.
    unsigned addGhostParticle(unsigned anchorBIndex, unsigned anchorCIndex, unsigned anchorRefIndex, float r,
                              float theta_deg, float delta_deg);

    // Redistributes every ghost particle's currently accumulated force
    // (from ordinary spring force computation, e.g. a dihedral ghost
    // spring between two ghost particles) onto its 3 anchors, then resets
    // the ghost's own force to zero (it is static, so it never reaches
    // updateParticlePositions's normal per-dynamic-particle resetForce()).
    // Called once per step, right after computeForces().
    void redistributeGhostForces();

    // Recomputes every ghost particle's position from its 3 anchors'
    // CURRENT positions. Called once per step, right after
    // updateParticlePositions() (i.e. after the anchors themselves have
    // been integrated).
    void updateGhostPositions();

    const std::vector<GhostParticleBinding> & getGhostParticles() const { return _ghostparticles; }

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

    // Returns each dihedral ghost-spring family's list, for NetCDF I/O.
    const std::vector<Spring> & getDihedralPhiSprings() const { return _dihedralphisprings; }
    const std::vector<Spring> & getDihedralPsiSprings() const { return _dihedralpsisprings; }
    const std::vector<Spring> & getDihedralOmegaSprings() const { return _dihedralomegasprings; }
    const std::vector<Spring> & getDihedralSidechainSprings() const { return _dihedralsidechainsprings; }
    const std::vector<Spring> & getDihedralPlanaritySprings() const { return _dihedralplanaritysprings; }
    const std::vector<Spring> & getStretchSprings() const { return _stretchsprings; }
    const std::vector<Spring> & getBendSprings() const { return _bendsprings; }

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

    // Steric force scale (also pushed into the ForceField via setStericScale,
    // see SpringNetwork::setup). Named getGridScale for historical reasons;
    // despite the name it has nothing to do with potentialgrid or
    // densitygrid, which have their own independent scale (see
    // ForceField::getForceFieldScale and getDensityGridScale below).
    float getGridScale() const { return _config.steric.gridscale; }

    float getDensityGridScale() const { return _config.densitygrid.scale; }

    float getViscosity() const { return _config.viscosity.value; }
    void setViscosity(float visc) { _config.viscosity.value = visc; } // Viscosity able to be updated during simulation

    float getStericCutoff() const { return _config.steric.cutoff; }
    float getElectrostaticCutoff() const { return _config.electrostatic.cutoff; }
    float getHydrophobicCutoff() const { return _config.hydrophobicity.cutoff; }
    float getHydrogenBondCutoff() const { return _config.hbond.cutoff; }
    float getNeighborSkin() const { return _config.sim.neighborskin; }

    // -1 if particle `index` is not currently engaged in an exclusive
    // hydrogen bond, otherwise the index of its current partner. See
    // _assignHydrogenBondPairs / Particle::addHydrogenBondCoreRepulsion.
    int getHydrogenBondPartner(size_t index) const
    {
        return index < _hydrogenBondPartner.size() ? _hydrogenBondPartner[index] : -1;
    }

    bool isSpringEnabled() const { return _config.spring.enable; }
    // Per-family runtime debug toggles for the bonded-force-field springs
    // -stretching/-bending/-dihedral* already decided, at build time, to
    // create (see Configuration.hpp's own comment on these 5 settings).
    bool isBendingEnabled() const { return _config.bending.enable; }
    bool isDihedralPhiEnabled() const { return _config.dihedralphi.enable; }
    bool isDihedralPsiEnabled() const { return _config.dihedralpsi.enable; }
    bool isDihedralOmegaEnabled() const { return _config.dihedralomega.enable; }
    bool isDihedralChiEnabled() const { return _config.dihedralchi.enable; }
    bool isViscosityEnabled() const { return _config.viscosity.enable; }
    bool isStericEnabled() const { return _config.steric.enable; }
    bool isElectrostaticEnabled() const { return _config.electrostatic.enable; }
    bool isElectrostaticCoulombEnabled() const { return _config.electrostatic.enable; }
    bool isElectrostaticFieldEnabled() const { return _config.potentialgrid.enable; }
    bool isIMPEnabled() const { return _config.imp.enable; }
    bool isDensityGridEnabled() const { return _config.densitygrid.enable; }
    bool isConstraintEnabled() const { return _constraintenabled; }
    bool isHydrophobicityEnabled() const { return _config.hydrophobicity.enable; }
    bool isHydrogenBondEnabled() const { return _config.hbond.enable; }
    bool isProbeEnabled() const { return _config.probe.enable; }
    bool isProbeElectrostaticEnabled() const { return _config.probe.enableelectrostatic; }
    bool isProbeStericEnabled() const { return _config.probe.enablesteric; }
    bool isProbeElectrostaticFieldEnabled() const { return false; }
    bool isProbeParticle(size_t index) const
    {
        return isProbeEnabled() && _probeparticule.getId() >= 0 &&
               index == static_cast<size_t>(_probeparticule.getId());
    }
    const Particle & getProbeParticle() const { return _probeparticule; }

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
    virtual void computeStretchForces();
    virtual void computeBendForces();
    virtual void computeDihedralForces();
    virtual void computeHydrogenBondForces();
    virtual void computeParticleForces();
    virtual void updateParticlePositions();

  private:
    // Shared parallel-compute/serial-accumulate loop behind
    // computeStretch/Bend/DihedralForces -- see its definition for the
    // pattern and why the accumulation stays serial. Returns the summed
    // energy of the collection.
    float _computeSpringCollectionForces(std::vector<Spring> & springs, bool ignoreDynamicState,
                                         bool subtractDcOffset);

  public:

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
    void _setupHydrogenBond();
    void _setupForceField();
    void _setupElectrostatic();
    void _setupDensityGrid();
    void _setupProbe();
    void _setupTrajectories();
    void _setupInsertionVector();
    void _setupSelections();
    void _setupConstraints();
    std::vector<size_t> _chargedParticleIndexes() const;
    std::vector<size_t> _hydrophobicParticleIndexes() const;
    std::vector<size_t> _donorAcceptorParticleIndexes() const;

    // Updates _hydrogenBondPartner: breaks any active pair that has drifted
    // beyond the hbond cutoff, then matches newly-free donors and acceptors
    // by mutual nearest neighbor (a pair forms only if each is the other's
    // closest still-free candidate within cutoff -- the same reciprocal-best-
    // hit criterion used to detect orthologs between two gene sets). Already-
    // engaged particles are excluded from this matching entirely.
    void _assignHydrogenBondPairs();

    void _excludeProbeFromNeighborSearch(NeighborSearch::Searcher & searcher);
    void _updateNeighborSearches();
    void _markNeighborSearchesDirty();
    void _syncProbeParticle();
    void _rebuildSpringNeighbors();

    // Resizes the nonbonded pair scratch buffers to the current number of
    // dynamic particles and clears their contents, reusing prior capacity.
    void _resizeNonbondedPairScratch();

    // Applies deferred nonbonded pair contributions to their target
    // particles and adds their energy to `energy`. Must run serially, after
    // the parallel region that filled `scratch`, since two buckets may defer
    // a contribution to the same target particle.
    void _applyNonbondedPairScratch(const std::vector<std::vector<spn::DeferredNonbondedContribution>> & scratch,
                                     float & energy);

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

    // The opaque viewer pointer is always present so enabling the optional
    // viewer never changes SpringNetwork's ABI or class layout.
    ::SpringNetworkViewer * _viewer;
    void setSpringNetworkViewer(::SpringNetworkViewer * viewer) { _viewer = viewer; }
    ::SpringNetworkViewer * getSpringNetworkViewer() const { return _viewer; }

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

    // Dihedral ghost springs, kept in their own arrays (rather than tagged
    // entries in _springs) so each family stays identifiable for NetCDF I/O
    // (see getDihedralPhiSprings etc.) without touching Spring/_springs at
    // all. Which families are BUILT is a build-time decision (see
    // -dihedral/--dihedral in pdb2spn-cli.cpp: PHI/PSI/OMEGA are always
    // built together); which of the built ones are actually APPLIED at
    // runtime is independently gated in computeDihedralForces by this
    // Configuration's dihedral.phi/psi/omega/chi settings, on top of the
    // same isSpringEnabled() master switch regular springs use. Always
    // fully iterated when enabled (no static/dynamic split): a ghost
    // spring connecting two fully static particles is an unusual, not a
    // performance-critical, case.
    std::vector<Spring> _dihedralphisprings;
    std::vector<Spring> _dihedralpsisprings;
    std::vector<Spring> _dihedralomegasprings;
    std::vector<Spring> _dihedralsidechainsprings;
    std::vector<Spring> _dihedralplanaritysprings;

    // STRETCH springs: a real 1-2 bond, retuned to real AMBER r0/k -- a new
    // spring between the SAME two real atoms as the (now-zeroed, kept only
    // for nonbonded exclusion) --rigidbody spring, in its own collection
    // rather than retuned in place, so it has its own energy channel
    // (_energies.stretch) separate from the plain ENM/rigidbody baseline
    // (_energies.spring) -- see Energies' own comment.
    std::vector<Spring> _stretchsprings;

    // BEND ghost-ghost springs, standing in for a real valence-angle
    // restraint (see topology::Topology::_bend_springs' own comment for
    // why they connect two vertex-anchored ghosts instead of the two real
    // 1-3 atoms directly). Own energy channel (_energies.bend), same
    // reasoning as _stretchsprings above.
    std::vector<Spring> _bendsprings;

    // Ghost (massless virtual-site) particle bindings -- see
    // GhostParticle.h. Each entry's own Particle lives in _particles like
    // any other (isStatic()=true, mass=0); this only records which 3
    // anchor particles (by index) drive its position/force.
    std::vector<GhostParticleBinding> _ghostparticles;

    // One force contribution per dynamic spring. Reused between steps to avoid
    // allocations in the simulation loop and to keep OpenMP writes disjoint.
    std::vector<Vector3f> _springForceScratch;

    // One anchor-force triple per ghost, filled in parallel by
    // redistributeGhostForces' Jacobian pass, then accumulated serially
    // (anchors are heavily shared -- up to 62 ghosts per anchor on
    // example/072 -- so they must not be written concurrently). Same
    // parallel-compute/serial-accumulate split as _springForceScratch.
    struct GhostForceContribution
    {
        Vector3f F_B;
        Vector3f F_C;
        Vector3f F_Ref;
    };
    std::vector<GhostForceContribution> _ghostForceScratch;

    // One bucket per dynamic-particle-loop index, filled while computing
    // nonbonded pair interactions in parallel: each pair is evaluated once,
    // and the contribution owed to the *other* particle of the pair is
    // recorded here rather than written directly (that particle may be
    // processed concurrently by another thread). Reused between steps to
    // avoid allocations. See computeParticleForces / _applyNonbondedPairScratch.
    std::vector<std::vector<spn::DeferredNonbondedContribution>> _stericPairScratch;
    std::vector<std::vector<spn::DeferredNonbondedContribution>> _electrostaticPairScratch;
    std::vector<std::vector<spn::DeferredNonbondedContribution>> _hydrophobicPairScratch;

    // For Particle::addHydrogenBondCoreRepulsion's always-on, short-range-
    // only repulsive floor -- a neighbor-summed interaction like the
    // buffers above (unlike the exclusive attractive mechanism below),
    // since a particle can be pushed on by several nearby donors/acceptors
    // at once regardless of its own engagement status.
    std::vector<std::vector<spn::DeferredNonbondedContribution>> _hydrogenBondCoreRepulsionPairScratch;

    // Hydrogen bonds are exclusive pairs, not a neighbor-summed interaction:
    // each particle (donor or acceptor) may be engaged with at most one
    // partner at a time (see _assignHydrogenBondPairs). _hydrogenBondPartner
    // is indexed by particle index and persists across steps: -1 means free,
    // otherwise the index of the current partner. A bonded particle is not
    // even considered as a candidate for re-matching until its bond breaks
    // (current distance exceeds the hbond cutoff), so bonds do not flicker
    // toward a momentarily closer alternative. One force contribution per
    // active pair, reused between steps like _springForceScratch.
    std::vector<int> _hydrogenBondPartner;
    std::vector<Vector3f> _hydrogenBondForceScratch;

    Energies _energies;
    NeighborSearch _nsearch;
    bool _neighborSearchesDirty;

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
    bool _isTimeToLogData() const { return static_cast<unsigned>(_nbiter) % (static_cast<unsigned>(getSampleRate())) == 0; }

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
