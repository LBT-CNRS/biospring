#ifndef _SPRINGNETWORK_H_
#define _SPRINGNETWORK_H_

#include <array>
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
        // The --rigidbody/ENM mesh's own energy. The dihedral ghost rings
        // keep a dedicated channel instead of being folded in here, since
        // "spring energy" is only a meaningful, self-contained quantity for
        // a plain ENM/rigid-body network.
        float spring = 0.0f;
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
          _dynamicsprings(), _dihedralsprings(),
          _ghostparticles(),
          _springForceScratch(), _stericPairScratch(), _electrostaticPairScratch(),
          _hydrophobicPairScratch(), _hydrogenBondCoreRepulsionPairScratch(), _hbDonorSlot(),
          _hbDonorOffset(), _hbAcceptorSlot(), _hbAcceptorOffset(),
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
    // The dihedral families, used as the index of _dihedralsprings and of
    // anything that has to walk them in a fixed order (NetCDF I/O, the
    // per-family runtime gates, Topology's own collections and
    // BondedForceFieldReader's DihedralFamily, which both alias this).
    // Keep DIHEDRAL_FAMILY_COUNT last.
    //
    // The first five name protein chemistry (phi/psi/omega, chi1-4, and the
    // improper that keeps aromatic rings planar); the nucleic ones are
    // separate families rather than reusing SIDECHAIN because a nucleotide's
    // axes do not map onto a residue's: alpha..zeta are six backbone
    // torsions where a protein has three, and NUCLEIC_SUGAR is split off
    // from the rest of the backbone because those four furanose ring bonds
    // are what governs the sugar pucker -- the single lever selecting the A
    // or B helical form -- so its energy has to be readable, and switchable,
    // on its own.
    enum DihedralFamilyIndex
    {
        DIHEDRAL_PHI = 0,
        DIHEDRAL_PSI,
        DIHEDRAL_OMEGA,
        DIHEDRAL_SIDECHAIN,
        DIHEDRAL_PLANARITY,
        DIHEDRAL_NUCLEIC_BACKBONE,
        DIHEDRAL_NUCLEIC_CHI,
        DIHEDRAL_NUCLEIC_SUGAR,
        DIHEDRAL_FAMILY_COUNT
    };

    // The .nc variable-group prefix of each family, in DihedralFamilyIndex
    // order. Stated once here so the writer and the reader cannot drift out
    // of step -- they used to spell these names out separately, four times
    // over.
    static constexpr const char * DIHEDRAL_FAMILY_NAMES[DIHEDRAL_FAMILY_COUNT] = {
        "dihedralphi",      "dihedralpsi",   "dihedralomega",           "dihedralsidechain", "dihedralplanarity",
        "dihedralnucleicbackbone", "dihedralnucleicchi", "dihedralnucleicsugar"};

    // Read-only view of one family's springs, for NetCDF I/O.
    const std::vector<Spring> & getDihedralSprings(unsigned family) const
    {
        return _dihedralsprings[family];
    }

    float getDihedralEnergy() const { return _energies.dihedral; }
    float getStericEnergy() const { return _energies.steric; }
    float getElectrostaticEnergy() const { return _energies.electrostatic; }
    float getIMPEnergy() const { return _energies.imp; }
    float getHydrophobicEnergy() const { return _energies.hydrophobic; }
    float getHydrogenBondEnergy() const { return _energies.hbond; }

    // How many hydrogen bonds are held right now. Walks the donor slots, so
    // each bond counts once. Worth reporting alongside the energy: the total
    // alone cannot distinguish many weak bonds from few strong ones, and the
    // count is what shows a base pair holding its two or three.
    size_t getHydrogenBondCount() const;

    // How many RESIDUE PAIRS hold exactly 1, 2, 3, or 4-or-more bonds right
    // now, index 0 being the 1-bond bucket. This is what shows a base pair
    // holding its two or its three: the total energy cannot tell many weak
    // bonds from few strong ones, and the bond count cannot say how they are
    // distributed. Pairs of ADJACENT residues in one chain are excluded --
    // in a helix their bases are stacked within reach and would crowd the
    // census while carrying almost none of the energy.
    std::array<size_t, 4> getHydrogenBondPairCensus() const;

    // Writes every bond held right now, one per line, to `path`:
    //
    //   donor_chain donor_resid donor_resname donor_atom
    //   acceptor_...  distance_A  angular_weight  energy_kJmol
    //
    // Appended, with a "# step N" header per sample. Written only when
    // hbond.log names a file. This exists because a total and a count cannot
    // settle a disagreement about WHICH bonds are held, and nothing else in
    // the output can.
    void dumpHydrogenBonds(const std::string & path, int step) const;

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
    void addDihedralSpring(unsigned family, unsigned id1, unsigned id2, float equilibrium, float stiffness,
                           float dcOffset = 0.0f);

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
    unsigned addGhostParticle(unsigned placement, unsigned anchorBIndex, unsigned anchorCIndex,
                              unsigned anchorRefIndex, float r,
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

    // One dihedral ghost-spring family's list is reached by index only, via
    // getDihedralSprings(family) above.

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

    // -1 if particle `index` currently holds no hydrogen bond at all,
    // otherwise the index of one of its partners. Used by the core-repulsion
    // term, which only needs to know whether two particles are already bound
    // to each other -- see Particle::addHydrogenBondCoreRepulsion.
    int getHydrogenBondPartner(size_t index) const { return _anyHydrogenBondPartner(index); }

    // True when `a` and `b` currently hold a bond together, in either
    // direction.
    bool areHydrogenBonded(size_t a, size_t b) const;

    bool isSpringEnabled() const { return _config.spring.enable; }
    // Per-family runtime debug toggles for the dihedral ghost springs that
    // -dihedral* already decided, at build time, to create (see
    // Configuration.hpp's own comment on these settings).
    bool isDihedralPhiEnabled() const { return _config.dihedralphi.enable; }
    bool isDihedralPsiEnabled() const { return _config.dihedralpsi.enable; }
    bool isDihedralOmegaEnabled() const { return _config.dihedralomega.enable; }
    bool isDihedralChiEnabled() const { return _config.dihedralchi.enable; }
    bool isDihedralPlanarityEnabled() const { return _config.dihedralplanarity.enable; }
    bool isDihedralNucleicBackboneEnabled() const { return _config.dihedralnucleicbackbone.enable; }
    bool isDihedralNucleicChiEnabled() const { return _config.dihedralnucleicchi.enable; }
    bool isDihedralNucleicSugarEnabled() const { return _config.dihedralnucleicsugar.enable; }
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
    virtual void computeDihedralForces();
    virtual void computeHydrogenBondForces();
    virtual void computeParticleForces();
    virtual void updateParticlePositions();

  private:
    // Shared parallel-compute/serial-accumulate loop behind
    // computeDihedralForces -- see its definition for the
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

    // Any partner of `index`, in either role, or -1. Cheap: capacities are
    // 1 or 2 in practice.
    int _anyHydrogenBondPartner(size_t index) const;

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
    // Indexed by family rather than held as parallel members, for the same
    // reason topology::Topology is: a new family used to mean repeating the
    // same declaration, clear, accumulate, add and clear-all lines here too.
    std::array<std::vector<Spring>, DIHEDRAL_FAMILY_COUNT> _dihedralsprings;

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
    // example 072 -- so they must not be written concurrently). Same
    // parallel-compute/serial-accumulate split as _springForceScratch.
    struct GhostForceContribution
    {
        Vector3f F_B;
        Vector3f F_C;
        Vector3f F_Ref;
    };
    std::vector<GhostForceContribution> _ghostForceScratch;

    // One entry per distinct (B, C) ghost axis. The reaction owed to the
    // two axis atoms is reconstructed once per axis from the running
    // force/torque totals of every ghost hanging off it -- see
    // GhostParticle::redistributeAxisReaction, which is only valid over
    // complete spring pairs, hence per axis rather than per ghost.
    struct GhostAxis
    {
        unsigned anchorBIndex;
        unsigned anchorCIndex;
        Vector3f sumGhostForces;
        Vector3f sumGhostTorquesAboutB;
        Vector3f sumAtomForces;
        Vector3f sumAtomTorquesAboutB;
    };
    std::vector<GhostAxis> _ghostaxes;

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

    // Hydrogen bonds are counted pairs, not a neighbor-summed interaction.
    // Each particle has as many donor slots as it has donatable hydrogens
    // and as many acceptor slots as it has lone pairs (see
    // ParticleProperties::donorCapacity) -- an amino nitrogen donates twice,
    // a carbonyl oxygen accepts twice, a hydroxyl does one of each. A bond
    // consumes one DONOR slot on one side and one ACCEPTOR slot on the
    // other, so iterating the donor slots walks every bond exactly once and
    // knows which side donated, which is what the angular weight needs.
    //
    // Both are CSR: slots for particle i live in
    // [offset[i], offset[i+1]), -1 meaning free. They persist across steps
    // -- an engaged slot is not offered again until its own bond breaks
    // (distance beyond the hbond cutoff), so bonds do not flicker toward a
    // momentarily closer alternative -- and a particle with one slot behaves
    // exactly as the old single-partner array did.
    std::vector<int> _hbDonorSlot;
    std::vector<size_t> _hbDonorOffset;
    std::vector<int> _hbAcceptorSlot;
    std::vector<size_t> _hbAcceptorOffset;
    // One force contribution per active bond, reused between steps like
    // _springForceScratch. Three entries per bond: antecedent, donor,
    // acceptor -- the angular weight makes the antecedent a third body.
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
