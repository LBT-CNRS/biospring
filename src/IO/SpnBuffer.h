//
// Defines structures that hold a SpringNetwork particle and spring data
// to help write it using external APIs such as NetCDF.
//

#include "Particle.h"
#include "SpringNetwork.h"
#include "Vector3f.h"

static const size_t CHAIN_NAME_LENGTH = 4;
static const size_t RESIDUE_NAME_LENGTH = 4;
static const size_t PARTICLE_NAME_LENGTH = 5;

struct SpringBuffer
{
    size_t number_of_springs;
    int (*springs)[2];
    float * springsstiffnesses;
    float * springsequilibriums;
    int * nbofspringsperparticle;

    ~SpringBuffer() { clear(); }

    SpringBuffer(const SpringBuffer &) = delete;
    SpringBuffer & operator=(const SpringBuffer &) = delete;

    void clear()
    {
        delete[] springs;
        delete[] springsstiffnesses;
        delete[] springsequilibriums;
        delete[] nbofspringsperparticle;
        springs = nullptr;
        springsstiffnesses = nullptr;
        springsequilibriums = nullptr;
        nbofspringsperparticle = nullptr;
        number_of_springs = 0;
    }

    //
    // Empty constructor.
    //
    SpringBuffer()
        : number_of_springs(0), springs(0), springsstiffnesses(0), springsequilibriums(0), nbofspringsperparticle(0)
    {
    }

    //
    // Constructor with number of particles.
    // Allocates memory.
    //
    SpringBuffer(size_t number_of_springs, size_t number_of_particles) : SpringBuffer()
    {
        initialize(number_of_springs, number_of_particles);
    }

    // Allocates memory for all buffers.
    void initialize(size_t nSprings, size_t nParticles)
    {
        clear();
        number_of_springs = nSprings;

        if (nSprings > 0)
        {
            springs = new int[nSprings][2]{};
            springsstiffnesses = new float[nSprings]{};
            springsequilibriums = new float[nSprings]{};
        }
        if (nParticles > 0)
            nbofspringsperparticle = new int[nParticles]{};
    }

    // Copies SpringNetwork particle data to buffers.
    void bufferize(const biospring::spn::SpringNetwork * const spn)
    {
        for (size_t i = 0; i < spn->getNumberOfSprings(); ++i)
        {
            const biospring::spn::Spring & s = spn->getSpring(i);
            springs[i][0] = s.getParticle1().getId();
            springs[i][1] = s.getParticle2().getId();
            springsequilibriums[i] = s.getEquilibrium();
            springsstiffnesses[i] = s.getStiffness();
        }

        for (size_t i = 0; i < spn->getNumberOfParticles(); ++i)
        {
            // nbofspringsperparticle mirrors the NetCDF "nbspringsperparticle"
            // variable's on-disk type (nc_INT).
            nbofspringsperparticle[i] = static_cast<int>(spn->getParticle(i).getNumberOfSprings());
        }
    }
};

// Buffers one dihedral ghost-spring family (phi/psi/omega/sidechain/planarity
// -- see SpringNetwork's _dihedral*springs) for NetCDF I/O. Deliberately
// simpler than SpringBuffer: no `nbofspringsperparticle` companion array,
// since ghost springs are never registered as spring-neighbours (see
// SpringNetwork::addDihedralPhiSpring's comment) and nothing reads
// that field back on the SpringBuffer side either (write-only, informational).
struct DihedralSpringBuffer
{
    size_t number_of_springs;
    int (*springs)[2];
    float * springsstiffnesses;
    float * springsequilibriums;
    // This spring's share of its axis's exact dihedral-energy correction
    // (see spn::Spring::getDcOffset). Zero for a family that never got this
    // far in a version of the .bi.ff written before the correction existed
    // (old .nc files stay readable, see readDihedralSpringGroup).
    float * springsdcoffsets;

    ~DihedralSpringBuffer() { clear(); }

    DihedralSpringBuffer(const DihedralSpringBuffer &) = delete;
    DihedralSpringBuffer & operator=(const DihedralSpringBuffer &) = delete;

    void clear()
    {
        delete[] springs;
        delete[] springsstiffnesses;
        delete[] springsequilibriums;
        delete[] springsdcoffsets;
        springs = nullptr;
        springsstiffnesses = nullptr;
        springsequilibriums = nullptr;
        springsdcoffsets = nullptr;
        number_of_springs = 0;
    }

    DihedralSpringBuffer()
        : number_of_springs(0), springs(0), springsstiffnesses(0), springsequilibriums(0), springsdcoffsets(0)
    {
    }

    DihedralSpringBuffer(size_t number_of_springs) : DihedralSpringBuffer() { initialize(number_of_springs); }

    void initialize(size_t nSprings)
    {
        clear();
        number_of_springs = nSprings;
        if (nSprings > 0)
        {
            springs = new int[nSprings][2]{};
            springsstiffnesses = new float[nSprings]{};
            springsequilibriums = new float[nSprings]{};
            springsdcoffsets = new float[nSprings]{};
        }
    }

    // Copies one dihedral spring family's data (from SpringNetwork's own
    // vector, e.g. getDihedralPhiSprings()) into these buffers.
    void bufferize(const std::vector<biospring::spn::Spring> & source)
    {
        for (size_t i = 0; i < source.size(); ++i)
        {
            const biospring::spn::Spring & s = source[i];
            springs[i][0] = s.getParticle1().getId();
            springs[i][1] = s.getParticle2().getId();
            springsequilibriums[i] = s.getEquilibrium();
            springsstiffnesses[i] = s.getStiffness();
            springsdcoffsets[i] = s.getDcOffset();
        }
    }
};

// Buffers ghost-particle anchor bindings (SpringNetwork's _ghostparticles /
// GhostParticleBinding, see GhostParticle.h) for NetCDF I/O. A ghost
// particle's own coordinates/mass/etc. are already covered by the regular
// ParticleBuffer above (it is a real, if massless and static, entry in
// SpringNetwork's particle list) -- this buffer only carries the extra
// binding info (which particle is a ghost, its 3 anchor particles, and its
// placement parameters) needed to reconstruct it on reload.
struct GhostParticleBuffer
{
    size_t number_of_ghostparticles;
    int * ownindices;
    int (*anchorindices)[3];
    float * rs;
    float * thetas;
    float * deltas;
    // spn::GhostPlacement: which construction places the ghost. Explicit
    // rather than inferred from the geometry -- the two modes read disjoint
    // parameters, so no value combination identifies them on its own.
    int * placements;

    ~GhostParticleBuffer() { clear(); }

    GhostParticleBuffer(const GhostParticleBuffer &) = delete;
    GhostParticleBuffer & operator=(const GhostParticleBuffer &) = delete;

    void clear()
    {
        delete[] ownindices;
        delete[] anchorindices;
        delete[] rs;
        delete[] thetas;
        delete[] deltas;
        delete[] placements;
        placements = nullptr;
        ownindices = nullptr;
        anchorindices = nullptr;
        rs = nullptr;
        thetas = nullptr;
        deltas = nullptr;
        number_of_ghostparticles = 0;
    }

    GhostParticleBuffer()
        : number_of_ghostparticles(0), ownindices(0), anchorindices(0), rs(0), thetas(0), deltas(0), placements(0)
    {
    }

    GhostParticleBuffer(size_t number_of_ghostparticles) : GhostParticleBuffer()
    {
        initialize(number_of_ghostparticles);
    }

    void initialize(size_t nGhosts)
    {
        clear();
        number_of_ghostparticles = nGhosts;
        if (nGhosts > 0)
        {
            ownindices = new int[nGhosts]{};
            anchorindices = new int[nGhosts][3]{};
            rs = new float[nGhosts]{};
            thetas = new float[nGhosts]{};
            deltas = new float[nGhosts]{};
            placements = new int[nGhosts]{};
        }
    }

    // Copies SpringNetwork ghost-particle bindings (getGhostParticles()) into these buffers.
    void bufferize(const std::vector<biospring::spn::GhostParticleBinding> & source)
    {
        for (size_t i = 0; i < source.size(); ++i)
        {
            const biospring::spn::GhostParticleBinding & g = source[i];
            ownindices[i] = static_cast<int>(g.ownIndex);
            anchorindices[i][0] = static_cast<int>(g.anchorBIndex);
            anchorindices[i][1] = static_cast<int>(g.anchorCIndex);
            anchorindices[i][2] = static_cast<int>(g.anchorRefIndex);
            rs[i] = g.r;
            thetas[i] = g.theta_deg;
            deltas[i] = g.delta_deg;
            placements[i] = static_cast<int>(g.placement);
        }
    }
};

struct ParticleBuffer
{
    size_t number_of_particles;

    float (*coordinates)[3];

    float * charges;
    float * radii;
    float * epsilons;
    float * masses;
    // Per-particle hydrophobicity, used by the pairwise hydrophobic force term
    // (Particle::getHydrophobicity()). Distinct from `hscales`, which carries
    // the IMPALA transfer energy per unit accessible surface -- the .nc file
    // calls that one `hydrophobicityscale`, and conflating the two is what lost
    // this field once already.
    float * hydrophobicities;
    float * hscales;
    float * surface_accessibilities;

    int * ids;
    int * resids;

    unsigned char * dynamic_states;
    char (*chainnames)[CHAIN_NAME_LENGTH];
    char (*particlenames)[PARTICLE_NAME_LENGTH];
    char (*resnames)[RESIDUE_NAME_LENGTH];

    // Allocates memory for all buffers.
    void initialize(size_t nParticles)
    {
        clear();
        number_of_particles = nParticles;

        if (nParticles > 0)
        {
            coordinates = new float[nParticles][3]{};
            charges = new float[nParticles]{};
            radii = new float[nParticles]{};
            epsilons = new float[nParticles]{};
            masses = new float[nParticles]{};
            hydrophobicities = new float[nParticles]{};
            hscales = new float[nParticles]{};
            surface_accessibilities = new float[nParticles]{};
            ids = new int[nParticles]{};
            resids = new int[nParticles]{};
            dynamic_states = new unsigned char[nParticles]{};
            chainnames = new char[nParticles][CHAIN_NAME_LENGTH]{};
            particlenames = new char[nParticles][PARTICLE_NAME_LENGTH]{};
            resnames = new char[nParticles][RESIDUE_NAME_LENGTH]{};
        }
    }

    // Copies SpringNetwork particle data to buffers.
    void bufferize(const biospring::spn::SpringNetwork * const spn)
    {
        for (size_t i = 0; i < spn->getNumberOfParticles(); ++i)
        {
            const biospring::spn::Particle & p = spn->getParticle(i);

            Vector3f c = p.getPosition();
            coordinates[i][0] = c.getX();
            coordinates[i][1] = c.getY();
            coordinates[i][2] = c.getZ();

            charges[i] = p.getCharge();
            radii[i] = p.getRadius();
            epsilons[i] = p.getEpsilon();
            masses[i] = p.getMass();
            hydrophobicities[i] = p.getHydrophobicity();
            hscales[i] = p.getTransferEnergyByAccessibleSurface();
            surface_accessibilities[i] = p.getSolventAccessibilitySurface();

            ids[i] = p.getId();
            // resids mirrors the NetCDF "resids" variable's on-disk type (nc_INT).
            resids[i] = static_cast<int>(p.getResId());
            dynamic_states[i] = p.isDynamic();

            strncpy(chainnames[i], p.getChainName().c_str(), CHAIN_NAME_LENGTH);
            strncpy(particlenames[i], p.getName().c_str(), PARTICLE_NAME_LENGTH);
            strncpy(resnames[i], p.getResName().c_str(), RESIDUE_NAME_LENGTH);

        }
    }

    //
    // Empty constructor.
    //
    ParticleBuffer()
        : number_of_particles(0), coordinates(0), charges(0), radii(0), epsilons(0), masses(0), hydrophobicities(0),
          hscales(0), surface_accessibilities(0), ids(0), resids(0), dynamic_states(0), chainnames(0),
          particlenames(0), resnames(0)
    {
    }

    //
    // Constructor with number of particles.
    // Allocates memory.
    //
    ParticleBuffer(size_t number_of_particles) : ParticleBuffer() { initialize(number_of_particles); }

    ParticleBuffer(const ParticleBuffer &) = delete;
    ParticleBuffer & operator=(const ParticleBuffer &) = delete;

    virtual ~ParticleBuffer() { clear(); }

    void clear()
    {
        delete[] coordinates;
        delete[] charges;
        delete[] radii;
        delete[] epsilons;
        delete[] masses;
        delete[] hydrophobicities;
        delete[] hscales;
        delete[] surface_accessibilities;
        delete[] ids;
        delete[] resids;
        delete[] dynamic_states;
        delete[] chainnames;
        delete[] particlenames;
        delete[] resnames;

        coordinates = nullptr;
        charges = nullptr;
        radii = nullptr;
        epsilons = nullptr;
        masses = nullptr;
        hydrophobicities = nullptr;
        hscales = nullptr;
        surface_accessibilities = nullptr;
        ids = nullptr;
        resids = nullptr;
        dynamic_states = nullptr;
        chainnames = nullptr;
        particlenames = nullptr;
        resnames = nullptr;
        number_of_particles = 0;
    }
};
