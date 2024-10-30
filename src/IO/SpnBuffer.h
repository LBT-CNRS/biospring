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

    ~SpringBuffer()
    {
        delete[] springs;
        delete springsstiffnesses;
        delete springsequilibriums;
        delete nbofspringsperparticle;
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
    SpringBuffer(size_t number_of_springs, size_t number_of_particles)
    {
        this->initialize(number_of_springs, number_of_particles);
    }

    // Allocates memory for all buffers.
    void initialize(size_t nSprings, size_t nParticles)
    {
        number_of_springs = nSprings;

        if (nSprings > 0)
        {
            springs = new int[nSprings][2];
            springsstiffnesses = new float[nSprings];
            springsequilibriums = new float[nSprings];
            nbofspringsperparticle = new int[nParticles];
        }
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
            nbofspringsperparticle[i] = spn->getParticle(i).getNumberOfSprings();
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
        number_of_particles = nParticles;

        if (nParticles > 0)
        {
            coordinates = new float[nParticles][3];
            charges = new float[nParticles];
            radii = new float[nParticles];
            epsilons = new float[nParticles];
            masses = new float[nParticles];
            hscales = new float[nParticles];
            surface_accessibilities = new float[nParticles];
            ids = new int[nParticles];
            resids = new int[nParticles];
            dynamic_states = new unsigned char[nParticles];
            chainnames = new char[nParticles][CHAIN_NAME_LENGTH];
            particlenames = new char[nParticles][PARTICLE_NAME_LENGTH];
            resnames = new char[nParticles][RESIDUE_NAME_LENGTH];
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
            hscales[i] = p.getTransferEnergyByAccessibleSurface();

            ids[i] = p.getId();
            resids[i] = p.getResId();
            dynamic_states[i] = p.isDynamic();

            strncpy(chainnames[i], p.getChainName().c_str(), CHAIN_NAME_LENGTH);
            strncpy(particlenames[i], p.getName().c_str(), PARTICLE_NAME_LENGTH);
            strncpy(resnames[i], p.getResName().c_str(), RESIDUE_NAME_LENGTH);

            // !!! WARNING: Surface accessibilities not saved to output file !!!
        }
    }

    //
    // Empty constructor.
    //
    ParticleBuffer()
        : number_of_particles(0), coordinates(0), charges(0), radii(0), epsilons(0), masses(0), hscales(0),
          surface_accessibilities(0), ids(0), resids(0), dynamic_states(0), chainnames(0), particlenames(0), resnames(0)
    {
    }

    //
    // Constructor with number of particles.
    // Allocates memory.
    //
    ParticleBuffer(size_t number_of_particles) { this->initialize(number_of_particles); }

    virtual ~ParticleBuffer()
    {
        delete[] coordinates;
        delete charges;
        delete radii;
        delete epsilons;
        delete masses;
        delete hscales;
        delete surface_accessibilities;
        delete ids;
        delete resids;
        delete dynamic_states;
        delete[] chainnames;
        delete[] particlenames;
        delete[] resnames;
    }
};
