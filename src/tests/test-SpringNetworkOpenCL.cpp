#include <gtest/gtest.h>

#include <algorithm>
#include <vector>

#include "Particle.h"
#include "SpringNetwork.h"
#include "SpringNetworkOpenCL.h"
#include "interactor/Interactor.h"
#include "Vector3f.h"
#include "configuration/Configuration.hpp"

using namespace biospring;

namespace
{

// The smallest system that exercises everything the two backends have to agree
// on: one spring, two particles of DIFFERENT mass, released away from
// equilibrium, with a spring scale that is not 1.
//
// Every one of those details is load-bearing, because each caught a real defect
// when this comparison was first done by hand:
//
//   - a scale of 25 rather than 1: the kernel used to be passed only the unit
//     conversion, never the force field's spring scale, so it pulled 25 times
//     too weakly. With scale = 1 the two paths agree and the bug is invisible.
//   - masses of 12 and 14 rather than 1: the integration kernel used to do
//     velocities += forces*timestep with no mass at all, which is the CPU's
//     answer only when everything weighs 1 Da.
//   - released at 3.5 A from an equilibrium of 2.0 A: a system at equilibrium
//     does not move, and two backends that both do nothing agree perfectly.
//
// It also covers the two structural fixes: the CSR spring offsets (the old
// -1 sentinel made a particle whose successor had no spring contribute
// nothing), and the copy of the device's results back into the Particle
// objects, without which the GPU computed correctly and reported the structure
// it started from.
void buildReleasedSpring(spn::SpringNetwork & network, configuration::Configuration & config)
{
    spn::Particle first;
    first.setPosition(Vector3f(0.0f, 0.0f, 0.0f));
    first.setMass(12.0f);
    network.addParticle(first);

    spn::Particle second;
    second.setPosition(Vector3f(3.5f, 0.0f, 0.0f));
    second.setMass(14.0f);
    network.addParticle(second);

    network.addSpring(0, 1, /*equilibrium=*/2.0f, /*stiffness=*/1.0f);

    config = configuration::defaultConfiguration();
    config.sim.nbsteps = 1000;
    config.sim.timestep = 0.1;
    config.spring.enable = true;
    config.spring.scale = 25.0;
    config.viscosity.enable = true;
    config.viscosity.value = 1.0;

    network.setup(config);
}

float separation(const spn::SpringNetwork & network)
{
    return (network.getParticle(1).getPosition() - network.getParticle(0).getPosition()).norm();
}

// Constructing SpringNetworkOpenCL calls exit() when no platform answers, so
// the availability check has to happen before that, not inside it.
bool hasOpenCLDevice()
{
    std::vector<cl::Platform> platforms;
    try
    {
        cl::Platform::get(&platforms);
    }
    catch (...)
    {
        return false;
    }
    for (auto & platform : platforms)
    {
        std::vector<cl::Device> devices;
        try
        {
            platform.getDevices(CL_DEVICE_TYPE_ALL, &devices);
        }
        catch (...)
        {
            continue;
        }
        if (!devices.empty())
            return true;
    }
    return false;
}

} // namespace

// The regression this exists to prevent is not a crash: it is the two backends
// quietly computing different physics. They are two implementations of one
// model -- the force law is shared source (forcefield/shared/spring_shared.h,
// compiled as C++ here and as OpenCL C into the kernel), but the memory layout,
// the traversal and the integration are written twice. Nothing else in the
// build fails when they drift.
TEST(SpringNetworkOpenCL, MatchesTheCPUOnAReleasedSpring)
{
    if (!hasOpenCLDevice())
        GTEST_SKIP() << "no OpenCL device available on this machine";

    spn::SpringNetwork cpu;
    configuration::Configuration cpuconfig;
    buildReleasedSpring(cpu, cpuconfig);
    cpu.run();

    SpringNetworkOpenCL gpu;
    configuration::Configuration gpuconfig;
    buildReleasedSpring(gpu, gpuconfig);
    gpu.run();

    const float cpuseparation = separation(cpu);
    const float gpuseparation = separation(gpu);

    // Guard against the test passing because nothing happened on either side:
    // two backends that both do nothing agree perfectly. Over 1000 steps the
    // spring pulls the pair from 3.5 A well past 3.1 A, so anything above that
    // means the run did not take. The long integration is also what makes the
    // agreement below worth asserting -- it gives any real difference a
    // thousand steps to show itself.
    EXPECT_LT(cpuseparation, 3.1f) << "the CPU run did not move the particles";
    EXPECT_LT(gpuseparation, 3.1f) << "the GPU run did not move the particles";

    // Both paths are float, and they sum in different orders, so they are not
    // expected to be bit-identical -- only to be the same physics. Measured
    // agreement on this case is within a thousandth of an Angstrom; the bar is
    // set a little wider, and still far below the 25x, mass-free and
    // never-copied-back failures it is there to catch.
    EXPECT_NEAR(cpuseparation, gpuseparation, 5.0e-3f)
        << "CPU and OpenCL disagree on the separation after " << cpuconfig.sim.nbsteps << " steps";

    // SpringNetworkOpenCL makes its getNumberOfParticles() override private, so
    // the count comes from the CPU side; both networks were built by the same
    // function, so they hold the same two particles.
    for (unsigned i = 0; i < cpu.getNumberOfParticles(); ++i)
    {
        const Vector3f delta = cpu.getParticle(i).getPosition() - gpu.getParticle(i).getPosition();
        EXPECT_LT(delta.norm(), 5.0e-3f) << "particle " << i << " ended up elsewhere on the GPU";
    }
}

namespace
{

// Same released spring, plus a third particle that is STATIC and carries a real
// mass, held to the first by a spring of its own.
//
// That combination is the one the mass-based filter got wrong. The CPU freezes
// a particle because it is absent from its dynamic list; its mass never enters
// into it. The kernel had no notion of the dynamic state and filtered on the
// mass alone, which agrees with the CPU only while every static particle is a
// massless ghost -- true of every example shipped today, and untrue the moment
// anyone builds a network with pdb2spn --static, which freezes particles
// without touching their masses.
void buildWithAnchoredStaticParticle(spn::SpringNetwork & network, configuration::Configuration & config)
{
    spn::Particle first;
    first.setPosition(Vector3f(0.0f, 0.0f, 0.0f));
    first.setMass(12.0f);
    network.addParticle(first);

    spn::Particle second;
    second.setPosition(Vector3f(3.5f, 0.0f, 0.0f));
    second.setMass(14.0f);
    network.addParticle(second);

    spn::Particle anchor;
    anchor.setPosition(Vector3f(0.0f, 3.5f, 0.0f));
    anchor.setMass(12.0f);
    anchor.setDynamic(false);
    network.addParticle(anchor);

    network.addSpring(0, 1, /*equilibrium=*/2.0f, /*stiffness=*/1.0f);
    network.addSpring(0, 2, /*equilibrium=*/2.0f, /*stiffness=*/1.0f);

    config = configuration::defaultConfiguration();
    config.sim.nbsteps = 1000;
    config.sim.timestep = 0.1;
    config.spring.enable = true;
    config.spring.scale = 25.0;
    config.viscosity.enable = true;
    config.viscosity.value = 1.0;

    network.setup(config);
}

} // namespace

TEST(SpringNetworkOpenCL, HoldsStaticParticlesTheWayTheCPUDoes)
{
    if (!hasOpenCLDevice())
        GTEST_SKIP() << "no OpenCL device available on this machine";

    const Vector3f anchorstart(0.0f, 3.5f, 0.0f);

    spn::SpringNetwork cpu;
    configuration::Configuration cpuconfig;
    buildWithAnchoredStaticParticle(cpu, cpuconfig);
    cpu.run();

    SpringNetworkOpenCL gpu;
    configuration::Configuration gpuconfig;
    buildWithAnchoredStaticParticle(gpu, gpuconfig);
    gpu.run();

    // The spring pulls hard on the anchor -- it starts 1.5 A from equilibrium --
    // so it staying put is a statement about the dynamic state being honoured,
    // not about there being nothing to move it.
    EXPECT_LT((cpu.getParticle(2).getPosition() - anchorstart).norm(), 1.0e-5f)
        << "the CPU moved a static particle";
    EXPECT_LT((gpu.getParticle(2).getPosition() - anchorstart).norm(), 1.0e-5f)
        << "the GPU moved a static particle; it is filtering on the mass instead of the dynamic state";

    // And the two dynamic particles, pulled by a spring whose other end is
    // pinned, must still end up in the same place on both paths.
    for (unsigned i = 0; i < 2; ++i)
    {
        const Vector3f delta = cpu.getParticle(i).getPosition() - gpu.getParticle(i).getPosition();
        EXPECT_LT(delta.norm(), 5.0e-3f) << "particle " << i << " ended up elsewhere on the GPU";
    }
}

// ---------------------------------------------------------------------------
// The cell list
// ---------------------------------------------------------------------------
//
// A neighbour structure that quietly misses pairs does not crash and does not
// look wrong. It makes every term built on it too weak, by an amount nothing
// reports -- so the only honest test is against the O(N^2) answer, which is
// what "near" means.
//
// What this exercises that a regular lattice would not: a cloud whose extent is
// not a whole number of cells, particles in the corner cells (whose stencil
// falls outside the grid and must be clipped rather than wrapped, since nothing
// here is periodic), and cells holding several particles at once -- which is
// where the atomic exchange in binParticles earns its keep.
namespace
{
void buildParticleCloud(spn::SpringNetwork & network, configuration::Configuration & config,
                        unsigned n, float extent)
{
    // A fixed sequence rather than a random one: a test that fails only on
    // some seeds is a test nobody can act on.
    unsigned state = 12345u;
    const auto next = [&state]() {
        state = state * 1103515245u + 12345u;
        return static_cast<float>((state >> 16) & 0x7fffu) / static_cast<float>(0x7fff);
    };

    for (unsigned i = 0; i < n; ++i)
    {
        spn::Particle p;
        // Deliberately off-origin and partly negative: cells are indexed from a
        // measured origin, and a version that assumed positive coordinates
        // would pass on a structure that happens to sit in the first octant.
        p.setPosition(Vector3f(next() * extent - 0.35f * extent,
                               next() * extent - 0.35f * extent,
                               next() * extent - 0.35f * extent));
        p.setMass(12.0f);
        p.setRadius(1.9f);
        // Half carry a charge and a third a hydrophobicity, so the Coulomb and
        // hydrophobic grids hold genuine SUBSETS rather than everyone. With
        // every particle in all three, a grid that ignored its mask would pass.
        if (i % 2 == 0)
            p.setCharge(0.4f);
        if (i % 3 == 0)
            p.setHydrophobicity(1.0f);
        network.addParticle(p);
    }

    config = configuration::defaultConfiguration();
    config.sim.nbsteps = 1;
    config.sim.timestep = 0.1;
    config.spring.enable = false;

    // Three terms, three DIFFERENT cutoffs, which is what makes the per-term
    // grids worth testing: three sets of bins, three cell widths, three
    // populations.
    config.steric.enable = true;
    config.steric.cutoff = 6.0;
    config.electrostatic.enable = true;
    config.electrostatic.cutoff = 14.0;
    config.hydrophobicity.enable = true;
    config.hydrophobicity.cutoff = 10.0;

    network.setup(config);
}

// The O(N^2) answer, which is what "near" means.
std::vector<unsigned> neighborsByBruteForce(const spn::SpringNetwork & network, unsigned i, float cutoff)
{
    std::vector<unsigned> expected;
    for (unsigned j = 0; j < network.getNumberOfParticles(); ++j)
    {
        if (j == i)
            continue;
        if ((network.getParticle(i).getPosition() - network.getParticle(j).getPosition()).norm() <= cutoff)
            expected.push_back(j);
    }
    return expected;
}
} // namespace

TEST(SpringNetworkOpenCL, EachTermGetsItsOwnGridAtItsOwnWidth)
{
    if (!hasOpenCLDevice())
        GTEST_SKIP() << "no OpenCL device available on this machine";

    const unsigned N = 900;
    const float EXTENT = 40.0f;

    SpringNetworkOpenCL gpu;
    configuration::Configuration config;
    buildParticleCloud(gpu, config, N, EXTENT);
    gpu.run();

    struct Term
    {
        const char * name;
        const SpringNetworkOpenCL::CellGrid * grid;
        float cutoff;
        bool (*belongs)(const spn::Particle &);
    };
    const Term terms[] = {
        {"steric", &gpu.cells(), 6.0f,
         [](const spn::Particle &) { return true; }},
        {"electrostatic", &gpu.chargedCells(), 14.0f,
         [](const spn::Particle & p) { return p.isCharged(); }},
        {"hydrophobic", &gpu.hydrophobicCells(), 10.0f,
         [](const spn::Particle & p) { return p.isHydrophobic(); }},
    };

    // Three widths, all different, because each is its term's own search
    // radius. Under the single shared grid these were equal by construction,
    // so this is the assertion that tells the two designs apart.
    EXPECT_NE(gpu.cells().width, gpu.chargedCells().width);
    EXPECT_NE(gpu.cells().width, gpu.hydrophobicCells().width);
    EXPECT_NE(gpu.chargedCells().width, gpu.hydrophobicCells().width);

    for (const Term & term : terms)
    {
        const SpringNetworkOpenCL::CellGrid & grid = *term.grid;
        ASSERT_GT(grid.ncellstotal, 0u) << term.name << " built no cell list at all";

        // A 40 A cloud is nowhere near the cell limit, so the cells built are
        // the cells asked for. Widening is the divergence path, not this one.
        EXPECT_FLOAT_EQ(grid.width, grid.requestedwidth) << term.name;

        // Cells the size of the search radius, so the walk is the 27 cells
        // around the particle's own and no more.
        EXPECT_EQ(biospring_stencil_radius(term.cutoff, grid.width), 1)
            << term.name << " is not walking its own cutoff in one step of cells";

        size_t pairs = 0;
        for (unsigned i = 0; i < N; ++i)
        {
            if (!term.belongs(gpu.getParticle(i)))
                continue;

            std::vector<unsigned> found = gpu.neighborsFromCellList(grid, i, term.cutoff);

            // Brute force over THIS term's population, not over everyone: a
            // grid that ignored its mask would find the uncharged particles
            // too, and this is what catches it.
            std::vector<unsigned> expected;
            for (unsigned j = 0; j < N; ++j)
            {
                if (j == i || !term.belongs(gpu.getParticle(j)))
                    continue;
                if ((gpu.getParticle(i).getPosition() - gpu.getParticle(j).getPosition()).norm()
                    <= term.cutoff)
                    expected.push_back(j);
            }

            std::sort(found.begin(), found.end());
            std::sort(expected.begin(), expected.end());

            // Duplicates matter as much as omissions: a particle counted twice
            // is a force applied twice, and the linked list makes that possible
            // if a stencil ever visits the same cell more than once.
            ASSERT_EQ(found, expected)
                << term.name << "'s walk over its own grid disagrees with the O(N^2) answer"
                << " for particle " << i;
            pairs += expected.size();
        }

        // Guards against passing because nothing is near anything.
        EXPECT_GT(pairs, 500u) << "the cloud is too sparse for " << term.name << " to mean anything";
    }
}

// What happens when the structure no longer fits the grid.
//
// Both backends face it the same way now -- an array over a measured box, which
// a structure can outgrow -- so both answer it the same way.
//
// Cells WIDER than asked for are still exact: every walk derives its stencil
// radius from the width it is handed, so it takes fewer and bigger steps and the
// distance test drops the surplus. The answer to a box too big is therefore
// wider cells, not no grid. This drives it with a cutoff of 1 A and one particle
// 5000 A away, which asks for more cells than the build will allocate and gets
// them widened instead. The neighbour sets must still be exactly the O(N^2) ones:
// degrading towards brute force has to stay correct, or it is not a degradation
// but a bug with a warning attached.
TEST(SpringNetworkOpenCL, TooBigABoxWidensTheCellsAndStaysExact)
{
    if (!hasOpenCLDevice())
        GTEST_SKIP() << "no OpenCL device available on this machine";

    const unsigned N = 200;
    const float CUTOFF = 1.0f;

    SpringNetworkOpenCL gpu;
    configuration::Configuration config;

    unsigned state = 7u;
    const auto next = [&state]() {
        state = state * 1103515245u + 12345u;
        return static_cast<float>((state >> 16) & 0x7fffu) / static_cast<float>(0x7fff);
    };
    for (unsigned i = 0; i < N; ++i)
    {
        spn::Particle p;
        p.setPosition(Vector3f(next() * 6.0f, next() * 6.0f, next() * 6.0f));
        p.setMass(12.0f);
        gpu.addParticle(p);
    }
    // The one that blows the box out.
    spn::Particle far;
    far.setPosition(Vector3f(5000.0f, 5000.0f, 5000.0f));
    far.setMass(12.0f);
    gpu.addParticle(far);

    config = configuration::defaultConfiguration();
    config.sim.nbsteps = 1;
    config.sim.timestep = 0.1;
    config.spring.enable = false;
    config.steric.enable = true;
    config.steric.cutoff = CUTOFF;
    config.electrostatic.enable = false;
    config.hydrophobicity.enable = false;
    gpu.setup(config);
    gpu.run();

    const SpringNetworkOpenCL::CellGrid & grid = gpu.cells();
    ASSERT_GT(grid.ncellstotal, 0u) << "the grid was refused outright instead of widening";
    EXPECT_GT(grid.width, grid.requestedwidth) << "the cells were not widened, so the box cannot have fitted";

    // And it is still the right answer.
    for (unsigned i = 0; i < N + 1; ++i)
    {
        std::vector<unsigned> found = gpu.neighborsFromCellList(grid, i, CUTOFF);
        std::vector<unsigned> expected = neighborsByBruteForce(gpu, i, CUTOFF);
        std::sort(found.begin(), found.end());
        std::sort(expected.begin(), expected.end());
        ASSERT_EQ(found, expected) << "widened cells gave the wrong neighbours for particle " << i;
    }
}

// Coulomb, on the device, against the same run on the CPU.
//
// The two reach it by different routes -- the CPU walks its unordered_map grid
// and evaluates each pair once, handing the other half over; the device walks
// its cell list and evaluates each pair from both ends -- so agreeing is a
// statement about the law and about the neighbour sets, not about shared code
// paths. The force module itself IS shared source
// (forcefield/shared/electrostatic_shared.h), which is what makes the residual
// difference float32 rounding rather than a second transcription.
//
// ALL CHARGES THE SAME SIGN, which is not a detail. The first version of this
// gave them both signs and no steric term, and two opposite charges fell into
// each other: 0.274 A apart after 200 steps, where the 1/r^2 force is some 1900
// times its value at the cutoff. From there the last bit of a float decides
// which way the pair flies, and the two backends ended up 15 A apart on a law
// they both computed correctly. Pure repulsion has no such singularity, so what
// is left to disagree about is the arithmetic. The guard below keeps it that
// way -- a test that quietly starts measuring a collapse again would go on
// passing for a while and then fail for the wrong reason.
namespace
{
void buildChargedCloud(spn::SpringNetwork & network, configuration::Configuration & config, unsigned n,
                       double skin = 0.0, float charge = 0.4f);

void buildChargedCloud(spn::SpringNetwork & network, configuration::Configuration & config, unsigned n, double skin,
                       float charge)
{
    unsigned state = 99u;
    const auto next = [&state]() {
        state = state * 1103515245u + 12345u;
        return static_cast<float>((state >> 16) & 0x7fffu) / static_cast<float>(0x7fff);
    };

    for (unsigned i = 0; i < n; ++i)
    {
        spn::Particle p;
        p.setPosition(Vector3f(next() * 25.0f - 9.0f, next() * 25.0f - 9.0f, next() * 25.0f - 9.0f));
        p.setMass(12.0f);
        p.setCharge(charge);
        network.addParticle(p);
    }
    // One sprung pair, held at 1.4 A. Coulomb would throw two like charges that
    // close apart at once, so the pair staying there is what says the exclusion
    // ran -- on both backends, by the same measurement.
    network.addSpring(0, 1, /*equilibrium=*/1.4f, /*stiffness=*/500.0f);

    config = configuration::defaultConfiguration();
    config.sim.nbsteps = 200;
    config.sim.timestep = 0.5;
    config.spring.enable = true;
    config.spring.scale = 1.0;
    config.viscosity.enable = true;
    config.viscosity.value = 1.0;
    config.steric.enable = false;
    config.hydrophobicity.enable = false;
    config.electrostatic.enable = true;
    config.electrostatic.scale = 1.0;
    config.electrostatic.cutoff = 12.0;
    config.sim.neighborskin = skin;

    network.setup(config);
}

// =====================================================================================
//
// The asymmetry of a pair: who is seen, and who looks.
//
// =====================================================================================

// A network of four particles, 3 A apart along x and 2 A in radius -- so every
// pair overlaps and the steric term actually pushes -- one of each kind:
//   0  dynamic, charged      -- looks, and is seen
//   1  STATIC,  charged      -- is seen by 0, does not look
//   2  dynamic, uncharged    -- looks and is seen for steric, neither for Coulomb
//   3  STATIC,  uncharged    -- seen for steric only, never looks
void buildFourKinds(spn::SpringNetwork & network, configuration::Configuration & config,
                    bool steric, bool coulomb, int nbsteps)
{
    const struct { float x; bool dynamic; float charge; } kinds[4] = {
        {0.0f, true, 0.5f}, {3.0f, false, 0.5f}, {6.0f, true, 0.0f}, {9.0f, false, 0.0f}};

    for (const auto & k : kinds)
    {
        spn::Particle p;
        p.setPosition(Vector3f(k.x, 0.0f, 0.0f));
        p.setMass(12.0f);
        p.setCharge(k.charge);
        p.setRadius(2.0f);
        p.setEpsilon(0.5f);
        p.setDynamic(k.dynamic);
        network.addParticle(p);
    }

    config = configuration::defaultConfiguration();
    // The step count is the caller's: looking at the lists wants one step, so
    // that they still describe the geometry they were built for, and comparing
    // the forces wants many, because one step moves the dynamic particles less
    // than the tolerance and a backend that never saw the static ones would
    // still agree.
    config.sim.nbsteps = nbsteps;
    config.sim.timestep = 2.0;
    config.spring.enable = false;
    config.steric.enable = steric;
    config.steric.cutoff = 20.0;
    config.electrostatic.enable = coulomb;
    config.electrostatic.cutoff = 20.0;
    config.hydrophobicity.enable = false;
    network.setup(config);
}
} // namespace

// The rule, stated once and checked on the device's own structure: a STATIC
// particle applies a force to its dynamic neighbours and must therefore be
// visible to them, but has no neighbours of its own, because no force will ever
// be computed on it. That asymmetry is the whole point of having two masks
// rather than one, and it holds for the steric term exactly as for Coulomb --
// the ONLY difference between them is that Coulomb additionally drops the
// particles that carry no charge, at both ends.
TEST(SpringNetworkOpenCL, StaticParticlesAreSeenButDoNotLook)
{
    if (!hasOpenCLDevice())
        GTEST_SKIP() << "no OpenCL device available on this machine";

    SpringNetworkOpenCL gpu;
    configuration::Configuration config;
    buildFourKinds(gpu, config, /*steric=*/true, /*coulomb=*/true, /*nbsteps=*/1);
    gpu.run();

    const SpringNetworkOpenCL::NeighbourList & st = gpu.stericList();
    const SpringNetworkOpenCL::NeighbourList & el = gpu.electrostaticList();
    ASSERT_TRUE(st.valid) << "no steric list was built";
    ASSERT_TRUE(el.valid) << "no electrostatic list was built";

    const auto list = [&](const SpringNetworkOpenCL::NeighbourList & l, unsigned i) {
        std::vector<unsigned> v = gpu.neighboursFromList(l, i);
        std::sort(v.begin(), v.end());
        return v;
    };

    // STERIC. Everything is a candidate, so a dynamic particle sees all three
    // others -- including both static ones, which is how they push it.
    EXPECT_EQ(list(st, 0), (std::vector<unsigned>{1, 2, 3})) << "the dynamic charged particle cannot see everyone";
    EXPECT_EQ(list(st, 2), (std::vector<unsigned>{0, 1, 3})) << "the dynamic uncharged particle cannot see everyone";
    // And neither static one looks back.
    EXPECT_TRUE(list(st, 1).empty()) << "a static particle was given a steric list it will never read";
    EXPECT_TRUE(list(st, 3).empty()) << "a static particle was given a steric list it will never read";

    // COULOMB. Same asymmetry, plus the charge filter at BOTH ends: 2 and 3
    // carry no charge, so they neither look nor are looked at.
    EXPECT_EQ(list(el, 0), (std::vector<unsigned>{1})) << "Coulomb should see only the other charged particle";
    EXPECT_TRUE(list(el, 1).empty()) << "a static particle was given a Coulomb list";
    EXPECT_TRUE(list(el, 2).empty()) << "an uncharged particle was given a Coulomb list";
    EXPECT_TRUE(list(el, 3).empty()) << "an uncharged static particle was given a Coulomb list";
}

// And the forces that come out of that have to be the CPU's, which reaches the
// same place by another route: it never asks a static particle for neighbours,
// because computeParticleForces only walks the dynamic ones.
//
// Run term by term, so that a mistake in one cannot be hidden by the other.
TEST(SpringNetworkOpenCL, TheAsymmetryGivesTheCPUsForces)
{
    if (!hasOpenCLDevice())
        GTEST_SKIP() << "no OpenCL device available on this machine";

    struct Case { const char * name; bool steric; bool coulomb; };
    for (const Case & c : {Case{"steric alone", true, false},
                           Case{"coulomb alone", false, true},
                           Case{"both", true, true}})
    {
        spn::SpringNetwork cpu;
        SpringNetworkOpenCL gpu;
        configuration::Configuration c1, c2;
        buildFourKinds(cpu, c1, c.steric, c.coulomb, /*nbsteps=*/200);
        buildFourKinds(gpu, c2, c.steric, c.coulomb, /*nbsteps=*/200);
        cpu.run();
        gpu.run();

        for (unsigned i = 0; i < 4; ++i)
        {
            const Vector3f a = cpu.getParticle(i).getPosition();
            const Vector3f b = gpu.getParticle(i).getPosition();
            EXPECT_LT((a - b).norm(), 1.0e-4f)
                << c.name << ": particle " << i << " ended " << (a - b).norm() << " A apart";
        }

        // The static ones must not have moved at all, on either backend.
        for (unsigned i : {1u, 3u})
        {
            EXPECT_FLOAT_EQ(cpu.getParticle(i).getPosition().getX(), i == 1 ? 3.0f : 9.0f)
                << c.name << ": the CPU moved static particle " << i;
            EXPECT_FLOAT_EQ(gpu.getParticle(i).getPosition().getX(), i == 1 ? 3.0f : 9.0f)
                << c.name << ": the device moved static particle " << i;
        }

        // And the dynamic ones must have moved, or none of the above means
        // anything: a run where nothing happens agrees perfectly.
        // The dynamic ones have to have moved FAR ENOUGH that a dropped pair
        // could not hide under the tolerance above.
        const float moved = std::max((cpu.getParticle(0).getPosition() - Vector3f(0, 0, 0)).norm(),
                                     (cpu.getParticle(2).getPosition() - Vector3f(6, 0, 0)).norm());
        EXPECT_GT(moved, 0.1f) << c.name << ": moved only " << moved << " A, which proves nothing";
    }
}

namespace
{
float closestPair(const spn::SpringNetwork & network)
{
    float closest = 1.0e9f;
    for (unsigned i = 0; i < network.getNumberOfParticles(); ++i)
        for (unsigned j = i + 1; j < network.getNumberOfParticles(); ++j)
            closest = std::min(closest,
                (network.getParticle(i).getPosition() - network.getParticle(j).getPosition()).norm());
    return closest;
}
} // namespace

TEST(SpringNetworkOpenCL, CoulombMatchesTheCPU)
{
    if (!hasOpenCLDevice())
        GTEST_SKIP() << "no OpenCL device available on this machine";

    const unsigned N = 300;

    spn::SpringNetwork cpu;
    configuration::Configuration cpuconfig;
    buildChargedCloud(cpu, cpuconfig, N);
    cpu.run();

    SpringNetworkOpenCL gpu;
    configuration::Configuration gpuconfig;
    buildChargedCloud(gpu, gpuconfig, N);
    gpu.run();

    // Two backends that both did nothing would agree perfectly and mean
    // nothing, so check the cloud actually blew apart.
    spn::SpringNetwork start;
    configuration::Configuration startconfig;
    buildChargedCloud(start, startconfig, N);
    float moved = 0.0f;
    for (unsigned i = 0; i < N; ++i)
        moved = std::max(moved, (cpu.getParticle(i).getPosition() - start.getParticle(i).getPosition()).norm());
    ASSERT_GT(moved, 1.0f) << "the cloud barely moved, so this comparison proves nothing";

    // The sprung pair is the exclusion's witness: at 1.4 A, two 0.4 e charges
    // repel far harder than a 500 kJ.mol-1.A-2 spring holds, so this pair
    // stays put only if Coulomb never saw it.
    for (const spn::SpringNetwork * net : {static_cast<const spn::SpringNetwork *>(&cpu),
                                           static_cast<const spn::SpringNetwork *>(&gpu)})
    {
        const float held = (net->getParticle(0).getPosition() - net->getParticle(1).getPosition()).norm();
        EXPECT_NEAR(held, 1.4f, 0.2f) << "the sprung pair was blown apart: Coulomb is not excluding it";
    }

    // And no OTHER pair may have collapsed, or the comparison below is
    // measuring a singularity rather than a force law.
    EXPECT_GT(closestPair(cpu), 1.0f) << "a pair collapsed on the CPU; this test is no longer about the law";
    EXPECT_GT(closestPair(gpu), 1.0f) << "a pair collapsed on the GPU; this test is no longer about the law";

    float worst = 0.0f;
    for (unsigned i = 0; i < N; ++i)
        worst = std::max(worst, (cpu.getParticle(i).getPosition() - gpu.getParticle(i).getPosition()).norm());

    // Loose enough for float32 against float64 and a different summation order
    // over 200 steps, tight enough that a missing exclusion or a sign error
    // cannot hide.
    EXPECT_LT(worst, 0.05f) << "the GPU's Coulomb ended up " << worst << " A from the CPU's";
}


// A skin must not change the answer, on either backend.
//
// That is the whole claim of a stored neighbour list: it is an optimisation and
// not an approximation, because what it holds beyond the cutoff is exactly what
// lets it stay right until the next rebuild. A list that were rebuilt too late
// would drop pairs, and on a cloud of like charges a dropped pair is a particle
// that quietly stops being pushed.
//
// Held against the SAME configuration without a skin, which is the reference
// this whole port has to reproduce -- separately for the CPU, whose list lives
// in nsearch.hpp, and for the device, whose list is two kernels and a packing
// pass, so that neither can be right by accident while the other is wrong.
TEST(SpringNetworkOpenCL, ASkinChangesNothingOnEitherBackend)
{
    if (!hasOpenCLDevice())
        GTEST_SKIP() << "no OpenCL device available on this machine";

    const unsigned N = 300;
    const double SKIN = 4.0;

    spn::SpringNetwork cpubare, cpuskinned;
    SpringNetworkOpenCL gpubare, gpuskinned;
    configuration::Configuration c1, c2, c3, c4;
    // A tenth of the usual charge, so the cloud creeps rather than explodes and
    // the list survives tens of steps: a list rebuilt every step would be right
    // whatever radius it was built at, and would prove nothing.
    const float CREEP = 0.04f;
    buildChargedCloud(cpubare, c1, N, 0.0, CREEP);
    buildChargedCloud(cpuskinned, c2, N, SKIN, CREEP);
    buildChargedCloud(gpubare, c3, N, 0.0, CREEP);
    buildChargedCloud(gpuskinned, c4, N, SKIN, CREEP);
    cpubare.run();
    cpuskinned.run();
    gpubare.run();
    gpuskinned.run();

    // The cloud has to have actually moved, or agreeing proves nothing.
    spn::SpringNetwork start;
    configuration::Configuration c0;
    buildChargedCloud(start, c0, N, 0.0, CREEP);
    float moved = 0.0f;
    for (unsigned i = 0; i < N; ++i)
        moved = std::max(moved, (cpubare.getParticle(i).getPosition() - start.getParticle(i).getPosition()).norm());
    ASSERT_GT(moved, 1.0f) << "the cloud barely moved, so this comparison proves nothing";

    // And the list has to have been REUSED, or the comparison is vacuous: a
    // list rebuilt at every step holds this step's neighbours whatever radius
    // it was built at, so it could be built at the bare cutoff and still agree.
    ASSERT_GT(gpuskinned.neighbourListRebuilds(), 0u) << "no list was ever built";
    ASSERT_LT(gpuskinned.neighbourListRebuilds(), 100u)
        << "the device rebuilt its list " << gpuskinned.neighbourListRebuilds()
        << " times in 200 steps, so nothing was ever reused";

    float worstcpu = 0.0f, worstgpu = 0.0f;
    for (unsigned i = 0; i < N; ++i)
    {
        worstcpu = std::max(worstcpu,
            (cpubare.getParticle(i).getPosition() - cpuskinned.getParticle(i).getPosition()).norm());
        worstgpu = std::max(worstgpu,
            (gpubare.getParticle(i).getPosition() - gpuskinned.getParticle(i).getPosition()).norm());
    }

    // The pairs are the same set either way, so only the summation order can
    // differ -- which over 200 steps is worth a little, and nowhere near what a
    // dropped pair would be worth.
    EXPECT_LT(worstcpu, 0.05f) << "the CPU's list moved a particle " << worstcpu << " A from the cell walk";
    EXPECT_LT(worstgpu, 0.05f) << "the device's list moved a particle " << worstgpu << " A from the cell walk";

    // And the structure itself, against the O(N^2) answer. Comparing positions
    // is far too blunt on its own: a pair dropped near the cutoff carries
    // almost no force, so a list built at the wrong radius passes the check
    // above and fails here.
    const SpringNetworkOpenCL::NeighbourList & list = gpuskinned.electrostaticList();
    ASSERT_TRUE(list.valid) << "the device built no list";
    EXPECT_NEAR(list.radius, 12.0f + static_cast<float>(SKIN), 1.0e-4f)
        << "the list was not built at cutoff + skin";

    // The list dates from its last rebuild, not from the positions it is being
    // read at, so it is not the neighbour set of this instant -- it is a
    // SUPERSET of it, and that is precisely the property that makes it correct:
    // every pair within the cutoff now has to be in it, whatever has moved
    // since, or the force on that pair was silently dropped.
    size_t held = 0, inside = 0;
    for (unsigned i = 0; i < N; ++i)
    {
        std::vector<unsigned> found = gpuskinned.neighboursFromList(list, i);
        std::vector<unsigned> within = neighborsByBruteForce(gpuskinned, i, 12.0f);
        std::sort(found.begin(), found.end());
        std::sort(within.begin(), within.end());
        ASSERT_TRUE(std::includes(found.begin(), found.end(), within.begin(), within.end()))
            << "particle " << i << " has neighbours inside the cutoff that its list never held";
        held += found.size();
        inside += within.size();
    }
    EXPECT_GT(inside, 1000u) << "the lists are too empty for the comparison to mean anything";
    EXPECT_GT(held, inside) << "the list holds nothing beyond the cutoff, so it has no margin to live on";
}

// Steric, on the device, against the CPU -- once per law.
//
// The four laws are the reason this term is not Coulomb with another formula:
// a .msp picks one by name, the host resolves it to an integer, and one kernel
// serves all four. Testing only the default would leave three untested, and one
// of those is what the rigid-body examples actually use.
//
// The particles sit on a LATTICE at 2.6 A, just outside the 2 A minimum the
// Good-Hope rule gives two 2 A radii, with a small jitter so the arrangement is
// not symmetric. That spacing is the whole difficulty of testing this term.
// The first version dropped 250 particles into a 14 A box, which is a mean
// spacing of 2.2 A -- inside the repulsive wall, where r^-13 is astronomical.
// The CPU alone reached 9.7e7 A in a hundred steps, and comparing two backends
// on an explosion measures which one rounded first. The linear law, the only
// bounded one, agreed to 4.8e-07 A throughout, which is what said the setup was
// at fault rather than the code.
//
// Neutral particles, so nothing but the steric law moves the cloud.
namespace
{
void buildLattice(spn::SpringNetwork & network, configuration::Configuration & config,
                  const std::string & mode)
{
    const int SIDE = 6;             // 216 particles
    const float SPACING = 2.6f;

    unsigned state = 4242u;
    const auto jitter = [&state]() {
        state = state * 1103515245u + 12345u;
        return (static_cast<float>((state >> 16) & 0x7fffu) / static_cast<float>(0x7fff) - 0.5f) * 1.2f;
    };

    for (int x = 0; x < SIDE; ++x)
        for (int y = 0; y < SIDE; ++y)
            for (int z = 0; z < SIDE; ++z)
            {
                spn::Particle p;
                p.setPosition(Vector3f(x * SPACING + jitter(), y * SPACING + jitter(),
                                       z * SPACING + jitter()));
                p.setMass(12.0f);
                p.setCharge(0.0f);
                p.setRadius(2.0f);
                p.setEpsilon(0.5f);
                network.addParticle(p);
            }

    // A sprung pair off to the side, placed AT its rest length: every one of
    // these laws would throw two particles 1.4 A apart violently apart, so the
    // pair staying there is what says the exclusion ran. Placing them at random
    // and trusting the spring to pull them together, as the first version did,
    // measures the spring rather than the exclusion.
    const float away = SIDE * SPACING + 20.0f;
    for (int i = 0; i < 2; ++i)
    {
        spn::Particle p;
        p.setPosition(Vector3f(away + i * 1.4f, away, away));
        p.setMass(12.0f);
        p.setCharge(0.0f);
        p.setRadius(2.0f);
        p.setEpsilon(0.5f);
        network.addParticle(p);
    }
    const unsigned first = SIDE * SIDE * SIDE;
    network.addSpring(first, first + 1, /*equilibrium=*/1.4f, /*stiffness=*/500.0f);

    config = configuration::defaultConfiguration();
    config.sim.nbsteps = 400;
    config.sim.timestep = 0.2;
    config.spring.enable = true;
    config.spring.scale = 1.0;
    config.viscosity.enable = true;
    config.viscosity.value = 0.2;
    config.electrostatic.enable = false;
    config.hydrophobicity.enable = false;
    config.steric.enable = true;
    config.steric.mode = mode;
    config.steric.gridscale = 1.0;
    config.steric.cutoff = 8.0;

    network.setup(config);
}
} // namespace

TEST(SpringNetworkOpenCL, StericMatchesTheCPUInEveryMode)
{
    if (!hasOpenCLDevice())
        GTEST_SKIP() << "no OpenCL device available on this machine";

    const char * modes[] = {"linear", "lennard-jones-12-6Amber", "lennard-jones-8-6Lewitt",
                            "lennard-jones-8-6Zacharias"};

    for (const char * mode : modes)
    {
        spn::SpringNetwork cpu;
        configuration::Configuration cpuconfig;
        buildLattice(cpu, cpuconfig, mode);
        cpu.run();

        SpringNetworkOpenCL gpu;
        configuration::Configuration gpuconfig;
        buildLattice(gpu, gpuconfig, mode);
        gpu.run();

        spn::SpringNetwork start;
        configuration::Configuration startconfig;
        buildLattice(start, startconfig, mode);

        const unsigned n = cpu.getNumberOfParticles();
        float moved = 0.0f, worst = 0.0f;
        for (unsigned i = 0; i < n; ++i)
        {
            moved = std::max(moved,
                (cpu.getParticle(i).getPosition() - start.getParticle(i).getPosition()).norm());
            worst = std::max(worst,
                (cpu.getParticle(i).getPosition() - gpu.getParticle(i).getPosition()).norm());
        }

        EXPECT_GT(moved, 0.05f) << mode << ": the lattice barely moved, so this comparison proves nothing";
        // And it must not have gone off: comparing two backends on an explosion
        // measures which one rounded first, not whether they agree on the law.
        EXPECT_LT(moved, 20.0f) << mode << ": the lattice blew up (" << moved
                                << " A); this is no longer a comparison of force laws";
        // 1e-3 A is 450 times the worst difference measured across the four
        // laws, which is 2.2e-06 A -- the two backends run the same text here,
        // so what is left is the order of a summation. Tight enough that a
        // wrong combination rule or a missing exclusion cannot hide in it.
        EXPECT_LT(worst, 1.0e-3f) << mode << ": the GPU ended up " << worst << " A from the CPU";

        const unsigned first = n - 2;
        const float held = (gpu.getParticle(first).getPosition()
                            - gpu.getParticle(first + 1).getPosition()).norm();
        EXPECT_NEAR(held, 1.4f, 0.2f) << mode << ": the sprung pair was blown apart; steric is not excluding it";
    }
}

// Hydrophobicity, on the device, against the CPU.
//
// The third pairwise term, and the simplest law of the three: -h1*h2*exp(-r),
// an attraction with no repulsive core at all. That is what makes its test
// setup delicate rather than its arithmetic -- two hydrophobic particles with
// nothing to hold them apart approach without limit, and comparing two
// backends on a collapse measures which one rounded first, as the Coulomb test
// found the hard way.
//
// So the cloud is held by a spring network: every particle is sprung to the
// next, at a rest length the attraction cannot pull them far below. Half the
// particles are hydrophobic and half are not, so the term has something to
// distinguish and the sum is not a single uniform contraction.
namespace
{
void buildHydrophobicChain(spn::SpringNetwork & network, configuration::Configuration & config,
                           unsigned n)
{
    unsigned state = 31337u;
    const auto next = [&state]() {
        state = state * 1103515245u + 12345u;
        return static_cast<float>((state >> 16) & 0x7fffu) / static_cast<float>(0x7fff);
    };

    // A self-avoiding-ish walk with a FIXED 3.8 A step, so the chain starts at
    // its springs' rest length. Dropping the particles at random instead, as
    // the first version did, makes the springs contract a tangle: the chain
    // then moved 10.4 A whatever the hydrophobicity was -- including at 0.3 --
    // which is the tell that the test was measuring the springs.
    const float STEP = 3.8f;
    Vector3f at(0.0f, 0.0f, 0.0f);
    for (unsigned i = 0; i < n; ++i)
    {
        spn::Particle p;
        p.setPosition(at);
        p.setMass(12.0f);
        p.setCharge(0.0f);
        // In blocks of two, not alternating. Alternating looks like the
        // better test -- neighbours differ, the sum is not uniform -- and
        // quietly makes the exclusion untestable: the law is a product, so
        // every sprung pair would have one hydrophobic end and one not, and
        // contribute nothing whether it is excluded or not. In blocks, half
        // the sprung pairs are hydrophobic at both ends.
        p.setHydrophobicity(((i / 2) % 2 == 0) ? 1.0f : 0.0f);
        network.addParticle(p);

        // Self-avoiding: a plain walk folds back on itself, and two
        // non-adjacent particles starting 0.3 A apart is a bad initial
        // structure rather than anything the law did.
        Vector3f candidate = at;
        for (int attempt = 0; attempt < 200; ++attempt)
        {
            Vector3f dir(next() - 0.5f, next() - 0.5f, next() - 0.5f);
            dir.normalize();
            candidate = at + dir * STEP;
            float closest = 1.0e9f;
            for (unsigned j = 0; j < network.getNumberOfParticles(); ++j)
                closest = std::min(closest, (candidate - network.getParticle(j).getPosition()).norm());
            if (closest > 3.0f)
                break;
        }
        at = candidate;
    }
    for (unsigned i = 0; i + 1 < n; ++i)
        network.addSpring(i, i + 1, /*equilibrium=*/STEP, /*stiffness=*/200.0f);

    config = configuration::defaultConfiguration();
    config.sim.nbsteps = 200;
    config.sim.timestep = 0.5;
    config.spring.enable = true;
    config.spring.scale = 1.0;
    config.viscosity.enable = true;
    config.viscosity.value = 1.0;
    config.steric.enable = false;
    config.electrostatic.enable = false;
    config.hydrophobicity.enable = true;
    // 1000 was the value while the law decayed over 1 A instead of 10, which
    // made it exp(3.8 - 0.38) = 30 times weaker at this chain's step. Left
    // there, the corrected law piles the chain onto itself (closest pair
    // 0.0045 A) and stops being a comparison of force laws.
    config.hydrophobicity.scale = 30.0;
    config.hydrophobicity.cutoff = 12.0;
    // decaylength is left at its default (10 A) on purpose, so that a change
    // to that default shows up here.

    network.setup(config);
}
} // namespace

TEST(SpringNetworkOpenCL, HydrophobicityMatchesTheCPU)
{
    if (!hasOpenCLDevice())
        GTEST_SKIP() << "no OpenCL device available on this machine";

    const unsigned N = 200;

    spn::SpringNetwork cpu;
    configuration::Configuration cpuconfig;
    buildHydrophobicChain(cpu, cpuconfig, N);
    cpu.run();

    SpringNetworkOpenCL gpu;
    configuration::Configuration gpuconfig;
    buildHydrophobicChain(gpu, gpuconfig, N);
    gpu.run();

    spn::SpringNetwork start;
    configuration::Configuration startconfig;
    buildHydrophobicChain(start, startconfig, N);

    float moved = 0.0f, worst = 0.0f;
    for (unsigned i = 0; i < N; ++i)
    {
        moved = std::max(moved,
            (cpu.getParticle(i).getPosition() - start.getParticle(i).getPosition()).norm());
        worst = std::max(worst,
            (cpu.getParticle(i).getPosition() - gpu.getParticle(i).getPosition()).norm());
    }

    // The chain contracts by about 0.58 A under the attraction, and the closest
    // any two particles come is 2.78 A, having started no closer than 3.0.
    EXPECT_GT(moved, 0.1f) << "the chain barely moved, so this comparison proves nothing";
    EXPECT_LT(moved, 20.0f) << "the chain ran away (" << moved
                            << " A); this is no longer a comparison of force laws";

    // This law is an ATTRACTION with no repulsive core at all, so it will pile
    // particles up given enough strength -- at scale 300 the closest pair
    // reaches 0.034 A. Unlike Coulomb's 1/r^2 the force stays bounded
    // (exp(-r/L) -> 1), so the collapse is slow rather than singular, which
    // makes it easy to miss: hence the guard, at a distance two atoms could
    // actually be. The usable window for this fixture is scale 10 to 100:
    // below it the chain barely moves, above it it piles up.
    EXPECT_GT(closestPair(cpu), 1.5f) << "the chain piled up; this is no longer a physical configuration";

    // Measured agreement is 5.4e-06 A. The bar is over two orders above it, and
    // still far below anything a wrong law or a missing exclusion would cause.
    EXPECT_LT(worst, 1.0e-3f) << "the GPU's hydrophobicity ended up " << worst << " A from the CPU's";
}

// An interactor is the only way a force from outside the model reaches a
// particle: MDDriver's pull arrives through Interactor::syncSystemStateData(),
// which calls SpringNetwork::setForce(), which adds to the particle's force for
// the step about to be computed. This is the same entry point, with the same
// call, so whatever it exercises is what an interactive pull exercises.
namespace
{

class ConstantPull : public Interactor
{
  public:
    explicit ConstantPull(float fx) : _fx(fx) {}
    bool continueInteractionThread() override { return false; }
    void stopInteractionThread() override {}
    void startInteractionThread() override {}   // no thread: this pulls in-line
    void syncParticleStateData(unsigned) override {}
    void setupInteraction() override {}
    void processInteractions() override {}

    void syncSystemStateData() override
    {
        float force[3] = {_fx, 0.0f, 0.0f};
        _springnetwork->setForce(0, force);
    }

  private:
    float _fx;
};

void buildFreeParticlePair(spn::SpringNetwork & network, configuration::Configuration & config)
{
    // Nothing but the integrator: no spring, no pairwise term, no viscosity.
    // Whatever the particle does is then entirely the external force's doing,
    // which is what this test is about.
    spn::Particle first;
    first.setPosition(Vector3f(0.0f, 0.0f, 0.0f));
    first.setMass(1.0f);
    first.setDynamic(true);
    network.addParticle(first);

    spn::Particle second;
    second.setPosition(Vector3f(50.0f, 0.0f, 0.0f));
    second.setMass(1.0f);
    second.setDynamic(true);
    network.addParticle(second);

    config = configuration::defaultConfiguration();
    config.sim.nbsteps = 100;
    config.sim.timestep = 1.0;
    config.spring.enable = false;
    config.steric.enable = false;
    config.electrostatic.enable = false;
    config.viscosity.enable = false;

    network.setup(config);
}

} // namespace

// A pull applied through an interactor has to move the particle on the device
// exactly as it does on the CPU.
//
// It did not. SpringNetworkOpenCL::setForce() wrote _particleforces, a host
// array that the end of every step overwrites with what the device just
// returned, and that nothing ever uploads. The buffer that WOULD have carried
// it, _inExternalForceBuffer, was fed from an array nobody wrote, so the device
// added zeros to its forces at every step and an interactive pull did nothing
// unless it went through the probe kernel.
TEST(SpringNetworkOpenCL, AnInteractorsPullReachesTheDevice)
{
    if (!hasOpenCLDevice())
        GTEST_SKIP() << "no OpenCL device available on this machine";

    const float PULL = 1.0f;

    spn::SpringNetwork cpu;
    configuration::Configuration cpuconfig;
    buildFreeParticlePair(cpu, cpuconfig);
    ConstantPull cpupull(PULL);
    cpupull.setSpringNetwork(&cpu);
    cpu.addInteractor(&cpupull);
    cpu.run();

    SpringNetworkOpenCL gpu;
    configuration::Configuration gpuconfig;
    buildFreeParticlePair(gpu, gpuconfig);
    ConstantPull gpupull(PULL);
    gpupull.setSpringNetwork(&gpu);
    gpu.addInteractor(&gpupull);
    gpu.run();

    const float cpux = cpu.getParticle(0).getPosition().getX();
    const float gpux = gpu.getParticle(0).getPosition().getX();

    // Guard against agreeing because neither moved, which is exactly how this
    // defect hid: a pull that does nothing on both sides compares perfectly.
    EXPECT_GT(cpux, 1.0f) << "the CPU did not apply the interactor's force at all";
    EXPECT_GT(gpux, 1.0f) << "the device never received the interactor's force";

    EXPECT_NEAR(cpux, gpux, 1.0e-3f)
        << "CPU and OpenCL disagree on where an external pull takes the particle";

    // The particle that was not pulled must not have moved on either side.
    EXPECT_NEAR(cpu.getParticle(1).getPosition().getX(), 50.0f, 1.0e-4f);
    EXPECT_NEAR(gpu.getParticle(1).getPosition().getX(), 50.0f, 1.0e-4f);
}

// The bounding box the cell grids are measured against, computed on the device
// instead of by a host loop over the positions.
//
// That host loop is one of the three reasons the positions have to come down
// at every step -- the others being the frame check and the list rebuild
// criterion -- so it is the first brick of getting them to follow the sample
// rate like everything else already does. It has to give the SAME box: a
// reduction that disagrees by a hair sizes every grid differently, which
// changes which pairs each walk finds, and nothing about that looks wrong from
// outside.
TEST(SpringNetworkOpenCL, TheDeviceMeasuresTheSameBoxAsTheHost)
{
    if (!hasOpenCLDevice())
        GTEST_SKIP() << "no OpenCL device available on this machine";

    SpringNetworkOpenCL gpu;
    configuration::Configuration config;
    // 2000 particles is eight blocks of 256 with a partial one at the end,
    // which is where a tree reduction goes wrong if it goes wrong at all.
    const unsigned N = 2000;
    buildParticleCloud(gpu, config, N, /*extent=*/60.0f);

    // The extremes are PLACED, and placed in the LAST block. A cloud alone
    // leaves them wherever chance puts them: a first version of this test
    // passed happily with the final block dropped from the reduction, because
    // the six extremes happened to sit earlier. Pinning them here is what
    // makes the test able to fail.
    gpu.getParticle(N - 1).setPosition(Vector3f(-999.0f, -888.0f, -777.0f));
    gpu.getParticle(N - 2).setPosition(Vector3f(555.0f, 666.0f, 777.0f));

    // Measured on the structure as given, before any step moves it, so the
    // expected answer is exactly the two corners above.
    gpu.initRun();

    float lo[3], hi[3];
    ASSERT_TRUE(gpu._measureBoundsOnDevice(lo, hi)) << "the device reported a non-finite coordinate";

    const float expectedlo[3] = {-999.0f, -888.0f, -777.0f};
    const float expectedhi[3] = {555.0f, 666.0f, 777.0f};

    // Exactly equal, not nearly: a minimum and a maximum are selections, not
    // sums, so no reordering can change them and there is no tolerance to
    // grant.
    for (int d = 0; d < 3; ++d)
    {
        EXPECT_FLOAT_EQ(lo[d], expectedlo[d]) << "lower bound on axis " << d;
        EXPECT_FLOAT_EQ(hi[d], expectedhi[d]) << "upper bound on axis " << d;
    }
}

// The frame check, asked of the device instead of walked on the host.
//
// A particle outside the grid's frame is binned nowhere and is invisible to
// every neighbour walk, so getting this wrong is silent: no crash, no warning,
// just a term quietly missing its pairs. The two answers have to agree.
TEST(SpringNetworkOpenCL, TheDeviceChecksTheFrameLikeTheHost)
{
    if (!hasOpenCLDevice())
        GTEST_SKIP() << "no OpenCL device available on this machine";

    const unsigned N = 2000;

    {   // Everyone inside: both must say so.
        SpringNetworkOpenCL gpu;
        configuration::Configuration config;
        buildParticleCloud(gpu, config, N, /*extent=*/60.0f);
        gpu.run();
        EXPECT_TRUE(gpu._frameStillHoldsOnDevice(gpu.stericCells()))
            << "a settled cloud should sit inside the frame measured around it";
    }

    {   // One particle far outside, and in the LAST block, where a tree
        // reduction drops it if it drops anything. Placed after the run so
        // nothing re-measures the frame around it.
        SpringNetworkOpenCL gpu;
        configuration::Configuration config;
        buildParticleCloud(gpu, config, N, /*extent=*/60.0f);
        gpu.run();
        ASSERT_TRUE(gpu._frameStillHoldsOnDevice(gpu.stericCells())) << "precondition";

        gpu.getParticle(N - 1).setPosition(Vector3f(1.0e6f, 0.0f, 0.0f));
        gpu.uploadPositionsForTesting();
        EXPECT_FALSE(gpu._frameStillHoldsOnDevice(gpu.stericCells()))
            << "a particle a million angstroms away is not inside anything";
    }
}

// The drift criterion, asked of the device instead of walked on the host.
//
// This is the one whose two failure modes are asymmetric. Too eager, and the
// lists are rebuilt for nothing, which only costs time. Too lazy, and a pair
// that has come within the cutoff is missing from a list nobody rebuilt, which
// costs correctness and says nothing.
TEST(SpringNetworkOpenCL, TheDeviceKnowsWhenTheListsHaveDrifted)
{
    if (!hasOpenCLDevice())
        GTEST_SKIP() << "no OpenCL device available on this machine";

    const unsigned N = 2000;
    SpringNetworkOpenCL gpu;
    configuration::Configuration config;
    buildParticleCloud(gpu, config, N, /*extent=*/60.0f);
    gpu.run();

    // With the reference taken from where the particles stand, nothing has
    // drifted by definition.
    gpu._snapshotListReferenceOnDevice();
    EXPECT_FALSE(gpu._listsNeedRebuildingOnDevice())
        << "a reference just taken cannot already be stale";

    // Half the skin is the threshold. Just under it is not a drift...
    const float half = 0.5f * gpu.getNeighborSkin();
    ASSERT_GT(half, 0.0f) << "this test needs a skin to measure against";
    Vector3f home = gpu.getParticle(N - 1).getPosition();
    gpu.getParticle(N - 1).setPosition(home + Vector3f(0.9f * half, 0.0f, 0.0f));
    gpu.uploadPositionsForTesting();
    EXPECT_FALSE(gpu._listsNeedRebuildingOnDevice())
        << "nine tenths of half a skin is still inside the margin";

    // ... and just over it is. The particle is the last one, so it lands in
    // the final block: a reduction that drops that block passes the case above
    // and fails here, which is the point of moving this one rather than any.
    gpu.getParticle(N - 1).setPosition(home + Vector3f(1.1f * half, 0.0f, 0.0f));
    gpu.uploadPositionsForTesting();
    EXPECT_TRUE(gpu._listsNeedRebuildingOnDevice())
        << "past half a skin the stored lists can be missing a pair";
}
