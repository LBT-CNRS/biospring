#include <gtest/gtest.h>

#include <vector>

#include "Particle.h"
#include "SpringNetwork.h"
#include "SpringNetworkOpenCL.h"
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
