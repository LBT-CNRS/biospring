
#include <gtest/gtest.h>

#include "Particle.h"
#include "SpringNetwork.h"
#include "configuration/Configuration.hpp"
#include "rigidbody/RigidBody.h"

using namespace biospring;

namespace
{

// A rigid body keeps its own per-particle arrays (_p0 and friends) indexed by
// position in _particulesIds. A particle also carries two unrelated numbers: its
// id in the spring network (its array index) and its extid, which is the PDB
// atom id -- or the residue id for a coarse-grained model. The extids here are
// deliberately far larger than the particle count, which is the normal case for
// any real structure, so indexing a body array with one is out of bounds.
struct TestRigidBody : public ::testing::Test
{
    configuration::Configuration config;
    spn::SpringNetwork spn;

    void SetUp() override
    {
        ::testing::Test::SetUp();
        config.sim.nbsteps = 1;
        config.sim.timestep = 0.01;

        for (int i = 0; i < 3; ++i)
        {
            spn::Particle p;
            p.setPosition(Vector3f(static_cast<float>(i), 0.0, 0.0));
            p.setMass(12.0);
            p.setExtid(static_cast<unsigned>(500 + i));
            spn.addParticle(p);
        }
        spn.setup(config);
    }
};

} // namespace

// The regression: local index must come from the particle's position in the
// body, never from its extid. With extids of 500..502 and a 3-particle body,
// the old _p0[getExtid()] read ~500 elements past the end (AddressSanitizer
// reports it as a heap-buffer-overflow). _p0 itself is private and only filled
// by init(), which is reached through paths that call logging::die() -- i.e.
// exit() -- so what is locked here is the index mapping the fix introduced:
// as long as initImpalaSampling() goes through localIndexOf(), it cannot go
// back to indexing a body array by a PDB atom id.
TEST_F(TestRigidBody, LocalIndexIsPositionInBodyNotExtid)
{
    rigidbody::RigidBody body(&spn, 0, {0, 1, 2});

    for (unsigned i = 0; i < 3; ++i)
    {
        const spn::Particle & p = spn.getParticle(i);
        ASSERT_NE(p.getExtid(), static_cast<unsigned>(p.getId())) << "test fixture must keep the numberings distinct";
        EXPECT_EQ(body.localIndexOf(p), static_cast<size_t>(i));
        EXPECT_LT(body.localIndexOf(p), body.getParticlesIds().size());
    }
}

// A body holding only part of the network must map its own members to 0..n-1,
// not to their network ids -- the case where "local index == particle id"
// silently stops being true.
TEST_F(TestRigidBody, LocalIndexIsRelativeToTheBodyNotTheNetwork)
{
    rigidbody::RigidBody body(&spn, 0, {1, 2});

    EXPECT_EQ(body.localIndexOf(spn.getParticle(1)), 0u);
    EXPECT_EQ(body.localIndexOf(spn.getParticle(2)), 1u);
}
