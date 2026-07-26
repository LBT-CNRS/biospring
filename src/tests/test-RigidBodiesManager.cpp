
#include <gtest/gtest.h>

#include "Particle.h"
#include "SpringNetwork.h"
#include "configuration/Configuration.hpp"
#include "rigidbody/RigidBodiesManager.h"

using namespace biospring;

namespace
{

// RigidBodiesManager keeps its collection in a static member, so it outlives any
// single test; clear it on both ends to keep cases independent. Note
// getCollection() hands back a copy of the pointer vector, so it can only be
// used to observe, never to mutate.
struct TestRigidBodiesManager : public ::testing::Test
{
    configuration::Configuration config;
    spn::SpringNetwork spn;

    void SetUp() override
    {
        ::testing::Test::SetUp();
        rigidbody::RigidBodiesManager::CleanRigidBodies();

        config.sim.nbsteps = 1;
        config.sim.timestep = 0.01;

        for (int i = 0; i < 4; ++i)
        {
            spn::Particle p;
            p.setPosition(Vector3f(static_cast<float>(i), 0.0, 0.0));
            p.setMass(12.0);
            spn.addParticle(p);
        }
        spn.setup(config);
    }

    void TearDown() override
    {
        rigidbody::RigidBodiesManager::CleanRigidBodies();
        ::testing::Test::TearDown();
    }
};

} // namespace

// The regression: a second group sharing no particle with the first is entirely
// new and must become its own body. The membership test was inverted, so
// "isAllNewParticle" was cleared as soon as any candidate was *absent* from an
// existing body -- which is exactly what being new means -- and every group
// after the first was silently dropped.
TEST_F(TestRigidBodiesManager, DisjointGroupCreatesASecondBody)
{
    rigidbody::RigidBodiesManager::InitRigidBodies(&spn, {0, 1});
    ASSERT_EQ(rigidbody::RigidBodiesManager::getCollection().size(), 1u);

    rigidbody::RigidBodiesManager::InitRigidBodies(&spn, {2, 3});
    EXPECT_EQ(rigidbody::RigidBodiesManager::getCollection().size(), 2u);
}

// The behaviour the flag is actually meant to provide: a group that reuses a
// particle already owned by a body must be rejected, not duplicated.
TEST_F(TestRigidBodiesManager, OverlappingGroupIsRejected)
{
    rigidbody::RigidBodiesManager::InitRigidBodies(&spn, {0, 1});
    ASSERT_EQ(rigidbody::RigidBodiesManager::getCollection().size(), 1u);

    rigidbody::RigidBodiesManager::InitRigidBodies(&spn, {1, 2});
    EXPECT_EQ(rigidbody::RigidBodiesManager::getCollection().size(), 1u);
}

// An empty request clears the collection rather than adding an empty body.
TEST_F(TestRigidBodiesManager, EmptyRequestClearsTheCollection)
{
    rigidbody::RigidBodiesManager::InitRigidBodies(&spn, {0, 1});
    ASSERT_EQ(rigidbody::RigidBodiesManager::getCollection().size(), 1u);

    rigidbody::RigidBodiesManager::InitRigidBodies(&spn, {});
    EXPECT_EQ(rigidbody::RigidBodiesManager::getCollection().size(), 0u);
}
