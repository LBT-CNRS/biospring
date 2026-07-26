
#include <gtest/gtest.h>

#include "Particle.h"
#include "SpringNetwork.h"
#include "configuration/Configuration.hpp"
#include "interactor/freesasa/InteractorFreeSASA.h"

using namespace biospring;

namespace
{

// The areas are produced by the interactor's worker thread. These tests drive
// the main-thread side on its own, which is exactly the window the guard
// protects: startInteractionThread() publishes _isRunning before the worker has
// computed anything, so syncSystemStateData can arrive first.
struct TestInteractorFreeSASA : public ::testing::Test
{
    configuration::Configuration config;
    spn::SpringNetwork spn;
    spn::Particle p1, p2;

    void SetUp() override
    {
        ::testing::Test::SetUp();
        config.sim.nbsteps = 1;
        config.sim.timestep = 0.01;

        p1.setPosition(Vector3f(0.0, 0.0, 0.0));
        p1.setRadius(1.5);
        p1.setMass(12.0);
        p1.setSolventAccessibilitySurface(42.0);

        p2.setPosition(Vector3f(3.0, 0.0, 0.0));
        p2.setRadius(1.5);
        p2.setMass(12.0);
        p2.setSolventAccessibilitySurface(43.0);

        spn.addParticle(p1);
        spn.addParticle(p2);
        spn.setup(config);
    }
};

} // namespace

// Regression test for the startup race: syncing before the worker has computed
// any area used to dereference a null _sasa and segfault (observed on a
// 1329-bead IMPALA system). It must now be a no-op.
TEST_F(TestInteractorFreeSASA, SyncBeforeFirstComputationDoesNotCrash)
{
    interactor::InteractorFreeSASA interactor;
    interactor.setSpringNetwork(&spn);

    EXPECT_NO_THROW(interactor.syncSystemStateData());
}

// The guard must not "recover" by publishing zeros. SASA is a purely
// multiplicative factor in the IMPALA energy and force, so a surface of 0 does
// not degrade the result, it removes the membrane term entirely -- and the run
// still completes, which is how a wrong result passes for a good one. Leaving
// the previous values untouched keeps that failure loud instead of silent.
TEST_F(TestInteractorFreeSASA, SyncBeforeFirstComputationLeavesSurfacesUntouched)
{
    interactor::InteractorFreeSASA interactor;
    interactor.setSpringNetwork(&spn);

    interactor.syncSystemStateData();

    EXPECT_FLOAT_EQ(spn.getParticle(0).getSolventAccessibilitySurface(), 42.0);
    EXPECT_FLOAT_EQ(spn.getParticle(1).getSolventAccessibilitySurface(), 43.0);
}

// Repeated syncs while the areas are still missing must stay harmless: idleRun
// calls this every step, not once.
TEST_F(TestInteractorFreeSASA, RepeatedSyncBeforeFirstComputationStaysSafe)
{
    interactor::InteractorFreeSASA interactor;
    interactor.setSpringNetwork(&spn);

    for (int i = 0; i < 10; ++i)
        EXPECT_NO_THROW(interactor.syncSystemStateData());

    EXPECT_FLOAT_EQ(spn.getParticle(0).getSolventAccessibilitySurface(), 42.0);
}
