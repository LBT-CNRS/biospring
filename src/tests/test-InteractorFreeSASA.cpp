
#include <gtest/gtest.h>

#include <chrono>
#include <thread>

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
        p1.setName("CA");
        p1.setResName("ALA");
        p1.setRadius(1.5);
        p1.setMass(12.0);
        p1.setSolventAccessibilitySurface(42.0);

        p2.setPosition(Vector3f(3.0, 0.0, 0.0));
        p2.setName("CB");
        p2.setResName("ALA");
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

// THE SURFACES MUST EXIST BEFORE THE RUN DOES.
//
// They used to be produced by the worker thread while the main thread was
// already taking steps. syncParticleStateData rightly refuses to publish a
// surface it has not got -- the surface MULTIPLIES the IMPALA energy and force,
// so publishing zero would switch the term off in silence -- so the first steps
// ran on whatever surface the particles already carried, and HOW MANY steps
// that was depended on how fast the worker happened to be.
//
// That made the trajectory of an IMPALA run depend on timing, and it had a
// reproducible trigger: registering the MDDriver interactor makes
// SpringNetwork::run() sleep 100 ms before the loop, which is ample for the
// worker to finish. Measured on 052.Trypsin at three steps of 10 fs, the spring
// energy was 1.82 kJ/mol with MDDriver registered and 0.57 without -- the same
// binary on the same input, a transport deciding the physics.
//
// Ten repetitions, because what is being guarded against was a race: a single
// pass could pass on broken code whenever the worker happened to win.
TEST_F(TestInteractorFreeSASA, TheSurfacesAreReadyBeforeTheRunStarts)
{
    for (int attempt = 0; attempt < 10; ++attempt)
    {
        spn.getParticle(0).setSolventAccessibilitySurface(42.0);
        spn.getParticle(1).setSolventAccessibilitySurface(43.0);

        interactor::InteractorFreeSASA interactor;
        interactor.setSpringNetwork(&spn);
        // The constructor's own default, "default_classifier", matches none of
        // the names setupFreesasaInteractions tests and makes it bail out; the
        // CLI always passes one, so a test has to as well.
        interactor.setRadiiClassifier("default");
        interactor.setDynamic(false);

        interactor.startInteractionThread();
        const float a = spn.getParticle(0).getSolventAccessibilitySurface();
        const float b = spn.getParticle(1).getSolventAccessibilitySurface();
        interactor.stopInteractionThread();
        interactor.waitForInteractionThread();

        EXPECT_NE(a, 42.0f) << "attempt " << attempt
                            << ": particle 0 still carries its placeholder, so no surface was "
                               "computed before the run could start";
        EXPECT_NE(b, 43.0f) << "attempt " << attempt << ": same for particle 1";
        EXPECT_GT(a, 0.0f) << "attempt " << attempt << ": a surface of " << a
                           << " A^2 would switch the IMPALA term off in silence";
        EXPECT_GT(b, 0.0f) << "attempt " << attempt;
    }
}

// And nothing must arrive AFTERWARDS. This is the property the trajectory
// depended on: if a later pass of the worker replaces the surfaces, then how
// many steps ran on the old ones decides the result, and the run is no longer
// reproducible. Without --sasa-dynamic there is nothing to refresh, so the
// values read right after startInteractionThread() must be the values still
// there once any worker would have had its turn.
TEST_F(TestInteractorFreeSASA, NothingReplacesTheSurfacesLater)
{
    interactor::InteractorFreeSASA interactor;
    interactor.setSpringNetwork(&spn);
    interactor.setRadiiClassifier("default");
    interactor.setDynamic(false);

    interactor.startInteractionThread();
    const float a0 = spn.getParticle(0).getSolventAccessibilitySurface();
    const float b0 = spn.getParticle(1).getSolventAccessibilitySurface();

    // Far longer than the worker's own 1500 microsecond cadence, and longer
    // than the 100 ms MDDriver's registration used to sleep -- which is exactly
    // the window that used to change the answer.
    std::this_thread::sleep_for(std::chrono::milliseconds(200));
    interactor.syncSystemStateData();

    const float a1 = spn.getParticle(0).getSolventAccessibilitySurface();
    const float b1 = spn.getParticle(1).getSolventAccessibilitySurface();
    interactor.stopInteractionThread();
    interactor.waitForInteractionThread();

    EXPECT_FLOAT_EQ(a0, a1) << "particle 0's surface changed from " << a0 << " to " << a1
                            << " A^2 after the run had started, so the result depends on timing";
    EXPECT_FLOAT_EQ(b0, b1) << "particle 1's surface changed from " << b0 << " to " << b1 << " A^2";
}
