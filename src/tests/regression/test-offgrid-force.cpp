
#include <gtest/gtest.h>

#include "Particle.h"
#include "SpringNetwork.h"
#include "configuration/Configuration.hpp"
#include "grid/PotentialGrid.hpp"

using namespace biospring;

// A steered or free-moving particle can leave the electrostatic/density grid.
// DenseGrid::at() throws std::out_of_range for an out-of-bounds cell, which used
// to abort the whole simulation (seen when steering prot-dna off-grid via IMD).
// The grids are shaped directly here rather than loaded from a .dx file, so the
// test exercises the force paths without needing an OpenDX fixture.
struct TestOffGridForce : public ::testing::Test
{
    configuration::Configuration config;
    spn::SpringNetwork spn;
    spn::Particle inside, outside;

    void SetUp() override
    {
        ::testing::Test::SetUp();
        config.sim.nbsteps = 1;
        config.sim.timestep = 0.01;

        // Well inside the grid built below.
        inside.setPosition(Vector3f(5.0, 5.0, 5.0));
        inside.setCharge(1.0);
        inside.setMass(12.0);

        // Far outside it, in every direction.
        outside.setPosition(Vector3f(1000.0, 1000.0, 1000.0));
        outside.setCharge(1.0);
        outside.setMass(12.0);

        spn.addParticle(inside);
        spn.addParticle(outside);
        spn.setup(config);

        for (grid::PotentialGrid * g : {&spn.getPotentialGrid(), &spn.getDensityGrid()})
        {
            g->reshape({0.0, 0.0, 0.0, 10.0, 10.0, 10.0}, {1.0, 1.0, 1.0});
            for (size_t i = 0; i < g->shape()[0]; ++i)
                for (size_t j = 0; j < g->shape()[1]; ++j)
                    for (size_t k = 0; k < g->shape()[2]; ++k)
                    {
                        grid::PotentialCell & cell =
                            g->at(grid::discrete_coordinates({static_cast<int>(i), static_cast<int>(j),
                                                              static_cast<int>(k)}));
                        cell.scalar = 1.0;
                        cell.vector = Vector3f(1.0, 0.0, 0.0);
                    }
        }
    }
};

// The regression itself: the off-grid lookup must not throw. Before the guard
// this propagated std::out_of_range out of the force loop and killed the run.
TEST_F(TestOffGridForce, OffGridElectrostaticFieldForceDoesNotThrow)
{
    EXPECT_NO_THROW(spn.getParticle(1).addElectrostaticFieldForce());
}

TEST_F(TestOffGridForce, OffGridDensityFieldForceDoesNotThrow)
{
    EXPECT_NO_THROW(spn.getParticle(1).addDensityFieldForce());
}

// Off-grid contributes exactly zero, matching the JAX port, rather than the
// value of the nearest edge cell -- clamping to the edge would give a constant
// force of unlimited range with a discontinuity at the boundary.
TEST_F(TestOffGridForce, OffGridParticleGetsZeroForceAndEnergy)
{
    spn::Particle & p = spn.getParticle(1);

    p.addElectrostaticFieldForce();
    p.addDensityFieldForce();

    EXPECT_FLOAT_EQ(p.getForce().getX(), 0.0);
    EXPECT_FLOAT_EQ(p.getForce().getY(), 0.0);
    EXPECT_FLOAT_EQ(p.getForce().getZ(), 0.0);
    EXPECT_FLOAT_EQ(p.getElectrostaticEnergy(), 0.0);
}

// The guard must only fire off-grid: a particle inside still feels the field.
TEST_F(TestOffGridForce, InGridParticleStillFeelsTheField)
{
    spn::Particle & p = spn.getParticle(0);

    p.addElectrostaticFieldForce();

    EXPECT_GT(p.getForce().getX(), 0.0);
}
