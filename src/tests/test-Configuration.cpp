#include <gtest/gtest.h>

#include "configuration/Configuration.hpp"

using namespace biospring::configuration;

TEST(Configuration, assignement_operator)
{
    Configuration source;

    // Sets every source parameter to non-default values.
    source.sim.nbsteps = 100;
    source.sim.timestep = 0.1;
    source.sim.samplerate = 20;

    source.steric.enable = true;
    source.steric.gridscale = 0.3;
    source.steric.cutoff = 0.4;

    source.spring.enable = true;
    source.spring.cutoff = 0.5;
    source.spring.scale = 0.6;

    source.hydrophobicity.enable = true;
    source.hydrophobicity.cutoff = 0.7;
    source.hydrophobicity.scale = 0.8;

    source.electrostatic.enable = true;
    source.electrostatic.cutoff = 0.9;
    source.electrostatic.scale = 1.8;
    source.electrostatic.dielectric = 1.1;

    source.ivector.enable = true;

    source.viscosity.enable = true;
    source.viscosity.value = 1.2;

    source.probe.enable = true;
    source.probe.enableelectrostatic = true;
    source.probe.enablesteric = true;
    source.probe.x = 1.3;
    source.probe.y = 1.4;
    source.probe.z = 1.5;
    source.probe.mass = 1.6;
    source.probe.epsilon = 1.7;
    source.probe.radius = 1.8;
    source.probe.charge = 1.9;

    // Assigns source to target.
    Configuration target;
    target = source;

    EXPECT_EQ(source.sim.nbsteps, target.sim.nbsteps);
    EXPECT_FLOAT_EQ(source.sim.timestep, target.sim.timestep);
    EXPECT_FLOAT_EQ(source.sim.samplerate, target.sim.samplerate);

    EXPECT_EQ(source.steric.enable, target.steric.enable);
    EXPECT_FLOAT_EQ(source.steric.gridscale, target.steric.gridscale);
    EXPECT_FLOAT_EQ(source.steric.cutoff, target.steric.cutoff);
    EXPECT_EQ(source.steric.mode, target.steric.mode);

    EXPECT_EQ(source.spring.enable, target.spring.enable);
    EXPECT_FLOAT_EQ(source.spring.cutoff, target.spring.cutoff);
    EXPECT_FLOAT_EQ(source.spring.scale, target.spring.scale);

    EXPECT_EQ(source.hydrophobicity.enable, target.hydrophobicity.enable);
    EXPECT_FLOAT_EQ(source.hydrophobicity.cutoff, target.hydrophobicity.cutoff);
    EXPECT_FLOAT_EQ(source.hydrophobicity.scale, target.hydrophobicity.scale);

    EXPECT_EQ(source.electrostatic.enable, target.electrostatic.enable);
    EXPECT_FLOAT_EQ(source.electrostatic.cutoff, target.electrostatic.cutoff);
    EXPECT_FLOAT_EQ(source.electrostatic.scale, target.electrostatic.scale);
    EXPECT_FLOAT_EQ(source.electrostatic.dielectric, target.electrostatic.dielectric);

    EXPECT_EQ(source.ivector.enable, target.ivector.enable);

    EXPECT_EQ(source.viscosity.enable, target.viscosity.enable);
    EXPECT_FLOAT_EQ(source.viscosity.value, target.viscosity.value);

    EXPECT_EQ(source.probe.enable, target.probe.enable);
    EXPECT_EQ(source.probe.enableelectrostatic, target.probe.enableelectrostatic);
    EXPECT_EQ(source.probe.enablesteric, target.probe.enablesteric);
    EXPECT_FLOAT_EQ(source.probe.x, target.probe.x);
    EXPECT_FLOAT_EQ(source.probe.y, target.probe.y);
    EXPECT_FLOAT_EQ(source.probe.z, target.probe.z);
    EXPECT_FLOAT_EQ(source.probe.mass, target.probe.mass);
    EXPECT_FLOAT_EQ(source.probe.epsilon, target.probe.epsilon);
    EXPECT_FLOAT_EQ(source.probe.radius, target.probe.radius);
}

// -- Main function  ----------------------------------------------------------
int main(int argc, char * argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
