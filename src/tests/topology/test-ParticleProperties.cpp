#include "topology/ParticleProperties.hpp"
#include <gtest/gtest.h>

// ======================================================================================
// IMPProperties
// ======================================================================================

TEST(TestIMPProperties, getters_and_setters)
{
    biospring::topology::IMPProperties p;
    EXPECT_FLOAT_EQ(p.transfert_energy_by_accessible_surface(), 0.0);
    EXPECT_FLOAT_EQ(p.solvent_accessible_surface(), 0.0);

    p.set_transfert_energy_by_accessible_surface(42.0);
    p.set_solvent_accessible_surface(17.0);

    EXPECT_FLOAT_EQ(p.transfert_energy_by_accessible_surface(), 42.0);
    EXPECT_FLOAT_EQ(p.solvent_accessible_surface(), 17.0);
}

TEST(TestIMPProperties, equality_operator)
{
    biospring::topology::IMPProperties lhs;
    biospring::topology::IMPProperties rhs;

    EXPECT_TRUE(lhs == rhs);

    lhs.set_transfert_energy_by_accessible_surface(42.0);
    EXPECT_FALSE(lhs == rhs);
}

TEST(TestIMPProperties, builder)
{
    biospring::topology::IMPProperties p = biospring::topology::IMPProperties::build()
                                               .transfert_energy_by_accessible_surface(42.0)
                                               .solvent_accessible_surface(17.0);

    EXPECT_FLOAT_EQ(p.transfert_energy_by_accessible_surface(), 42.0);
    EXPECT_FLOAT_EQ(p.solvent_accessible_surface(), 17.0);
}

// ======================================================================================
// ParticleProperties
// ======================================================================================

TEST(ParticleProperties, builder)
{
    biospring::topology::ParticleProperties p = biospring::topology::ParticleProperties::build()
                                                    .mass(42.0)
                                                    .charge(17.0)
                                                    .radius(24.0)
                                                    .epsilon(30.0)
                                                    .temperature_factor(-2.0)
                                                    .occupancy(0.7)
                                                    .hydrophobicity(12.9)
                                                    .burying(1.394)
                                                    .imp(biospring::topology::IMPProperties::build()
                                                             .transfert_energy_by_accessible_surface(42.0)
                                                             .solvent_accessible_surface(17.0));

    EXPECT_FLOAT_EQ(p.mass(), 42.0);
    EXPECT_FLOAT_EQ(p.charge(), 17.0);
    EXPECT_FLOAT_EQ(p.radius(), 24.0);
    EXPECT_FLOAT_EQ(p.epsilon(), 30.0);
    EXPECT_FLOAT_EQ(p.temperature_factor(), -2.0);
    EXPECT_FLOAT_EQ(p.occupancy(), 0.7);
    EXPECT_FLOAT_EQ(p.hydrophobicity(), 12.9);
    EXPECT_FLOAT_EQ(p.burying(), 1.394);
    EXPECT_FLOAT_EQ(p.imp().transfert_energy_by_accessible_surface(), 42.0);
    EXPECT_FLOAT_EQ(p.imp().solvent_accessible_surface(), 17.0);

    p = p.build().mass(1.0).charge(0.0).radius(1.0);

    EXPECT_FLOAT_EQ(p.mass(), 1.0);
    EXPECT_FLOAT_EQ(p.charge(), 0.0);
    EXPECT_FLOAT_EQ(p.radius(), 1.0);
}

TEST(ParticleProperties, equality_operator)
{
    biospring::topology::ParticleProperties lhs;
    biospring::topology::ParticleProperties rhs;
    EXPECT_TRUE(lhs == rhs);

    lhs.set_mass(42.0);
    EXPECT_FALSE(lhs == rhs);
}

TEST(ParticleProperties, is_charged)
{
    biospring::topology::ParticleProperties p;
    EXPECT_FALSE(p.is_charged());

    p.set_charge(1.0);
    EXPECT_TRUE(p.is_charged());
}

TEST(ParticleProperties, is_hydrophogic)
{
    biospring::topology::ParticleProperties p;
    EXPECT_FALSE(p.is_hydrophobic());

    p.set_hydrophobicity(1.0);
    EXPECT_TRUE(p.is_hydrophobic());
}

TEST(ParticleProperties, mass_initialized_to_one)
{
    biospring::topology::ParticleProperties p;
    EXPECT_FLOAT_EQ(p.mass(), 1.0);
}

TEST(ParticleProperties, radius_initialized_to_one)
{
    biospring::topology::ParticleProperties p;
    EXPECT_FLOAT_EQ(p.radius(), 1.0);
}

TEST(ParticleProperties, burying_initialized_to_one)
{
    biospring::topology::ParticleProperties p;
    EXPECT_FLOAT_EQ(p.burying(), 1.0);
}

TEST(ParticleProperties, set_position)
{
    biospring::topology::ParticleProperties p;
    Vector3f position = p.position();
    EXPECT_FLOAT_EQ(position.getX(), 0.0);
    EXPECT_FLOAT_EQ(position.getY(), 0.0);
    EXPECT_FLOAT_EQ(position.getZ(), 0.0);

    position = {1, 2, 3};
    p.set_position(position);
    EXPECT_FLOAT_EQ(p.position().getX(), 1.0);
    EXPECT_FLOAT_EQ(p.position().getY(), 2.0);
    EXPECT_FLOAT_EQ(p.position().getZ(), 3.0);

    p.set_position({4, 5, 6});
    EXPECT_FLOAT_EQ(p.position().getX(), 4.0);
    EXPECT_FLOAT_EQ(p.position().getY(), 5.0);
    EXPECT_FLOAT_EQ(p.position().getZ(), 6.0);

    EXPECT_FLOAT_EQ(p.previous_position().getX(), 1.0);
    EXPECT_FLOAT_EQ(p.previous_position().getY(), 2.0);
    EXPECT_FLOAT_EQ(p.previous_position().getZ(), 3.0);
}

// -- Main function  ----------------------------------------------------------
int main(int argc, char * argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
