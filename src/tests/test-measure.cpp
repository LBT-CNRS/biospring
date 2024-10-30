#include <gtest/gtest.h>

#include "Particle.h"
#include "measure.hpp"

#include <functional>
#include <vector>

using Particle = biospring::spn::Particle;

// =====================================================================================
// measure::box
TEST(TestMeasure, box)
{
    std::vector<Particle> particles(4);

    particles[0].setPosition({0.0, 0.0, 0.0});
    particles[1].setPosition({1.0, 0.0, 0.0});
    particles[2].setPosition({0.0, 1.0, 0.0});
    particles[3].setPosition({0.0, 0.0, 1.0});

    auto box = biospring::measure::box(particles);
    EXPECT_FLOAT_EQ(0.0, box.boundaries()[0]);
    EXPECT_FLOAT_EQ(0.0, box.boundaries()[1]);
    EXPECT_FLOAT_EQ(0.0, box.boundaries()[2]);

    EXPECT_FLOAT_EQ(1.0, box.boundaries()[3]);
    EXPECT_FLOAT_EQ(1.0, box.boundaries()[4]);
    EXPECT_FLOAT_EQ(1.0, box.boundaries()[5]);

    EXPECT_FLOAT_EQ(1.0, box.length()[0]);
    EXPECT_FLOAT_EQ(1.0, box.length()[1]);
    EXPECT_FLOAT_EQ(1.0, box.length()[2]);
}

// =====================================================================================
// measure::rmsd
TEST(TestMeasure, rmsd_with_empty_arrays)
{
    std::array<Vector3f, 0> lhs, rhs;
    EXPECT_FLOAT_EQ(0.0, biospring::measure::rmsd(lhs, rhs));
}

TEST(TestMeasure, rmsd_with_same_coordinates)
{
    std::vector<Particle> lhs(4), rhs(4);

    lhs[0].setPosition({0.0, 0.0, 0.0});
    lhs[1].setPosition({0.0, 0.0, 1.0});
    lhs[2].setPosition({0.0, 0.0, 2.0});
    lhs[3].setPosition({0.0, 0.0, 3.0});

    rhs[0].setPosition({0.0, 0.0, 0.0});
    rhs[1].setPosition({0.0, 0.0, 1.0});
    rhs[2].setPosition({0.0, 0.0, 2.0});
    rhs[3].setPosition({0.0, 0.0, 3.0});

    EXPECT_FLOAT_EQ(0.0, biospring::measure::rmsd(lhs, rhs));
}

TEST(TestMeasure, rmsd_with_different_coordinates)
{
    std::vector<Particle> lhs(4), rhs(4);

    lhs[0].setPosition({0.0, 0.0, 0.0});
    lhs[1].setPosition({0.0, 0.0, 0.0});
    lhs[2].setPosition({0.0, 0.0, 0.0});
    lhs[3].setPosition({0.0, 0.0, 0.0});

    rhs[0].setPosition({0.0, 0.0, 1.0});
    rhs[1].setPosition({0.0, 0.0, 1.0});
    rhs[2].setPosition({0.0, 0.0, 1.0});
    rhs[3].setPosition({0.0, 0.0, 1.0});

    EXPECT_FLOAT_EQ(1.0, biospring::measure::rmsd(lhs, rhs));
}

TEST(TestMeasure, rmsd_with_container_of_pointers)
{
    std::vector<Particle *> lhs, rhs;

    lhs.push_back(new Particle());
    lhs.push_back(new Particle());
    lhs.push_back(new Particle());
    lhs.push_back(new Particle());

    rhs.push_back(new Particle());
    rhs.push_back(new Particle());
    rhs.push_back(new Particle());
    rhs.push_back(new Particle());

    lhs[0]->setPosition({0.0, 0.0, 0.0});
    lhs[1]->setPosition({0.0, 0.0, 0.0});
    lhs[2]->setPosition({0.0, 0.0, 0.0});
    lhs[3]->setPosition({0.0, 0.0, 0.0});

    rhs[0]->setPosition({0.0, 0.0, 1.0});
    rhs[1]->setPosition({0.0, 0.0, 1.0});
    rhs[2]->setPosition({0.0, 0.0, 1.0});
    rhs[3]->setPosition({0.0, 0.0, 1.0});

    EXPECT_FLOAT_EQ(1.0, biospring::measure::rmsd(lhs, rhs));

    for (auto & p : lhs)
        delete p;

    for (auto & p : rhs)
        delete p;
}

// =====================================================================================
// measure::closest_to
TEST(TestMeasure, closest_to_with_empty_array)
{
    std::array<Vector3f, 0> particles;
    auto c = biospring::measure::closest_to(particles, {0.0, 0.0, 0.0});
    EXPECT_FLOAT_EQ(0.0, c[0]);
    EXPECT_FLOAT_EQ(0.0, c[1]);
    EXPECT_FLOAT_EQ(0.0, c[2]);
}

TEST(TestMeasure, closest_to_with_array_of_particles)
{
    std::array<Particle, 4> particles;
    particles[0].setPosition({0.0, 0.0, 0.0});
    particles[1].setPosition({0.0, 0.0, 1.0});
    particles[2].setPosition({0.0, 0.0, 2.0});
    particles[3].setPosition({0.0, 0.0, 3.0});

    auto result = biospring::measure::closest_to(particles, {0.0, 0.0, 0.0});
    EXPECT_FLOAT_EQ(0.0, result[0]);
    EXPECT_FLOAT_EQ(0.0, result[1]);
    EXPECT_FLOAT_EQ(0.0, result[2]);

    result = biospring::measure::closest_to(particles, {0.0, 0.0, 1.2});
    EXPECT_FLOAT_EQ(0.0, result[0]);
    EXPECT_FLOAT_EQ(0.0, result[1]);
    EXPECT_FLOAT_EQ(1.0, result[2]);

    result = biospring::measure::closest_to(particles, {0.0, 0.0, 2.4});
    EXPECT_FLOAT_EQ(0.0, result[0]);
    EXPECT_FLOAT_EQ(0.0, result[1]);
    EXPECT_FLOAT_EQ(2.0, result[2]);

    result = biospring::measure::closest_to(particles, {0.0, 0.0, -3.0});
    EXPECT_FLOAT_EQ(0.0, result[0]);
    EXPECT_FLOAT_EQ(0.0, result[1]);
    EXPECT_FLOAT_EQ(0.0, result[2]);
}

TEST(TestMeasure, closest_to_origin)
{
    std::vector<Particle> particles(4);

    particles[0].setPosition({0.0, 0.0, -1.0});
    particles[1].setPosition({0.0, 0.0, 1.0});
    particles[2].setPosition({0.0, 0.0, 2.0});
    particles[3].setPosition({0.0, 0.0, 3.0});

    auto result = biospring::measure::closest_to_origin(particles);
    EXPECT_FLOAT_EQ(0.0, result[0]);
    EXPECT_FLOAT_EQ(0.0, result[1]);
    EXPECT_FLOAT_EQ(-1.0, result[2]);
}

TEST(TestMeasure, closest_to_centroid)
{
    std::vector<Particle> particles(4);

    particles[0].setPosition({0.0, 0.0, 0.0});
    particles[1].setPosition({0.0, 0.0, 1.0});
    particles[2].setPosition({0.0, 0.0, 2.0});
    particles[3].setPosition({0.0, 0.0, 3.0});

    auto result = biospring::measure::closest_to_centroid(particles);

    auto centroid = biospring::measure::centroid(particles);
    ASSERT_FLOAT_EQ(0.0, centroid[0]);
    ASSERT_FLOAT_EQ(0.0, centroid[1]);
    ASSERT_FLOAT_EQ(1.5, centroid[2]);

    EXPECT_FLOAT_EQ(0.0, result[0]);
    EXPECT_FLOAT_EQ(0.0, result[1]);
    EXPECT_FLOAT_EQ(1.0, result[2]);
}

// =====================================================================================
// measure::radius
TEST(TestMeasure, radius)
{
    std::vector<Particle> particles(4);

    particles[0].setPosition({0.0, 0.0, 0.0});
    particles[1].setPosition({0.0, 0.0, 1.0});
    particles[2].setPosition({0.0, 0.0, 2.0});
    particles[3].setPosition({0.0, 0.0, 3.0});

    EXPECT_FLOAT_EQ(1.5, biospring::measure::radius(particles));
}

// =====================================================================================
// measure::centroid
TEST(TestMeasure, centroid_with_empty_array)
{
    std::array<Vector3f, 0> particles;
    auto c = biospring::measure::centroid(particles);
    EXPECT_FLOAT_EQ(0.0, c[0]);
    EXPECT_FLOAT_EQ(0.0, c[1]);
    EXPECT_FLOAT_EQ(0.0, c[2]);
}

TEST(TestMeasure, centroid_with_array_of_one_particle)
{
    std::array<Vector3f, 1> particles;
    particles[0] = {0.0, 0.0, 0.0};
    auto c = biospring::measure::centroid(particles);
    EXPECT_FLOAT_EQ(0.0, c[0]);
    EXPECT_FLOAT_EQ(0.0, c[1]);
    EXPECT_FLOAT_EQ(0.0, c[2]);
}

TEST(TestMeasure, centroid_with_array_of_two_particles)
{
    std::array<Vector3f, 2> particles;
    particles[0] = {0.0, 0.0, 0.0};
    particles[1] = {0.0, 0.0, 1.0};
    auto c = biospring::measure::centroid(particles);
    EXPECT_FLOAT_EQ(0.0, c[0]);
    EXPECT_FLOAT_EQ(0.0, c[1]);
    EXPECT_FLOAT_EQ(0.5, c[2]);
}

TEST(TestMeasure, centroid_with_array_of_two_particles_pointers)
{
    std::array<Particle *, 2> particles;
    particles[0] = new Particle();
    particles[1] = new Particle();

    particles[0]->setPosition({0.0, 0.0, 0.0});
    particles[1]->setPosition({0.0, 0.0, 1.0});

    auto c = biospring::measure::centroid(particles);
    EXPECT_FLOAT_EQ(0.0, c[0]);
    EXPECT_FLOAT_EQ(0.0, c[1]);
    EXPECT_FLOAT_EQ(0.5, c[2]);

    delete particles[0];
    delete particles[1];
}

TEST(TestMeasure, centroid_with_vector_of_references)
{
    std::vector<Particle> topology(2);
    topology[0].setPosition({0.0, 0.0, 0.0});
    topology[1].setPosition({0.0, 0.0, 1.0});

    std::vector<std::reference_wrapper<Particle>> particles;
    particles.push_back(topology[0]);
    particles.push_back(topology[1]);

    auto c = biospring::measure::centroid(particles);
    EXPECT_FLOAT_EQ(0.0, c[0]);
    EXPECT_FLOAT_EQ(0.0, c[1]);
    EXPECT_FLOAT_EQ(0.5, c[2]);
}

// =====================================================================================
// measure::distance
TEST(TestMeasure, distance_with_particles)
{
    Particle lhs, rhs;
    lhs.setPosition({0.0, 0.0, 0.0});

    rhs.setPosition({0.0, 0.0, 0.0});
    EXPECT_FLOAT_EQ(0.0, biospring::measure::distance(lhs, rhs));

    rhs.setPosition({0.0, 0.0, 1.0});
    EXPECT_FLOAT_EQ(1.0, biospring::measure::distance(lhs, rhs));

    rhs.setPosition({0.0, 0.0, 2.0});
    EXPECT_FLOAT_EQ(2.0, biospring::measure::distance(lhs, rhs));
}

TEST(TestMeasure, distance_with_particle_pointers)
{
    Particle lhs, rhs;
    lhs.setPosition({0.0, 0.0, 0.0});

    rhs.setPosition({0.0, 0.0, 0.0});
    EXPECT_FLOAT_EQ(0.0, biospring::measure::distance(&lhs, &rhs));

    rhs.setPosition({0.0, 0.0, 1.0});
    EXPECT_FLOAT_EQ(1.0, biospring::measure::distance(&lhs, &rhs));

    rhs.setPosition({0.0, 0.0, 2.0});
    EXPECT_FLOAT_EQ(2.0, biospring::measure::distance(&lhs, &rhs));
}

TEST(TestMeasure, distance_with_arrays)
{
    std::array<double, 3> lhs, rhs;
    lhs = {0.0, 0.0, 0.0};

    rhs = {0.0, 0.0, 0.0};
    EXPECT_FLOAT_EQ(0.0, biospring::measure::distance(lhs, rhs));

    rhs = {0.0, 0.0, 1.0};
    EXPECT_FLOAT_EQ(1.0, biospring::measure::distance(lhs, rhs));

    rhs = {0.0, 0.0, 2.0};
    EXPECT_FLOAT_EQ(2.0, biospring::measure::distance(lhs, rhs));
}

TEST(TestMeasure, distance_with_mixed_types)
{
    Particle lhs;
    std::array<double, 3> rhs;
    lhs.setPosition({0.0, 0.0, 0.0});

    // =======================================
    // distance(Particle, array)
    rhs = {0.0, 0.0, 0.0};
    EXPECT_FLOAT_EQ(0.0, biospring::measure::distance(lhs, rhs));

    rhs = {0.0, 0.0, 1.0};
    EXPECT_FLOAT_EQ(1.0, biospring::measure::distance(lhs, rhs));

    rhs = {0.0, 0.0, 2.0};
    EXPECT_FLOAT_EQ(2.0, biospring::measure::distance(lhs, rhs));

    // =======================================
    // distance(Particle, Vector3f)
    Vector3f rhs2;
    rhs2 = {0.0, 0.0, 0.0};
    EXPECT_FLOAT_EQ(0.0, biospring::measure::distance(lhs, rhs2));

    rhs2 = {0.0, 0.0, 1.0};
    EXPECT_FLOAT_EQ(1.0, biospring::measure::distance(lhs, rhs2));

    rhs2 = {0.0, 0.0, 2.0};
    EXPECT_FLOAT_EQ(2.0, biospring::measure::distance(lhs, rhs2));
}

// =====================================================================================
// measure::norm
TEST(TestMeasure, norm)
{
    Particle lhs;
    lhs.setPosition({0.0, 0.0, 0.0});
    EXPECT_FLOAT_EQ(0.0, biospring::measure::norm(lhs));

    lhs.setPosition({0.0, 0.0, 1.0});
    EXPECT_FLOAT_EQ(1.0, biospring::measure::norm(lhs));

    lhs.setPosition({0.0, 0.0, 2.0});
    EXPECT_FLOAT_EQ(2.0, biospring::measure::norm(lhs));
}

TEST(TestMeasure, norm_with_pointers)
{
    Particle lhs;
    lhs.setPosition({0.0, 0.0, 0.0});
    EXPECT_FLOAT_EQ(0.0, biospring::measure::norm(&lhs));

    lhs.setPosition({0.0, 0.0, 1.0});
    EXPECT_FLOAT_EQ(1.0, biospring::measure::norm(&lhs));

    lhs.setPosition({0.0, 0.0, 2.0});
    EXPECT_FLOAT_EQ(2.0, biospring::measure::norm(&lhs));
}

// -- Main function  ----------------------------------------------------------
int main(int argc, char * argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
