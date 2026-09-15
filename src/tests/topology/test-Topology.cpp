#include <gtest/gtest.h>

#include <vector>

#include "IO/io.h"
#include "SpringNetwork.h"
#include "forcefield/constants.hpp"
#include "measure.hpp"
#include "topology/Spring.hpp"
#include "topology/Topology.hpp"

using namespace biospring;

// ====================================================================================
//
// Helper functions.
//
// ====================================================================================

static auto generate_random_particles(size_t n);
// Returns a random 3D-position.
// Each coordinate is drawn from a uniform distribution on the interval [-100, 100].
static std::array<double, 3> random_position()
{
    std::array<double, 3> pos;
    for (double & x : pos)
        x = 200.0 * (double)rand() / RAND_MAX - 100.0;
    return pos;
}

// Generates a collection of particles with random positions.
static auto generate_random_particles(size_t n)
{
    std::vector<topology::Particle> particles(n);
    for (topology::Particle & p : particles)
        p.set_position(random_position());
    return particles;
}

// ====================================================================================
//
// Test for I/O functions.
//
// ====================================================================================
TEST(Topology, readTopology)
{
    std::string path = "../data/model.pdb";
    topology::Topology top = io::readTopology(path);
    EXPECT_EQ(top.number_of_particles(), 2814);
    EXPECT_EQ(top.number_of_springs(), 0);
    EXPECT_EQ(top.particles().size(), top.springs().particles().size());
}

// ====================================================================================
//
// Merge topologies.
//
// ====================================================================================

TEST(Topology, merge_static)
{
    topology::Topology top1, top2;
    top1.add_particles(generate_random_particles(100));
    top2.add_particles(generate_random_particles(100));

    top1.add_springs_from_cutoff(10.0);
    top2.add_springs_from_cutoff(10.0);

    ASSERT_GE(top1.number_of_springs(), 1);
    ASSERT_GE(top2.number_of_springs(), 1);

    size_t expected_number_of_particles = top1.number_of_particles() + top2.number_of_particles();
    size_t expected_number_of_springs = top1.number_of_springs() + top2.number_of_springs();

    topology::Topology merged = topology::Topology::merge(top1, top2);
    EXPECT_EQ(merged.number_of_particles(), expected_number_of_particles);
    EXPECT_EQ(merged.number_of_springs(), expected_number_of_springs);
}

TEST(Topology, merge)
{
    topology::Topology top1, top2;
    top1.add_particles(generate_random_particles(100));
    top2.add_particles(generate_random_particles(100));

    top1.add_springs_from_cutoff(10.0);
    top2.add_springs_from_cutoff(10.0);

    ASSERT_GE(top1.number_of_springs(), 1);
    ASSERT_GE(top2.number_of_springs(), 1);

    size_t expected_number_of_particles = top1.number_of_particles() + top2.number_of_particles();
    size_t expected_number_of_springs = top1.number_of_springs() + top2.number_of_springs();

    topology::Topology merged = top1.merge(top2);
    EXPECT_EQ(merged.number_of_particles(), expected_number_of_particles);
    EXPECT_EQ(merged.number_of_springs(), expected_number_of_springs);
}

// ====================================================================================
//
// Add particles
//
// ====================================================================================

TEST(Topology, add_particle_one_by_one)
{
    topology::Topology top;
    topology::Particle p1, p2;

    top.particles().push_back(p1);
    EXPECT_EQ(top.number_of_particles(), 1);

    top.add_particle(p2);
    EXPECT_EQ(top.number_of_particles(), 2);
}

TEST(Topology, add_particle_container)
{
    topology::Topology top;
    std::vector<topology::Particle> particles(2);

    top.add_particles(particles);
    EXPECT_EQ(top.number_of_particles(), 2);
}

TEST(Topology, add_particle_initializer_list)
{
    topology::Topology top;
    topology::Particle p1, p2;

    top.add_particles({p1, p2});
    EXPECT_EQ(top.number_of_particles(), 2);
}

TEST(Topology, add_spring)
{
    topology::Topology top;
    topology::Particle p1, p2;
    top.add_particles({p1, p2});

    top.add_spring(top.get_particle(0), top.get_particle(1), 1.0, 1.0);

    EXPECT_EQ(top.number_of_springs(), 1);
    EXPECT_FLOAT_EQ(top.get_spring(0).equilibrium(), 1.0);
    EXPECT_FLOAT_EQ(top.get_spring(0).stiffness(), 1.0);
}

// Checks that springs are added between all particles within a cutoff distance.
TEST(Topology, add_springs_from_cutoff)
{
    topology::Topology top;
    top.add_particles(generate_random_particles(1000));
    top.add_springs_from_cutoff(10.0);
    ASSERT_GE(top.number_of_springs(), 1);

    for (size_t i = 0; i < top.number_of_particles(); ++i)
    {
        for (size_t j = i + 1; j < top.number_of_particles(); ++j)
        {
            const topology::Particle & p1 = top.get_particle(i);
            const topology::Particle & p2 = top.get_particle(j);
            double distance = measure::distance(p1, p2);

            if (distance < 10.0)
                EXPECT_TRUE(top.has_spring_between(p1, p2));

            // else
            //     EXPECT_FALSE(top.has_spring_between(p1, p2));
        }
    }
}

TEST(Topology, to_spring_network)
{
    topology::Topology top;
    top.add_particles(generate_random_particles(1000));
    top.add_springs_from_cutoff(10.0);
    ASSERT_GE(top.number_of_springs(), 1);

    spn::SpringNetwork spn;
    top.to_spring_network(spn);

    EXPECT_EQ(spn.getNumberOfParticles(), top.number_of_particles());
    EXPECT_EQ(spn.getNumberOfSprings(), top.number_of_springs());
}

// -- Main function  ----------------------------------------------------------
// The hydrogen bond's angular weight makes the antecedent a third body, and
// its gradient is easy to transcribe with a sign flipped -- which is exactly
// what happened once and what nothing but this check caught (the energy is
// unaffected by it, and so is every pairing count).
TEST(Topology, hydrogen_bond_forces_match_energy_gradient_by_finite_differences)
{
    topology::Topology top;

    // Antecedent, donor, acceptor. Deliberately off-axis and off-equilibrium
    // so the radial and the angular term are both doing real work, and the
    // angle is nowhere near 0 or 90 degrees where one of them vanishes.
    const std::array<std::array<double, 3>, 3> pos = {{{0.0, 0.0, 0.0},
                                                       {1.35, 0.22, -0.11},
                                                       {3.4, 1.9, 0.35}}};
    for (size_t i = 0; i < pos.size(); ++i)
    {
        topology::ParticleProperties p;
        p.set_position(Vector3f(static_cast<float>(pos[i][0]), static_cast<float>(pos[i][1]),
                                static_cast<float>(pos[i][2])));
        p.set_mass(12.0f);
        // Distinct residues: the pairing rule skips a candidate inside the
        // donor's own residue (a backbone N and O sit at covalent distance
        // and would otherwise always win).
        p.set_residue_id(static_cast<int>(i) + 1);
        top.add_particle(topology::Particle(p));
    }

    spn::SpringNetwork spn;
    top.to_spring_network(spn);
    spn.getParticle(1).setDonorCapacity(1);
    spn.getParticle(1).setAntecedentIndex(0);
    spn.getParticle(2).setAcceptorCapacity(1);

    configuration::Configuration conf = configuration::defaultConfiguration();
    conf.hbond.enable = true;
    conf.hbond.cutoff = 7.0;
    spn.setup(conf);

    auto energy = [&]() {
        for (size_t i = 0; i < spn.getNumberOfParticles(); ++i)
            spn.getParticle(i).resetForce();
        spn.computeForces();
        return spn.getHydrogenBondEnergy();
    };

    const float e0 = energy();
    std::array<Vector3f, 3> analytic;
    for (size_t i = 0; i < 3; ++i)
        analytic[i] = spn.getParticle(i).getForce();

    // The bond must actually exist, and the angular weight must be strictly
    // between 0 and 1 -- otherwise this test would pass on a two-body force.
    ASSERT_LT(e0, -1.0f);
    ASSERT_GT(analytic[0].norm(), 1e-6f) << "the antecedent carries no force: the angular term is inert";

    const float h = 1e-2f;
    for (size_t i = 0; i < 3; ++i)
    {
        for (int dim = 0; dim < 3; ++dim)
        {
            const Vector3f saved = spn.getParticle(i).getPosition();
            Vector3f displaced = saved;
            auto set_dim = [&](Vector3f & v, float value) {
                if (dim == 0) v.setX(value);
                else if (dim == 1) v.setY(value);
                else v.setZ(value);
            };
            auto get_dim = [&](const Vector3f & v) { return dim == 0 ? v.getX() : (dim == 1 ? v.getY() : v.getZ()); };

            set_dim(displaced, get_dim(saved) + h);
            spn.getParticle(i).setPosition(displaced);
            const float e_plus = energy();

            set_dim(displaced, get_dim(saved) - h);
            spn.getParticle(i).setPosition(displaced);
            const float e_minus = energy();

            spn.getParticle(i).setPosition(saved);

            const float fd = -(e_plus - e_minus) / (2.0f * h);
            const float an = get_dim(analytic[i]) /
                             static_cast<float>(biospring::forcefield::GLOBAL_SPRING_FORCE_CONVERT);
            EXPECT_NEAR(fd, an, 3e-2f * std::max(1.0f, std::abs(an)))
                << "particle " << i << " dim " << dim;
        }
    }
}

TEST(Topology, hydrogen_bond_bisector_forces_match_energy_gradient)
{
    topology::Topology top;

    // TWO antecedents, donor, acceptor. With two, the hydrogen's direction is
    // the sum of the two away-from-antecedent unit vectors, so the angular
    // term's gradient reaches four atoms instead of three, and each new term
    // carries a projector that is easy to transcribe with a sign or a
    // normalisation wrong. Nothing but a finite-difference check would see
    // it: the energy, the pairing and the bond count are all unaffected.
    //
    // Placed so neither antecedent lies along the donor-acceptor line and the
    // two are not symmetric about it, which would let a wrong split cancel.
    const std::array<std::array<double, 3>, 4> pos = {{{0.0, 0.0, 0.0},
                                                       {0.7, -1.15, 0.31},
                                                       {1.35, 0.22, -0.11},
                                                       {3.4, 1.9, 0.35}}};
    for (size_t i = 0; i < pos.size(); ++i)
    {
        topology::ParticleProperties p;
        p.set_position(Vector3f(static_cast<float>(pos[i][0]), static_cast<float>(pos[i][1]),
                                static_cast<float>(pos[i][2])));
        p.set_mass(12.0f);
        // Distinct residues: the pairing rule skips a candidate inside the
        // donor's own residue (a backbone N and O sit at covalent distance
        // and would otherwise always win).
        p.set_residue_id(static_cast<int>(i) + 1);
        top.add_particle(topology::Particle(p));
    }

    spn::SpringNetwork spn;
    top.to_spring_network(spn);
    spn.getParticle(2).setDonorCapacity(1);
    spn.getParticle(2).setAntecedentIndex(0);
    spn.getParticle(2).setAntecedentIndex2(1);
    spn.getParticle(3).setAcceptorCapacity(1);

    configuration::Configuration conf = configuration::defaultConfiguration();
    conf.hbond.enable = true;
    conf.hbond.cutoff = 7.0;
    spn.setup(conf);

    auto energy = [&]() {
        for (size_t i = 0; i < spn.getNumberOfParticles(); ++i)
            spn.getParticle(i).resetForce();
        spn.computeForces();
        return spn.getHydrogenBondEnergy();
    };

    const float e0 = energy();
    std::array<Vector3f, 4> analytic;
    for (size_t i = 0; i < 4; ++i)
        analytic[i] = spn.getParticle(i).getForce();

    // The bond must actually exist, and the angular weight must be strictly
    // between 0 and 1 -- otherwise this test would pass on a two-body force.
    ASSERT_LT(e0, -1.0f);
    ASSERT_GT(analytic[0].norm(), 1e-6f) << "the first antecedent carries no force";
    ASSERT_GT(analytic[1].norm(), 1e-6f) << "the second carries none: the bisector is not in play";

    const float h = 1e-2f;
    for (size_t i = 0; i < 4; ++i)
    {
        for (int dim = 0; dim < 3; ++dim)
        {
            const Vector3f saved = spn.getParticle(i).getPosition();
            Vector3f displaced = saved;
            auto set_dim = [&](Vector3f & v, float value) {
                if (dim == 0) v.setX(value);
                else if (dim == 1) v.setY(value);
                else v.setZ(value);
            };
            auto get_dim = [&](const Vector3f & v) { return dim == 0 ? v.getX() : (dim == 1 ? v.getY() : v.getZ()); };

            set_dim(displaced, get_dim(saved) + h);
            spn.getParticle(i).setPosition(displaced);
            const float e_plus = energy();

            set_dim(displaced, get_dim(saved) - h);
            spn.getParticle(i).setPosition(displaced);
            const float e_minus = energy();

            spn.getParticle(i).setPosition(saved);

            const float fd = -(e_plus - e_minus) / (2.0f * h);
            const float an = get_dim(analytic[i]) /
                             static_cast<float>(biospring::forcefield::GLOBAL_SPRING_FORCE_CONVERT);
            EXPECT_NEAR(fd, an, 3e-2f * std::max(1.0f, std::abs(an)))
                << "particle " << i << " dim " << dim;
        }
    }
}

TEST(Topology, hydrogen_bond_acceptor_angle_forces_match_energy_gradient)
{
    topology::Topology top;

    // FIVE atoms: two antecedents on the donor, the donor, the acceptor, and
    // the acceptor's own antecedent. The angular weight is now a PRODUCT of
    // two factors, one per side, so each side's gradient carries the other
    // side's weight as a scale -- drop that and the forces are wrong by a
    // factor that varies bond by bond while every energy stays right.
    //
    // Original note, still true of the donor side: with two antecedents the
    // hydrogen's direction is
    // the sum of the two away-from-antecedent unit vectors, so the angular
    // term's gradient reaches four atoms instead of three, and each new term
    // carries a projector that is easy to transcribe with a sign or a
    // normalisation wrong. Nothing but a finite-difference check would see
    // it: the energy, the pairing and the bond count are all unaffected.
    //
    // Placed so neither antecedent lies along the donor-acceptor line and the
    // two are not symmetric about it, which would let a wrong split cancel.
    const std::array<std::array<double, 3>, 5> pos = {{{0.0, 0.0, 0.0},
                                                       {0.7, -1.15, 0.31},
                                                       {1.35, 0.22, -0.11},
                                                       {3.4, 1.9, 0.35},
                                                       {4.55, 2.7, -0.4}}};
    for (size_t i = 0; i < pos.size(); ++i)
    {
        topology::ParticleProperties p;
        p.set_position(Vector3f(static_cast<float>(pos[i][0]), static_cast<float>(pos[i][1]),
                                static_cast<float>(pos[i][2])));
        p.set_mass(12.0f);
        // Distinct residues: the pairing rule skips a candidate inside the
        // donor's own residue (a backbone N and O sit at covalent distance
        // and would otherwise always win).
        p.set_residue_id(static_cast<int>(i) + 1);
        top.add_particle(topology::Particle(p));
    }

    spn::SpringNetwork spn;
    top.to_spring_network(spn);
    spn.getParticle(2).setDonorCapacity(1);
    spn.getParticle(2).setAntecedentIndex(0);
    spn.getParticle(2).setAntecedentIndex2(1);
    spn.getParticle(3).setAcceptorCapacity(1);
    spn.getParticle(3).setAntecedentIndex(4);

    configuration::Configuration conf = configuration::defaultConfiguration();
    conf.hbond.enable = true;
    conf.hbond.cutoff = 7.0;
    spn.setup(conf);

    auto energy = [&]() {
        for (size_t i = 0; i < spn.getNumberOfParticles(); ++i)
            spn.getParticle(i).resetForce();
        spn.computeForces();
        return spn.getHydrogenBondEnergy();
    };

    const float e0 = energy();
    std::array<Vector3f, 5> analytic;
    for (size_t i = 0; i < 5; ++i)
        analytic[i] = spn.getParticle(i).getForce();

    // The bond must actually exist, and the angular weight must be strictly
    // between 0 and 1 -- otherwise this test would pass on a two-body force.
    ASSERT_LT(e0, -1.0f);
    ASSERT_GT(analytic[0].norm(), 1e-6f) << "the first antecedent carries no force";
    ASSERT_GT(analytic[1].norm(), 1e-6f) << "the second carries none: the bisector is not in play";
    ASSERT_GT(analytic[4].norm(), 1e-6f) << "the acceptor's antecedent carries none: its angular term is inert";

    const float h = 1e-2f;
    for (size_t i = 0; i < 5; ++i)
    {
        for (int dim = 0; dim < 3; ++dim)
        {
            const Vector3f saved = spn.getParticle(i).getPosition();
            Vector3f displaced = saved;
            auto set_dim = [&](Vector3f & v, float value) {
                if (dim == 0) v.setX(value);
                else if (dim == 1) v.setY(value);
                else v.setZ(value);
            };
            auto get_dim = [&](const Vector3f & v) { return dim == 0 ? v.getX() : (dim == 1 ? v.getY() : v.getZ()); };

            set_dim(displaced, get_dim(saved) + h);
            spn.getParticle(i).setPosition(displaced);
            const float e_plus = energy();

            set_dim(displaced, get_dim(saved) - h);
            spn.getParticle(i).setPosition(displaced);
            const float e_minus = energy();

            spn.getParticle(i).setPosition(saved);

            const float fd = -(e_plus - e_minus) / (2.0f * h);
            const float an = get_dim(analytic[i]) /
                             static_cast<float>(biospring::forcefield::GLOBAL_SPRING_FORCE_CONVERT);
            EXPECT_NEAR(fd, an, 3e-2f * std::max(1.0f, std::abs(an)))
                << "particle " << i << " dim " << dim;
        }
    }
}

int main(int argc, char * argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
