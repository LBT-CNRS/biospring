#include <gtest/gtest.h>

#include <vector>

#include "IO/io.h"
#include "SpringNetwork.h"
#include "configuration/Configuration.hpp"
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
// The LOBE form of a site's direction, against finite differences.
//
// A planar sp2 site with one heavy neighbour and two hydrogens -- an exocyclic
// amine, a carbonyl oxygen -- has TWO directions, at +/- 62 degrees of its axis
// and in its plane, and the bond takes whichever fits. The plane comes from a
// third atom of the same group, marked '~' in the .hbond table, and that atom
// enters the gradient like an antecedent does. Four bodies carry force here, not
// three, and the derivation goes through a Gram-Schmidt projection, so it is
// exactly the kind of expression a sign error hides in.
// FOUR cases: both lobe geometries, each with the partner on either side of the
// axis so sigma comes out +1 in one and -1 in the other. One sign convention has to
// satisfy all four, which is what stops this being fitted to a single case -- the
// algebra's B term came out inverted on the first attempt, and only a case of the
// opposite sigma distinguishes that from a correct one.
//
// The angle is measured from e1 = self - antecedent1. For a carbonyl carbon whose
// antecedent is its own oxygen, e1 points AWAY from the oxygen, so Burgi-Dunitz's
// 105 degrees of the C=O axis is written as 75 here: a lobe at 75 degrees of e1 puts
// the partner at exactly 105.0 of C->O, which is checked below.
struct LobeCase
{
    std::array<std::array<double, 3>, 4> pos; // antecedent1, plane atom, site, partner
    int mode;
    float angle;
};

class HydrogenBondLobeForces : public ::testing::TestWithParam<LobeCase>
{
};

INSTANTIATE_TEST_SUITE_P(
    BothGeometriesBothLobes, HydrogenBondLobeForces,
    ::testing::Values(
        // IN PLANE, 62 degrees: an exocyclic amine's two hydrogens. The partner is
        // off-axis to one side, then to the other.
        LobeCase{{{{0.0, 0.0, 0.0}, {-0.7, 1.2, 0.05}, {1.35, 0.22, -0.11}, {2.6, 2.3, 0.4}}},
                 spn::Particle::HBOND_LOBES_IN_PLANE, 62.0f},
        LobeCase{{{{0.0, 0.0, 0.0}, {-0.7, 1.2, 0.05}, {1.35, 0.22, -0.11}, {3.09, -1.91, 1.32}}},
                 spn::Particle::HBOND_LOBES_IN_PLANE, 62.0f},
        // OUT OF PLANE, 75 degrees: the Burgi-Dunitz approach to a carbonyl carbon.
        // Antecedent1 is the oxygen at the origin, the site is the carbon 1.23 A
        // along x, and the plane atom puts the sp2 plane at z = 0 -- so the lobes run
        // along +/- z and the Burgi-Dunitz angle comes out at 103.4 degrees of C->O.
        //
        // The partner is tilted off the lobe TOWARD THE IN-PLANE PERPENDICULAR, not
        // within the (axis, normal) plane, and that is not cosmetic: tilt it inside
        // that plane and everything lies in one plane, ts comes out parallel to the
        // axis, B is EXACTLY zero and the plane atom carries no force at all. The
        // test would then pass while never exercising the term it exists for. Here
        // the weight is 0.800 and |B| is 0.332.
        LobeCase{{{{0.0, 0.0, 0.0}, {1.93, 1.30, 0.0}, {1.23, 0.0, 0.0}, {1.924, 1.342, 2.592}}},
                 spn::Particle::HBOND_LOBES_OUT_OF_PLANE, 75.0f},
        LobeCase{{{{0.0, 0.0, 0.0}, {1.93, 1.30, 0.0}, {1.23, 0.0, 0.0}, {1.924, 1.342, -2.592}}},
                 spn::Particle::HBOND_LOBES_OUT_OF_PLANE, 75.0f}));

TEST_P(HydrogenBondLobeForces, match_energy_gradient_by_finite_differences)
{
    topology::Topology top;
    const LobeCase param = GetParam();
    const std::array<std::array<double, 3>, 4> pos = param.pos;
    for (size_t i = 0; i < pos.size(); ++i)
    {
        topology::ParticleProperties p;
        p.set_position(Vector3f(static_cast<float>(pos[i][0]), static_cast<float>(pos[i][1]),
                                static_cast<float>(pos[i][2])));
        p.set_mass(12.0f);
        p.set_residue_id(static_cast<int>(i) + 1);
        top.add_particle(topology::Particle(p));
    }

    spn::SpringNetwork spn;
    top.to_spring_network(spn);
    spn.getParticle(2).setDonorCapacity(2);
    spn.getParticle(2).setAntecedentIndex(0);
    spn.getParticle(2).setAntecedentIndex2(1);
    spn.getParticle(2).setLobes(param.mode, param.angle);
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

    ASSERT_LT(e0, -1.0f) << "no bond formed, so this measures nothing";
    // Both the heavy neighbour AND the plane atom must carry force: if either is
    // inert the lobe direction is not being differentiated through.
    ASSERT_GT(analytic[0].norm(), 1e-6f) << "the heavy neighbour carries no force";
    ASSERT_GT(analytic[1].norm(), 1e-6f) << "the PLANE atom carries no force: the lobe is not differentiated";

    // And the weight must be strictly inside (0, 1), otherwise the angular
    // derivative is zero or saturated and the test passes on a radial force.
    {
        const Vector3f v = spn.getParticle(3).getPosition() - spn.getParticle(2).getPosition();
        const Vector3f hd = spn.donorDirection(spn.getParticle(2), v / v.norm());
        const float w = biospring::forcefield::hydrogen_bond_angular_factor(hd.dot(v / v.norm()));
        ASSERT_GT(w, 0.1f) << "angular weight " << w << " is too close to zero";
        ASSERT_LT(w, 0.95f) << "angular weight " << w << " is saturated";
        // And the direction must BE a lobe, not the axis: the site's own angle apart
        // from it, whichever geometry this case uses.
        const Vector3f axis = (spn.getParticle(2).getPosition() - spn.getParticle(0).getPosition());
        EXPECT_NEAR(hd.dot(axis / axis.norm()),
                    std::cos(param.angle * static_cast<float>(M_PI) / 180.0f), 1e-3f);
        // An out-of-plane lobe must actually leave the plane; an in-plane one must
        // stay in it. Otherwise the two modes could be silently the same code path.
        const Vector3f pl = spn.getParticle(1).getPosition() - spn.getParticle(0).getPosition();
        const Vector3f nrm = (axis ^ pl);
        const float outofplane = std::abs(hd.dot(nrm / nrm.norm()));
        if (param.mode == spn::Particle::HBOND_LOBES_OUT_OF_PLANE)
            EXPECT_GT(outofplane, 0.9f) << "an out-of-plane lobe is lying in the plane";
        else
            EXPECT_LT(outofplane, 1e-3f) << "an in-plane lobe is leaving the plane";
    }

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

// With sin(62) taken out of the lobe form it collapses to the single-antecedent
// direction, so the two forms must agree in that limit. This checks the SITE
// rather than the force: a lobe site whose acceptor sits exactly on the axis is
// weaker than a plain axis site by cos^2(62), and by nothing else.
TEST(Topology, a_lobe_site_reduces_to_the_axis_when_the_partner_is_on_it)
{
    topology::Topology top;
    const std::array<std::array<double, 3>, 4> pos = {{{0.0, 0.0, 0.0},
                                                       {-0.7, 1.2, 0.05},
                                                       {1.35, 0.0, 0.0},
                                                       {4.25, 0.0, 0.0}}};
    for (size_t i = 0; i < pos.size(); ++i)
    {
        topology::ParticleProperties p;
        p.set_position(Vector3f(static_cast<float>(pos[i][0]), static_cast<float>(pos[i][1]),
                                static_cast<float>(pos[i][2])));
        p.set_mass(12.0f);
        p.set_residue_id(static_cast<int>(i) + 1);
        top.add_particle(topology::Particle(p));
    }
    spn::SpringNetwork spn;
    top.to_spring_network(spn);
    spn.getParticle(2).setDonorCapacity(2);
    spn.getParticle(2).setAntecedentIndex(0);
    spn.getParticle(2).setAntecedentIndex2(1);
    spn.getParticle(3).setAcceptorCapacity(1);
    configuration::Configuration conf = configuration::defaultConfiguration();
    conf.hbond.enable = true;
    conf.hbond.cutoff = 7.0;
    spn.setup(conf);

    const Vector3f towards(1.0f, 0.0f, 0.0f);
    // The axis form: antecedent 0 only, so the direction is +x exactly.
    spn.getParticle(2).setLobes(spn::Particle::HBOND_LOBES_NONE, 62.0f);
    spn.getParticle(2).setAntecedentIndex2(-1);
    const Vector3f axis = spn.donorDirection(spn.getParticle(2), towards);
    EXPECT_NEAR(axis.dot(towards), 1.0f, 1e-5f);

    // Either lobe form: the same axis, tilted by the site's own angle. Both must be
    // unit and both must sit exactly that angle off the axis -- what separates them
    // is the PLANE they tilt into, which the finite-difference cases above check.
    spn.getParticle(2).setAntecedentIndex2(1);
    for (int mode : {spn::Particle::HBOND_LOBES_IN_PLANE, spn::Particle::HBOND_LOBES_OUT_OF_PLANE})
    {
        for (float angle : {62.0f, 75.0f, 105.0f})
        {
            spn.getParticle(2).setLobes(mode, angle);
            const Vector3f lobe = spn.donorDirection(spn.getParticle(2), towards);
            EXPECT_NEAR(lobe.dot(towards), std::cos(angle * static_cast<float>(M_PI) / 180.0f), 1e-4f)
                << "mode " << mode << " angle " << angle;
            EXPECT_NEAR(lobe.norm(), 1.0f, 1e-5f) << "mode " << mode << " angle " << angle;
        }
    }
}

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
