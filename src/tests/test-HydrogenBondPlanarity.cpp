#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <vector>

#include "configuration/Configuration.hpp"
#include "forcefield/constants.hpp"
#include "spn/SpringNetwork.h"

using namespace biospring;

namespace
{

// Two bases hydrogen bonded through one donor and one acceptor, each hanging
// off a PLANE ATOM of its own base -- which is the real arrangement (guanine's
// N1 hangs off C2 and C6, both ring atoms) and which also puts a plane atom in
// the role of an antecedent, where a double-counted gradient would show.
//
//   0,1,2  plane of residue 1        3  its donor, antecedent 2
//   4,5,6  plane of residue 2        7  its acceptor, antecedent 4
//
// The two planes are deliberately SKEWED. Computed independently, this geometry
// gives out-of-plane distances of 1.45/1.90/1.15 A one way and 2.31/2.69/2.03 A
// the other, so at delta = 2.5 A
//
//   Fi = 0.6950   Fj = 0.4180   Pl = 0.2905
//
// with all six per-atom derivatives live and of comparable size, a
// donor-acceptor distance of 3.054 A (off the 2.90 A Morse minimum, so the
// radial term is live too) and an angular weight of 0.8677. The closest pair of
// atoms is 1.336 A, a real covalent distance.
const std::array<Vector3f, 8> SKEWED = {{
    Vector3f(0.00f, 0.00f, 0.00f), Vector3f(2.40f, 0.00f, 0.00f), Vector3f(1.00f, 2.00f, 0.00f),
    Vector3f(1.30f, 3.30f, 0.20f), Vector3f(3.60f, 6.80f, 1.45f), Vector3f(5.85f, 6.45f, 1.90f),
    Vector3f(4.40f, 8.70f, 1.15f), Vector3f(2.75f, 5.85f, 1.05f),
}};

// Both planes in z = 0 and the pair edge to edge: a perfect Watson-Crick
// arrangement, where the factor must be exactly 1.
const std::array<Vector3f, 8> COPLANAR = {{
    Vector3f(0.00f, 0.00f, 0.00f), Vector3f(2.40f, 0.00f, 0.00f), Vector3f(1.00f, 2.00f, 0.00f),
    Vector3f(1.30f, 3.30f, 0.00f), Vector3f(2.10f, 6.70f, 0.00f), Vector3f(4.50f, 6.70f, 0.00f),
    Vector3f(3.10f, 8.70f, 0.00f), Vector3f(1.60f, 6.00f, 0.00f),
}};

void place(spn::SpringNetwork & network, const std::string & name, const std::string & resname, unsigned resid,
           const Vector3f & position)
{
    spn::Particle p;
    p.setName(name);
    p.setResName(resname);
    p.setChainName("A");
    p.setResId(resid);
    p.setPosition(position);
    p.setMass(14.0f);
    p.setRadius(1.8f);
    p.setEpsilon(0.5f);
    network.addParticle(p);
}

// `frames` false leaves the network with no base plane at all, which is what
// every protein hydrogen bond looks like to this factor.
void build(spn::SpringNetwork & spn, const std::array<Vector3f, 8> & pos, bool frames = true)
{
    static const char * const names[8] = {"C8", "C6", "N3", "N1", "C6", "C4", "C2", "O2"};
    for (unsigned i = 0; i < 8; ++i)
        place(spn, names[i], i < 4 ? "G" : "C", i < 4 ? 1u : 2u, pos[i]);

    spn.getParticle(3).setDonorCapacity(1);
    spn.getParticle(3).setAntecedentIndex(2);
    spn.getParticle(7).setAcceptorCapacity(1);
    spn.getParticle(7).setAntecedentIndex(4);

    if (!frames)
        return;
    std::vector<spn::SpringNetwork::BaseFrame> f(2);
    f[0].atoms[0] = 0;
    f[0].atoms[1] = 1;
    f[0].atoms[2] = 2;
    f[1].atoms[0] = 4;
    f[1].atoms[1] = 5;
    f[1].atoms[2] = 6;
    spn.setBaseFrames(f);
    spn.setBaseFrameOfParticle({0, 0, 0, 0, 1, 1, 1, 1});
}

configuration::Configuration hbondOnly(bool planarity, double delta = 2.5)
{
    configuration::Configuration conf = configuration::defaultConfiguration();
    conf.hbond.enable = true;
    conf.hbond.cutoff = 7.0;
    conf.hbond.planarity = planarity;
    conf.hbond.planaritywidth = delta;
    return conf;
}

float energyOf(spn::SpringNetwork & spn)
{
    for (size_t i = 0; i < spn.getNumberOfParticles(); ++i)
        spn.getParticle(i).resetForce();
    spn.computeForces();
    return spn.getHydrogenBondEnergy();
}

} // namespace

TEST(HydrogenBondPlanarity, the_factor_matches_an_independent_calculation)
{
    spn::SpringNetwork spn;
    build(spn, SKEWED);
    spn.setup(hbondOnly(true));

    const spn::SpringNetwork::HydrogenBondPlanarity pl = spn.hydrogenBondPlanarity(3, 7);
    EXPECT_TRUE(pl.active);
    EXPECT_NEAR(pl.factor, 0.2905f, 1e-3f);
    // Symmetric in its two arguments: it is a property of the PAIR of planes.
    EXPECT_NEAR(spn.hydrogenBondPlanarity(7, 3).factor, pl.factor, 1e-6f);
}

// THE test for the factor: the force against the numerical gradient of the
// energy it multiplies. Eight atoms carry force -- the donor, the acceptor and
// the six plane atoms, two of which are also antecedents -- and the factor's
// own gradient runs through a centroid, a cross product and a normalisation on
// each side.
TEST(HydrogenBondPlanarity, forces_match_energy_gradient_by_finite_differences)
{
    spn::SpringNetwork spn;
    build(spn, SKEWED);
    spn.setup(hbondOnly(true));

    const float e0 = energyOf(spn);
    std::array<Vector3f, 8> analytic;
    for (size_t i = 0; i < 8; ++i)
        analytic[i] = spn.getParticle(i).getForce();

    // The bond must exist and the factor must be strictly between 0 and 1, or
    // this would pass on an inert term.
    ASSERT_LT(e0, -0.5f);
    const float factor = spn.hydrogenBondPlanarity(3, 7).factor;
    ASSERT_GT(factor, 0.05f);
    ASSERT_LT(factor, 0.95f);
    for (size_t i = 0; i < 8; ++i)
        ASSERT_GT(analytic[i].norm(), 1e-7f) << "atom " << i << " carries no force";

    const float h = 1e-3f;
    for (size_t i = 0; i < 8; ++i)
        for (int dim = 0; dim < 3; ++dim)
        {
            auto set_dim = [&](Vector3f & v, float value) {
                if (dim == 0)
                    v.setX(value);
                else if (dim == 1)
                    v.setY(value);
                else
                    v.setZ(value);
            };
            auto get_dim = [&](const Vector3f & v) {
                return dim == 0 ? v.getX() : (dim == 1 ? v.getY() : v.getZ());
            };

            const Vector3f saved = spn.getParticle(i).getPosition();
            Vector3f displaced = saved;

            set_dim(displaced, get_dim(saved) + h);
            spn.getParticle(i).setPosition(displaced);
            const float e_plus = energyOf(spn);

            set_dim(displaced, get_dim(saved) - h);
            spn.getParticle(i).setPosition(displaced);
            const float e_minus = energyOf(spn);

            spn.getParticle(i).setPosition(saved);

            const float fd = -(e_plus - e_minus) / (2.0f * h);
            const float an =
                get_dim(analytic[i]) / static_cast<float>(biospring::forcefield::GLOBAL_SPRING_FORCE_CONVERT);
            EXPECT_NEAR(fd, an, 3e-2f * std::max(1.0f, std::abs(an))) << "atom " << i << " dim " << dim;
        }
}

// The factor is normalised to be exactly 1 on a perfectly coplanar pair, which
// is what keeps an already calibrated well depth untouched: a real base pair
// pays nothing for this term, only a skewed one does.
TEST(HydrogenBondPlanarity, a_coplanar_pair_keeps_its_full_well)
{
    spn::SpringNetwork on, off;
    build(on, COPLANAR);
    build(off, COPLANAR);
    on.setup(hbondOnly(true));
    off.setup(hbondOnly(false));

    EXPECT_NEAR(on.hydrogenBondPlanarity(3, 7).factor, 1.0f, 1e-5f);
    const float e_on = energyOf(on), e_off = energyOf(off);
    ASSERT_LT(e_off, -1.0f);
    EXPECT_NEAR(e_on, e_off, 1e-4f * std::abs(e_off));
}

// And the containment that matters most: a hydrogen bond whose two ends carry
// no base plane -- every protein backbone bond -- is untouched to the BIT, not
// merely to a tolerance. Nothing about a protein can change because this
// setting was turned on.
TEST(HydrogenBondPlanarity, a_bond_without_base_planes_is_untouched_to_the_bit)
{
    spn::SpringNetwork on, off;
    build(on, SKEWED, false);
    build(off, SKEWED, false);
    on.setup(hbondOnly(true));
    off.setup(hbondOnly(false));

    EXPECT_FALSE(on.hydrogenBondPlanarity(3, 7).active);
    EXPECT_FLOAT_EQ(on.hydrogenBondPlanarity(3, 7).factor, 1.0f);

    const float e_on = energyOf(on), e_off = energyOf(off);
    ASSERT_LT(e_off, -1.0f);
    EXPECT_EQ(e_on, e_off) << "the planarity setting changed an energy it must not touch";
    for (size_t i = 0; i < 8; ++i)
    {
        EXPECT_EQ(on.getParticle(i).getForce().getX(), off.getParticle(i).getForce().getX()) << "atom " << i;
        EXPECT_EQ(on.getParticle(i).getForce().getY(), off.getParticle(i).getForce().getY()) << "atom " << i;
        EXPECT_EQ(on.getParticle(i).getForce().getZ(), off.getParticle(i).getForce().getZ()) << "atom " << i;
    }
}

// What the term is FOR: a skewed contact loses most of its energy while the
// coplanar pair above loses none.
TEST(HydrogenBondPlanarity, a_skewed_contact_loses_most_of_its_energy)
{
    spn::SpringNetwork on, off;
    build(on, SKEWED);
    build(off, SKEWED);
    on.setup(hbondOnly(true));
    off.setup(hbondOnly(false));

    const float e_on = energyOf(on), e_off = energyOf(off);
    ASSERT_LT(e_off, -1.0f);
    EXPECT_LT(e_on / e_off, 0.35f) << "a skewed pair kept " << 100.0f * e_on / e_off << " % of its well";
    EXPECT_GT(e_on / e_off, 0.25f);
}

// Narrowing delta tightens the requirement monotonically, and no width may ever
// make the factor exceed 1.
TEST(HydrogenBondPlanarity, a_narrower_width_is_always_stricter_and_never_rewards)
{
    float previous = 1.1f;
    for (const double delta : {4.0, 3.0, 2.5, 2.0, 1.5, 1.0})
    {
        spn::SpringNetwork spn;
        build(spn, SKEWED);
        spn.setup(hbondOnly(true, delta));
        const float factor = spn.hydrogenBondPlanarity(3, 7).factor;
        EXPECT_LE(factor, 1.0f) << "delta " << delta;
        EXPECT_LT(factor, previous) << "delta " << delta << " was not stricter than the width before it";
        previous = factor;
    }
}

// The angular exponent generalises max(cos theta, 0)^2 to ^n, which is the one
// part of HiRE-RNA's hydrogen bond treatment that says nothing about bases and
// so can reach a protein. Its derivative is n cos^(n-1), a second place to get
// an exponent wrong, so the gradient is checked at every exponent the setting
// is likely to carry -- including the default, which must not have moved.
class AngularExponent : public ::testing::TestWithParam<int>
{
};

TEST_P(AngularExponent, forces_match_energy_gradient_by_finite_differences)
{
    spn::SpringNetwork spn;
    build(spn, SKEWED);
    configuration::Configuration conf = hbondOnly(true);
    conf.hbond.angularexponent = GetParam();
    spn.setup(conf);

    const float e0 = energyOf(spn);
    std::array<Vector3f, 8> analytic;
    for (size_t i = 0; i < 8; ++i)
        analytic[i] = spn.getParticle(i).getForce();
    ASSERT_LT(e0, -0.05f) << "exponent " << GetParam() << " left no bond to differentiate";

    const float h = 1e-3f;
    for (size_t i = 0; i < 8; ++i)
        for (int dim = 0; dim < 3; ++dim)
        {
            auto set_dim = [&](Vector3f & v, float value) {
                if (dim == 0)
                    v.setX(value);
                else if (dim == 1)
                    v.setY(value);
                else
                    v.setZ(value);
            };
            auto get_dim = [&](const Vector3f & v) {
                return dim == 0 ? v.getX() : (dim == 1 ? v.getY() : v.getZ());
            };
            const Vector3f saved = spn.getParticle(i).getPosition();
            Vector3f displaced = saved;
            set_dim(displaced, get_dim(saved) + h);
            spn.getParticle(i).setPosition(displaced);
            const float e_plus = energyOf(spn);
            set_dim(displaced, get_dim(saved) - h);
            spn.getParticle(i).setPosition(displaced);
            const float e_minus = energyOf(spn);
            spn.getParticle(i).setPosition(saved);
            const float fd = -(e_plus - e_minus) / (2.0f * h);
            const float an =
                get_dim(analytic[i]) / static_cast<float>(biospring::forcefield::GLOBAL_SPRING_FORCE_CONVERT);
            EXPECT_NEAR(fd, an, 3e-2f * std::max(1.0f, std::abs(an)))
                << "exponent " << GetParam() << " atom " << i << " dim " << dim;
        }
}

INSTANTIATE_TEST_SUITE_P(TwoFourAndSix, AngularExponent, ::testing::Values(1, 2, 3, 4, 6, 8));

// A higher exponent is strictly stricter off-axis and leaves a perfectly
// aligned bond alone, and the default must still be two.
TEST(HydrogenBondPlanarity, a_higher_angular_exponent_is_stricter_off_axis)
{
    EXPECT_FLOAT_EQ(forcefield::hydrogen_bond_angular_factor(1.0f, 2),
                    forcefield::hydrogen_bond_angular_factor(1.0f, 6));
    EXPECT_LT(forcefield::hydrogen_bond_angular_factor(0.8f, 6),
              forcefield::hydrogen_bond_angular_factor(0.8f, 2));
    // Behind the donor is not a weak bond, it is not one at all, at any
    // exponent -- including the odd ones, where the clamp is what prevents a
    // NEGATIVE weight rather than a revived one.
    for (const int n : {1, 2, 3, 6})
    {
        EXPECT_FLOAT_EQ(forcefield::hydrogen_bond_angular_factor(-0.5f, n), 0.0f) << "exponent " << n;
        EXPECT_FLOAT_EQ(forcefield::hydrogen_bond_angular_derivative(-0.5f, n), 0.0f) << "exponent " << n;
    }
    EXPECT_FLOAT_EQ(forcefield::hydrogen_bond_angular_factor(0.7f), forcefield::hydrogen_bond_angular_factor(0.7f, 2));
}
