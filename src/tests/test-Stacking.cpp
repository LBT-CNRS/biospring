#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>

#include "IO/BasePlaneRuleReader.h"
#include "configuration/Configuration.hpp"
#include "forcefield/constants.hpp"
#include "spn/SpringNetwork.h"

using namespace biospring;

namespace
{

// A .stack table that lives as long as the test does.
class Table
{
  public:
    explicit Table(const std::string & body)
    {
        _path = std::string(std::tmpnam(nullptr)) + ".stack";
        std::ofstream(_path) << body;
    }
    ~Table() { std::remove(_path.c_str()); }
    const std::string & path() const { return _path; }

  private:
    std::string _path;
};

unsigned place(spn::SpringNetwork & network, const std::string & name, const std::string & resname, unsigned resid,
               const Vector3f & position)
{
    spn::Particle p;
    p.setName(name);
    p.setResName(resname);
    p.setChainName("A");
    p.setResId(resid);
    p.setPosition(position);
    p.setMass(12.0f);
    p.setRadius(1.8f);
    p.setEpsilon(0.5f);
    network.addParticle(p);
    return network.getNumberOfParticles() - 1;
}

// Two base planes, three atoms each, as two residues of one chain. The caller
// gives the six positions; the frames are installed directly rather than
// through a table, so a geometry test does not also depend on name resolution.
void build(spn::SpringNetwork & network, const std::array<Vector3f, 6> & pos)
{
    static const char * const names[6] = {"C6", "C4", "C2", "C6", "C4", "C2"};
    for (unsigned i = 0; i < 6; ++i)
        place(network, names[i], "C", i < 3 ? 1u : 2u, pos[i]);
    std::vector<spn::SpringNetwork::BaseFrame> frames(2);
    frames[0].atoms[0] = 0;
    frames[0].atoms[1] = 1;
    frames[0].atoms[2] = 2;
    frames[1].atoms[0] = 3;
    frames[1].atoms[1] = 4;
    frames[1].atoms[2] = 5;
    network.setBaseFrames(frames);
}

configuration::Configuration stackingOnly(double epsilon = 8.0, double equilibrium = 3.92, double width = 0.85)
{
    configuration::Configuration conf = configuration::defaultConfiguration();
    conf.stacking.enable = true;
    conf.stacking.cutoff = 7.0;
    conf.stacking.epsilon = epsilon;
    conf.stacking.equilibrium = equilibrium;
    conf.stacking.width = width;
    return conf;
}

// Deliberately generic: the two planes sit off the equilibrium distance, their
// normals are neither parallel nor perpendicular to the axis joining them, and
// the two alignment factors differ. Computed independently, this geometry gives
//
//   d = 3.3187 A      G = 0.6535   dG = 0.9246
//   n_i.n_j = 0.8323  P = 0.6927
//   u_i = 0.8236      A_i = 0.8965  dA_i = 1.0597
//   u_j = 0.6447      A_j = 0.6586  dA_j = 1.5069
//   E = -2.1384 kJ.mol-1
//
// so every one of the four factors AND every one of their derivatives is doing
// real work. A geometry at the bottom of the well, or with parallel normals, or
// face on, would let a wrong term hide behind a vanishing derivative.
const std::array<Vector3f, 6> GENERIC = {{Vector3f(0.00f, 0.00f, 0.00f), Vector3f(2.40f, 0.00f, 0.00f),
                                          Vector3f(1.00f, 2.00f, 0.00f), Vector3f(1.55f, 1.05f, 2.60f),
                                          Vector3f(3.60f, 1.60f, 3.35f), Vector3f(2.60f, 2.95f, 2.25f)}};

} // namespace

// THE test for a new force: the analytic force against the numerical gradient
// of the energy it claims to be the gradient of. Six atoms carry force here,
// because a frame's normal and its centre are both built from all three of its
// atoms, and the normal's gradient goes through a cross product and a
// normalisation -- two places where a sign or a projector is easy to transcribe
// wrong without changing the energy at all.
TEST(Stacking, forces_match_energy_gradient_by_finite_differences)
{
    spn::SpringNetwork spn;
    build(spn, GENERIC);
    spn.setup(stackingOnly());

    auto energy = [&]() {
        for (size_t i = 0; i < spn.getNumberOfParticles(); ++i)
            spn.getParticle(i).resetForce();
        spn.computeForces();
        return spn.getStackingEnergy();
    };

    const float e0 = energy();
    std::array<Vector3f, 6> analytic;
    for (size_t i = 0; i < 6; ++i)
        analytic[i] = spn.getParticle(i).getForce();

    // The term must be well off both ends of its range, or this would pass on
    // something inert.
    ASSERT_NEAR(e0, -2.1384f, 2e-3f);
    for (size_t i = 0; i < 6; ++i)
        ASSERT_GT(analytic[i].norm(), 1e-6f) << "atom " << i << " carries no force";

    const float h = 1e-3f;
    for (size_t i = 0; i < 6; ++i)
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
            const float e_plus = energy();

            set_dim(displaced, get_dim(saved) - h);
            spn.getParticle(i).setPosition(displaced);
            const float e_minus = energy();

            spn.getParticle(i).setPosition(saved);

            const float fd = -(e_plus - e_minus) / (2.0f * h);
            const float an =
                get_dim(analytic[i]) / static_cast<float>(biospring::forcefield::GLOBAL_SPRING_FORCE_CONVERT);
            EXPECT_NEAR(fd, an, 2e-2f * std::max(1.0f, std::abs(an))) << "atom " << i << " dim " << dim;
        }
}

// Translation invariance, which the force distribution gives exactly rather
// than by balancing: each frame's three normal-forces sum to zero by
// construction, and the centre's share is split in three equal parts of
// opposite sign.
TEST(Stacking, the_six_forces_sum_to_zero)
{
    spn::SpringNetwork spn;
    build(spn, GENERIC);
    spn.setup(stackingOnly());
    spn.computeForces();

    Vector3f total;
    float largest = 0.0f;
    for (size_t i = 0; i < 6; ++i)
    {
        total = total + spn.getParticle(i).getForce();
        largest = std::max(largest, spn.getParticle(i).getForce().norm());
    }
    ASSERT_GT(largest, 1e-6f);
    EXPECT_LT(total.norm(), 1e-5f * largest);
}

// A base plane has no preferred face for stacking, so naming its three atoms in
// another order -- which flips the normal -- must change nothing. Both of the
// term's orientation factors are even in the normal for exactly this reason,
// and this is what checks that neither of them lost its square.
TEST(Stacking, the_order_of_the_three_plane_atoms_changes_nothing)
{
    float reference = 0.0f;
    // Every permutation of the second frame's three atoms. Three of them flip
    // the normal's sign and three do not; the centroid is the same for all six.
    static const unsigned perms[6][3] = {{3, 4, 5}, {3, 5, 4}, {4, 3, 5}, {4, 5, 3}, {5, 3, 4}, {5, 4, 3}};
    for (unsigned p = 0; p < 6; ++p)
    {
        spn::SpringNetwork spn;
        build(spn, GENERIC);
        std::vector<spn::SpringNetwork::BaseFrame> frames = spn.getBaseFrames();
        for (unsigned k = 0; k < 3; ++k)
            frames[1].atoms[k] = perms[p][k];
        spn.setBaseFrames(frames);
        spn.setup(stackingOnly());
        spn.computeForces();

        if (p == 0)
        {
            reference = spn.getStackingEnergy();
            ASSERT_LT(reference, -1.0f);
        }
        else
            EXPECT_NEAR(spn.getStackingEnergy(), reference, 1e-4f * std::abs(reference))
                << "permutation " << p;
    }
}

// The one property that keeps this term from double counting the hydrogen
// bond: two bases SIDE BY SIDE in a shared plane -- which is what a base pair
// is -- collect no stacking energy at all, because the axis joining them lies
// in both planes and the alignment factor is exactly zero there.
TEST(Stacking, a_coplanar_pair_collects_nothing)
{
    const std::array<Vector3f, 6> coplanar = {{Vector3f(0.00f, 0.00f, 0.0f), Vector3f(2.40f, 0.00f, 0.0f),
                                               Vector3f(1.00f, 2.00f, 0.0f), Vector3f(3.92f, 0.00f, 0.0f),
                                               Vector3f(6.32f, 0.00f, 0.0f), Vector3f(4.92f, 2.00f, 0.0f)}};
    spn::SpringNetwork spn;
    build(spn, coplanar);
    spn.setup(stackingOnly());
    spn.computeForces();

    EXPECT_FLOAT_EQ(spn.getStackingEnergy(), 0.0f);
    for (size_t i = 0; i < 6; ++i)
        EXPECT_LT(spn.getParticle(i).getForce().norm(), 1e-9f) << "atom " << i;
}

// And the other end of the same scale: one base directly above the other, at
// the equilibrium distance, with parallel normals sits at exactly the well
// depth. Fixes the sign and the scale of the whole term in one number.
TEST(Stacking, a_perfect_stack_sits_at_the_well_depth)
{
    const float h = 3.92f;
    const std::array<Vector3f, 6> stacked = {{Vector3f(0.00f, 0.00f, 0.0f), Vector3f(2.40f, 0.00f, 0.0f),
                                              Vector3f(1.00f, 2.00f, 0.0f), Vector3f(0.00f, 0.00f, h),
                                              Vector3f(2.40f, 0.00f, h), Vector3f(1.00f, 2.00f, h)}};
    spn::SpringNetwork spn;
    build(spn, stacked);
    spn.setup(stackingOnly(8.0, 3.92, 0.85));
    spn.computeForces();

    EXPECT_NEAR(spn.getStackingEnergy(), -8.0f, 1e-4f);
    // At the bottom of the well and perfectly aligned, every derivative
    // vanishes: no force at all.
    for (size_t i = 0; i < 6; ++i)
        EXPECT_LT(spn.getParticle(i).getForce().norm(), 1e-4f) << "atom " << i;
}

// A residue the table does not name, and one missing a plane atom, take no part
// in the term -- which is how a protein stays untouched by it.
TEST(Stacking, the_table_installs_one_frame_per_named_residue)
{
    Table table("# comment\nC  C6 C4 C2\nA  C8 C6 N3\n");
    io::BasePlaneRuleReader reader(table.path());
    reader.read();
    EXPECT_EQ(reader.rules().size(), 2u);

    spn::SpringNetwork spn;
    // A complete cytosine, a complete adenine, a cytosine missing C2, and an
    // alanine the table never names.
    place(spn, "C6", "C", 1, Vector3f(0.0f, 0.0f, 0.0f));
    place(spn, "C4", "C", 1, Vector3f(2.4f, 0.0f, 0.0f));
    place(spn, "C2", "C", 1, Vector3f(1.0f, 2.0f, 0.0f));
    place(spn, "C8", "A", 2, Vector3f(0.0f, 0.0f, 4.0f));
    place(spn, "C6", "A", 2, Vector3f(2.4f, 0.0f, 4.0f));
    place(spn, "N3", "A", 2, Vector3f(1.0f, 2.0f, 4.0f));
    place(spn, "C6", "C", 3, Vector3f(0.0f, 0.0f, 8.0f));
    place(spn, "C4", "C", 3, Vector3f(2.4f, 0.0f, 8.0f));
    place(spn, "CA", "ALA", 4, Vector3f(0.0f, 0.0f, 12.0f));

    reader.buildFrames(spn);
    ASSERT_EQ(spn.getBaseFrames().size(), 2u);
    EXPECT_EQ(spn.getBaseFrames()[0].atoms[0], 0u);
    EXPECT_EQ(spn.getBaseFrames()[1].atoms[0], 3u);
}

TEST(Stacking, a_table_naming_one_atom_twice_is_refused)
{
    Table table("C  C6 C4 C6\n");
    io::BasePlaneRuleReader reader(table.path());
    EXPECT_DEATH(reader.read(), "names 'C6' twice");
}
