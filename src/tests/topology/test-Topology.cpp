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

TEST(Topology, to_spring_network_with_ghost_particle)
{
    topology::Topology top;

    topology::ParticleProperties propsB, propsC, propsRef;
    propsB.set_name("B");
    propsB.set_position(Vector3f(0.0, 0.0, 0.0));
    propsC.set_name("C");
    propsC.set_position(Vector3f(0.0, 0.0, 2.0));
    propsRef.set_name("REF");
    propsRef.set_position(Vector3f(1.0, 0.0, 0.0));

    top.add_particle(topology::Particle(propsB));
    top.add_particle(topology::Particle(propsC));
    top.add_particle(topology::Particle(propsRef));

    topology::ParticleProperties ghost_props;
    ghost_props.set_name("GHOST");
    ghost_props.set_static(true);
    ghost_props.set_mass(0.0f);

    // A ghost is the image of its reference ATOM under a rotation of delta
    // about the B->C axis (see spn::GhostParticle): r/theta are no longer
    // read, the radius and axial position come from that atom itself.
    // Axis is +z here and REF sits at (1,0,0), so delta=90deg must land the
    // ghost on (0,1,0).
    const float r = 1.5f, theta_deg = 90.0f, delta_deg = 90.0f;
    top.add_ghost_particle(topology::Particle(ghost_props), top.get_particle(0), top.get_particle(1),
                           top.get_particle(2), r, theta_deg, delta_deg, static_cast<unsigned>(spn::GhostPlacement::AxisRotation));

    ASSERT_EQ(top.number_of_particles(), 4);
    EXPECT_TRUE(top.is_ghost_particle(top.get_particle(3).unique_id()));
    EXPECT_FALSE(top.is_ghost_particle(top.get_particle(0).unique_id()));

    spn::SpringNetwork spn;
    top.to_spring_network(spn);

    ASSERT_EQ(spn.getNumberOfParticles(), 4);
    ASSERT_EQ(spn.getGhostParticles().size(), 1);

    const spn::GhostParticleBinding & binding = spn.getGhostParticles()[0];
    EXPECT_EQ(binding.ownIndex, 3u);
    EXPECT_EQ(binding.anchorBIndex, 0u);
    EXPECT_EQ(binding.anchorCIndex, 1u);
    EXPECT_EQ(binding.anchorRefIndex, 2u);

    const spn::Particle & ghost = spn.getParticle(3);
    EXPECT_TRUE(ghost.isStatic());
    EXPECT_FLOAT_EQ(ghost.getMass(), 0.0f);
    // REF=(1,0,0) rotated 90deg about the +z axis lands on (0,1,0), at the
    // same distance from the axis and the same height along it. Placed here
    // through the full Topology -> SpringNetwork conversion instead of
    // calling the placement formula directly.
    EXPECT_NEAR(ghost.getPosition().getX(), 0.0f, 1e-4);
    EXPECT_NEAR(ghost.getPosition().getY(), 1.0f, 1e-4);
    EXPECT_NEAR(ghost.getPosition().getZ(), 0.0f, 1e-4);

    // Force redistribution TRANSFERS what acted on the ghost: the three
    // anchors must end up carrying exactly it, and the ghost's own force
    // must be cleared.
    const Vector3f applied(1.0f, 2.0f, 3.0f);
    spn.getParticle(3).addForce(applied);
    spn.redistributeGhostForces();

    EXPECT_FLOAT_EQ(spn.getParticle(3).getForce().getX(), 0.0f);
    EXPECT_FLOAT_EQ(spn.getParticle(3).getForce().getY(), 0.0f);
    EXPECT_FLOAT_EQ(spn.getParticle(3).getForce().getZ(), 0.0f);

    Vector3f total_redistributed =
        spn.getParticle(0).getForce() + spn.getParticle(1).getForce() + spn.getParticle(2).getForce();
    EXPECT_NEAR(total_redistributed.getX(), applied.getX(), 1e-3);
    EXPECT_NEAR(total_redistributed.getY(), applied.getY(), 1e-3);
    EXPECT_NEAR(total_redistributed.getZ(), applied.getZ(), 1e-3);

    // ... and the torque about B too, which is the part a naive split gets
    // wrong (see doc/BondedForceFieldSprings.md).
    const Vector3f B = spn.getParticle(0).getPosition();
    Vector3f torque_out = ((spn.getParticle(1).getPosition() - B) ^ spn.getParticle(1).getForce()) +
                          ((spn.getParticle(2).getPosition() - B) ^ spn.getParticle(2).getForce());
    Vector3f torque_in = (ghost.getPosition() - B) ^ applied;
    EXPECT_NEAR((torque_out - torque_in).norm(), 0.0f, 1e-3);

    // Position tracking: the ghost follows its reference ATOM, not a fixed
    // offset -- move REF out to twice its radius and the ghost doubles too.
    spn.getParticle(2).setPosition(Vector3f(2.0f, 0.0f, 0.0f));
    spn.updateGhostPositions();
    EXPECT_NEAR(spn.getParticle(3).getPosition().getX(), 0.0f, 1e-4);
    EXPECT_NEAR(spn.getParticle(3).getPosition().getY(), 2.0f, 1e-4);
    EXPECT_NEAR(spn.getParticle(3).getPosition().getZ(), 0.0f, 1e-4);
}

// Regression test for a bug found while adding NetCDF ghost-particle I/O:
// Topology's copy constructor/assignment operator only ever copied
// particles and springs, never _ghost_particles -- silently dropping every
// ghost binding whenever a Topology was copied (e.g. io::readTopology's
// by-value return from NetCDFReader::getTopology()). Also exercises the
// uid-remapping this required: ParticleCollection::operator= mints a fresh
// unique id per particle (see its own comment), so a naive `_ghost_particles
// = other._ghost_particles` would carry over stale, non-matching keys.
TEST(Topology, copy_preserves_ghost_particle)
{
    topology::Topology top;

    topology::ParticleProperties propsB, propsC, propsRef;
    propsB.set_name("B");
    propsB.set_position(Vector3f(0.0, 0.0, 0.0));
    propsC.set_name("C");
    propsC.set_position(Vector3f(0.0, 0.0, 2.0));
    propsRef.set_name("REF");
    propsRef.set_position(Vector3f(1.0, 0.0, 0.0));

    top.add_particle(topology::Particle(propsB));
    top.add_particle(topology::Particle(propsC));
    top.add_particle(topology::Particle(propsRef));

    topology::ParticleProperties ghost_props;
    ghost_props.set_name("GHOST");
    ghost_props.set_static(true);
    ghost_props.set_mass(0.0f);

    const float r = 1.5f, theta_deg = 90.0f, delta_deg = 0.0f;
    top.add_ghost_particle(topology::Particle(ghost_props), top.get_particle(0), top.get_particle(1),
                           top.get_particle(2), r, theta_deg, delta_deg, static_cast<unsigned>(spn::GhostPlacement::AxisRotation));

    // Copy constructor.
    topology::Topology copy(top);
    ASSERT_EQ(copy.number_of_particles(), 4);
    EXPECT_TRUE(copy.is_ghost_particle(copy.get_particle(3).unique_id()));
    EXPECT_FALSE(copy.is_ghost_particle(copy.get_particle(0).unique_id()));
    // The copy's uids are freshly minted, distinct from the original's.
    EXPECT_NE(copy.get_particle(3).unique_id(), top.get_particle(3).unique_id());

    const topology::Topology::GhostParticleInfo & info = copy.get_ghost_particle_info(copy.get_particle(3).unique_id());
    EXPECT_EQ(info.anchor_B_uid, copy.get_particle(0).unique_id());
    EXPECT_EQ(info.anchor_C_uid, copy.get_particle(1).unique_id());
    EXPECT_EQ(info.anchor_ref_uid, copy.get_particle(2).unique_id());
    EXPECT_FLOAT_EQ(static_cast<float>(info.r), r);

    spn::SpringNetwork spn;
    copy.to_spring_network(spn);
    ASSERT_EQ(spn.getGhostParticles().size(), 1);

    // Assignment operator.
    topology::Topology assigned;
    assigned = top;
    ASSERT_EQ(assigned.number_of_particles(), 4);
    EXPECT_TRUE(assigned.is_ghost_particle(assigned.get_particle(3).unique_id()));

    spn::SpringNetwork spn2;
    assigned.to_spring_network(spn2);
    ASSERT_EQ(spn2.getGhostParticles().size(), 1);
}

// Regression test for a bug found while running the first real dynamics
// validation with the ghost-particle dihedral mechanism active: every
// single dihedral ghost-ghost spring connects two particles that are BOTH
// static (spn::GhostParticle virtual sites are always massless/static by
// design). spn::Spring::computeForce had a long-standing early exit
// ("skip if both endpoints are non-dynamic" -- a valid optimization for two
// genuinely frozen real atoms) that silently zeroed out force *and* energy
// for every dihedral ghost spring, with no error or warning: a whole real
// dynamics run showed exactly 0.00 dihedral energy for 5000 steps before
// this was caught. None of the existing ghost-particle tests exercised
// this path, since they call GhostParticle::redistributeForce/
// SpringNetwork::redistributeGhostForces directly with a manually-pushed
// force, bypassing Spring::computeForce entirely.
TEST(Topology, dihedral_ghost_spring_applies_force)
{
    topology::Topology top;

    // First ghost's anchors.
    topology::ParticleProperties propsB1, propsC1, propsRef1;
    propsB1.set_position(Vector3f(0.0, 0.0, 0.0));
    propsC1.set_position(Vector3f(0.0, 0.0, 2.0));
    propsRef1.set_position(Vector3f(1.0, 0.0, 0.0));
    top.add_particle(topology::Particle(propsB1));
    top.add_particle(topology::Particle(propsC1));
    top.add_particle(topology::Particle(propsRef1));

    // Second ghost's anchors, far away, so the two ghosts end up clearly
    // separated (a large, unambiguous distance mismatch against the small
    // d0 chosen below, regardless of the exact placement geometry).
    topology::ParticleProperties propsB2, propsC2, propsRef2;
    propsB2.set_position(Vector3f(20.0, 0.0, 0.0));
    propsC2.set_position(Vector3f(20.0, 0.0, 2.0));
    propsRef2.set_position(Vector3f(21.0, 0.0, 0.0));
    top.add_particle(topology::Particle(propsB2));
    top.add_particle(topology::Particle(propsC2));
    top.add_particle(topology::Particle(propsRef2));

    topology::ParticleProperties ghost_props;
    ghost_props.set_static(true);
    ghost_props.set_mass(0.0f);

    topology::Particle & ghost1 = top.add_ghost_particle(topology::Particle(ghost_props), top.get_particle(0),
                                                          top.get_particle(1), top.get_particle(2), 1.0f, 90.0f, 0.0f, static_cast<unsigned>(spn::GhostPlacement::AxisRotation));
    topology::Particle & ghost2 = top.add_ghost_particle(topology::Particle(ghost_props), top.get_particle(3),
                                                          top.get_particle(4), top.get_particle(5), 1.0f, 90.0f, 0.0f, static_cast<unsigned>(spn::GhostPlacement::AxisRotation));

    // d0 = 0.1: guaranteed far from the ~20 A actual distance between the
    // two ghosts, so (d - d0)^2 is unambiguously large.
    top.add_dihedral_spring(spn::SpringNetwork::DIHEDRAL_SIDECHAIN, ghost1, ghost2, 0.1, 10.0);

    spn::SpringNetwork spn;
    top.to_spring_network(spn);
    ASSERT_EQ(spn.getGhostParticles().size(), 2);
    // computeDihedralForces() needs a force field (SpringNetwork::_ff) to
    // compute the harmonic term through -- only setup() constructs it.
    spn.setup(configuration::defaultConfiguration());

    spn.computeDihedralForces();

    // Before this bug's fix, both would be exactly 0.0f: computeForce's
    // "both endpoints non-dynamic" early exit silently skipped all real
    // work for this spring.
    EXPECT_GT(spn.getDihedralEnergy(), 1.0f);

    bool any_ghost_force_nonzero = false;
    for (const spn::GhostParticleBinding & binding : spn.getGhostParticles())
    {
        if (spn.getParticle(binding.ownIndex).getForce().norm() > 1e-6f)
            any_ghost_force_nonzero = true;
    }
    EXPECT_TRUE(any_ghost_force_nonzero);

    // Redistributing should move that force onto the (real, dynamic)
    // anchors instead, and reset each ghost's own force back to zero.
    spn.redistributeGhostForces();
    for (const spn::GhostParticleBinding & binding : spn.getGhostParticles())
        EXPECT_FLOAT_EQ(spn.getParticle(binding.ownIndex).getForce().norm(), 0.0f);

    float total_anchor_force = 0.0f;
    for (size_t i = 0; i < 6; ++i)
        total_anchor_force += spn.getParticle(i).getForce().norm();
    EXPECT_GT(total_anchor_force, 0.0f);
}

// Global gradient/energy consistency check by central finite differences,
// over a synthetic system exercising every spring type at once: a regular
// spring, a STRETCH spring, a BEND ghost-ghost spring, and a DIHEDRAL
// ghost-ghost spring with a nonzero dc_offset. For every dynamic particle
// and every coordinate, the force the network applies (after ghost-force
// redistribution) must equal -dE/dx of the energy it reports, with the
// ghosts treated as implicit functions of their anchors (repositioned from
// the displaced anchors before each energy evaluation, exactly as
// computeStep's update order produces on the next step). This is the test
// that catches the deep-consistency bug classes hit before: a wrong
// Jacobian-transpose in GhostParticle::redistributeForce (would break
// anchors' FD match), a dc_offset leaking into forces (dc shifts E by a
// constant, so FD is blind to it ONLY if forces ignore it too), or a
// spring type silently skipped on one side but not the other.
TEST(Topology, forces_match_energy_gradient_by_finite_differences)
{
    topology::Topology top;

    // Two anchor triplets, deliberately asymmetric positions (no axis
    // aligned with a coordinate plane, so every FD component is nonzero).
    const std::array<std::array<double, 3>, 6> pos = {{{0.1, -0.2, 0.05},
                                                       {0.4, 0.3, 2.1},
                                                       {1.2, 0.15, 0.3},
                                                       {5.6, 0.4, 0.2},
                                                       {5.2, -0.3, 2.3},
                                                       {6.5, 0.25, 0.6}}};
    for (size_t i = 0; i < pos.size(); ++i)
    {
        topology::ParticleProperties p;
        p.set_position(Vector3f(static_cast<float>(pos[i][0]), static_cast<float>(pos[i][1]),
                                static_cast<float>(pos[i][2])));
        top.add_particle(topology::Particle(p));
    }

    topology::ParticleProperties ghost_props;
    ghost_props.set_static(true);
    ghost_props.set_mass(0.0f);

    // Non-trivial placement (theta != 90, delta != 0) so the placement
    // Jacobian has every term active.
    topology::Particle & ghost1 = top.add_ghost_particle(topology::Particle(ghost_props), top.get_particle(0),
                                                          top.get_particle(1), top.get_particle(2), 1.3f, 110.0f, 35.0f, static_cast<unsigned>(spn::GhostPlacement::AxisRotation));
    topology::Particle & ghost2 = top.add_ghost_particle(topology::Particle(ghost_props), top.get_particle(3),
                                                          top.get_particle(4), top.get_particle(5), 1.1f, 75.0f, -20.0f, static_cast<unsigned>(spn::GhostPlacement::AxisRotation));

    top.add_spring(top.get_particle(2), top.get_particle(5), 2.0, 4.0);
    top.add_stretch_spring(top.get_particle(0), top.get_particle(1), 1.5, 6.0);
    top.add_bend_spring(ghost1, ghost2, 3.0, 2.5);
    // d0 far from the actual ghost-ghost distance so the dihedral force is
    // large; dc_offset nonzero to prove it shifts reported energy only.
    top.add_dihedral_spring(spn::SpringNetwork::DIHEDRAL_SIDECHAIN, ghost1, ghost2, 0.5, 7.0).set_dc_offset(3.21);

    spn::SpringNetwork spn;
    top.to_spring_network(spn);
    ASSERT_EQ(spn.getGhostParticles().size(), 2);
    // Default config has spring.enable = false; computeForces() honours that
    // master switch (unlike the direct computeDihedralForces() call the
    // other ghost test uses), so it must be on for anything to happen.
    configuration::Configuration conf = configuration::defaultConfiguration();
    conf.spring.enable = true;
    spn.setup(conf);

    auto reset_forces = [&]() {
        for (size_t i = 0; i < spn.getNumberOfParticles(); ++i)
            spn.getParticle(i).resetForce();
    };
    auto total_energy = [&]() {
        reset_forces();
        spn.updateGhostPositions();   // ghosts are functions of the anchors
        spn.computeForces();
        return spn.getSpringEnergy() + spn.getStretchEnergy() + spn.getBendEnergy() + spn.getDihedralEnergy();
    };

    // Analytic forces at the base configuration.
    const float e0 = total_energy();
    spn.redistributeGhostForces();
    std::array<Vector3f, 6> analytic;
    for (size_t i = 0; i < 6; ++i)
        analytic[i] = spn.getParticle(i).getForce();
    // Sanity: the dihedral spring is doing real work, and its reported
    // energy carries the dc correction (energy would differ by exactly
    // 3.21 without it -- checked implicitly by the FD match below, since a
    // constant shift cannot change any gradient).
    EXPECT_GT(std::abs(e0), 1.0f);

    const float h = 1e-2f;   // A; large enough to beat float32 noise on E
    for (size_t i = 0; i < 6; ++i)
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
            const float e_plus = total_energy();

            set_dim(displaced, get_dim(saved) - h);
            spn.getParticle(i).setPosition(displaced);
            const float e_minus = total_energy();

            spn.getParticle(i).setPosition(saved);

            // Forces on particles are stored in integrator units
            // (Da.A.fs-2, see GLOBAL_SPRING_FORCE_CONVERT); energies are
            // reported in kJ/mol. Convert the analytic force back to
            // kJ/mol/A to compare against -dE/dx.
            const float fd = -(e_plus - e_minus) / (2.0f * h);
            const float an = get_dim(analytic[i]) /
                             static_cast<float>(biospring::forcefield::GLOBAL_SPRING_FORCE_CONVERT);
            EXPECT_NEAR(fd, an, 2e-2f * std::max(1.0f, std::abs(an)))
                << "particle " << i << " dim " << dim << " FD " << fd << " vs analytic " << an;
        }
    }
}

// -- Main function  ----------------------------------------------------------
int main(int argc, char * argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
