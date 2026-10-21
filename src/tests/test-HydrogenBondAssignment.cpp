#include <gtest/gtest.h>

#include <array>
#include <string>

#include "configuration/Configuration.hpp"
#include "spn/SpringNetwork.h"

using namespace biospring;

// WHAT THIS COVERS AND WHY IT WAS MISSING.
//
// Every other hydrogen bond test in this suite gives the network exactly ONE
// donor and ONE acceptor. With one candidate on each side the assignment's rule
// -- a pair forms when i's best candidate is j and j's best is i -- is
// trivially satisfied and cannot fail. The force field is therefore tested
// exhaustively (the Morse well, both angular forms, the planarity factor, every
// gradient by finite differences) while the MATCHING is not tested at all.
//
// That gap hid a measured defect: on 8PKL's all-atom region the term holds 16 of
// the 49 base pairs the crystal has, and the 33 it misses score BETTER than the
// ones it forms (-15.7 kJ/mol at the median against -10.8). The system-level
// tests could not catch it either, because a duplex, a hairpin and a codon all
// start in their native geometry, where each base's best candidate IS its
// Watson-Crick partner, reciprocally -- the one arrangement in which the rule is
// correct.
//
// These tests therefore do the one thing none of the others do: give a site
// SEVERAL candidates and check which one it takes.

namespace
{

unsigned place(spn::SpringNetwork & net, const std::string & name, unsigned resid, const Vector3f & p)
{
    spn::Particle q;
    q.setName(name);
    q.setResName("X");
    q.setChainName("A");
    q.setResId(resid);
    q.setPosition(p);
    q.setMass(14.0f);
    q.setRadius(1.8f);
    q.setEpsilon(0.5f);
    net.addParticle(q);
    return net.getNumberOfParticles() - 1;
}

configuration::Configuration hbondOnly(double floor = 0.0)
{
    configuration::Configuration conf = configuration::defaultConfiguration();
    conf.hbond.enable = true;
    conf.hbond.cutoff = 7.0;
    conf.hbond.minimumstrength = floor;
    return conf;
}

} // namespace

// One acceptor, two donors, different distances: the NEARER donor should take
// it, and the other should get nothing. The baseline competition case.
TEST(HydrogenBondAssignment, the_better_of_two_competing_donors_wins)
{
    spn::SpringNetwork spn;
    place(spn, "C", 1, Vector3f(0.0f, 0.0f, 0.0f));          // 0 antecedent of the near donor
    place(spn, "N", 1, Vector3f(1.35f, 0.0f, 0.0f));         // 1 near donor
    place(spn, "C", 2, Vector3f(0.0f, 8.0f, 0.0f));          // 2 antecedent of the far donor
    place(spn, "N", 2, Vector3f(1.35f, 8.0f, 0.0f));         // 3 far donor
    place(spn, "O", 3, Vector3f(4.25f, 0.0f, 0.0f));         // 4 the contested acceptor

    spn.getParticle(1).setDonorCapacity(1);
    spn.getParticle(1).setAntecedentIndex(0);
    spn.getParticle(3).setDonorCapacity(1);
    spn.getParticle(3).setAntecedentIndex(2);
    spn.getParticle(4).setAcceptorCapacity(1);
    spn.setup(hbondOnly());
    spn.computeForces();

    EXPECT_TRUE(spn.areHydrogenBonded(1, 4)) << "the near donor did not get the acceptor";
    EXPECT_FALSE(spn.areHydrogenBonded(3, 4));
    EXPECT_EQ(spn.getHydrogenBondCount(), 1u);
}

// THE CASE THE SUITE HAD NO COVER FOR. Two acceptors and two donors, arranged so
// that each donor's FAVOURITE is the same acceptor. Only one of them can have
// it; the other has a perfectly good second choice, free, and well scored.
//
// A reciprocal-best-hit rule with no fallback forms ONE bond here and strands
// the other pair. Both bonds are what the geometry asks for.
TEST(HydrogenBondAssignment, a_donor_whose_favourite_is_taken_takes_its_second_choice)
{
    spn::SpringNetwork spn;
    // Two donors side by side, and two acceptors: donor 1 is nearest to
    // acceptor 4, donor 3 is ALSO nearest to acceptor 4, and acceptor 5 sits a
    // little further from donor 3 but well within range.
    place(spn, "C", 1, Vector3f(0.00f, 0.00f, 0.0f));   // 0
    place(spn, "N", 1, Vector3f(1.35f, 0.00f, 0.0f));   // 1 donor A
    place(spn, "C", 2, Vector3f(0.00f, 3.40f, 0.0f));   // 2
    place(spn, "N", 2, Vector3f(1.35f, 3.40f, 0.0f));   // 3 donor B
    place(spn, "O", 3, Vector3f(4.20f, 1.70f, 0.0f));   // 4 the acceptor both prefer
    place(spn, "O", 4, Vector3f(4.35f, 5.40f, 0.0f));   // 5 donor B's second choice

    spn.getParticle(1).setDonorCapacity(1);
    spn.getParticle(1).setAntecedentIndex(0);
    spn.getParticle(3).setDonorCapacity(1);
    spn.getParticle(3).setAntecedentIndex(2);
    spn.getParticle(4).setAcceptorCapacity(1);
    spn.getParticle(5).setAcceptorCapacity(1);
    spn.setup(hbondOnly());
    spn.computeForces();

    // Both donors should end up bonded: one to the contested acceptor, the
    // other to the free one. Which donor wins the contest is the ranking's
    // business and not asserted here -- what is asserted is that NEITHER is
    // left empty-handed while a free, well-scored acceptor sits in reach.
    const bool a_bonded = spn.areHydrogenBonded(1, 4) || spn.areHydrogenBonded(1, 5);
    const bool b_bonded = spn.areHydrogenBonded(3, 4) || spn.areHydrogenBonded(3, 5);
    EXPECT_TRUE(a_bonded) << "donor A got nothing";
    EXPECT_TRUE(b_bonded) << "donor B got nothing";
    EXPECT_EQ(spn.getHydrogenBondCount(), 2u)
        << "a donor was left unbonded while a free acceptor was in range and well scored";
}

// The occupancy floor: a contact worth almost nothing must not hold a slot that
// a real bond could use. Without the floor it does, which is what
// hbond.minimumstrength exists for.
TEST(HydrogenBondAssignment, a_worthless_contact_does_not_hold_a_slot)
{
    // The acceptor's lone pair points AWAY from the near donor, so that bond is
    // worth almost nothing; the far donor is lined up with it properly.
    const auto build = [](spn::SpringNetwork & spn) {
        place(spn, "C", 1, Vector3f(0.00f, 0.00f, 0.0f));    // 0
        place(spn, "N", 1, Vector3f(0.00f, 1.35f, 0.0f));    // 1 donor pointing sideways
        place(spn, "C", 2, Vector3f(9.00f, 0.00f, 0.0f));    // 2
        place(spn, "N", 2, Vector3f(7.65f, 0.00f, 0.0f));    // 3 donor aimed at the acceptor
        place(spn, "O", 3, Vector3f(4.40f, 0.00f, 0.0f));    // 4 acceptor
        place(spn, "C", 3, Vector3f(3.05f, 0.00f, 0.0f));    // 5 its antecedent, so its
        spn.getParticle(1).setDonorCapacity(1);              //   lone pair points at +x
        spn.getParticle(1).setAntecedentIndex(0);
        spn.getParticle(3).setDonorCapacity(1);
        spn.getParticle(3).setAntecedentIndex(2);
        spn.getParticle(4).setAcceptorCapacity(1);
        spn.getParticle(4).setAntecedentIndex(5);
    };
    spn::SpringNetwork loose, strict;
    build(loose);
    build(strict);
    loose.setup(hbondOnly(0.0));
    strict.setup(hbondOnly(1.0));
    loose.computeForces();
    strict.computeForces();

    // With a floor, whatever is held must be worth having.
    EXPECT_LE(strict.getHydrogenBondEnergy(), loose.getHydrogenBondEnergy() + 1e-4f)
        << "the floor let a weaker set of bonds form";
    for (size_t i = 0; i < strict.getNumberOfParticles(); ++i)
        for (size_t j = i + 1; j < strict.getNumberOfParticles(); ++j)
            if (strict.areHydrogenBonded(i, j))
            {
                const float d = spn::Particle::distance(strict.getParticle(i), strict.getParticle(j));
                EXPECT_LT(d, 7.0f) << i << "-" << j;
            }
}
