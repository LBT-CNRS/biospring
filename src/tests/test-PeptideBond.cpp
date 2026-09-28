#include <gtest/gtest.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>

#include "configuration/Configuration.hpp"
#include "spn/PeptideBond.h"
#include "spn/SpringNetwork.h"

using namespace biospring;

namespace
{

// A .rbody holding only what this is about: the peptide plane ahead of an
// alanine. Written to a temporary file rather than pointing at the shipped
// data/reducerules/ProteinAtomRigidGroups.rbody, so the test says in one place
// what it expects the rule to be and does not fail the day that file gains an
// unrelated group.
class RulesFile
{
  public:
    RulesFile()
    {
        _path = std::string(std::tmpnam(nullptr)) + ".rbody";
        std::ofstream out(_path);
        out << "# name resname atoms...\n";
        out << "A_CA ALA N CA C CB\n";
        out << "A_PSI ALA CA C O +N +H +CA\n";
        out.close();
    }
    ~RulesFile() { std::remove(_path.c_str()); }
    const std::string & path() const { return _path; }

  private:
    std::string _path;
};

unsigned place(spn::SpringNetwork & network, const std::string & name, unsigned resid, const Vector3f & position)
{
    spn::Particle p;
    p.setName(name);
    p.setResName("ALA");
    p.setChainName("A");
    p.setResId(resid);
    p.setPosition(position);
    p.setMass(12.0f);
    p.setRadius(1.7f);
    p.setEpsilon(0.5f);
    network.addParticle(p);
    return network.getNumberOfParticles() - 1;
}

// Two alanines: one bringing a carbonyl, the other a free amine. The carbonyl
// sits in the xy plane, so its sp2 NORMAL is z, and the amine is placed on that
// normal -- which is the Burgi-Dunitz approach, and the only place an
// out-of-plane lobe site scores.
//
// `attack` is how far up the normal the nitrogen sits, in A. Large puts it out
// of range; `inPlane` instead puts it in the carbonyl's own plane at the same
// distance, where the angular factor must kill it however close it gets.
struct Pair
{
    unsigned ca1, c, o, oxt, n, h, ca2;
};

Pair build(spn::SpringNetwork & network, configuration::Configuration & config, const std::string & rules,
           float attack, bool inPlane, unsigned dwell, unsigned ramp, unsigned steps)
{
    Pair idx;
    // Residue 1: the carbonyl, in the xy plane. C at the origin, O along -x-ish
    // and CA on the other side, the ideal trans geometry this file also ramps
    // towards.
    idx.ca1 = place(network, "CA", 1, {-0.6725f, 1.3687f, 0.0f});
    idx.c = place(network, "C", 1, {0.0f, 0.0f, 0.0f});
    idx.o = place(network, "O", 1, {-0.6705f, -1.0324f, 0.0f});
    // The leaving group: a carboxyl's second oxygen, bound to the same carbon.
    // Not named by the _PSI rule, so it takes no part in the peptide plane --
    // it is there to be let go of.
    idx.oxt = place(network, "OXT", 1, {-1.2800f, 0.7400f, 0.9000f});

    // Residue 2: the amine, on the carbonyl's normal (or in its plane). Its own
    // two atoms are placed RELATIVE to the approach axis, not in fixed
    // directions, so the nitrogen's lone pair points back at the carbon by the
    // same angle in both cases -- otherwise the in-plane test would be refused
    // by the nitrogen's angle and would say nothing about the carbon's.
    const Vector3f axis = inPlane ? Vector3f(1.0f, 0.0f, 0.0f) : Vector3f(0.0f, 0.0f, 1.0f);
    const Vector3f side = inPlane ? Vector3f(0.0f, 0.0f, 1.0f) : Vector3f(1.0f, 0.0f, 0.0f);
    const Vector3f other = Vector3f(0.0f, 1.0f, 0.0f);
    const Vector3f approach = axis * attack;
    idx.n = place(network, "N", 2, approach);
    idx.h = place(network, "H", 2, approach + other * 0.7f + axis * 0.5f);
    idx.ca2 = place(network, "CA", 2, approach + side * 0.8f + axis * 1.2f);

    // The carbonyl carbon's Burgi-Dunitz site: axis from its own oxygen, plane
    // fixed by CA, lobes OUT of that plane at 75 degrees -- 105 degrees of the
    // C=O axis, because e1 points away from O.
    network.getParticle(idx.c).setAntecedentIndex(static_cast<int>(idx.o));
    network.getParticle(idx.c).setAntecedentIndex2(static_cast<int>(idx.ca1));
    network.getParticle(idx.c).setLobes(spn::Particle::HBOND_LOBES_OUT_OF_PLANE, 75.0f);

    // The nitrogen's lone pair: one antecedent, so the direction is the CA->N
    // axis extended, and CA2 is placed so that it points back at the carbon.
    network.getParticle(idx.n).setAntecedentIndex(static_cast<int>(idx.ca2));

    // The two residues' own internal springs, so each is a body rather than
    // three free points. Nothing crosses between them: that is what the bond
    // has to create.
    const auto tie = [&](unsigned a, unsigned b)
    { network.addSpring(a, b, (network.getParticle(a).getPosition() - network.getParticle(b).getPosition()).norm(), 650.0f); };
    tie(idx.c, idx.oxt);
    tie(idx.ca1, idx.oxt);
    tie(idx.ca1, idx.c);
    tie(idx.c, idx.o);
    tie(idx.ca1, idx.o);
    tie(idx.n, idx.h);
    tie(idx.n, idx.ca2);
    tie(idx.h, idx.ca2);

    config = configuration::defaultConfiguration();
    config.sim.nbsteps = steps;
    config.sim.timestep = 1.0;
    config.spring.enable = true;
    config.viscosity.enable = true;
    config.viscosity.value = 0.5;
    config.peptidebond.enable = true;
    config.peptidebond.path = rules;
    config.peptidebond.group = "PSI";
    config.peptidebond.bond = "C:N";
    config.peptidebond.stiffness = 650.0;
    config.peptidebond.distance = 3.2;
    config.peptidebond.weight = 0.5;
    config.peptidebond.dwell = dwell;
    config.peptidebond.ramp = ramp;
    network.setup(config);
    return idx;
}

} // namespace

// The whole point of stage two: the bond is made, and it is made of the SAME
// springs the model would have had all along.
TEST(PeptideBond, FormsTheFifteenSpringsOfItsOwnPsiRule)
{
    RulesFile rules;
    spn::SpringNetwork network;
    configuration::Configuration config;
    const Pair idx = build(network, config, rules.path(), 2.9f, false, 1, 0, 2);

    const unsigned before = network.getNumberOfSprings();
    network.run();

    ASSERT_EQ(network.getPeptideBonds().getNumberOfBonds(), 1u)
        << "no bond formed at 2.9 A on the normal, which IS the attack conformation; closest approach "
        << network.getPeptideBonds().getClosestApproach() << " A, best weight "
        << network.getPeptideBonds().getBestWeight();

    // Six atoms pairwise is fifteen springs, of which six already existed (the
    // two residues' own), so nine are new.
    EXPECT_EQ(network.getNumberOfSprings() - before, 9u)
        << "the _PSI group is six atoms sprung pairwise; six of the fifteen were already there";

    // And they are on the right atoms: every pair among the six, and nothing
    // reaching outside them.
    const unsigned six[6] = {idx.ca1, idx.c, idx.o, idx.n, idx.h, idx.ca2};
    for (unsigned a : six)
        for (unsigned b : six)
            if (a != b)
                EXPECT_TRUE(network.getParticle(a).isInSpringNeighbors(b))
                    << "particles " << a << " and " << b << " of the peptide plane are not sprung together";
}

// A bond that forms because two atoms happened to be close, whatever the angle,
// would be a distance criterion wearing a chemistry name.
TEST(PeptideBond, RefusesAnApproachInTheCarbonylPlane)
{
    RulesFile rules;
    spn::SpringNetwork network;
    configuration::Configuration config;
    build(network, config, rules.path(), 2.9f, true, 1, 0, 2);

    const unsigned before = network.getNumberOfSprings();
    network.run();

    EXPECT_EQ(network.getPeptideBonds().getNumberOfBonds(), 0u)
        << "bonded from inside the sp2 plane, where a nucleophile cannot reach the pi*";
    EXPECT_EQ(network.getNumberOfSprings(), before) << "springs were created for a bond that must not have formed";
    // The distance condition WAS met, so this test is about the angle and not
    // about the pair being out of range.
    EXPECT_LT(network.getPeptideBonds().getClosestApproach(), config.peptidebond.distance)
        << "the pair was out of range, so the angular condition was never the thing being tested";
    EXPECT_LT(network.getPeptideBonds().getBestWeight(), config.peptidebond.weight);
}

// Too far is too far, however good the angle.
TEST(PeptideBond, RefusesAnApproachOutOfRange)
{
    RulesFile rules;
    spn::SpringNetwork network;
    configuration::Configuration config;
    build(network, config, rules.path(), 6.0f, false, 1, 0, 2);

    network.run();
    EXPECT_EQ(network.getPeptideBonds().getNumberOfBonds(), 0u);
    EXPECT_GT(network.getPeptideBonds().getClosestApproach(), config.peptidebond.distance);
}

// The dwell is what separates an approach from a fly-by, so it has to be a
// count of steps and not a formality.
TEST(PeptideBond, WaitsForTheDwellAndThenFormsExactlyOnce)
{
    RulesFile rules;
    const unsigned DWELL = 40;

    {
        spn::SpringNetwork tooshort;
        configuration::Configuration config;
        build(tooshort, config, rules.path(), 2.9f, false, DWELL, 0, DWELL - 1);
        tooshort.run();
        EXPECT_EQ(tooshort.getPeptideBonds().getNumberOfBonds(), 0u)
            << "formed after " << DWELL - 1 << " steps although the dwell is " << DWELL;
    }
    {
        spn::SpringNetwork longenough;
        configuration::Configuration config;
        build(longenough, config, rules.path(), 2.9f, false, DWELL, 0, 4 * DWELL);
        longenough.run();
        // Exactly one: the pair must leave both candidate sets the moment it
        // bonds, or the next step makes the same bond again -- fifteen more
        // springs on the same six atoms, at whatever length they then have.
        EXPECT_EQ(longenough.getPeptideBonds().getNumberOfBonds(), 1u)
            << "a run four times the dwell made " << longenough.getPeptideBonds().getNumberOfBonds() << " bonds";
        ASSERT_FALSE(longenough.getPeptideBonds().getBonds().empty());
        EXPECT_GE(longenough.getPeptideBonds().getBonds().front().step, DWELL)
            << "formed before the dwell had elapsed";
    }
}

// Born where the atoms already are, so that making a bond does not create
// energy: at 650 kJ/mol/A2 the 1.6 A between an attack distance and a peptide
// bond would be 830 kJ/mol appearing inside one timestep.
TEST(PeptideBond, CostsNoEnergyAtTheInstantItForms)
{
    RulesFile rules;
    spn::SpringNetwork network;
    configuration::Configuration config;
    const Pair idx = build(network, config, rules.path(), 2.9f, false, 1, 0, 2);
    network.run();

    ASSERT_EQ(network.getPeptideBonds().getNumberOfBonds(), 1u);
    const auto & bond = network.getPeptideBonds().getBonds().front();
    for (size_t k = 0; k < bond.springs.size(); ++k)
    {
        const spn::Spring & s = network.getSpring(bond.springs[k]);
        EXPECT_NEAR(s.getEquilibrium(), bond.born[k], 1e-4f)
            << "spring " << k << " was not born at the length its two atoms already had";
    }
    // And the new C-N spring is the attack distance, not a peptide bond yet.
    EXPECT_NEAR((network.getParticle(idx.c).getPosition() - network.getParticle(idx.n).getPosition()).norm(), 2.9f,
                0.2f);
}

// The thermostatted path is a DIFFERENT function: computeStep() hands over to
// computeStepBAOAB() and returns, so anything called only at the end of the
// first never runs under a thermostat. It went unnoticed because every other
// test here runs without one -- and a thermostat is exactly how one asks for
// the case this whole term is about, two molecules warm enough to meet.
TEST(PeptideBond, FormsUnderAThermostatToo)
{
    RulesFile rules;
    spn::SpringNetwork network;
    configuration::Configuration config;
    build(network, config, rules.path(), 2.9f, false, 1, 0, 50);
    config.thermostat.enable = true;
    config.thermostat.temperature = 50.0;   // warm enough to move, cold enough
    network.setup(config);                  // to leave the attack geometry alone

    network.run();

    EXPECT_EQ(network.getPeptideBonds().getNumberOfBonds(), 1u)
        << "no bond under a thermostat, where the same fixture makes one without; "
        << "closest " << network.getPeptideBonds().getClosestApproach() << " A, weight "
        << network.getPeptideBonds().getBestWeight();
}

// A real carboxyl does not gain a fourth partner: it loses its hydroxyl. What
// distinguishes a bonded atom from a free one here is its springs, so leaving
// means those springs break -- and break BOTH ways, because a spring is also
// what keeps a pair out of the steric and Coulomb terms. A leaving group still
// excluded from them would sit on top of the carbon for ever.
TEST(PeptideBond, LetsTheLeavingGroupGo)
{
    RulesFile rules;
    spn::SpringNetwork network;
    configuration::Configuration config;
    const unsigned RAMP = 40;
    const Pair idx = build(network, config, rules.path(), 2.9f, false, 1, RAMP, RAMP + 10);
    config.peptidebond.leaving = "OXT";
    network.setup(config);

    ASSERT_TRUE(network.getParticle(idx.oxt).isInSpringNeighbors(idx.c))
        << "the fixture did not bind the leaving group in the first place";
    const unsigned springs = network.getNumberOfSprings();

    network.run();

    ASSERT_EQ(network.getPeptideBonds().getNumberOfBonds(), 1u);
    EXPECT_FALSE(network.getParticle(idx.oxt).isInSpringNeighbors(idx.c))
        << "the leaving group is still bound to the carbon it left";
    // And it left from FAR ENOUGH OUT. Released in place it would sit 1.2 A
    // from the carbon with the exclusion gone and the full steric wall on: on a
    // concentrated soup that took the temperature to 3.8e+18 K in one step.
    EXPECT_GT((network.getParticle(idx.oxt).getPosition() - network.getParticle(idx.c).getPosition()).norm(), 2.5f)
        << "the leaving group was let go while still on top of its carbon";
    EXPECT_FALSE(network.getParticle(idx.oxt).isInSpringNeighbors(idx.ca1))
        << "the leaving group is still bound to the alpha carbon";
    // Broken, not erased: the slot stays so that nothing else has to be
    // renumbered, and it is the stiffness that says the bond is gone.
    EXPECT_GE(network.getNumberOfSprings(), springs);
    for (unsigned id = 0; id < network.getNumberOfSprings(); ++id)
    {
        const spn::Spring & s = network.getSpring(id);
        const unsigned a = s.getParticle1().getId(), b = s.getParticle2().getId();
        if (a == idx.oxt || b == idx.oxt)
            EXPECT_FLOAT_EQ(s.getStiffness(), 0.0f)
                << "spring " << id << " still pulls on the leaving group";
    }
}

// And nothing leaves unless it was named: the default is a carbon that simply
// gains a partner, which is what a tRNA-esterified one does.
TEST(PeptideBond, KeepsEverythingWhenNoLeavingGroupIsNamed)
{
    RulesFile rules;
    spn::SpringNetwork network;
    configuration::Configuration config;
    const Pair idx = build(network, config, rules.path(), 2.9f, false, 1, 0, 2);
    network.run();

    ASSERT_EQ(network.getPeptideBonds().getNumberOfBonds(), 1u);
    EXPECT_TRUE(network.getParticle(idx.oxt).isInSpringNeighbors(idx.c));
}

// A hydrogen bond held on the pair must be RELEASED when the pair becomes
// covalent, or the two terms fight over the same two atoms.
//
// The assignment refuses to propose a sprung pair, which was enough while
// nothing changed the springs during a run. It does not re-examine a bond it
// is already holding, so one made before the topology changed survives it --
// and at the peptide C-N length of 1.33 A the Morse well is deep into its
// repulsive wall. Measured on two alanines over MDDriver before the fix: a
// "hydrogen bond energy" of +4.31 kJ/mol, pushing apart exactly the pair the
// new springs hold together.
TEST(PeptideBond, ReleasesTheHydrogenBondItHasJustMadeCovalent)
{
    RulesFile rules;
    spn::SpringNetwork network;
    configuration::Configuration config;
    const Pair idx = build(network, config, rules.path(), 2.9f, false, 1, 20, 200);

    // The same two atoms the bond will be made between, as a donor/acceptor
    // couple: at 2.9 A with both sites aimed at each other this is exactly what
    // the term is for, and it holds the pair from the first step.
    network.getParticle(idx.n).setDonorCapacity(1);
    network.getParticle(idx.c).setAcceptorCapacity(1);
    config.hbond.enable = true;
    config.hbond.cutoff = 7.0;
    config.hbond.scale = 1.0;
    network.setup(config);

    network.run();

    ASSERT_EQ(network.getPeptideBonds().getNumberOfBonds(), 1u);
    EXPECT_TRUE(network.getParticle(idx.c).isInSpringNeighbors(idx.n));
    EXPECT_FALSE(network.areHydrogenBonded(idx.c, idx.n))
        << "the pair is covalently bonded and still holds a hydrogen bond slot on each other";
    EXPECT_GE(network.getHydrogenBondEnergy(), -1e-6f)
        << "a hydrogen bond is still being evaluated on a pair that is now a peptide bond";
}

// What makes the product a peptide rather than two residues held at arm's
// length: the rest lengths are drawn to an ideal trans peptide plane.
TEST(PeptideBond, TheRampPullsTheBondToItsIdealLength)
{
    RulesFile rules;
    const unsigned RAMP = 50;
    spn::SpringNetwork network;
    configuration::Configuration config;
    const Pair idx = build(network, config, rules.path(), 2.9f, false, 1, RAMP, RAMP + 20);
    network.run();

    ASSERT_EQ(network.getPeptideBonds().getNumberOfBonds(), 1u);
    const auto & bond = network.getPeptideBonds().getBonds().front();

    // The C-N spring specifically: 1.329 A is Engh & Huber's peptide bond, and
    // it is the one distance a reader will check.
    int cn = -1;
    for (size_t k = 0; k < bond.springs.size(); ++k)
    {
        const spn::Spring & s = network.getSpring(bond.springs[k]);
        const unsigned a = s.getParticle1().getId();
        const unsigned b = s.getParticle2().getId();
        if ((a == idx.c && b == idx.n) || (a == idx.n && b == idx.c))
            cn = static_cast<int>(k);
    }
    ASSERT_GE(cn, 0) << "the rule created no spring between the two atoms whose bond this is";
    EXPECT_NEAR(network.getSpring(bond.springs[cn]).getEquilibrium(), 1.329f, 1e-3f)
        << "the ramp finished somewhere other than a peptide bond";

    // And the atoms followed their springs, which is the only thing that makes
    // the rest length mean anything.
    EXPECT_LT((network.getParticle(idx.c).getPosition() - network.getParticle(idx.n).getPosition()).norm(), 2.0f)
        << "the rest length reached 1.329 A but the two atoms stayed apart";
}
