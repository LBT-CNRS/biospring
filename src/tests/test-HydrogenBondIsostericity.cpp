#include <gtest/gtest.h>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>

#include "IO/DonorAcceptorRuleReader.h"
#include "configuration/Configuration.hpp"
#include "forcefield/constants.hpp"
#include "spn/SpringNetwork.h"

#include "TempPath.h"

using namespace biospring;

namespace
{

class Table
{
  public:
    explicit Table(const std::string & body)
    {
        _path = biospring::tests::temp_path(".hbond");
        std::ofstream(_path) << body;
    }
    ~Table() { std::remove(_path.c_str()); }
    const std::string & path() const { return _path; }

  private:
    std::string _path;
};

unsigned place(spn::SpringNetwork & network, const std::string & name, const std::string & resname, unsigned resid,
               const std::string & chain, const Vector3f & position)
{
    spn::Particle p;
    p.setName(name);
    p.setResName(resname);
    p.setChainName(chain);
    p.setResId(resid);
    p.setPosition(position);
    p.setMass(14.0f);
    p.setRadius(1.8f);
    p.setEpsilon(0.5f);
    network.addParticle(p);
    return network.getNumberOfParticles() - 1;
}

// Vector3f exposes named accessors rather than an index, and the tests sweep
// three axes, so one helper keeps the loops readable.
float component(const Vector3f & v, size_t axis) { return axis == 0 ? v.getX() : axis == 1 ? v.getY() : v.getZ(); }

Vector3f axisStep(size_t axis, float h)
{
    return Vector3f(axis == 0 ? h : 0.0f, axis == 1 ? h : 0.0f, axis == 2 ? h : 0.0f);
}

const float DEG = 0.01745329252f;
const float D0 = 10.46f, L0 = 54.9f, SD = 1.5f, SL = 18.0f;

// One nucleotide's four atoms, laid out so that its glycosidic bond makes a
// chosen angle with the C1'-C1' line, which is the x axis by construction.
//
//   C1' at (x, 0, 0); the glycosidic N one bond length away at `lambda` from
//   the direction `sense` * x; the base nitrogen that carries the bond put
//   beyond it, and one antecedent for that nitrogen so the site is directional.
struct Nucleotide
{
    unsigned c1, n9, site, antecedent;
};

Nucleotide placeNucleotide(spn::SpringNetwork & network, const std::string & resname, const std::string & glyco,
                           unsigned resid, const std::string & chain, float x, float sense, float lambda_deg,
                           float azimuth_deg, const Vector3f & site_position, const Vector3f & antecedent_position)
{
    const float a = lambda_deg * DEG, z = azimuth_deg * DEG;
    const Vector3f c1(x, 0.0f, 0.0f);
    // 1.47 A, a glycosidic bond. The azimuth turns the bond about the C1'-C1'
    // line, which leaves lambda untouched -- exactly the freedom the factor is
    // blind to and the pairing term is blind to as well.
    const Vector3f dir(sense * std::cos(a), std::sin(a) * std::cos(z), std::sin(a) * std::sin(z));
    Nucleotide nt{};
    nt.c1 = place(network, "C1'", resname, resid, chain, c1);
    nt.n9 = place(network, glyco, resname, resid, chain, c1 + dir * 1.47f);
    nt.site = place(network, "N3", resname, resid, chain, site_position);
    nt.antecedent = place(network, "C2", resname, resid, chain, antecedent_position);
    return nt;
}

// computeForces() ACCUMULATES into each particle: the reset lives in the
// integration step, not in the force evaluation, so a test that calls it twice
// reads the sum of both. Measured the hard way -- the force came out 9 and 16
// times the energy derivative, which is exactly how many evaluations had piled
// up by the time it was read.
float forcesFrom(spn::SpringNetwork & network)
{
    for (unsigned k = 0; k < network.getNumberOfParticles(); ++k)
        network.getParticle(k).resetForce();
    network.computeForces();
    return network.getHydrogenBondEnergy();
}

float gaussian(float v, float c, float w)
{
    const float z = (v - c) / w;
    return std::exp(-z * z);
}

} // namespace

// THE FRAME COMES FROM THE TOPOLOGY, not from a table. A purine hangs its base
// from N9 and a pyrimidine from N1; a residue carrying both an N9 and an N1 is
// a purine and must take N9, because N1 is a ring atom on the FAR side of the
// base and pointing lambda at it would measure a different angle entirely.
TEST(HydrogenBondIsostericity, TheFrameIsN9OnAPurineAndN1OnAPyrimidine)
{
    Table table("G N3 0 1 C2\n");
    spn::SpringNetwork network;
    // A guanine, with the N1 a purine also has.
    const unsigned g_c1 = place(network, "C1'", "G", 1, "A", {0.0f, 0.0f, 0.0f});
    const unsigned g_n9 = place(network, "N9", "G", 1, "A", {1.47f, 0.0f, 0.0f});
    place(network, "N1", "G", 1, "A", {4.0f, 2.0f, 0.0f});
    // A cytosine, which has an N1 and no N9.
    const unsigned c_c1 = place(network, "C1'", "C", 2, "A", {0.0f, 8.0f, 0.0f});
    const unsigned c_n1 = place(network, "N1", "C", 2, "A", {1.47f, 8.0f, 0.0f});
    // And a residue that is not a nucleotide at all.
    const unsigned ala = place(network, "CA", "ALA", 3, "A", {0.0f, 16.0f, 0.0f});

    configuration::Configuration config = configuration::defaultConfiguration();
    config.hbond.enable = true;
    config.sim.nbsteps = 1;
    // Tagged BEFORE setup, which is the order biospring-cli uses and the only
    // one that works: _setupHydrogenBond builds its slots from whoever is
    // already a donor or an acceptor, so a table read afterwards tags
    // particles that have no slot to hold a bond in.
    io::DonorAcceptorRuleReader reader(table.path());
    reader.read();
    reader.tagParticles(network);
    network.setup(config);

    EXPECT_TRUE(network.getParticle(g_c1).hasGlycosidicFrame());
    EXPECT_EQ(network.getParticle(g_c1).glycosidicC1Index(), static_cast<int>(g_c1));
    EXPECT_EQ(network.getParticle(g_c1).glycosidicNIndex(), static_cast<int>(g_n9))
        << "a purine took its N1 instead of its N9";
    EXPECT_EQ(network.getParticle(c_c1).glycosidicNIndex(), static_cast<int>(c_n1));
    EXPECT_FALSE(network.getParticle(ala).hasGlycosidicFrame())
        << "an amino acid was given a glycosidic frame";
}

// The factor is 1 AT the Watson-Crick geometry and nowhere near it elsewhere,
// and the numbers are checked against the formula computed independently here
// rather than against a remembered constant.
TEST(HydrogenBondIsostericity, TheFactorIsOneAtTheWatsonCrickGeometry)
{
    Table table("G N3 0 1 C2\nC N3 1 0 C2\n");
    spn::SpringNetwork network;
    const Nucleotide i =
        placeNucleotide(network, "G", "N9", 1, "A", 0.0f, +1.0f, L0, 0.0f, {4.0f, 1.0f, 0.0f}, {3.0f, 0.5f, 0.0f});
    const Nucleotide j = placeNucleotide(network, "C", "N1", 2, "B", D0, -1.0f, L0, 0.0f, {6.9f, 1.0f, 0.0f},
                                         {7.9f, 0.5f, 0.0f});

    configuration::Configuration config = configuration::defaultConfiguration();
    config.hbond.enable = true;
    config.hbond.isostericity = true;
    config.sim.nbsteps = 1;
    // Tagged BEFORE setup, which is the order biospring-cli uses and the only
    // one that works: _setupHydrogenBond builds its slots from whoever is
    // already a donor or an acceptor, so a table read afterwards tags
    // particles that have no slot to hold a bond in.
    io::DonorAcceptorRuleReader reader(table.path());
    reader.read();
    reader.tagParticles(network);
    network.setup(config);

    const auto iso =
        network.isostericityFactor(network.getParticle(i.site), network.getParticle(j.site), false);
    ASSERT_TRUE(iso.applies);
    EXPECT_NEAR(iso.value, 1.0f, 1e-5f) << "the centres are not where the geometry was built";

    // Now move one nucleotide's glycosidic angle to the measured MIRROR mean,
    // 36.6 degrees, leaving everything else alone. The expected value is the
    // formula, not a remembered number.
    const float mirror_lambda = 36.6f;
    const float a = mirror_lambda * DEG;
    network.getParticle(i.n9).setPosition(Vector3f(std::cos(a), std::sin(a), 0.0f) * 1.47f);
    const auto flipped =
        network.isostericityFactor(network.getParticle(i.site), network.getParticle(j.site), false);
    EXPECT_NEAR(flipped.value, gaussian(mirror_lambda, L0, SL), 1e-4f);
    EXPECT_LT(flipped.value, 0.4f);

    // The azimuth about the C1'-C1' line is exactly the freedom the factor does
    // not constrain, and saying so here is what keeps a future "improvement"
    // from quietly adding a fourth coordinate nobody measured.
    network.getParticle(i.n9).setPosition(
        Vector3f(std::cos(a), std::sin(a) * std::cos(1.0f), std::sin(a) * std::sin(1.0f)) * 1.47f);
    const auto turned =
        network.isostericityFactor(network.getParticle(i.site), network.getParticle(j.site), false);
    EXPECT_NEAR(turned.value, flipped.value, 1e-5f);
}

// THE GRADIENT, against finite differences of the factor itself. Four atoms,
// three directions each: the C1' of each nucleotide appears in all three
// coordinates of the factor -- it is the vertex of its own angle, one ray of
// the other's, and one end of the distance -- so it is the one whose gradient a
// missing term shows up in.
TEST(HydrogenBondIsostericity, TheGradientMatchesFiniteDifferences)
{
    Table table("G N3 0 1 C2\nC N3 1 0 C2\n");
    spn::SpringNetwork network;
    // Deliberately OFF the centres and off the plane, so every Gaussian has a
    // non-zero derivative and no symmetry can hide a sign.
    const Nucleotide i =
        placeNucleotide(network, "G", "N9", 1, "A", 0.3f, +1.0f, 48.0f, 23.0f, {4.0f, 1.0f, 0.4f}, {3.0f, 0.5f, 0.2f});
    const Nucleotide j = placeNucleotide(network, "C", "N1", 2, "B", 11.1f, -1.0f, 61.0f, -40.0f, {6.9f, 1.2f, -0.3f},
                                         {7.9f, 0.6f, -0.1f});

    configuration::Configuration config = configuration::defaultConfiguration();
    config.hbond.enable = true;
    config.hbond.isostericity = true;
    config.sim.nbsteps = 1;
    // Tagged BEFORE setup, which is the order biospring-cli uses and the only
    // one that works: _setupHydrogenBond builds its slots from whoever is
    // already a donor or an acceptor, so a table read afterwards tags
    // particles that have no slot to hold a bond in.
    io::DonorAcceptorRuleReader reader(table.path());
    reader.read();
    reader.tagParticles(network);
    network.setup(config);

    const auto analytic =
        network.isostericityFactor(network.getParticle(i.site), network.getParticle(j.site), true);
    ASSERT_TRUE(analytic.applies);
    ASSERT_GT(analytic.value, 1e-3f) << "the test geometry scores so low that nothing is being checked";

    const float h = 1e-4f;
    for (size_t k = 0; k < 4; ++k)
    {
        const size_t atom = static_cast<size_t>(analytic.atom[k]);
        const Vector3f origin = network.getParticle(atom).getPosition();
        for (size_t axis = 0; axis < 3; ++axis)
        {
            const Vector3f step = axisStep(axis, h);
            network.getParticle(atom).setPosition(origin + step);
            const float plus =
                network.isostericityFactor(network.getParticle(i.site), network.getParticle(j.site), false).value;
            network.getParticle(atom).setPosition(origin - step);
            const float minus =
                network.isostericityFactor(network.getParticle(i.site), network.getParticle(j.site), false).value;
            network.getParticle(atom).setPosition(origin);
            const float numerical = (plus - minus) / (2.0f * h);
            EXPECT_NEAR(component(analytic.grad[k], axis), numerical, 2e-3f)
                << "atom " << k << " (" << network.getParticle(atom).getName() << "), axis " << axis;
        }
    }
}

// AND THROUGH THE REAL FORCE PATH. The factor's force lands on four atoms that
// are not in the bond, and the donor and the acceptor do not move here, so the
// whole of the energy's change is the factor's: comparing the force on a frame
// atom against the finite difference of the TOTAL hydrogen bond energy checks
// the chain rule, the conversion and the accumulation at once.
//
// This is the shape that a stride mismatch broke silently on the GPU once
// before -- every energy stayed correct while third-party forces landed on the
// wrong atoms -- which is why it is checked here and not only at the factor.
TEST(HydrogenBondIsostericity, TheForceOnTheFrameMatchesTheEnergyDerivative)
{
    Table table("G N3 0 1 C2\nC N3 1 0 C2\n");
    spn::SpringNetwork network;
    const Nucleotide i =
        placeNucleotide(network, "G", "N9", 1, "A", 0.3f, +1.0f, 48.0f, 23.0f, {4.0f, 1.0f, 0.4f}, {3.0f, 0.5f, 0.2f});
    const Nucleotide j = placeNucleotide(network, "C", "N1", 2, "B", 11.1f, -1.0f, 61.0f, -40.0f, {6.9f, 1.2f, -0.3f},
                                         {7.9f, 0.6f, -0.1f});

    configuration::Configuration config = configuration::defaultConfiguration();
    config.hbond.enable = true;
    config.hbond.isostericity = true;
    config.hbond.cutoff = 7.0;
    config.spring.enable = false;
    config.steric.enable = false;
    config.electrostatic.enable = false;
    config.sim.nbsteps = 1;
    // Tagged BEFORE setup, which is the order biospring-cli uses and the only
    // one that works: _setupHydrogenBond builds its slots from whoever is
    // already a donor or an acceptor, so a table read afterwards tags
    // particles that have no slot to hold a bond in.
    io::DonorAcceptorRuleReader reader(table.path());
    reader.read();
    reader.tagParticles(network);
    network.setup(config);

    const auto frame = network.isostericityFactor(network.getParticle(i.site), network.getParticle(j.site), false);
    ASSERT_TRUE(frame.applies);

    auto energy = [&]() { return forcesFrom(network); };
    ASSERT_LT(energy(), -1e-3f) << "no bond formed, so there is no force to check";

    const float conv = static_cast<float>(forcefield::GLOBAL_SPRING_FORCE_CONVERT);
    const float h = 1e-3f;
    for (size_t k = 0; k < 4; ++k)
    {
        const size_t atom = static_cast<size_t>(frame.atom[k]);
        const Vector3f origin = network.getParticle(atom).getPosition();
        forcesFrom(network);
        const Vector3f force = network.getParticle(atom).getForce();
        for (size_t axis = 0; axis < 3; ++axis)
        {
            const Vector3f step = axisStep(axis, h);
            network.getParticle(atom).setPosition(origin + step);
            const float plus = energy();
            network.getParticle(atom).setPosition(origin - step);
            const float minus = energy();
            network.getParticle(atom).setPosition(origin);
            const float expected = -conv * (plus - minus) / (2.0f * h);
            EXPECT_NEAR(component(force, axis), expected, 1e-3f * std::max(1.0f, std::fabs(expected)))
                << "frame atom " << k << " (" << network.getParticle(atom).getName() << "), axis " << axis;
        }
    }
}

// A BOND WITH NO FRAME IS UNTOUCHED TO THE BIT, not within a tolerance. Every
// protein hydrogen bond in every existing .msp goes through this path once the
// term is compiled in, and "almost the same" is how a trajectory drifts without
// anyone being able to say when it started.
TEST(HydrogenBondIsostericity, AProteinBondIsBitIdenticalWithTheTermOnOrOff)
{
    Table table("ALA N 1 0 CA\nALA O 0 2 C\n");
    auto build = [&](bool isostericity, float & energy, Vector3f * forces) {
        spn::SpringNetwork network;
        // The antecedent sits BEHIND its own site: the donor's hydrogen points
        // from CA through N towards the acceptor, and the acceptor's lone pair
        // from C through O back towards the donor. Placed the other way round
        // the angular weight is exactly zero and no bond forms at all, which is
        // how this test first read as "nothing to compare".
        place(network, "N", "ALA", 1, "A", {0.0f, 0.0f, 0.0f});
        place(network, "CA", "ALA", 1, "A", {-0.2f, -1.46f, -0.1f});
        place(network, "O", "ALA", 2, "A", {0.1f, 2.9f, 0.3f});
        place(network, "C", "ALA", 2, "A", {0.2f, 4.1f, 0.4f});

        configuration::Configuration config = configuration::defaultConfiguration();
        config.hbond.enable = true;
        config.hbond.cutoff = 7.0;
        config.hbond.isostericity = isostericity;
        config.spring.enable = false;
        config.steric.enable = false;
        config.electrostatic.enable = false;
        config.sim.nbsteps = 1;
        // Tagged BEFORE setup, which is the order biospring-cli uses and the only
        // one that works: _setupHydrogenBond builds its slots from whoever is
        // already a donor or an acceptor, so a table read afterwards tags
        // particles that have no slot to hold a bond in.
        io::DonorAcceptorRuleReader reader(table.path());
        reader.read();
        reader.tagParticles(network);
        network.setup(config);
        energy = forcesFrom(network);
        for (unsigned k = 0; k < 4; ++k)
            forces[k] = network.getParticle(k).getForce();
    };

    float off = 0.0f, on = 0.0f;
    Vector3f foff[4], fon[4];
    build(false, off, foff);
    build(true, on, fon);

    ASSERT_LT(off, -1e-3f) << "no bond formed, so nothing is being compared";
    EXPECT_EQ(std::memcmp(&off, &on, sizeof(float)), 0) << "the energy of a protein bond moved: " << off << " -> " << on;
    for (unsigned k = 0; k < 4; ++k)
        EXPECT_EQ(std::memcmp(&foff[k], &fon[k], sizeof(Vector3f)), 0)
            << "the force on particle " << k << " moved";
}

// THE GLYCOSIDIC PARTNER IS CHOSEN BY BOND LENGTH, which is what makes
// pseudouridine work. Psi is a C-GLYCOSIDE: its sugar is bonded to the base's
// C5, and its N1 -- which a name-ordered rule picks -- is a ring atom 3.98 A
// away, measured on 8PKL's A-site tRNA residue 55. A tRNA carries several, so
// the anticodon arm's own geometry depended on this.
TEST(HydrogenBondIsostericity, PseudouridineTakesItsCGlycosidicBondAndNotItsDistantN1)
{
    Table table("PSU N3 1 0 C2\n");
    spn::SpringNetwork network;
    const unsigned c1 = place(network, "C1'", "PSU", 55, "A", {0.0f, 0.0f, 0.0f});
    // The real distances from 8PKL: C5 bonded at 1.56 A, N1 a ring atom at 3.98.
    const unsigned c5 = place(network, "C5", "PSU", 55, "A", {1.56f, 0.0f, 0.0f});
    place(network, "N1", "PSU", 55, "A", {3.60f, 1.70f, 0.0f});
    place(network, "N3", "PSU", 55, "A", {4.60f, 0.40f, 0.0f});
    place(network, "C2", "PSU", 55, "A", {5.60f, 0.90f, 0.0f});

    configuration::Configuration config = configuration::defaultConfiguration();
    config.hbond.enable = true;
    config.hbond.isostericity = true;
    config.sim.nbsteps = 1;
    io::DonorAcceptorRuleReader reader(table.path());
    reader.read();
    reader.tagParticles(network);
    network.setup(config);

    EXPECT_TRUE(network.getParticle(c1).hasGlycosidicFrame());
    EXPECT_EQ(network.getParticle(c1).glycosidicNIndex(), static_cast<int>(c5))
        << "pseudouridine was given its N1, which its sugar is not bonded to";
}

// A SUGAR OR PHOSPHATE SITE IS NOT A BASE SITE. The shipped tables declare
// rOP1, rOP2, rO5, rO4 and rO3 as acceptors and rO2 as both, and the glycosidic
// frame is written onto every particle of the residue, so without this test the
// factor would score a backbone contact against a centre measured only for
// base-base pairs. Measured on 074.DNADuplex before the test existed: 24 of 77
// bonds were between sequence neighbours, mostly through sugar oxygens, and the
// factor cut them to between 3e-11 and 6e-6.
TEST(HydrogenBondIsostericity, ABackboneContactIsNotScoredAgainstABasePairCentre)
{
    Table table("G N3 0 1 C2\nC O2' 1 0 C2'\n");
    spn::SpringNetwork network;
    const Nucleotide i =
        placeNucleotide(network, "G", "N9", 1, "A", 0.0f, +1.0f, L0, 0.0f, {4.0f, 1.0f, 0.0f}, {3.0f, 0.5f, 0.0f});
    // A second nucleotide whose SITE is its 2'-OH oxygen, not a base atom.
    const float a = L0 * DEG;
    place(network, "C1'", "C", 2, "B", {D0, 0.0f, 0.0f});
    place(network, "N1", "C", 2, "B",
          {D0 - std::cos(a) * 1.47f, -std::sin(a) * 1.47f, 0.0f});
    const unsigned o2 = place(network, "O2'", "C", 2, "B", {6.9f, 1.0f, 0.0f});
    place(network, "C2'", "C", 2, "B", {7.9f, 0.5f, 0.0f});

    configuration::Configuration config = configuration::defaultConfiguration();
    config.hbond.enable = true;
    config.hbond.isostericity = true;
    config.sim.nbsteps = 1;
    io::DonorAcceptorRuleReader reader(table.path());
    reader.read();
    reader.tagParticles(network);
    network.setup(config);

    // Both residues HAVE a frame -- that is the point; it is the ATOM that
    // disqualifies the bond, not the residue.
    ASSERT_TRUE(network.getParticle(o2).hasGlycosidicFrame());
    const auto iso =
        network.isostericityFactor(network.getParticle(i.site), network.getParticle(o2), false);
    EXPECT_FALSE(iso.applies) << "a 2'-OH contact was scored as a base pair";
    EXPECT_FLOAT_EQ(iso.value, 1.0f);
}

// A FACTOR BELOW THE FLOOR IS EXACTLY ZERO, so the candidate ranking can refuse
// the bond. The ranking tests strength < 0 against a nearest-so-far of 0, and a
// Gaussian is never exactly zero: before the floor, a bond the term had cut to
// 1e-11 still won a free donor slot and held it. Measured: 8 such bonds on
// 074.DNADuplex, occupying slots for -2.3e-6 kJ.mol-1 in total.
TEST(HydrogenBondIsostericity, AFactorBelowTheFloorIsExactlyZeroSoNoSlotIsTaken)
{
    Table table("G N3 0 1 C2\nC N3 1 0 C2\n");
    spn::SpringNetwork network;
    // lambda_i driven far off the centre: at 150 degrees against a centre of
    // 54.9 and a width of 18 the acceptance is e^(-27.9).
    const Nucleotide i =
        placeNucleotide(network, "G", "N9", 1, "A", 0.0f, +1.0f, 150.0f, 0.0f, {4.0f, 1.0f, 0.0f}, {3.0f, 0.5f, 0.0f});
    const Nucleotide j = placeNucleotide(network, "C", "N1", 2, "B", D0, -1.0f, L0, 0.0f, {6.9f, 1.0f, 0.0f},
                                         {7.9f, 0.5f, 0.0f});

    configuration::Configuration config = configuration::defaultConfiguration();
    config.hbond.enable = true;
    config.hbond.isostericity = true;
    config.sim.nbsteps = 1;
    io::DonorAcceptorRuleReader reader(table.path());
    reader.read();
    reader.tagParticles(network);
    network.setup(config);

    const auto iso =
        network.isostericityFactor(network.getParticle(i.site), network.getParticle(j.site), true);
    ASSERT_TRUE(iso.applies) << "the factor must still claim the bond, it is the VALUE that is zero";
    EXPECT_EQ(iso.value, 0.0f) << "a factor under the floor must be exactly zero, not merely small";
    for (size_t k = 0; k < 4; ++k)
    {
        EXPECT_FLOAT_EQ(iso.grad[k].getX(), 0.0f);
        EXPECT_FLOAT_EQ(iso.grad[k].getY(), 0.0f);
        EXPECT_FLOAT_EQ(iso.grad[k].getZ(), 0.0f);
    }
}

// THE DISTANCE ACCEPTANCE IS A FLOOR, NOT A WINDOW, and this is the test that
// says why. Two bases on their way to pairing are FURTHER apart than the
// Watson-Crick 10.46 A, so a Gaussian on both sides switches the well off
// exactly during the approach -- 0.157 of it left at a C1'-C1' of 12.5 A, 7e-4
// at 14.5, and zero at the 18.6 A a free tRNA starts from. The term could then
// only judge a pair already formed and could never help one form.
//
// Costing nothing, because the parasites it exists to reject are SHORTER than
// a real pair: the +1 cross-strand diagonal measures 7.85 to 8.83 A against
// the native 10.5, and the 180-degree mirror is caught by lambda, not by
// distance. Measured on those populations the discrimination goes 7.3x -> 7.0x
// while the approach goes 5.3 % -> 100 %.
TEST(HydrogenBondIsostericity, BeyondTheCentreTheDistanceCostsNothing)
{
    Table table("G N3 0 1 C2\nC N3 1 0 C2\n");

    auto factorAt = [&](float separation, bool twosided) {
        spn::SpringNetwork network;
        const Nucleotide i = placeNucleotide(network, "G", "N9", 1, "A", 0.0f, +1.0f, L0, 0.0f,
                                             {4.0f, 1.0f, 0.0f}, {3.0f, 0.5f, 0.0f});
        const Nucleotide j = placeNucleotide(network, "C", "N1", 2, "B", separation, -1.0f, L0, 0.0f,
                                             {separation - 3.5f, 1.0f, 0.0f}, {separation - 2.5f, 0.5f, 0.0f});
        configuration::Configuration config = configuration::defaultConfiguration();
        config.hbond.enable = true;
        config.hbond.isostericity = true;
        config.hbond.isostericitytwosided = twosided;
        config.sim.nbsteps = 1;
        io::DonorAcceptorRuleReader reader(table.path());
        reader.read();
        reader.tagParticles(network);
        network.setup(config);
        const auto iso = network.isostericityFactor(network.getParticle(i.site), network.getParticle(j.site), true);
        EXPECT_TRUE(iso.applies);
        return iso;
    };

    // At the centre both forms agree, and both are 1.
    EXPECT_NEAR(factorAt(D0, false).value, 1.0f, 1e-5f);
    EXPECT_NEAR(factorAt(D0, true).value, 1.0f, 1e-5f);

    // Beyond it the floor keeps the full weight where the window throws it away.
    for (float d : {11.5f, 12.5f, 14.5f, 18.0f})
    {
        const auto floored = factorAt(d, false);
        const auto window = factorAt(d, true);
        EXPECT_NEAR(floored.value, 1.0f, 1e-5f) << "the floor taxed an approach at " << d << " A";
        EXPECT_LT(window.value, floored.value) << "the window did not tax it at " << d << " A";
        // And the floor's gradient vanishes there too: a flat acceptance must
        // exert no force, or it would pull on the frame for no reason.
        // And the DISTANCE part of the gradient vanishes there, so a flat
        // acceptance exerts no pull. What is left is not exactly zero but
        // float noise: this geometry puts lambda exactly at its centre, where
        // the lambda derivative is zero, and acos plus the normalisations
        // leave it a few 1e-5 degrees off. Measured, the residue runs 1e-7 to
        // 6e-6 against a factor of 1, which is why the bound is 1e-5 and not 0.
        for (size_t k = 0; k < 4; ++k)
            EXPECT_LT(floored.grad[k].norm(), 1e-5f)
                << "a flat acceptance still pushed, at " << d << " A, atom " << k;
    }

    // And SHORTER than the centre both forms penalise, which is what rejects
    // the cross-strand diagonal.
    EXPECT_LT(factorAt(8.5f, false).value, 0.3f) << "a too-short pair was not penalised";
    EXPECT_NEAR(factorAt(8.5f, false).value, factorAt(8.5f, true).value, 1e-5f)
        << "the two forms must agree below the centre";
}
