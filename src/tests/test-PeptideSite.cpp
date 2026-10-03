#include <gtest/gtest.h>

#include <cstdio>
#include <fstream>
#include <string>

#include "IO/PeptideSiteRuleReader.h"
#include "spn/Particle.h"

// HBOND_LOBES_* are static const ints with no out-of-line definition, so they
// are usable by value but not by reference -- and EXPECT_EQ binds both of its
// arguments to a const reference. Hence the static_cast on each one.

using namespace biospring;

namespace
{

// A .psite file that lives as long as the test does.
class SiteFile
{
  public:
    explicit SiteFile(const std::string & body)
    {
        _path = std::string(std::tmpnam(nullptr)) + ".psite";
        std::ofstream(_path) << body;
    }
    ~SiteFile() { std::remove(_path.c_str()); }
    const std::string & path() const { return _path; }

  private:
    std::string _path;
};

} // namespace

// The three columns a reactive site needs, and the fourth that splits it in two.
TEST(PeptideSite, ReadsTheAntecedentAndTheOutOfPlaneAngle)
{
    SiteFile file("# a comment, and a blank line\n"
                  "\n"
                  "ALA  C  O  ^CA:75\n"
                  "GLY  C  O  ^CA:105\n");
    io::PeptideSiteRuleReader reader(file.path());
    reader.read();

    ASSERT_EQ(reader.rules().size(), 2u);
    const io::PeptideSite & ala = reader.rules().at({"ALA", "C"});
    EXPECT_EQ(ala.antecedent, "O");
    EXPECT_EQ(ala.antecedent2, "CA");
    EXPECT_EQ(ala.lobeMode, static_cast<int>(spn::Particle::HBOND_LOBES_OUT_OF_PLANE));
    EXPECT_FLOAT_EQ(ala.lobeAngle, 75.0f);
    EXPECT_FLOAT_EQ(reader.rules().at({"GLY", "C"}).lobeAngle, 105.0f);
}

// Three columns is a complete site: a direction with no faces. Used on its own
// it means "anywhere around this axis", which is right for a site that is not
// sp2.
TEST(PeptideSite, AcceptsASiteWithNoPlaneAtAll)
{
    SiteFile file("ALA  C  O\n");
    io::PeptideSiteRuleReader reader(file.path());
    reader.read();
    const io::PeptideSite & site = reader.rules().at({"ALA", "C"});
    EXPECT_EQ(site.antecedent, "O");
    EXPECT_TRUE(site.antecedent2.empty());
    EXPECT_EQ(site.lobeMode, static_cast<int>(spn::Particle::HBOND_LOBES_NONE));
}

// A plain fourth column is a second BOND, not a plane: the direction becomes
// the opposite of the two antecedents' bisector and there are no lobes.
TEST(PeptideSite, APlainFourthColumnIsASecondAntecedentAndNotAPlane)
{
    SiteFile file("ALA  C  O  CA\n");
    io::PeptideSiteRuleReader reader(file.path());
    reader.read();
    const io::PeptideSite & site = reader.rules().at({"ALA", "C"});
    EXPECT_EQ(site.antecedent2, "CA");
    EXPECT_EQ(site.lobeMode, static_cast<int>(spn::Particle::HBOND_LOBES_NONE));
}

// '~' is the in-plane form and keeps 62 degrees, so a table written before the
// angle existed reads the same. This is the same code the .hbond table uses --
// io::parse_lobe_spec -- which is why the two files cannot drift apart on it.
TEST(PeptideSite, TheInPlaneFormKeepsItsHistoricalDefaultAngle)
{
    SiteFile file("ASN  ND2  CG  ~OD1\n");
    io::PeptideSiteRuleReader reader(file.path());
    reader.read();
    const io::PeptideSite & site = reader.rules().at({"ASN", "ND2"});
    EXPECT_EQ(site.lobeMode, static_cast<int>(spn::Particle::HBOND_LOBES_IN_PLANE));
    EXPECT_FLOAT_EQ(site.lobeAngle, 62.0f);
}

// An out-of-plane approach has no canonical angle, so the file has to say.
// Silently defaulting would put the lobes at 62 degrees -- 118 off the C=O axis
// instead of 105 -- and nothing would report it.
TEST(PeptideSite, RefusesAnOutOfPlaneSiteWithNoAngle)
{
    SiteFile file("ALA  C  O  ^CA\n");
    io::PeptideSiteRuleReader reader(file.path());
    EXPECT_DEATH(reader.read(), "needs an explicit angle");
}

TEST(PeptideSite, RefusesAnAngleOutsideTheOpenInterval)
{
    SiteFile file("ALA  C  O  ^CA:180\n");
    io::PeptideSiteRuleReader reader(file.path());
    EXPECT_DEATH(reader.read(), "outside .0, 180.");
}

// The antecedent is not optional here, unlike in the .hbond table where an
// undirected acceptor still has a distance-only meaning. A reactive site with
// no direction cannot be attacked from anywhere in particular.
TEST(PeptideSite, RefusesALineWithNoAntecedent)
{
    SiteFile file("ALA  C\n");
    io::PeptideSiteRuleReader reader(file.path());
    EXPECT_DEATH(reader.read(), "invalid number of tokens");
}

// Two lines for the same atom would leave which one wins to map ordering.
TEST(PeptideSite, RefusesTheSameSiteDeclaredTwice)
{
    SiteFile file("ALA  C  O  ^CA:75\n"
                  "ALA  C  O  ^CA:105\n");
    io::PeptideSiteRuleReader reader(file.path());
    EXPECT_DEATH(reader.read(), "declared twice");
}

// A whole protein table in one read, in the shape the shipped
// data/reducerules/ProteinPeptideSite.psite has. Written out here rather than
// read from data/, following test-PeptideBond.cpp: a test says in one place
// what it depends on. What guards the SHIPPED file is the run itself, which
// reports "N of M C atoms carry no out-of-plane lobe site" -- a missing entry
// does not make the model stricter, it removes the angular factor altogether.
TEST(PeptideSite, ReadsAWholeProteinTableAndKeepsEveryResidueDistinct)
{
    std::string body = "# twenty residues, one site each\n";
    const char * residues[] = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
                               "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"};
    for (const char * res : residues)
        body += std::string(res) + "  C  O  ^CA:75\n";

    SiteFile file(body);
    io::PeptideSiteRuleReader reader(file.path());
    reader.read();

    ASSERT_EQ(reader.rules().size(), 20u);
    for (const char * res : residues)
    {
        const auto it = reader.rules().find({res, "C"});
        ASSERT_NE(it, reader.rules().end()) << res << " has no reactive site";
        EXPECT_EQ(it->second.antecedent, "O") << res;
        EXPECT_EQ(it->second.antecedent2, "CA") << res;
        EXPECT_EQ(it->second.lobeMode, static_cast<int>(spn::Particle::HBOND_LOBES_OUT_OF_PLANE)) << res;
        EXPECT_FLOAT_EQ(it->second.lobeAngle, 75.0f) << res;
    }
}

