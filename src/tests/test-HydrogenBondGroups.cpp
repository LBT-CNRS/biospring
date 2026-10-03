#include <gtest/gtest.h>

#include <cstdio>
#include <fstream>
#include <string>

#include "IO/DonorAcceptorRuleReader.h"
#include "configuration/Configuration.hpp"
#include "spn/SpringNetwork.h"

using namespace biospring;

namespace
{

// A donor/acceptor table that declares one group and puts a couple in it.
class Table
{
  public:
    explicit Table(const std::string & body)
    {
        _path = std::string(std::tmpnam(nullptr)) + ".hbond";
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
    p.setMass(14.0f);
    p.setRadius(1.8f);
    p.setEpsilon(0.5f);
    network.addParticle(p);
    return network.getNumberOfParticles() - 1;
}

} // namespace

// A PAIR carries the parameters, and only when BOTH ends opt in.
//
// That condition is the whole design and not a convenience. A peptide bond's
// well is two orders of magnitude deeper than a hydrogen bond's and sits at
// 1.33 A rather than 2.90; a carbonyl carbon carrying it alone would drag every
// amine and every hydroxyl within reach down to a covalent distance. It is the
// couple that reacts, so it is the couple that has to say so.
TEST(HydrogenBondGroups, OnlyAPairWhoseTwoEndsAgreeLeavesTheForceFieldsWell)
{
    Table table("GROUP peptide 400.0 1.329 2.0\n"
                "ALA N  1 0 CA @peptide\n"
                "ALA C  0 1 O  @peptide\n"
                "GLY N  1 0 CA\n"
                "GLY O  0 2 C\n");

    spn::SpringNetwork network;
    const unsigned n_pep = place(network, "N", "ALA", 1, {0.0f, 0.0f, 0.0f});
    place(network, "CA", "ALA", 1, {1.46f, 0.0f, 0.0f});
    const unsigned c_pep = place(network, "C", "ALA", 2, {0.0f, 3.0f, 0.0f});
    place(network, "O", "ALA", 2, {0.0f, 4.2f, 0.0f});
    const unsigned n_ord = place(network, "N", "GLY", 3, {8.0f, 0.0f, 0.0f});
    place(network, "CA", "GLY", 3, {9.46f, 0.0f, 0.0f});
    const unsigned o_ord = place(network, "O", "GLY", 4, {8.0f, 3.0f, 0.0f});
    place(network, "C", "GLY", 4, {8.0f, 4.2f, 0.0f});

    configuration::Configuration config = configuration::defaultConfiguration();
    config.hbond.enable = true;
    config.hbond.cutoff = 7.0;
    config.sim.nbsteps = 1;
    network.setup(config);

    io::DonorAcceptorRuleReader reader(table.path());
    reader.read();
    reader.tagParticles(network);

    ASSERT_EQ(reader.groups().size(), 1u);
    EXPECT_EQ(reader.groups().front().name, "peptide");

    const auto pep = network.hydrogenBondParameters(n_pep, c_pep);
    EXPECT_FLOAT_EQ(pep.welldepth, 400.0f);
    EXPECT_FLOAT_EQ(pep.equilibrium, 1.329f);
    EXPECT_FLOAT_EQ(pep.width, 2.0f);

    // Neither end in a group: the force field's own, whatever they are.
    const auto ordinary = network.hydrogenBondParameters(n_ord, o_ord);
    EXPECT_FLOAT_EQ(ordinary.welldepth, network.getForceField()->getHydrogenBondWellDepth());
    EXPECT_FLOAT_EQ(ordinary.equilibrium, network.getForceField()->getHydrogenBondEquilibrium());

    // AND THE CASE THAT MATTERS: one end in the group, the other not. If this
    // took the group's well, the peptide site would pull every ordinary donor
    // in the structure to 1.33 A.
    const auto mixed = network.hydrogenBondParameters(n_ord, c_pep);
    EXPECT_FLOAT_EQ(mixed.welldepth, network.getForceField()->getHydrogenBondWellDepth())
        << "a site dragged a partner that never opted in into its own well";
    EXPECT_FLOAT_EQ(mixed.equilibrium, network.getForceField()->getHydrogenBondEquilibrium());
}

// A table written before groups existed must read exactly the same.
TEST(HydrogenBondGroups, ATableWithNoGroupIsUnchanged)
{
    Table table("ALA N  1 0 CA\n"
                "ALA O  0 2 C\n");
    spn::SpringNetwork network;
    const unsigned a = place(network, "N", "ALA", 1, {0.0f, 0.0f, 0.0f});
    place(network, "CA", "ALA", 1, {1.46f, 0.0f, 0.0f});
    const unsigned b = place(network, "O", "ALA", 2, {0.0f, 3.0f, 0.0f});
    place(network, "C", "ALA", 2, {0.0f, 4.2f, 0.0f});

    configuration::Configuration config = configuration::defaultConfiguration();
    config.hbond.enable = true;
    config.hbond.cutoff = 7.0;
    config.sim.nbsteps = 1;
    network.setup(config);

    io::DonorAcceptorRuleReader reader(table.path());
    reader.read();
    reader.tagParticles(network);

    EXPECT_TRUE(reader.groups().empty());
    EXPECT_EQ(network.getParticle(a).hydrogenBondGroup(), 0);
    const auto p = network.hydrogenBondParameters(a, b);
    EXPECT_FLOAT_EQ(p.welldepth, network.getForceField()->getHydrogenBondWellDepth());
}

// A table in plain PDB names resolves against a RENAMED topology with no
// run-time setting at all, because each particle carries the name it had
// BEFORE the .grp renamed it -- recorded by pdb2spn, which is the only code
// that had the .grp in front of it.
//
// Asking for the .grp again at run time was opaque: the user had to remember
// which one a topology was reduced with, and remembering the wrong one matched
// nothing and reported nothing. Guessing instead -- dropping the first
// character of the type -- was worse, and silently wrong for any prefix that
// is not one letter.
TEST(HydrogenBondGroups, APlainTableMatchesARenamedTopologyThroughTheRecordedName)
{
    Table table("ALA  N  1 0 CA\n");

    spn::SpringNetwork network;
    const unsigned n = place(network, "AN", "ALA", 1, Vector3f(0.0f, 0.0f, 0.0f));
    const unsigned ca = place(network, "ACA", "ALA", 1, Vector3f(1.5f, 0.0f, 0.0f));
    network.getParticle(n).setOriginalName("N");
    network.getParticle(ca).setOriginalName("CA");

    io::DonorAcceptorRuleReader reader(table.path());
    reader.read();
    reader.tagParticles(network);

    EXPECT_EQ(network.getParticle(n).donorCapacity(), 1u) << "'N' has to find the particle named 'AN'";
    EXPECT_EQ(network.getParticle(n).antecedentIndex(), static_cast<int>(ca))
        << "and 'CA' has to find the particle named 'ACA'";
}

// And with nothing recorded, nothing is guessed. The old code dropped the
// first character of every type, which made this pass by accident; it now
// fails to match, which is the honest answer and is reported.
TEST(HydrogenBondGroups, WithNoRecordedNameAPlainTableMatchesNothing)
{
    Table table("ALA  N  1 0 CA\n");

    spn::SpringNetwork network;
    const unsigned n = place(network, "AN", "ALA", 1, Vector3f(0.0f, 0.0f, 0.0f));
    place(network, "ACA", "ALA", 1, Vector3f(1.5f, 0.0f, 0.0f));

    io::DonorAcceptorRuleReader reader(table.path());
    reader.read();
    reader.tagParticles(network);

    EXPECT_EQ(network.getParticle(n).donorCapacity(), 0u)
        << "no recorded name, so 'N' cannot be known to mean 'AN'";
}
