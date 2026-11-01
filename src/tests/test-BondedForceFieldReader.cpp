#include <gtest/gtest.h>

#include "IO/BondedForceFieldReader.h"
#include "topology/Topology.hpp"

#include <algorithm>
#include <fstream>
#include <string>
#include <vector>

#include "TempPath.h"

using namespace biospring;
using biospring::rigidbodygroup::BondedForceFieldReader;

namespace
{

// Two tables and one record, the smallest file the format admits. bins = 2
// means 3 samples, so 6 interleaved energy/torque numbers.
const char * const VALID = "TORSIONTABLE 0 2 0.0 0.0 1.0 0.0 0.0 0.0\n"
                           "TORSIONTABLE 1 2 0.0 0.0 2.0 0.0 0.0 0.0\n"
                           "TORSION\tARG\tSIDECHAIN\tNE\tCZ\tNH1\tHH11\t1\n";

std::string write(const std::string & content)
{
    const std::string path = biospring::tests::temp_path(".bi.ff");
    std::ofstream out(path);
    out << content;
    out.close();
    return path;
}

} // namespace

TEST(BondedForceFieldReader, DistinctTableIdsParse)
{
    BondedForceFieldReader reader(write(VALID));
    reader.read();  // no death
}

// The defect this guards. Each generated .bi.ff numbers its tables from 0, so
// concatenating the protein file (32 tables, 0..31) with the RNA one (39) gives
// 32 collisions, and before this the later entry silently replaced the earlier
// -- which on 078/079 made all 894 protein torsions read a nucleic table.
TEST(BondedForceFieldReader, ADuplicateTableIdIsRefused)
{
    const std::string content = std::string(VALID) + "TORSIONTABLE 1 2 0.0 0.0 9.0 0.0 0.0 0.0\n";
    BondedForceFieldReader reader(write(content));
    EXPECT_DEATH(reader.read(), "TORSIONTABLE 1 is declared twice");
}

// The id does not have to be the one just read, and the colliding file does not
// have to be adjacent: a concatenation collides on whatever ids the shorter
// file used, wherever they fall.
TEST(BondedForceFieldReader, ADuplicateIsCaughtOutOfOrder)
{
    const std::string content = std::string(VALID) + "TORSIONTABLE 5 2 0.0 0.0 3.0 0.0 0.0 0.0\n"
                                                     "TORSIONTABLE 0 2 0.0 0.0 9.0 0.0 0.0 0.0\n";
    BondedForceFieldReader reader(write(content));
    EXPECT_DEATH(reader.read(), "TORSIONTABLE 0 is declared twice");
}

// And the fix the message asks for has to work: a union whose second half is
// renumbered is a legitimate file.
TEST(BondedForceFieldReader, ARenumberedUnionParses)
{
    const std::string content = std::string(VALID) + "TORSIONTABLE 2 2 0.0 0.0 3.0 0.0 0.0 0.0\n"
                                                     "TORSIONTABLE 3 2 0.0 0.0 4.0 0.0 0.0 0.0\n"
                                                     "TORSION\tA\tNUCLEIC_CHI\tO4'\tC1'\tN9\tC4\t3\n";
    BondedForceFieldReader reader(write(content));
    reader.read();  // no death
}

// Sparse ids are not a collision: the reader grows the vector and leaves the
// holes at bins = 0, which is how it tells "never declared" from "declared".
TEST(BondedForceFieldReader, SparseIdsAreNotADuplicate)
{
    const std::string content = "TORSIONTABLE 7 2 0.0 0.0 1.0 0.0 0.0 0.0\n"
                                "TORSIONTABLE 3 2 0.0 0.0 1.0 0.0 0.0 0.0\n"
                                "TORSION\tARG\tSIDECHAIN\tNE\tCZ\tNH1\tHH11\t7\n";
    BondedForceFieldReader reader(write(content));
    reader.read();  // no death
}

// Several files into one reader: this is how a mixed system is meant to be
// parametrised, and the whole point is that each file keeps numbering its own
// tables from 0.
namespace
{

// The resolved tables are what matters and they are protected, so expose them.
struct Probe : public BondedForceFieldReader
{
    size_t tables() const { return _torsiontables.size(); }
    size_t torsions() const { return _torsion.size(); }
    unsigned tableOf(size_t i) const { return _torsion[i].table; }
    float peakOf(size_t i) const
    {
        const auto & e = _torsiontables[_torsion[i].table].energy;
        return *std::max_element(e.begin(), e.end());
    }
};

// Same ids as VALID (0 and 1), different peaks, so a record that ends up on
// the wrong file's table is visible in the energy rather than only in an index.
const char * const SECOND = "TORSIONTABLE 0 2 0.0 0.0 7.0 0.0 0.0 0.0\n"
                            "TORSIONTABLE 1 2 0.0 0.0 8.0 0.0 0.0 0.0\n"
                            "TORSION\tA\tNUCLEIC_CHI\tO4'\tC1'\tN9\tC4\t1\n";

} // namespace

TEST(BondedForceFieldReader, TwoFilesKeepTheirOwnTables)
{
    Probe reader;
    reader.setFileName(write(VALID));
    reader.read();
    reader.setFileName(write(SECOND));
    reader.read();

    EXPECT_EQ(reader.tables(), 4u);     // 2 + 2, nothing overwritten
    EXPECT_EQ(reader.torsions(), 2u);
    EXPECT_EQ(reader.tableOf(0), 1u);   // first file's table 1, unshifted
    EXPECT_EQ(reader.tableOf(1), 3u);   // second file's table 1, shifted by 2
    EXPECT_FLOAT_EQ(reader.peakOf(0), 2.0f);
    EXPECT_FLOAT_EQ(reader.peakOf(1), 8.0f);
}

// Reading the files the other way round must give each record the same table
// contents -- otherwise the order on the command line would change the physics.
TEST(BondedForceFieldReader, TheOrderOfTheFilesDoesNotChangeTheParameters)
{
    Probe reader;
    reader.setFileName(write(SECOND));
    reader.read();
    reader.setFileName(write(VALID));
    reader.read();

    EXPECT_FLOAT_EQ(reader.peakOf(0), 8.0f);   // the nucleic record, still 8
    EXPECT_FLOAT_EQ(reader.peakOf(1), 2.0f);   // the protein record, still 2
}

// A single file is unchanged: the shift is 0 for the first one read.
TEST(BondedForceFieldReader, OneFileIsUnshifted)
{
    Probe reader;
    reader.setFileName(write(VALID));
    reader.read();
    EXPECT_EQ(reader.tableOf(0), 1u);
}

// With a shift in play, a record naming a table its file never declared would
// reach into a neighbouring file's tables instead of merely being empty.
TEST(BondedForceFieldReader, ARecordNamingAnUndeclaredTableIsRefused)
{
    const std::string content = "TORSIONTABLE 0 2 0.0 0.0 1.0 0.0 0.0 0.0\n"
                                "TORSION\tARG\tSIDECHAIN\tNE\tCZ\tNH1\tHH11\t4\n";
    BondedForceFieldReader reader(write(content));
    EXPECT_DEATH(reader.read(), "names table 4, which the file never declares");
}

// ---------------------------------------------------------------------------
// Terminal qualifiers.
//
// AMBER types a free chain end from its own residue template, so the records
// that apply there are not the interior ones: amber99sb gives N-CT-C-N three
// harmonics and N3-CT-C-N none, and the ammonium rotor and the carboxylate
// improper exist only at a terminus. A terminal residue therefore takes the
// :NTER/:CTER records INSTEAD of the unqualified ones.
//
// What is NOT testable from here, because the reader cannot see it: a torsion
// spanning two residues can be spelled from either of them, so the interior
// psi written as "-N -CA -C N" under the NEXT residue would re-create exactly
// the quadruplet the :NTER set omits. Nothing in a record says what its
// neighbour must be, so that copy must not exist in the first place --
// scripts/amber_tabulate.py's owners_for() drops it, and asserts afterwards
// that no interior record reaches for a terminal-sensitive neighbour atom.
namespace
{

// bins = 2 again, and one peak per table so a torsion that ends up on the
// wrong one is visible in the energy and not only in an index.
const char * const TERMINAL_TABLES = "TORSIONTABLE 0 2 0.0 0.0 1.0 0.0 0.0 0.0\n"  // interior
                                     "TORSIONTABLE 1 2 0.0 0.0 2.0 0.0 0.0 0.0\n"  // N-terminal
                                     "TORSIONTABLE 2 2 0.0 0.0 3.0 0.0 0.0 0.0\n"; // C-terminal

const char * const INTERIOR_PSI = "TORSION\tALA\tPSI\tN\tCA\tC\t+N\t0\n";
const char * const NTER_ROTOR = "TORSION\tALA:NTER\tPHI\tH1\tN\tCA\tC\t1\n";
const char * const CTER_CARBOXYLATE = "TORSION\tALA:CTER\tPLANARITY\tCA\tO\tC\tOXT\t2\n";

size_t place(topology::Topology & top, const std::string & name, const std::string & chain, int resid)
{
    topology::ParticleProperties properties;
    properties.set_name(name);
    properties.set_residue_name("ALA");
    properties.set_chain_name(chain);
    properties.set_residue_id(resid);
    top.add_particle(topology::Particle(properties));
    return top.number_of_particles() - 1;
}

// `nresidues` alanines in one chain, numbered from 1, in particle order. The
// first carries the ammonium's H1 and the last the carboxylate's OXT, so the
// terminal records above have atoms to resolve; every residue carries
// N/CA/C/O.
void chain(topology::Topology & top, const std::string & name, int nresidues)
{
    for (int r = 1; r <= nresidues; ++r)
    {
        place(top, "N", name, r);
        if (r == 1)
            place(top, "H1", name, r);
        place(top, "CA", name, r);
        place(top, "C", name, r);
        place(top, "O", name, r);
        if (r == nresidues)
            place(top, "OXT", name, r);
    }
}

// The residue id of a torsion's first atom, so a test can say WHERE a record
// fired and not only how often.
int first_residue_of(const topology::Topology & top, const topology::Topology::Torsion & t)
{
    for (size_t i = 0; i < top.number_of_particles(); ++i)
        if (top.get_particle(i).unique_id() == t.atoms[0])
            return top.get_particle(i).properties().residue_id();
    return -1;
}

std::vector<unsigned> families(const topology::Topology & top)
{
    std::vector<unsigned> out;
    for (const auto & t : top.torsions())
        out.push_back(t.family);
    std::sort(out.begin(), out.end());
    return out;
}

} // namespace

TEST(BondedForceFieldReaderTerminal, AQualifierIsParsed)
{
    Probe reader;
    reader.setFileName(write(std::string(TERMINAL_TABLES) + INTERIOR_PSI + NTER_ROTOR + CTER_CARBOXYLATE));
    reader.read();
    EXPECT_EQ(reader.torsions(), 3u);
}

TEST(BondedForceFieldReaderTerminal, AnUnknownQualifierIsRefused)
{
    const std::string content = std::string(TERMINAL_TABLES) + "TORSION\tALA:MIDDLE\tPSI\tN\tCA\tC\t+N\t0\n";
    BondedForceFieldReader reader(write(content));
    EXPECT_DEATH(reader.read(), "unknown terminal qualifier 'MIDDLE'");
}

TEST(BondedForceFieldReaderTerminal, TwoQualifiersAreRefused)
{
    const std::string content = std::string(TERMINAL_TABLES) + "TORSION\tALA:NTER:CTER\tPSI\tN\tCA\tC\t+N\t0\n";
    BondedForceFieldReader reader(write(content));
    EXPECT_DEATH(reader.read(), "at most one terminal qualifier");
}

// The defect this guards, measured on 071 before the fix: the interior psi
// record fired on the chain's first residue, where amber99sb assigns no term at
// all -- a well 17.8 kJ.mol-1 deep and 20.3 kJ.mol-1.rad-1 of torque on a
// torsion the force field leaves free.
TEST(BondedForceFieldReaderTerminal, TheInteriorSetDoesNotApplyAtATerminus)
{
    BondedForceFieldReader reader(write(std::string(TERMINAL_TABLES) + INTERIOR_PSI + NTER_ROTOR));
    reader.read();

    topology::Topology top;
    chain(top, "A", 3);
    reader.buildSprings(top, nullptr, true, true, true);

    // Residue 1 is the N-terminus: it gets the rotor and NOT the psi, although
    // "N CA C +N" resolves there perfectly well. Residue 2 gets the psi.
    // Residue 3 has no following residue, so its psi does not exist.
    ASSERT_EQ(top.torsions().size(), 2u);
    for (const auto & t : top.torsions())
    {
        if (t.family == spn::SpringNetwork::DIHEDRAL_PHI)
        {
            EXPECT_EQ(first_residue_of(top, t), 1); // the rotor, on the terminus
            EXPECT_EQ(t.table, 1u);                 // and from the N-terminal table
        }
        else
        {
            EXPECT_EQ(t.family, static_cast<unsigned>(spn::SpringNetwork::DIHEDRAL_PSI));
            EXPECT_EQ(first_residue_of(top, t), 2); // the interior psi, interior only
            EXPECT_EQ(t.table, 0u);
        }
    }
}

// And the terms AMBER adds only at a terminus do get built there.
TEST(BondedForceFieldReaderTerminal, ACTerminalRecordAppliesAtTheLastResidueOnly)
{
    BondedForceFieldReader reader(write(std::string(TERMINAL_TABLES) + CTER_CARBOXYLATE));
    reader.read();

    topology::Topology top;
    chain(top, "A", 3);
    reader.buildSprings(top, nullptr, true, true, true);

    ASSERT_EQ(top.torsions().size(), 1u);
    EXPECT_EQ(first_residue_of(top, top.torsions()[0]), 3);
    EXPECT_EQ(top.torsions()[0].table, 2u);
}

// A terminus is read off the chain, so every chain has its own pair of ends.
TEST(BondedForceFieldReaderTerminal, EachChainHasItsOwnTermini)
{
    BondedForceFieldReader reader(write(std::string(TERMINAL_TABLES) + INTERIOR_PSI + NTER_ROTOR));
    reader.read();

    topology::Topology top;
    chain(top, "A", 3);
    chain(top, "B", 3);
    reader.buildSprings(top, nullptr, true, true, true);

    // Two rotors (one per chain's first residue) and two interior psis (one per
    // chain's second). Chain A's residue 3 is a C-terminus, so no psi joins it
    // to chain B's residue 1 -- that is _resolve_atom's own chain check, and
    // what makes the position reading consistent with it.
    const std::vector<unsigned> expected{spn::SpringNetwork::DIHEDRAL_PHI, spn::SpringNetwork::DIHEDRAL_PHI,
                                         spn::SpringNetwork::DIHEDRAL_PSI, spn::SpringNetwork::DIHEDRAL_PSI};
    EXPECT_EQ(families(top), expected);
}

// Backward compatibility, and not a nicety: the DNA and RNA files carry no
// qualifier at all, so without this fallback every nucleic chain end would
// silently lose all of its torsions instead of keeping today's parameters.
TEST(BondedForceFieldReaderTerminal, AFileWithNoQualifierKeepsTheInteriorSetAtATerminus)
{
    BondedForceFieldReader reader(write(std::string(TERMINAL_TABLES) + INTERIOR_PSI));
    reader.read();

    topology::Topology top;
    chain(top, "A", 3);
    reader.buildSprings(top, nullptr, true, true, true);

    // Residues 1 and 2 both get it now: nothing in the file says the terminus
    // should be treated differently.
    ASSERT_EQ(top.torsions().size(), 2u);
    EXPECT_EQ(first_residue_of(top, top.torsions()[0]), 1);
    EXPECT_EQ(first_residue_of(top, top.torsions()[1]), 2);
}

// The fallback is per residue name, not per file: a file that qualifies ALA
// says nothing about GLY, which keeps its interior records at a terminus.
TEST(BondedForceFieldReaderTerminal, TheFallbackIsPerResidueName)
{
    BondedForceFieldReader reader(write(std::string(TERMINAL_TABLES) + INTERIOR_PSI + NTER_ROTOR +
                                        "TORSION\tGLY\tPSI\tN\tCA\tC\t+N\t0\n"));
    reader.read();

    topology::Topology top;
    chain(top, "A", 2);
    top.get_particle(0).properties().set_residue_name("GLY");
    top.get_particle(1).properties().set_residue_name("GLY");
    top.get_particle(2).properties().set_residue_name("GLY");
    top.get_particle(3).properties().set_residue_name("GLY");
    top.get_particle(4).properties().set_residue_name("GLY");
    reader.buildSprings(top, nullptr, true, true, true);

    // The glycine N-terminus keeps its psi, where the alanine one would lose it.
    ASSERT_EQ(top.torsions().size(), 1u);
    EXPECT_EQ(first_residue_of(top, top.torsions()[0]), 1);
    EXPECT_EQ(top.torsions()[0].family, static_cast<unsigned>(spn::SpringNetwork::DIHEDRAL_PSI));
}

// A lone residue is both ends at once. AMBER has no template for one, so it
// takes the N-terminal set -- see _wanted_position.
TEST(BondedForceFieldReaderTerminal, ALoneResidueTakesTheNTerminalSet)
{
    BondedForceFieldReader reader(write(std::string(TERMINAL_TABLES) + NTER_ROTOR + CTER_CARBOXYLATE));
    reader.read();

    topology::Topology top;
    chain(top, "A", 1);
    reader.buildSprings(top, nullptr, true, true, true);

    ASSERT_EQ(top.torsions().size(), 1u);
    EXPECT_EQ(top.torsions()[0].family, static_cast<unsigned>(spn::SpringNetwork::DIHEDRAL_PHI));
    EXPECT_EQ(top.torsions()[0].table, 1u);
}
