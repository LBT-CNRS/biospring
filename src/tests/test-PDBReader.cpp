
#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <string>

#include "IO/PDBReader.h"
#include "topology/Topology.hpp"

using namespace biospring;

namespace
{

// Three atoms, and CONECT records that also name serials 1 and 2 (valid) plus
// 99 (absent from the file). Looking 99 up with std::map::operator[] would
// default-insert index 0 and bond the atom to the first particle instead.
const std::string PDB_WITH_DANGLING_CONECT = "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
                                             "ATOM      2  CA  ALA A   1       1.500   0.000   0.000  1.00  0.00           C\n"
                                             "ATOM      3  C   ALA A   1       3.000   0.000   0.000  1.00  0.00           C\n"
                                             "CONECT    1    2\n"
                                             "CONECT    3   99\n";

std::string write_temp_pdb(const std::string & name, const std::string & content)
{
    const std::string path = (std::filesystem::temp_directory_path() / name).string();
    std::ofstream out(path);
    out << content;
    out.close();
    return path;
}

} // namespace

// A CONECT record naming an atom that is not in the model (filtered out, or
// simply absent) must not silently create a bond to particle 0.
TEST(TestPDBReader, ConectToUnknownAtomDoesNotBondToParticleZero)
{
    const std::string path = write_temp_pdb("conect-dangling.pdb", PDB_WITH_DANGLING_CONECT);

    PDBReader reader(path);
    reader.read();

    EXPECT_EQ(reader.getTopology().number_of_particles(), 3u);
    // Only the 1-2 bond is resolvable; the 3-99 pair must be dropped, not
    // silently turned into a 3-0 bond.
    EXPECT_EQ(reader.getTopology().number_of_springs(), 1u);
}

// The same hazard through the atom filter: the CONECT targets exist in the file
// but were excluded from the model, so their ids are still absent from the map.
TEST(TestPDBReader, ConectToFilteredOutAtomDoesNotBondToParticleZero)
{
    const std::string path = write_temp_pdb("conect-filtered.pdb", PDB_WITH_DANGLING_CONECT);

    PDBReader reader(path);
    reader.addAtomfilter("CA"); // keeps only atom 2
    reader.read();

    EXPECT_EQ(reader.getTopology().number_of_particles(), 1u);
    EXPECT_EQ(reader.getTopology().number_of_springs(), 0u);
}

// A CONECT line whose own serial field is not a number used to throw
// std::invalid_argument out of parseConectLine and abort the read; the
// connected serials were already parsed defensively.
TEST(TestPDBReader, MalformedConectSerialIsSkippedNotThrown)
{
    const std::string path = write_temp_pdb("conect-malformed.pdb",
                                            "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
                                            "ATOM      2  CA  ALA A   1       1.500   0.000   0.000  1.00  0.00           C\n"
                                            "CONECT  ****    2\n");

    PDBReader reader(path);

    EXPECT_NO_THROW(reader.read());
    EXPECT_EQ(reader.getTopology().number_of_particles(), 2u);
    EXPECT_EQ(reader.getTopology().number_of_springs(), 0u);
}

// Well-formed records must still produce their bonds.
TEST(TestPDBReader, ValidConectRecordsStillCreateSprings)
{
    const std::string path = write_temp_pdb("conect-valid.pdb",
                                            "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
                                            "ATOM      2  CA  ALA A   1       1.500   0.000   0.000  1.00  0.00           C\n"
                                            "ATOM      3  C   ALA A   1       3.000   0.000   0.000  1.00  0.00           C\n"
                                            "CONECT    1    2    3\n");

    PDBReader reader(path);
    reader.read();

    EXPECT_EQ(reader.getTopology().number_of_particles(), 3u);
    EXPECT_EQ(reader.getTopology().number_of_springs(), 2u);
}

// SSBOND is the PDB's own record for a disulfide bridge, and it names its
// atoms by residue rather than by serial -- and sits BEFORE them in the file.
// GKinase.1S4Q.pdb declares its bridge this way and nothing else; before this
// was read, such a file came in with no bond at all.
TEST(TestPDBReader, SSBondRecordCreatesTheBridge)
{
    const std::string content =
        "SSBOND   1 CYS A   22    CYS A   45                          1555   1555  2.03\n"
        "ATOM      1  SG  CYS A  22       0.000   0.000   0.000  1.00  0.00           S\n"
        "ATOM      2  CB  CYS A  22       1.800   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  SG  CYS A  45       2.040   0.000   0.000  1.00  0.00           S\n";
    const std::string path = write_temp_pdb("ssbond.pdb", content);

    PDBReader reader(path);
    reader.read();

    EXPECT_EQ(reader.getTopology().number_of_particles(), 3u);
    ASSERT_EQ(reader.getTopology().number_of_springs(), 1u);
    const auto & spring = reader.getTopology().get_spring(0);
    EXPECT_EQ(spring.first().properties().residue_id(), 22);
    EXPECT_EQ(spring.second().properties().residue_id(), 45);
}

// The same bridge declared twice, the standard way and by hand, is one bond.
TEST(TestPDBReader, SSBondAlreadyDeclaredByConectIsNotDoubled)
{
    const std::string content =
        "SSBOND   1 CYS A   22    CYS A   45\n"
        "ATOM      1  SG  CYS A  22       0.000   0.000   0.000  1.00  0.00           S\n"
        "ATOM      2  SG  CYS A  45       2.040   0.000   0.000  1.00  0.00           S\n"
        "CONECT    1    2\n";
    const std::string path = write_temp_pdb("ssbond-and-conect.pdb", content);

    PDBReader reader(path);
    reader.read();

    EXPECT_EQ(reader.getTopology().number_of_springs(), 1u);
}

// LINK carries any other inter-residue bond, and names both atoms explicitly.
TEST(TestPDBReader, LinkRecordCreatesTheBond)
{
    const std::string content =
        "LINK         O   ALA A   1                 N   GLY A  10     1555   1555  2.90\n"
        "ATOM      1  O   ALA A   1       0.000   0.000   0.000  1.00  0.00           O\n"
        "ATOM      2  N   GLY A  10       2.900   0.000   0.000  1.00  0.00           N\n";
    const std::string path = write_temp_pdb("link.pdb", content);

    PDBReader reader(path);
    reader.read();

    ASSERT_EQ(reader.getTopology().number_of_springs(), 1u);
}

// A record naming an atom the model does not carry is routine on a real file
// (a HETATM the filter dropped): skipped, never bonded to particle 0.
TEST(TestPDBReader, SSBondToAbsentAtomCreatesNoBond)
{
    const std::string content =
        "SSBOND   1 CYS A   22    CYS A   99\n"
        "ATOM      1  SG  CYS A  22       0.000   0.000   0.000  1.00  0.00           S\n"
        "ATOM      2  SG  CYS A  45       2.040   0.000   0.000  1.00  0.00           S\n";
    const std::string path = write_temp_pdb("ssbond-dangling.pdb", content);

    PDBReader reader(path);
    reader.read();

    EXPECT_EQ(reader.getTopology().number_of_springs(), 0u);
}
