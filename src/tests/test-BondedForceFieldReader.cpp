#include <gtest/gtest.h>

#include "IO/BondedForceFieldReader.h"

#include <algorithm>
#include <fstream>

#include "TempPath.h"

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
