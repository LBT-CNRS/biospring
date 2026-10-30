#include <gtest/gtest.h>

#include "IO/BondedForceFieldReader.h"

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
