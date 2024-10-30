#include <gtest/gtest.h>

#include "IO/ForceFieldReader.h"
#include "ParticleProperty.h"

#include <fstream>

struct TestForceFieldReader : public ::testing::Test
{
    std::string content;
    std::string path;
    ForceFieldReader reader;

    void SetUp() override
    {
        ::testing::Test::SetUp();
        content = "#type	charge(e)	radius(A)	epsilon(kJ.mol-1) 	mass(Da)	transferIMP(kJ.mol-1.A-2)	"
                  "Hydrophobicity\n"
                  "ACA	12.00	2.88	123.00000	89.00   0.011407 	42.0\n"
                  "RCA	1.00	3.47	1.00000	174.00  -0.017504	17.0\n"
                  "NCA	0.00	3.13	1.00000	132.00  -0.017252	17.0\n"
                  "DCA	-1.00	3.00	1.00000	133.00  -0.022949	17.0\n"
                  "CCA	0.00	3.03	1.00000	121.00  0.047259	17.0\n"
                  "QCA	0.00	3.34	1.00000	146.00  -0.005096	17.0\n"
                  "ECA	-1.00	3.22	456.00000	147.00  -0.015356	17.0\n";
        path = "/tmp/dummy.ff";
        reader.setFileName(path);
    }

    void writeFF()
    {
        std::ofstream outstream;
        outstream.open(path);
        outstream << content;
        outstream.close();
    }
};

TEST_F(TestForceFieldReader, ReadProperties)
{
    writeFF();
    reader.read();

    const biospring::forcefield::ForceField & ff = reader.getForceField();
    EXPECT_EQ(7, ff.getNumberOfProperties());
    EXPECT_TRUE(ff.hasProperty("ACA"));
    EXPECT_TRUE(ff.hasProperty("RCA"));
    EXPECT_TRUE(ff.hasProperty("NCA"));
    EXPECT_TRUE(ff.hasProperty("DCA"));
    EXPECT_TRUE(ff.hasProperty("CCA"));
    EXPECT_TRUE(ff.hasProperty("QCA"));
    EXPECT_TRUE(ff.hasProperty("ECA"));

    biospring::spn::ParticleProperty pp = ff.getPropertiesFromName("ACA");
    EXPECT_FLOAT_EQ(12.00, pp.getCharge());
    EXPECT_FLOAT_EQ(2.88, pp.getRadius());
    EXPECT_FLOAT_EQ(123.00, pp.getEpsilon());
    EXPECT_FLOAT_EQ(89.00, pp.getMass());
    EXPECT_FLOAT_EQ(0.011407, pp.getTransferEnergyByAccessibleSurface());
    EXPECT_FLOAT_EQ(42.00, pp.getHydrophobicity());
}

TEST_F(TestForceFieldReader, LineCanHave6or7Tokens)
{
    content = "#type	charge(e)	radius(A)	epsilon(kJ.mol-1) 	mass(g.mol-1)	transferIMP(kJ.mol-1.A-2)\n"
              "1 2 3 4 5 6 7 8";
    writeFF();

    EXPECT_DEATH(reader.read(),
                 "!! ERROR: ForcefieldReader: line 2: invalid number of tokens .expected 6 or 7, found 8.");

    content = "#type	charge(e)	radius(A)	epsilon(kJ.mol-1) 	mass(g.mol-1)	transferIMP(kJ.mol-1.A-2)\n"
              "1 2 3 4 5";
    writeFF();
    EXPECT_DEATH(reader.read(),
                 "!! ERROR: ForcefieldReader: line 2: invalid number of tokens .expected 6 or 7, found 5.");
}

TEST_F(TestForceFieldReader, LineCanHave6Tokens)
{
    content = "#type	charge(e)	radius(A)	epsilon(kJ.mol-1) 	mass(g.mol-1)	transferIMP(kJ.mol-1.A-2)\n"
              "ACA	12.00	2.88	123.00000	89.00   0.011407";
    writeFF();

    try
    {
        reader.read();
    }
    catch (const std::exception & e)
    {
        FAIL() << "ForceFieldReader should be able to handle lines with 6 tokens: \"" << e.what() << "\"";
    }
}

TEST_F(TestForceFieldReader, LineCanHave7Tokens)
{
    content = "#type	charge(e)	radius(A)	epsilon(kJ.mol-1) 	mass(g.mol-1)	transferIMP(kJ.mol-1.A-2)	"
              "Hydrophobicity\n"
              "ACA	12.00	2.88	123.00000	89.00   0.011407 	42.0";
    writeFF();

    try
    {
        reader.read();
    }
    catch (const std::exception & e)
    {
        FAIL() << "ForceFieldReader should be able to handle lines with 7 tokens: \"" << e.what() << "\"";
    }
}

TEST_F(TestForceFieldReader, SupportsEmptyLines)
{
    content = "#\n\n\n";
    writeFF();

    try
    {
        reader.read();
    }
    catch (const std::exception & e)
    {
        FAIL() << "ForceFieldReader should be able to handle empty lines: \"" << e.what() << "\"";
    }
}

TEST_F(TestForceFieldReader, SupportsCommentLines)
{
    content = "# this is a comment line\n"
              "# this one as well";
    writeFF();

    try
    {
        reader.read();
    }
    catch (const std::exception & e)
    {
        FAIL() << "ForceFieldReader should be able to handle comment lines: \"" << e.what() << "\"";
    }
}

// == Tests file not found =============================================================

TEST_F(TestForceFieldReader, FileNotFoundFailure)
{
    path = "/not/a/valid/path.msp";
    reader.setFileName(path);
    EXPECT_DEATH(reader.read(), "!! ERROR: can't open file: '" + path + "'");
}

TEST_F(TestForceFieldReader, FileNameNotSetFailure)
{
    reader.setFileName("");
    EXPECT_DEATH(reader.read();, "!! ERROR: openread: empty file name");
}

// -- Main function  ----------------------------------------------------------
int main(int argc, char * argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
