
#include <gtest/gtest.h>

#include "IO/OpenDXReader.h"
#include "SpringNetwork.h"

struct TestOpenDXReader : public ::testing::Test
{
    std::string path;
    OpenDXReader reader;

    void SetUp() override
    {
        ::testing::Test::SetUp();
        path = "data/sample.dx";
        reader.setFileName(path);
        reader.read();
    }
};

TEST_F(TestOpenDXReader, GetScalingFactors)
{
    EXPECT_FLOAT_EQ(reader.getGrid().cell_size()[0], 1);
    EXPECT_FLOAT_EQ(reader.getGrid().cell_size()[1], 1);
    EXPECT_FLOAT_EQ(reader.getGrid().cell_size()[2], 1);
}

TEST_F(TestOpenDXReader, GetOrigin)
{
    EXPECT_FLOAT_EQ(reader.getGrid().origin()[0], 1.356501);
    EXPECT_FLOAT_EQ(reader.getGrid().origin()[1], -1.401500);
    EXPECT_FLOAT_EQ(reader.getGrid().origin()[2], 0.500499);
}

TEST_F(TestOpenDXReader, GetDimensions)
{
    EXPECT_EQ(reader.getGrid().shape()[0], 58);
    EXPECT_EQ(reader.getGrid().shape()[1], 54);
    EXPECT_EQ(reader.getGrid().shape()[2], 47);
}

TEST(OpenDXReader, FileNotFoundFailure)
{
    std::string path = "/not/a/valid/path.dx";
    OpenDXReader reader;
    reader.setFileName(path);
    EXPECT_DEATH(reader.read(), "!! ERROR: can't open file: '" + path + "'");
}

TEST(OpenDXReader, FileNameNotSetFailure)
{
    OpenDXReader reader;
    EXPECT_DEATH(reader.read();, "!! ERROR: openread: empty file name");
}

// -- Main function  ----------------------------------------------------------
int main(int argc, char * argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
