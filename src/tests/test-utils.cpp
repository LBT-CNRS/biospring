
#include "../utils.hpp"

#include <fstream>
#include <gtest/gtest.h>
#include <string>


class StringCompareTestFixture : public ::testing::TestWithParam<std::tuple<std::string, std::string>> {};

using namespace biospring::utils;


////////////////////////////////////////////////////////////////////////////////////////
//
// Path functions
//
////////////////////////////////////////////////////////////////////////////////////////


// -- has/hasNoExtension functions ---------------------------------------------

TEST(HasNoExtension, NoExtension)
{
    EXPECT_TRUE(path::hasNoExtension("foo"));
}

TEST(HasNoExtension, Extension)
{
    EXPECT_FALSE(path::hasNoExtension("foo.pdb"));
    EXPECT_FALSE(path::hasNoExtension("foo.pdb.io"));
}

TEST(HasExtension, NoExtension)
{
    EXPECT_FALSE(path::hasExtension("foo"));
}

TEST(HasExtension, Extension)
{
    EXPECT_TRUE(path::hasExtension("foo.pdb"));
    EXPECT_TRUE(path::hasExtension("foo.pdb.io"));
}


// -- splitExtension function --------------------------------------------------
TEST(splitExtension, SimpleCase)
{
    auto result = path::splitExtension("foo.bar");
    EXPECT_EQ("foo", result.first);
    EXPECT_EQ("bar", result.second);
}

TEST(splitExtension, NoExtension)
{
    auto result = path::splitExtension("foo");
    EXPECT_EQ("foo", result.first);
    EXPECT_EQ("", result.second);
}

TEST(splitExtension, MultipleExtensions)
{
    auto result = path::splitExtension("foo.bar.baz");
    EXPECT_EQ("foo.bar", result.first);
    EXPECT_EQ("baz", result.second);
}

// -- getExtension function ----------------------------------------------------
class getExtensionTest : public StringCompareTestFixture {};

TEST_P(getExtensionTest, ReturnsActualExtension)
{
    const std::string expected = std::get<0>(GetParam());
    const std::string str = std::get<1>(GetParam());
    EXPECT_EQ(expected, path::getExtension(str));
}

INSTANTIATE_TEST_SUITE_P(
    getExtensionTestSuite, getExtensionTest,
    ::testing::Values(
        std::make_tuple("", ""),
        std::make_tuple("", "foo"),
        std::make_tuple("bar", "foo.bar"),
        std::make_tuple("baz", "foo.bar.baz")));


////////////////////////////////////////////////////////////////////////////////////////
//
// String functions
//
////////////////////////////////////////////////////////////////////////////////////////

// -- String join method --------------------------------------------------------------
TEST(join, JoinWithSeparator)
{
    std::vector<std::string> tokens = { "Hello", "World!", "It", "is", "sunny", "today!" };
    const std::string expected = "Hello-World!-It-is-sunny-today!";
    EXPECT_EQ(expected, string::join(tokens, "-"));
}

TEST(join, JoinWithSeparatorCarriageReturn)
{
    std::vector<std::string> tokens = { "Hello", "World!", "It", "is", "sunny", "today!" };
    const std::string expected = "Hello\nWorld!\nIt\nis\nsunny\ntoday!";
    EXPECT_EQ(expected, string::join(tokens, "\n"));
}


TEST(join, JoinWithNoSeparator)
{
    std::vector<std::string> tokens = { "Hello", "World!", "It", "is", "sunny", "today!" };
    const std::string expected = "HelloWorld!Itissunnytoday!";
    EXPECT_EQ(expected, string::join(tokens));
}


// -- String split method -------------------------------------------------------------

TEST(split, SplitUsingSeparator)
{
    std::vector<std::string> expected = {"Hello", "World!"};
    EXPECT_EQ(expected, string::split("Hello,World!", ","));

    expected = {"Hello", "World!", "Foo", "Bar"};
    EXPECT_EQ(expected, string::split("Hello,World!,Foo,Bar", ","));
}

TEST(split, SplitUsingSpace)
{
    std::vector<std::string> expected = {"Hello", "World!"};
    EXPECT_EQ(expected, string::split("Hello World!"));

    expected = {"Hello", "World!", "Foo", "Bar"};
    EXPECT_EQ(expected, string::split("Hello World! Foo Bar"));
}

// -- toupper/lower functions ---------------------------------------------------------
TEST(toupper, EmptyString) { EXPECT_EQ("", string::toupper("")); }

TEST(toupper, ActualString) { EXPECT_EQ("FOO BAR BAZ 1 2 3", string::toupper("foo bar baz 1 2 3")); }

TEST(tolower, EmptyString) { EXPECT_EQ("", string::tolower("")); }

TEST(tolower, ActualString) { EXPECT_EQ("foo bar baz 1 2 3", string::tolower("FOO BAR BAZ 1 2 3")); }



// ----------------------------------------------------------------------------------
//
// from_string function
//
// ----------------------------------------------------------------------------------


TEST(from_string, Fails)
{
    // When failing, expects from_string to return false and to stay
    // input value unchanged.

    int n = 0;
    EXPECT_FALSE(string::from_string(n, ""));
    EXPECT_EQ(0, n);

    n = 1;
    EXPECT_FALSE(string::from_string(n, ""));
    EXPECT_EQ(1, n);
}


TEST(from_string, StringToInteger)
{
    int n;
    ASSERT_TRUE(string::from_string(n, "2"));
    EXPECT_EQ(2, n);

    ASSERT_TRUE(string::from_string(n, "+2"));
    EXPECT_EQ(2, n);

    ASSERT_TRUE(string::from_string(n, "-2"));
    EXPECT_EQ(-2, n);
}

TEST(from_string, StringToFloat)
{
    float n;
    ASSERT_TRUE(string::from_string(n, "2.6"));
    EXPECT_FLOAT_EQ(2.6, n);

    ASSERT_TRUE(string::from_string(n, "+2.6"));
    EXPECT_FLOAT_EQ(2.6, n);

    ASSERT_TRUE(string::from_string(n, "-2.6"));
    EXPECT_FLOAT_EQ(-2.6, n);
}


// -- from_string<bool> specialization ------------------------------------------------------
class FromStringToBoolParameterizedTestFixture : public ::testing::TestWithParam<std::tuple<std::string, bool>> {};

TEST_P(FromStringToBoolParameterizedTestFixture, FromStringToBool)
{
    const std::string str = std::get<0>(GetParam());
    const bool expected = std::get<1>(GetParam());
    bool b = !expected;
    ASSERT_TRUE(string::from_string(b, str));
    EXPECT_EQ(expected, b);
}

INSTANTIATE_TEST_SUITE_P(
    FromStringTrueTests,
    FromStringToBoolParameterizedTestFixture,
    ::testing::Values(
        std::make_tuple("true", true),
        std::make_tuple("TRUE", true),
        std::make_tuple("1", true),
        std::make_tuple("on", true),
        std::make_tuple("false", false),
        std::make_tuple("FALSE", false),
        std::make_tuple("0", false),
        std::make_tuple("off", false)));


TEST(from_string, BoolFails)
{
    bool b;
    EXPECT_FALSE(string::from_string(b, "not a bool"));
    EXPECT_FALSE(string::from_string(b, "42"));
}


// ----------------------------------------------------------------------------------
//
// Trim functions
//
// ----------------------------------------------------------------------------------

class LTrimTest : public StringCompareTestFixture {};
class RTrimTest : public StringCompareTestFixture {};
class TrimTest : public StringCompareTestFixture {};


// -- Left Trim function ------------------------------------------------------
TEST_P(LTrimTest, ActuallyTrimsLeft)
{
    const std::string str = std::get<1>(GetParam());
    const std::string expected = std::get<0>(GetParam());
    EXPECT_EQ(expected, string::ltrim(str));
}

INSTANTIATE_TEST_SUITE_P(
    LTrimTestSuite, LTrimTest,
    ::testing::Values(
        std::make_tuple("Hello World!", "Hello World!"),
        std::make_tuple("Hello World!", "    Hello World!"),
        std::make_tuple("Hello World!    ", "Hello World!    "),
        std::make_tuple("Hello World!    ", "    Hello World!    ")));


// -- Right Trim function  ----------------------------------------------------
TEST_P(RTrimTest, ActuallyTrimsRight)
{
    const std::string str = std::get<1>(GetParam());
    const std::string expected = std::get<0>(GetParam());
    EXPECT_EQ(expected, string::rtrim(str));
}

INSTANTIATE_TEST_SUITE_P(
    RTrimTestSuite, RTrimTest,
    ::testing::Values(
        std::make_tuple("Hello World!", "Hello World!"),
        std::make_tuple("    Hello World!", "    Hello World!"),
        std::make_tuple("Hello World!", "Hello World!    "),
        std::make_tuple("    Hello World!", "    Hello World!    ")));


// -- Left & Right Trim functions  ---------------------------------------------
TEST_P(TrimTest, ActuallyTrimsBothSides)
{
    const std::string str = std::get<1>(GetParam());
    const std::string expected = std::get<0>(GetParam());
    EXPECT_EQ(expected, string::trim(str));
}

INSTANTIATE_TEST_SUITE_P(
    TrimTestSuite, TrimTest,
    ::testing::Values(
        std::make_tuple("Hello World!", "Hello World!"),
        std::make_tuple("Hello World!", "    Hello World!"),
        std::make_tuple("Hello World!", "Hello World!    "),
        std::make_tuple("Hello World!", "    Hello World!    ")));


////////////////////////////////////////////////////////////////////////////////////////
//
// File functions
//
////////////////////////////////////////////////////////////////////////////////////////

TEST(File, FileDoesNotExists)
{
    const std::string fname = "/tmp/foo";
    EXPECT_FALSE(file::exists(fname));
}

TEST(File, FileExists)
{
    const std::string fname = "/tmp/foo";
    std::ofstream outfile(fname, std::ofstream::out);
    EXPECT_TRUE(outfile);
    EXPECT_TRUE(file::exists(fname));
    outfile.close();
    std::remove(fname.c_str());
}

// -- Open file for writing ---------------------------------------------------
TEST(File, OpenwriteSuccess)
{
    const std::string fname = "/tmp/foo";
    std::ofstream outfile;
    file::openwrite(fname, outfile);
    EXPECT_TRUE(outfile);
    outfile.close();
    std::remove(fname.c_str());
}

TEST(File, OpenwriteFailure)
{
    const std::string fname = "/not-a-directory/foo"; // cannot write to a non-existing directory
    std::ofstream outfile;
    EXPECT_DEATH(file::openwrite(fname, outfile), "!! ERROR: can't open file: '" + fname + "'");
}

// -- Open file for reading ---------------------------------------------------
TEST(File, OpenreadSuccess)
{
    const std::string fname = "dummy";

    // Creates a temporary file.
    std::ofstream cfgout(fname);
    EXPECT_TRUE(cfgout.is_open());
    cfgout.close();

    // Tries opening back this file (should not fail).
    std::ifstream infile;
    file::openread(fname, infile);
    infile.close();
    std::remove(fname.c_str());
}

TEST(File, OpenreadFailure)
{
    // Tries opening a file that does not exist (should fail).
    const std::string fname = "foo";
    std::ifstream infile;
    EXPECT_DEATH(file::openread(fname, infile), "!! ERROR: can't open file: '" + fname + "'");
}

TEST(File, OpenreadEmptyPath)
{
    // Tries opening a file that does not exist (should fail).
    const std::string fname = "";
    std::ifstream infile;
    EXPECT_DEATH(file::openread(fname, infile), "!! ERROR: openread: empty file name");
}

// -- Main function  ----------------------------------------------------------
int main(int argc, char * argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
