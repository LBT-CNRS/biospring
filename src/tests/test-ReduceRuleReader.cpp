#include <gtest/gtest.h>

#include "IO/ReduceRuleReader.h"

#include <fstream>

struct TestReduceRuleReader : public ::testing::Test
{
    std::string content;
    std::string path;
    biospring::reduce::ReduceRuleReader reader;

    void SetUp() override
    {
        ::testing::Test::SetUp();
        content = "#Reduce rules for creation of the coarse-grain model by Zaccharias\n"
                  "ACA ALA CA\n"
                  "ACB ALA CA CB\n"
                  "RCA ARG CA\n"
                  "RCB ARG CG\n"
                  "RCG ARG NE CZ\n";
        path = "/tmp/dummy.grp";
        reader.setFileName(path);
    }

    void writeGRP()
    {
        std::ofstream outstream;
        outstream.open(path);
        outstream << content;
        outstream.close();
    }
};

TEST_F(TestReduceRuleReader, read_properties)
{
    writeGRP();
    reader.read();

    const auto & rules = reader.rules();
    ASSERT_EQ(5, rules.size());

    ASSERT_TRUE(rules.contains("ACA"));
    ASSERT_TRUE(rules.contains("ACB"));
    ASSERT_TRUE(rules.contains("RCA"));
    ASSERT_TRUE(rules.contains("RCB"));
    ASSERT_TRUE(rules.contains("RCG"));

    auto rule = rules["ACA"];
    EXPECT_EQ("ACA", rule.getName());
    EXPECT_EQ("ALA", rule.getResidueName());
    EXPECT_EQ(1, rule.getNumberOfAtoms());
    EXPECT_TRUE(rule.hasAtomNamed("CA"));

    rule = rules["ACB"];
    EXPECT_EQ("ACB", rule.getName());
    EXPECT_EQ("ALA", rule.getResidueName());
    EXPECT_EQ(2, rule.getNumberOfAtoms());
    EXPECT_TRUE(rule.hasAtomNamed("CA"));
    EXPECT_TRUE(rule.hasAtomNamed("CB"));

    rule = rules["RCA"];
    EXPECT_EQ("RCA", rule.getName());
    EXPECT_EQ("ARG", rule.getResidueName());
    EXPECT_EQ(1, rule.getNumberOfAtoms());
    EXPECT_TRUE(rule.hasAtomNamed("CA"));

    rule = rules["RCB"];
    EXPECT_EQ("RCB", rule.getName());
    EXPECT_EQ("ARG", rule.getResidueName());
    EXPECT_EQ(1, rule.getNumberOfAtoms());
    EXPECT_TRUE(rule.hasAtomNamed("CG"));

    rule = rules["RCG"];
    EXPECT_EQ("RCG", rule.getName());
    EXPECT_EQ("ARG", rule.getResidueName());
    EXPECT_EQ(2, rule.getNumberOfAtoms());
    EXPECT_TRUE(rule.hasAtomNamed("NE"));
    EXPECT_TRUE(rule.hasAtomNamed("CZ"));
}

TEST_F(TestReduceRuleReader, fails_if_not_at_least_3_tokens)
{
    content = "1 2";
    writeGRP();

    EXPECT_DEATH(reader.read(),
                 "!! ERROR: ReduceRuleReader: line 1: invalid number of tokens .expected at least 3, found 2.");
}

TEST_F(TestReduceRuleReader, supports_empty_lines)
{
    content = "\n\n\n";
    writeGRP();

    try
    {
        reader.read();
    }
    catch (const std::exception & e)
    {
        FAIL() << "ReduceRuleReader should be able to handle empty lines: \"" << e.what() << "\"";
    }
}

TEST_F(TestReduceRuleReader, supports_comment_lines)
{
    content = "# this is a comment line\n"
              "# this one as well";
    writeGRP();

    try
    {
        reader.read();
    }
    catch (const std::exception & e)
    {
        FAIL() << "ReduceRuleReader should be able to handle comment lines: \"" << e.what() << "\"";
    }
}

// == Tests file not found =============================================================

TEST_F(TestReduceRuleReader, file_not_found_failure)
{
    path = "/not/a/valid/path.msp";
    reader.setFileName(path);
    EXPECT_DEATH(reader.read(), "!! ERROR: can't open file: '" + path + "'");
}

TEST_F(TestReduceRuleReader, file_name_not_set_failure)
{
    reader.setFileName("");
    EXPECT_DEATH(reader.read();, "!! ERROR: openread: empty file name");
}

// ======================================================================================
//
// Tests for legacy ReduceRuleReader.
//
// ======================================================================================

struct TestReduceRuleReaderLegacy : public ::testing::Test
{
    std::string content;
    std::string path;
    biospring::reduce::legacy::ReduceRuleReader reader;

    void SetUp() override
    {
        ::testing::Test::SetUp();
        content = "#Reduce rules for creation of the coarse-grain model by Zaccharias\n"
                  "ACA ALA CA\n"
                  "ACB ALA CA CB\n"
                  "RCA ARG CA\n"
                  "RCB ARG CG\n"
                  "RCG ARG NE CZ\n";
        path = "/tmp/dummy.grp";
        reader.setFileName(path);
    }

    void writeGRP()
    {
        std::ofstream outstream;
        outstream.open(path);
        outstream << content;
        outstream.close();
    }
};

TEST_F(TestReduceRuleReaderLegacy, ReadProperties)
{
    writeGRP();
    reader.read();

    const biospring::reduce::legacy::Reduce * const reduce = reader.getReduce();
    EXPECT_EQ(5, reduce->getNumberOfRules());
    EXPECT_TRUE(reduce->hasRuleNamed("ACA"));
    EXPECT_TRUE(reduce->hasRuleNamed("ACB"));
    EXPECT_TRUE(reduce->hasRuleNamed("RCA"));
    EXPECT_TRUE(reduce->hasRuleNamed("RCB"));
    EXPECT_TRUE(reduce->hasRuleNamed("RCG"));

    biospring::reduce::ReduceRule rule = reduce->getReduceRule("ACA");
    EXPECT_EQ("ACA", rule.getName());
    EXPECT_EQ("ALA", rule.getResidueName());
    EXPECT_EQ(1, rule.getNumberOfAtoms());
    EXPECT_TRUE(rule.hasAtomNamed("CA"));

    rule = reduce->getReduceRule("ACB");
    EXPECT_EQ("ACB", rule.getName());
    EXPECT_EQ("ALA", rule.getResidueName());
    EXPECT_EQ(2, rule.getNumberOfAtoms());
    EXPECT_TRUE(rule.hasAtomNamed("CA"));
    EXPECT_TRUE(rule.hasAtomNamed("CB"));

    rule = reduce->getReduceRule("RCA");
    EXPECT_EQ("RCA", rule.getName());
    EXPECT_EQ("ARG", rule.getResidueName());
    EXPECT_EQ(1, rule.getNumberOfAtoms());
    EXPECT_TRUE(rule.hasAtomNamed("CA"));

    rule = reduce->getReduceRule("RCB");
    EXPECT_EQ("RCB", rule.getName());
    EXPECT_EQ("ARG", rule.getResidueName());
    EXPECT_EQ(1, rule.getNumberOfAtoms());
    EXPECT_TRUE(rule.hasAtomNamed("CG"));

    rule = reduce->getReduceRule("RCG");
    EXPECT_EQ("RCG", rule.getName());
    EXPECT_EQ("ARG", rule.getResidueName());
    EXPECT_EQ(2, rule.getNumberOfAtoms());
    EXPECT_TRUE(rule.hasAtomNamed("NE"));
    EXPECT_TRUE(rule.hasAtomNamed("CZ"));
}

TEST_F(TestReduceRuleReaderLegacy, FailsIfNotAtLeast3Tokens)
{
    content = "1 2";
    writeGRP();

    EXPECT_DEATH(reader.read(),
                 "!! ERROR: ReduceRuleReader: line 1: invalid number of tokens .expected at least 3, found 2.");
}

TEST_F(TestReduceRuleReaderLegacy, SupportsEmptyLines)
{
    content = "\n\n\n";
    writeGRP();

    try
    {
        reader.read();
    }
    catch (const std::exception & e)
    {
        FAIL() << "ReduceRuleReader should be able to handle empty lines: \"" << e.what() << "\"";
    }
}

TEST_F(TestReduceRuleReaderLegacy, SupportsCommentLines)
{
    content = "# this is a comment line\n"
              "# this one as well";
    writeGRP();

    try
    {
        reader.read();
    }
    catch (const std::exception & e)
    {
        FAIL() << "ReduceRuleReader should be able to handle comment lines: \"" << e.what() << "\"";
    }
}

// == Tests file not found =============================================================

TEST_F(TestReduceRuleReaderLegacy, FileNotFoundFailure)
{
    path = "/not/a/valid/path.msp";
    reader.setFileName(path);
    EXPECT_DEATH(reader.read(), "!! ERROR: can't open file: '" + path + "'");
}

TEST_F(TestReduceRuleReaderLegacy, FileNameNotSetFailure)
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
