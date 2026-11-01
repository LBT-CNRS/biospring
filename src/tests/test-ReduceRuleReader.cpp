#include <gtest/gtest.h>

#include "IO/ReduceRuleReader.h"

#include <filesystem>
#include <fstream>

#include "TempPath.h"

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
        path = biospring::tests::temp_path(".grp");
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

// ======================================================================================
// Terminal qualifier
// ======================================================================================
//
// The residue column may carry ":NTER" or ":CTER", the same qualifier the
// .bi.ff carries on the same column, read by the same helper
// (topology::parse_residue_field) so the two files cannot drift apart. Here:
// that the reader keeps AMBER's own residue name while recording the position,
// and refuses anything else.

TEST_F(TestReduceRuleReader, a_residue_name_without_a_qualifier_is_interior)
{
    writeGRP();
    reader.read();

    EXPECT_EQ(biospring::topology::POSITION_INTERIOR, reader.rules()["ACA"].position());
}

TEST_F(TestReduceRuleReader, a_qualifier_sets_the_position_and_leaves_the_residue_name_alone)
{
    content = "AN ALA:NTER N\n"
              "AOXT ALA:CTER OXT\n";
    writeGRP();
    reader.read();

    const auto & rules = reader.rules();
    ASSERT_EQ(2, rules.size());

    EXPECT_EQ("AN", rules[0].getName());
    EXPECT_EQ("ALA", rules[0].getResidueName());
    EXPECT_EQ(biospring::topology::POSITION_NTER, rules[0].position());
    EXPECT_TRUE(rules[0].hasAtomNamed("N"));

    EXPECT_EQ("AOXT", rules[1].getName());
    EXPECT_EQ("ALA", rules[1].getResidueName());
    EXPECT_EQ(biospring::topology::POSITION_CTER, rules[1].position());
}

// The grain name is deliberately the same as the interior rule's: it is the
// particle's identity for the .hbond, .psite and .rbody lookups, and only the
// PARAMETERS change at a terminus. So the two rules must coexist, and be told
// apart by position rather than by name.
TEST_F(TestReduceRuleReader, a_terminal_rule_may_reuse_the_interior_rule_s_name)
{
    content = "AN ALA N\n"
              "AN ALA:NTER N\n";
    writeGRP();
    reader.read();

    const auto & rules = reader.rules();
    ASSERT_EQ(2, rules.size());
    EXPECT_EQ(1, rules.get_rules_for_residue("ALA", biospring::topology::POSITION_INTERIOR).size());
    EXPECT_EQ(1, rules.get_rules_for_residue("ALA", biospring::topology::POSITION_NTER).size());
    EXPECT_EQ(0, rules.get_rules_for_residue("ALA", biospring::topology::POSITION_CTER).size());
    // And the position-blind overload, which is what resolves a .rbody or
    // .bi.ff atom name onto a renamed particle, still sees both.
    EXPECT_EQ(2, rules.get_rules_for_residue("ALA").size());
}

// The force-field key is where the qualifier shows up, because that is what the
// position actually changes.
TEST_F(TestReduceRuleReader, the_force_field_key_carries_the_qualifier)
{
    content = "AN ALA N\n"
              "AN ALA:NTER N\n"
              "AOXT ALA:CTER OXT\n";
    writeGRP();
    reader.read();

    const auto & rules = reader.rules();
    EXPECT_EQ("AN", rules[0].forcefield_key());
    EXPECT_EQ("AN:NTER", rules[1].forcefield_key());
    EXPECT_EQ("AOXT:CTER", rules[2].forcefield_key());
}

TEST_F(TestReduceRuleReader, an_unknown_qualifier_is_refused)
{
    content = "AN ALA:MIDDLE N\n";
    writeGRP();

    EXPECT_DEATH(reader.read(),
                 "ReduceRuleReader: line 1: 'ALA:MIDDLE' is not a residue name with at most one terminal qualifier");
}

TEST_F(TestReduceRuleReader, two_qualifiers_are_refused)
{
    content = "AN ALA:NTER:CTER N\n";
    writeGRP();

    EXPECT_DEATH(reader.read(), "at most one terminal qualifier");
}

TEST_F(TestReduceRuleReader, a_qualifier_with_no_residue_name_is_refused)
{
    content = "AN :NTER N\n";
    writeGRP();

    EXPECT_DEATH(reader.read(), "at most one terminal qualifier");
}
