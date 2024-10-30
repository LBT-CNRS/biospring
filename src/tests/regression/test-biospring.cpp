////////////////////////////////////////////////////////////////////////////////////////
//
// Regression tests for biospring.
//
// Makes sure biospring with various commands still give the same results as expected.
//
////////////////////////////////////////////////////////////////////////////////////////

#include <gtest/gtest.h>

#include "CommandLine.h"
#include "biospring-cli.h"
#include "version.h"


TEST(biospring, MissingConfiguration)
{
    CommandLine cmd({"biospring", "--nc", "foo.nc"});
    EXPECT_DEATH(biospring::biospringcli::main(cmd.size(), cmd.argv), biospring::logging::ERROR_PREFIX + "Required option '-c' has not been provided");
}

TEST(biospring, MissingTopology)
{
    CommandLine cmd({"biospring", "--msp", "foo.msp"});
    EXPECT_DEATH(biospring::biospringcli::main(cmd.size(), cmd.argv), biospring::logging::ERROR_PREFIX + "Required option '-s' has not been provided");
}

TEST(biospring, ConfigurationRequiresAValue)
{
    CommandLine cmd({"biospring", "-s", "foo.nc", "--msp"});
    EXPECT_DEATH(biospring::biospringcli::main(cmd.size(), cmd.argv), biospring::logging::ERROR_PREFIX + "Option '-c' requires an argument");
}

TEST(biospring, TopologyRequiresAValue)
{
    CommandLine cmd({"biospring", "--nc"});
    EXPECT_DEATH(biospring::biospringcli::main(cmd.size(), cmd.argv), biospring::logging::ERROR_PREFIX + "Option '-s' requires an argument");
}

TEST(biospring, HasVersionOption)
{
    // Test succeeds if program exists with status 0 and prints the expected description string.
    CommandLine cmd({"biospring", "--version"});
    ASSERT_EXIT(biospring::biospringcli::main(cmd.size(), cmd.argv), ::testing::ExitedWithCode(0), biospring::biospringcli::PROGRAM_VERSION);
}


TEST(biospring, HasHelpOption)
{
    // Test succeeds if program exists with status 0 and prints the expected description string.
    CommandLine cmd({"biospring", "--help"});
    ASSERT_EXIT(biospring::biospringcli::main(cmd.size(), cmd.argv), ::testing::ExitedWithCode(0), ".*"); // doesn't check prints description string
}

TEST(biospring, NoArgumentFailure)
{
    CommandLine cmd({"biospring"});
    EXPECT_DEATH(biospring::biospringcli::main(cmd.size(), cmd.argv),
                 biospring::logging::ERROR_PREFIX + "Required option '-s' has not been provided");
}

// -- Main function  ----------------------------------------------------------
int main(int argc, char * argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
