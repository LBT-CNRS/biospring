
////////////////////////////////////////////////////////////////////////////////////////
//
// Regression tests for pdb2spn.
//
// Makes sure pdb2spn with various commands still give the same results as expected.
//
// TODO:
//   - write tests for the --ignore-duplicate option
//
////////////////////////////////////////////////////////////////////////////////////////

#include <array>
#include <cstdio>
#include <fstream>

#include <gtest/gtest.h>

#include "CommandLine.h"
#include "IO/io.h"
#include "pdb2spn-cli.h"
#include "topology.hpp"
#include "utils.hpp"

const int DEFAULT_TOPOLOGY_PARTICLE_COUNT = 10;
const int DEFAULT_TOPOLOGY_SPRING_COUNT = 18;

// Dummy PDB string for a default topology with CONECT lines.
// Particle properties are easy to check and hard-coded.
const std::array<std::string, 20> DEFAULT_TOPOLOGY_PDB = {
    "ATOM      1  N   LYS A   1       1.000   2.000   3.000  2.00  3.0           N  ",
    "ATOM      2  N   LYS A   1       1.000   2.000   3.000  2.00  3.0           N  ",
    "ATOM      3  N   LYS A   1       1.000   2.000   3.000  2.00  3.0           N  ",
    "ATOM      4  N   LYS A   1       1.000   2.000   3.000  2.00  3.0           N  ",
    "ATOM      5  N   LYS A   1       1.000   2.000   3.000  2.00  3.0           N  ",
    "ATOM      6  N   LYS A   1       1.000   2.000   3.000  2.00  3.0           N  ",
    "ATOM      7  N   LYS A   1       1.000   2.000   3.000  2.00  3.0           N  ",
    "ATOM      8  N   LYS A   1       1.000   2.000   3.000  2.00  3.0           N  ",
    "ATOM      9  N   LYS A   1       1.000   2.000   3.000  2.00  3.0           N  ",
    "ATOM     10  N   LYS A   1       1.000   2.000   3.000  2.00  3.0           N  ",
    "CONECT    1    2",
    "CONECT    2    3    4",
    "CONECT    3    4    5",
    "CONECT    4    5    6",
    "CONECT    5    6    7",
    "CONECT    6    7    8",
    "CONECT    7    8    9",
    "CONECT    8    9   10",
    "CONECT    9    6   10",
    "CONECT   10    7",
};

// Run fixture.
//
// By default, the command-line argument is setup with the path to the topology file.
// When user calls run(), output files are automatically added to the command-line.
// Upon destruction, output files are deleted.
//
// Parameters:
//
//   - args, list of command-line arguments
//   - pathTopology, path to the input topology file
//   - pathOutput, list of the output files that will be created
//   - spn, SpringNetwork instance to check writing and reading back goes as expected
//
struct Run : public ::testing::Test
{
    std::vector<std::string> args;
    std::string pathTopology;
    std::vector<std::string> pathOutput;
    biospring::topology::Topology topology; // used to read output back

    void SetUp() override
    {
        ::testing::Test::SetUp();
        pathTopology = "model.pdb"; // default input file
    }

    void TearDown() override
    {
        // Delete output files.
        for (const auto & path : pathOutput)
            std::remove(path.c_str());
        ::testing::Test::TearDown();
    }

    void run()
    {
        std::vector<std::string> command = {"pdb2spn", "-s", pathTopology};
        if (pathOutput.size() > 0)
        {
            command.push_back("-o");
            for (const auto & path : pathOutput)
                command.push_back(path);
        }

        for (const auto & arg : args)
        {
            command.push_back(arg);
        }
        CommandLine cmd(command);
        std::cout << cmd << std::endl;
        biospring::pdb2spn::main(cmd.size(), cmd.argv);
    }

    void readOutput() { topology = biospring::io::readTopology(pathOutput[0]); }
};

// Run with default input and output.
//
// At setup-time, a dummy topology is created (see DEFAULT_TOPOLOGY_PDB).
struct RunDefaultTopology : virtual public Run
{
    void SetUp() override
    {
        Run::SetUp();
        writeDefaultTopology();
    }

    void TearDown() override
    {
        // Delete input file that has been created during setup.
        std::remove(pathTopology.c_str());
        Run::TearDown();
    }

    void writeDefaultTopology()
    {
        std::ofstream outfile(pathTopology);
        for (const std::string & line : DEFAULT_TOPOLOGY_PDB)
        {
            outfile << line << std::endl;
        }
        outfile.close();
    }
};

// Run where output is biospring's default binary NetCDF file (nc).
// Calling run(), runs pdb2spn and reads output file as well.
struct RunDefaultOutput : virtual public Run
{
    void SetUp() override
    {
        Run::SetUp();
        pathOutput.push_back("foo.nc"); // default output file
    }

    void run()
    {
        Run::run();
        readOutput();
    }
};

// Run with dummy topology and default output file.
// Calling run(), runs pdb2spn and reads output file as well.
struct RunDefaultTopologyDefaultOutput : virtual public RunDefaultOutput, RunDefaultTopology
{
    void SetUp() override
    {
        RunDefaultOutput::SetUp();
        RunDefaultTopology::SetUp();
    }

    void TearDown() override
    {
        RunDefaultOutput::TearDown();
        RunDefaultTopology::TearDown();
    }
};

// Run pdb2spn on the GK example.
struct RunGK : public RunDefaultOutput
{
    void SetUp() override
    {
        RunDefaultOutput::SetUp();
        pathTopology = "../data/model.pdb"; // default input file
    }
};

// Run pdb2spn on the GK example with reduction from all-atom to coarse grain.
struct RunGKReduce : public RunGK
{
    void SetUp() override
    {
        RunGK::SetUp();
        args.push_back("--grp");
        args.push_back("../data/model.grp");
        args.push_back("--ff");
        args.push_back("../data/model.ff");
    }
};

////////////////////////////////////////////////////////////////////////////////////////
//
// Tests on GK
//
////////////////////////////////////////////////////////////////////////////////////////

// TEST_F(RunGKReduce, Spring)
// {
//     args.push_back("--cutoff");
//     args.push_back("9.0");
//     args.push_back("--ignore-missing"); // gives the same number of spring as in GK example directory
//     run();
//     ASSERT_EQ(spn.getNumberOfSprings(), 1106);
// }

TEST_F(RunGKReduce, BasicReduce)
{
    run();
    ASSERT_EQ(topology.number_of_particles(), 174);
}

TEST_F(RunGKReduce, BasicReduceIgnoreMissing)
{
    args.push_back("--ignore-missing");
    run();
    ASSERT_EQ(topology.number_of_particles(), 183);
}

TEST_F(RunGK, NumberOfParticles)
{
    run();
    ASSERT_EQ(topology.number_of_particles(), 2814);
}

////////////////////////////////////////////////////////////////////////////////////////
//
// Tests on dummy topologies.
//
////////////////////////////////////////////////////////////////////////////////////////

TEST_F(RunDefaultTopology, createsPQRFileWithCharge)
{
    pathOutput.push_back("output.pqr");
    args.push_back("--charge");
    args.push_back("42");
    run();

    // Reads output.
    topology = biospring::io::readTopology(pathOutput[0]);

    // Checks results.
    ASSERT_EQ(topology.number_of_particles(), DEFAULT_TOPOLOGY_PARTICLE_COUNT);
    for (size_t i = 0; i < DEFAULT_TOPOLOGY_PARTICLE_COUNT; i++)
    {
        ASSERT_EQ(topology.get_particle(i).properties().charge(), 42.0);
    }
}

TEST_F(RunDefaultTopology, createsPQRFile)
{
    pathOutput.push_back("output.pqr");
    run();
    readOutput();

    // Checks results.
    ASSERT_EQ(topology.number_of_particles(), DEFAULT_TOPOLOGY_PARTICLE_COUNT);
    for (size_t i = 0; i < DEFAULT_TOPOLOGY_PARTICLE_COUNT; i++)
    {
        const biospring::topology::Particle & p = topology.get_particle(i);
        ASSERT_EQ(p.properties().charge(), 0.0);
        ASSERT_EQ(p.properties().radius(), 1.0);
    }
}

TEST_F(RunDefaultTopologyDefaultOutput, setCharges)
{
    args.push_back("--charge");
    args.push_back("42");
    run();

    ASSERT_EQ(topology.number_of_particles(), DEFAULT_TOPOLOGY_PARTICLE_COUNT);
    for (size_t i = 0; i < DEFAULT_TOPOLOGY_PARTICLE_COUNT; i++)
    {
        ASSERT_EQ(topology.get_particle(i).properties().charge(), 42);
    }
}

TEST_F(RunDefaultTopologyDefaultOutput, setStatic)
{
    args.push_back("--static");
    run();
    ASSERT_EQ(topology.number_of_particles(), DEFAULT_TOPOLOGY_PARTICLE_COUNT);
    for (size_t i = 0; i < DEFAULT_TOPOLOGY_PARTICLE_COUNT; i++)
    {
        ASSERT_FALSE(topology.get_particle(i).properties().is_dynamic());
    }
}

TEST_F(RunDefaultTopologyDefaultOutput, cutoffCreateSprings)
{
    args.push_back("--cutoff");
    args.push_back("100");               // creates a spring between all pairs of particle
    size_t expectedNumberOfSprings = 45; // combination(nparticles, 2)

    run();

    ASSERT_EQ(topology.number_of_particles(), DEFAULT_TOPOLOGY_PARTICLE_COUNT);
    ASSERT_EQ(topology.number_of_springs(), expectedNumberOfSprings);
    for (size_t i = 0; i < expectedNumberOfSprings; i++)
    {
        const biospring::topology::Spring & s = topology.get_spring(i);
        ASSERT_FLOAT_EQ(s.stiffness(), 1.0);
    }
}

TEST_F(RunDefaultTopologyDefaultOutput, setStiffness)
{
    // Creates springs with 42 as stiffness.
    args.push_back("--cutoff");
    args.push_back("100"); // creates a spring between all pairs of particle
    args.push_back("--stiffness");
    args.push_back("42");

    run();

    // Checks results
    size_t expectedNumberOfSprings = 45; // combination(nparticles, 2)

    ASSERT_EQ(topology.number_of_particles(), DEFAULT_TOPOLOGY_PARTICLE_COUNT);
    ASSERT_EQ(topology.number_of_springs(), expectedNumberOfSprings);

    for (size_t i = 0; i < expectedNumberOfSprings; i++)
    {
        const biospring::topology::Spring & s = topology.get_spring(i);
        ASSERT_FLOAT_EQ(s.stiffness(), 42);
    }
}

TEST_F(RunDefaultTopologyDefaultOutput, conectCreateSprings)
{
    run();

    ASSERT_EQ(topology.number_of_particles(), DEFAULT_TOPOLOGY_PARTICLE_COUNT);
    ASSERT_EQ(topology.number_of_springs(), DEFAULT_TOPOLOGY_SPRING_COUNT);
}

TEST_F(RunDefaultTopologyDefaultOutput, CreatesRightNumberOfParticles)
{
    run();

    ASSERT_EQ(topology.number_of_particles(), DEFAULT_TOPOLOGY_PARTICLE_COUNT);
    for (size_t i = 0; i < DEFAULT_TOPOLOGY_PARTICLE_COUNT; i++)
    {
        const biospring::topology::Particle & p = topology.get_particle(i);
        ASSERT_EQ(p.properties().atom_id(), i); // PDB id not written to output; internal id is.
        ASSERT_EQ(p.properties().name(), "N");
        ASSERT_EQ(p.properties().residue_name(), "LYS");
        ASSERT_EQ(p.properties().chain_name(), "A");
        ASSERT_EQ(p.properties().residue_id(), 1);
        ASSERT_FLOAT_EQ(p.getX(), 1.0);
        ASSERT_FLOAT_EQ(p.getY(), 2.0);
        ASSERT_FLOAT_EQ(p.getZ(), 3.0);
        // Occupancies and temperature factors are not written to nc file.
        ASSERT_TRUE(p.properties().is_dynamic());
    }
}

TEST_F(RunDefaultTopology, CanCreateMultipleOutputFiles)
{
    pathOutput.clear();
    pathOutput.push_back("foo.cdl");
    pathOutput.push_back("foo.nc");
    pathOutput.push_back("foo.pdb");

    run();

    ASSERT_TRUE(biospring::utils::file::exists("foo.cdl"));
    ASSERT_TRUE(biospring::utils::file::exists("foo.nc"));
    ASSERT_TRUE(biospring::utils::file::exists("foo.pdb"));
}

////////////////////////////////////////////////////////////////////////////////////////
//
// Asserts pdb2spn dies upon invalid command-line.
//
////////////////////////////////////////////////////////////////////////////////////////
TEST(pdb2spn, UnrecognizedArgument)
{
    CommandLine cmd({"pdb2spn", "-s", "topology.foo", "--static", "yes"});
    EXPECT_DEATH(biospring::pdb2spn::main(cmd.size(), cmd.argv),
                 biospring::logging::ERROR_PREFIX + "Unknown positional argument 'yes'");
}

TEST(pdb2spn, ChargeInvalidFloatValue)
{
    CommandLine cmd({"pdb2spn", "-s", "topology.foo", "--charge", "not-a-float"});
    EXPECT_DEATH(biospring::pdb2spn::main(cmd.size(), cmd.argv),
                 biospring::logging::ERROR_PREFIX +
                     "Invalid argument: --charge must be a float \\(got 'not-a-float'\\)");
}

TEST(pdb2spn, StiffnessInvalidFloatValue)
{
    CommandLine cmd({"pdb2spn", "-s", "topology.foo", "--stiffness", "not-a-float"});
    EXPECT_DEATH(biospring::pdb2spn::main(cmd.size(), cmd.argv),
                 biospring::logging::ERROR_PREFIX +
                     "Invalid argument: --stiffness must be a float \\(got 'not-a-float'\\)");
}

TEST(pdb2spn, CutoffInvalidFloatValue)
{
    CommandLine cmd({"pdb2spn", "-s", "topology.foo", "--cutoff", "not-a-float"});
    EXPECT_DEATH(biospring::pdb2spn::main(cmd.size(), cmd.argv),
                 biospring::logging::ERROR_PREFIX +
                     "Invalid argument: --cutoff must be a float \\(got 'not-a-float'\\)");
}

TEST(pdb2spn, MissingGRPWhenFFProvided)
{
    CommandLine cmd({"pdb2spn", "-s", "topology.foo", "--ff", "file.ff"});
    EXPECT_DEATH(biospring::pdb2spn::main(cmd.size(), cmd.argv),
                 biospring::logging::ERROR_PREFIX + "--grp <file> is mandatory when --ff is provided");
}

TEST(pdb2spn, MissingFFWhenGRPProvided)
{
    CommandLine cmd({"pdb2spn", "-s", "topology.foo", "--grp", "grpfile.grp"});
    EXPECT_DEATH(biospring::pdb2spn::main(cmd.size(), cmd.argv),
                 biospring::logging::ERROR_PREFIX + "--ff <file> is mandatory when --grp is provided");
}

TEST(pdb2spn, UnknownTopologyFormat)
{
    std::string path = "topology.foo";
    CommandLine cmd({"pdb2spn", "-s", path});
    EXPECT_DEATH(biospring::pdb2spn::main(cmd.size(), cmd.argv), "!! ERROR: Unrecognized topology format 'foo'");
}

TEST(pdb2spn, TopologyNotFound)
{
    std::string path = "/not/a/topology.pdb";
    CommandLine cmd({"pdb2spn", "-s", path});
    EXPECT_DEATH(biospring::pdb2spn::main(cmd.size(), cmd.argv), "!! ERROR: can't open file: '" + path + "'");
}

TEST(pdb2spn, OutputRequiresAValue)
{
    CommandLine cmd({"pdb2spn", "-s", "foo.pdb", "-o"});
    EXPECT_DEATH(biospring::pdb2spn::main(cmd.size(), cmd.argv),
                 biospring::logging::ERROR_PREFIX + "Option '-o' requires at least 1 argument.");
}

TEST(pdb2spn, ForceFieldRequiresAValue)
{
    CommandLine cmd({"pdb2spn", "-s", "foo.pdb", "--grp", "grp.grp", "--ff"});
    EXPECT_DEATH(biospring::pdb2spn::main(cmd.size(), cmd.argv),
                 biospring::logging::ERROR_PREFIX + "Option '--ff' requires an argument.");
}

TEST(pdb2spn, GroupRequiresAValue)
{
    CommandLine cmd({"pdb2spn", "-s", "foo.pdb", "--ff", "file.ff", "--grp"});
    EXPECT_DEATH(biospring::pdb2spn::main(cmd.size(), cmd.argv),
                 biospring::logging::ERROR_PREFIX + "Option '--grp' requires an argument.");
}

TEST(pdb2spn, TopologyRequiresAValue)
{
    CommandLine cmd({"pdb2spn", "-s"});
    EXPECT_DEATH(biospring::pdb2spn::main(cmd.size(), cmd.argv),
                 biospring::logging::ERROR_PREFIX + "Option '-s' requires an argument");
}

TEST(pdb2spn, HasHelpOption)
{
    // Test succeeds if program exists with status 0 and prints the expected description string.
    CommandLine cmd({"pdb2spn", "--help"});
    ASSERT_EXIT(biospring::pdb2spn::main(cmd.size(), cmd.argv), ::testing::ExitedWithCode(0),
                ".*"); // doesn't check prints description string
}

TEST(pdb2spn, NoArgumentFailure)
{
    CommandLine cmd({"pdb2spn"});
    EXPECT_DEATH(biospring::pdb2spn::main(cmd.size(), cmd.argv),
                 biospring::logging::ERROR_PREFIX + "Required option '-s' has not been provided");
}

// -- Main function  ----------------------------------------------------------
int main(int argc, char * argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
