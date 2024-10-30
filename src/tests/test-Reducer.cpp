#include <gtest/gtest.h>

#include <filesystem>
#include <string>

#include "IO/io.h"
#include "reduce/Reducer.h"
#include "topology.hpp"

static const std::string PATH_FORCEFIELD = "data/model.ff";
static const std::string PATH_REDUCE_RULES = "data/model.grp";
static const std::string PATH_TOPOLOGY = "data/model.pdb";

struct TestReducer : public ::testing::Test
{
    std::string path_topology;
    std::string path_forcefield;
    std::string path_reduce_rules;
    biospring::topology::Topology topology;

    void SetUp() override
    {
        ::testing::Test::SetUp();
        path_topology = PATH_TOPOLOGY;
        path_forcefield = PATH_FORCEFIELD;
        path_reduce_rules = PATH_REDUCE_RULES;

        ASSERT_TRUE(std::filesystem::exists(path_topology));
        ASSERT_TRUE(std::filesystem::exists(path_forcefield));
        ASSERT_TRUE(std::filesystem::exists(path_reduce_rules));

        topology = biospring::io::readTopology(path_topology);
        ASSERT_EQ(topology.number_of_particles(), 2814);
    }
};

// ===========================================================================
// Reduce tests
// ===========================================================================

TEST_F(TestReducer, reduce)
{
    biospring::reduce::Reducer reducer(topology);
    reducer.initialize_forcefield(path_forcefield);
    reducer.initialize_rules(path_reduce_rules);
    reducer.reduce();

    EXPECT_EQ(reducer.target_topology().number_of_particles(), 174);
}

TEST_F(TestReducer, reduce_ignore_missing)
{
    biospring::reduce::Reducer reducer(topology);
    reducer.initialize_forcefield(path_forcefield);
    reducer.initialize_rules(path_reduce_rules);
    reducer.set_ignore_missing_particles(true);
    reducer.reduce();

    EXPECT_EQ(reducer.target_topology().number_of_particles(), 183);
}

// ===========================================================================
// Basic tests
// ===========================================================================

TEST_F(TestReducer, initialize_forcefield)
{
    biospring::reduce::Reducer reducer(topology);
    reducer.initialize_forcefield(path_forcefield);
    EXPECT_EQ(reducer.forcefield().getNumberOfProperties(), 20);
}

TEST_F(TestReducer, initialize_rules)
{
    biospring::reduce::Reducer reducer(topology);
    reducer.initialize_rules(path_reduce_rules);
    EXPECT_EQ(reducer.rules().size(), 20);
}

// -- Main function  ----------------------------------------------------------
int main(int argc, char * argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
