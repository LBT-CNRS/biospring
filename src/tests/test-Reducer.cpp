#include <gtest/gtest.h>

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "IO/io.h"
#include "IO/ReduceRuleReader.h"
#include "reduce/Reducer.h"
#include "TempPath.h"
#include "measure.hpp"
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

    // model.grp's SER rule used to require a nonexistent "CO" atom, so every SER
    // grain was always incomplete and dropped by default; now that it correctly
    // requires "OG", the 9 SER residues form complete grains and are kept too.
    EXPECT_EQ(reducer.target_topology().number_of_particles(), 183);
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

// A declared bond -- a PDB CONECT record -- must survive reduction AND the
// assignment pdb2spn makes immediately after it.
//
// This is the test that was missing, and its absence let a real bug through.
// A Spring holds a Particle&, not an index, so any spring that outlives a
// replacement of the particle list points at whatever landed at its old
// position. The failure is SILENT: the spring keeps its equilibrium length and
// simply joins two different atoms. Measured on 074's duplex, 48 base-pair
// springs with a correct r0 of 2.7-3.3 A ended up on atoms 5 to 12 A apart,
// which at mesh stiffness threw the structure to 1e8 A.
//
// The invariant that catches it: a carried spring's equilibrium is the length
// it was declared with, so it must still equal the distance between the two
// particles it joins. Every other check -- spring count, equilibrium value,
// even the endpoints' residue ids -- passes while the bond is wrong.
TEST_F(TestReducer, declared_bonds_survive_reduction_and_the_assignment)
{
    // Learn first which grains this model actually keeps: most residues here
    // lose one, and a bond whose end lands in a dropped grain is correctly
    // discarded -- which would make this test pass for the wrong reason.
    std::set<std::pair<int, std::string>> kept;
    {
        biospring::reduce::Reducer probe(topology);
        probe.initialize_forcefield(path_forcefield);
        probe.initialize_rules(path_reduce_rules);
        probe.reduce();
        const auto & reduced = probe.target_topology();
        for (size_t i = 0; i < reduced.number_of_particles(); ++i)
        {
            const auto & p = reduced.get_particle(i).properties();
            kept.insert({p.residue_id(), p.name()});
        }
    }
    ASSERT_FALSE(kept.empty());

    // Two atoms whose grains both survive, in residues far enough apart that
    // they cannot merge into one grain.
    biospring::reduce::ReduceRuleReader rules_reader(path_reduce_rules);
    rules_reader.read();
    const auto all_rules = rules_reader.rules();
    auto grain_of = [&](const biospring::topology::Particle & particle) -> std::string
    {
        const auto & p = particle.properties();
        const auto rules = all_rules.get_rules_for_residue(p.residue_name());
        for (size_t r = 0; r < rules.size(); ++r)
            if (rules[r].hasAtomNamed(p.name()))
                return rules[r].getName();
        return "";
    };

    size_t first = topology.number_of_particles(), second = topology.number_of_particles();
    int first_residue = -1;
    for (size_t i = 0; i < topology.number_of_particles(); ++i)
    {
        const auto & p = topology.get_particle(i).properties();
        const std::string grain = grain_of(topology.get_particle(i));
        if (grain.empty() || kept.find({p.residue_id(), grain}) == kept.end())
            continue;
        if (first == topology.number_of_particles())
        {
            first = i;
            first_residue = p.residue_id();
        }
        else if (p.residue_id() > first_residue + 4)
        {
            second = i;
            break;
        }
    }
    ASSERT_LT(first, topology.number_of_particles());
    ASSERT_LT(second, topology.number_of_particles());
    topology.add_spring(topology.get_particle(first), topology.get_particle(second));
    ASSERT_EQ(topology.number_of_springs(), 1u);

    biospring::reduce::Reducer reducer(topology);
    reducer.initialize_forcefield(path_forcefield);
    reducer.initialize_rules(path_reduce_rules);
    reducer.reduce();

    // Exactly what pdb2spn does, assignment included -- and the assignment is
    // where this used to break. Topology::operator= kept its own springs
    // whenever the two particle counts happened to match, while
    // _copy_particles replaced, in place, the very list those springs hold
    // references into. Equal counts say nothing about equal ordering, and the
    // reducer emits grains in rule order.
    topology = reducer.target_topology();

    {
        std::set<biospring::topology::pid_t> uids;
        size_t dup = 0;
        for (size_t i = 0; i < topology.number_of_particles(); ++i)
            if (!uids.insert(topology.get_particle(i).unique_id()).second)
                ++dup;
        EXPECT_EQ(dup, 0u) << "the reduced topology reuses particle unique ids, which is what the .nc "
                              "writer resolves spring endpoints through";
    }
    ASSERT_EQ(topology.number_of_springs(), 1u) << "the declared bond was lost by reduction";
    const auto & spring = topology.get_spring(0);
    EXPECT_NEAR(spring.equilibrium(), biospring::measure::distance(spring.first(), spring.second()), 1e-3)
        << "the spring kept its length but moved to different atoms";
}

// -- Main function  ----------------------------------------------------------
int main(int argc, char * argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

// ===========================================================================
// Rules that match nothing
// ===========================================================================
//
// amber.grp declares the N-terminal ammonium hydrogens and the C-terminal OXT
// for all 20 residues, so most of those rules necessarily match nothing in any
// given model. Reported per residue, as it was, that came to 42% of the whole
// example suite's build output and 10334 of 076's 10434 lines, which buried the
// four warnings that example actually had. These two tests pin down the shape of
// the replacement: a rule matching nowhere is COUNTED, costs no grain, and a
// rule matching somewhere is not counted at all.

static std::string grp_plus(const std::string & path, const std::string & extra)
{
    std::ifstream in(path);
    std::string out((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    out += extra;
    const std::string tmp = biospring::tests::temp_path(".grp");
    std::ofstream(tmp) << out;
    return tmp;
}

TEST_F(TestReducer, a_rule_matching_no_atom_anywhere_is_counted_and_costs_no_grain)
{
    // "ZZZ" is in no residue of model.pdb, so this rule can never build.
    const std::string grp = grp_plus(path_reduce_rules, "AOXT ALA ZZZ\n");

    biospring::reduce::Reducer reducer(topology);
    reducer.initialize_forcefield(path_forcefield);
    reducer.initialize_rules(grp);
    reducer.reduce();

    EXPECT_EQ(reducer.grains_never_built(), 1u);
    // The same model as without the rule: an unmatched rule removes nothing.
    EXPECT_EQ(reducer.target_topology().number_of_particles(), 183);
}

TEST_F(TestReducer, a_rule_matching_somewhere_is_not_counted_however_many_residues_skip_it)
{
    // CB exists in every residue of model.pdb but glycine, so this rule builds
    // for one residue name and is skipped by none of ITS OWN -- what must not
    // happen is for it to be reported as never built.
    const std::string grp = grp_plus(path_reduce_rules, "AXTRA ALA CB\n");

    biospring::reduce::Reducer reducer(topology);
    reducer.initialize_forcefield(path_forcefield);
    reducer.initialize_rules(grp);
    reducer.reduce();

    EXPECT_EQ(reducer.grains_never_built(), 0u);
}

// A grain name is keyed by residue: one residue type satisfying its own
// terminal rule must not silence another type's. ALA's extra rule matches, ARG's
// cannot, and the count has to be exactly one.
TEST_F(TestReducer, the_count_is_per_residue_not_per_grain_name)
{
    const std::string grp = grp_plus(path_reduce_rules, "XX ALA CB\nXX ARG ZZZ\n");

    biospring::reduce::Reducer reducer(topology);
    reducer.initialize_forcefield(path_forcefield);
    reducer.initialize_rules(grp);
    reducer.reduce();

    EXPECT_EQ(reducer.grains_never_built(), 1u);
}

// ===========================================================================
// Terminal rules
// ===========================================================================
//
// A free chain end is AMBER's own residue -- NALA, not ALA with three hydrogens
// bolted on -- and the retyping reaches every atom: 9 of the 10 atoms NALA
// shares with ALA carry a different charge. So the rules qualified for a
// position REPLACE the unqualified ones, and what the qualifier changes is the
// force-field type the grain is parameterised from, not the grain's name (see
// ReduceRule::forcefield_key).
//
// Which residue is at a free terminus is read off its own atoms, not just off
// its place in the chain: a carved-out fragment has neither the ammonium nor
// the carboxylate, and must keep the interior parameters (see
// topology::is_nterminal_only_atom). These tests are built on a three-residue
// alanine chain written here, with and without those atoms, rather than on
// model.pdb -- whose own chain ends are fragments, which is the point.

struct TestTerminalRules : public ::testing::Test
{
    // ALA 1-2-3, chain A. The ammonium hydrogens and the second carboxylate
    // oxygen are present only when asked for; without them the two ends are
    // carved-out fragments.
    static std::string write_chain(bool nterminal, bool cterminal)
    {
        const std::string path = biospring::tests::temp_path(".pdb");
        std::ofstream out(path);
        int serial = 0;
        for (int residue = 1; residue <= 3; residue++)
        {
            std::vector<std::string> atoms = {"N", "CA", "CB", "C", "O"};
            if (residue == 1 && nterminal)
                atoms.insert(atoms.begin() + 1, {"H1", "H2", "H3"});
            if (residue == 3 && cterminal)
                atoms.push_back("OXT");

            for (size_t i = 0; i < atoms.size(); i++)
            {
                char line[82];
                std::snprintf(line, sizeof(line), "ATOM  %5d %-4s ALA A%4d    %8.3f%8.3f%8.3f", ++serial,
                              atoms[i].c_str(), residue, 3.8 * residue, 0.5 * i, 0.0);
                out << line << "\n";
            }
        }
        return path;
    }

    // One rule per atom, interior and both ends, exactly as amber.grp carries
    // them: the terminal rules reuse the interior grain names.
    static std::string write_rules()
    {
        const std::string path = biospring::tests::temp_path(".grp");
        std::ofstream(path) << "AN ALA N\nAH ALA H\nACA ALA CA\nACB ALA CB\nAC ALA C\nAO ALA O\n"
                               "AN ALA:NTER N\nAH1 ALA:NTER H1\nAH2 ALA:NTER H2\nAH3 ALA:NTER H3\n"
                               "ACA ALA:NTER CA\nACB ALA:NTER CB\nAC ALA:NTER C\nAO ALA:NTER O\n"
                               "AN ALA:CTER N\nACA ALA:CTER CA\nACB ALA:CTER CB\nAC ALA:CTER C\n"
                               "AO ALA:CTER O\nAOXT ALA:CTER OXT\n";
        return path;
    }

    // Charges chosen so that a grain's own value says which rule built it: the
    // interior types are 0, the N-terminal ones +1 and the C-terminal ones -1.
    // `qualified` false leaves the qualified types out, which is the force field
    // an old .ff pairs with a new .grp.
    static std::string write_forcefield(bool qualified)
    {
        const std::string path = biospring::tests::temp_path(".ff");
        std::ofstream out(path);
        out << "#type\tcharge\tradius\tepsilon\tmass\ttransfer\thydrophobicity\n";
        for (const std::string & name : {"AN", "AH", "ACA", "ACB", "AC", "AO", "AH1", "AH2", "AH3", "AOXT"})
            out << name << "\t0.0\t1.908\t0.457730\t12.01\t0.0\t0.0\n";
        if (qualified)
        {
            for (const std::string & name : {"AN", "AH1", "AH2", "AH3", "ACA", "ACB", "AC", "AO"})
                out << name << ":NTER\t1.0\t1.908\t0.457730\t12.01\t0.0\t0.0\n";
            for (const std::string & name : {"AN", "ACA", "ACB", "AC", "AO", "AOXT"})
                out << name << ":CTER\t-1.0\t1.908\t0.457730\t12.01\t0.0\t0.0\n";
        }
        return path;
    }

    // Reduces `pdb` and returns the net charge of each residue, in order.
    static std::vector<float> net_charges(const std::string & pdb, const std::string & grp, const std::string & ff)
    {
        biospring::topology::Topology topology = biospring::io::readTopology(pdb);
        biospring::reduce::Reducer reducer(topology);
        reducer.initialize_forcefield(ff);
        reducer.initialize_rules(grp);
        reducer.reduce();

        const auto & reduced = reducer.target_topology();
        std::vector<float> net(3, 0.0f);
        for (size_t i = 0; i < reduced.number_of_particles(); i++)
        {
            const auto & properties = reduced.get_particle(i).properties();
            net.at(properties.residue_id() - 1) += properties.charge();
        }
        return net;
    }
};

TEST_F(TestTerminalRules, a_free_terminus_takes_the_qualified_rules_and_their_parameters)
{
    const std::vector<float> net =
        net_charges(write_chain(true, true), write_rules(), write_forcefield(true));

    // The N-terminal residue is built from 8 N-terminal types at +1 each, the
    // interior one from 5 interior types at 0, the C-terminal one from 6
    // C-terminal types at -1.
    EXPECT_NEAR(net[0], 8.0, 1e-5);
    EXPECT_NEAR(net[1], 0.0, 1e-5);
    EXPECT_NEAR(net[2], -6.0, 1e-5);
}

// What tells the position apart from the chemistry: the same chain, the same
// rules, only the terminal atoms removed.
TEST_F(TestTerminalRules, a_carved_out_fragment_keeps_the_interior_rules)
{
    const std::vector<float> net =
        net_charges(write_chain(false, false), write_rules(), write_forcefield(true));

    EXPECT_NEAR(net[0], 0.0, 1e-5);
    EXPECT_NEAR(net[1], 0.0, 1e-5);
    EXPECT_NEAR(net[2], 0.0, 1e-5);
}

// One end free and the other carved: each is decided on its own.
TEST_F(TestTerminalRules, the_two_ends_are_decided_separately)
{
    const std::vector<float> net =
        net_charges(write_chain(true, false), write_rules(), write_forcefield(true));

    EXPECT_NEAR(net[0], 8.0, 1e-5);
    EXPECT_NEAR(net[2], 0.0, 1e-5);
}

// The grain NAME must not change at a terminus: it is what the .hbond, .psite
// and .rbody lookups are keyed on, and renaming the terminal grains would drop
// every terminal residue out of all three, silently.
TEST_F(TestTerminalRules, the_grain_names_are_the_interior_ones)
{
    const std::string pdb = write_chain(true, true);
    biospring::topology::Topology topology = biospring::io::readTopology(pdb);
    biospring::reduce::Reducer reducer(topology);
    reducer.initialize_forcefield(write_forcefield(true));
    reducer.initialize_rules(write_rules());
    reducer.reduce();

    const auto & reduced = reducer.target_topology();
    std::set<std::string> nterminal, cterminal;
    for (size_t i = 0; i < reduced.number_of_particles(); i++)
    {
        const auto & properties = reduced.get_particle(i).properties();
        if (properties.residue_id() == 1)
            nterminal.insert(properties.name());
        if (properties.residue_id() == 3)
            cterminal.insert(properties.name());
    }

    EXPECT_EQ(nterminal, (std::set<std::string>{"AN", "AH1", "AH2", "AH3", "ACA", "ACB", "AC", "AO"}));
    EXPECT_EQ(cterminal, (std::set<std::string>{"AN", "ACA", "ACB", "AC", "AO", "AOXT"}));
}

// A .grp whose terminal rules are newer than the force field it is paired with
// must degrade to the interior parameters, not to particles with no mass --
// which is what a force-field type nobody declares gives.
TEST_F(TestTerminalRules, a_missing_qualified_force_field_type_keeps_the_interior_parameters)
{
    const std::string pdb = write_chain(true, true);
    biospring::topology::Topology topology = biospring::io::readTopology(pdb);
    biospring::reduce::Reducer reducer(topology);
    reducer.initialize_forcefield(write_forcefield(false));
    reducer.initialize_rules(write_rules());
    reducer.reduce();

    const auto & reduced = reducer.target_topology();
    ASSERT_EQ(reduced.number_of_particles(), 8u + 5u + 6u);
    for (size_t i = 0; i < reduced.number_of_particles(); i++)
    {
        EXPECT_NEAR(reduced.get_particle(i).properties().charge(), 0.0, 1e-5);
        EXPECT_NEAR(reduced.get_particle(i).properties().mass(), 12.01, 1e-5);
    }
}

// And a .grp that declares no qualifier at all -- every .grp here but amber.grp
// -- behaves exactly as before, free terminus or not.
TEST_F(TestReducer, a_terminal_rule_for_a_type_that_sits_nowhere_terminal_changes_nothing)
{
    // No ALA sits at either end of model.pdb's chain, so this rule can never
    // apply -- which is the normal case for a positional rule, and is counted
    // rather than warned about (see _report_grains_never_built).
    const std::string grp = grp_plus(path_reduce_rules, "ACA ALA:NTER N CA C O CB\n");

    biospring::reduce::Reducer reducer(topology);
    reducer.initialize_forcefield(path_forcefield);
    reducer.initialize_rules(grp);
    reducer.reduce();

    EXPECT_EQ(reducer.target_topology().number_of_particles(), 183);
}
