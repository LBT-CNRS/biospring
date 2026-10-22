
#include "mergespn-cli.h"
#include "IO/DonorAcceptorRuleReader.h"
#include "IO/io.h"
#include "topology.hpp"

namespace biospring
{
namespace mergespn
{

const std::string PROGRAM_VERSION = "0.1.0";
const biospring::argparse::description_t PROGRAM_DESCRIPTION = {
    "mergespn merges two or more spring-network topologies.",
    "",
    "Input topology formats are selected from the -s/--topology filename extension:",
    "  .nc  : binary NetCDF spring-network file",
    "  .pdb : Protein Data Bank file",
    "  .pqr : PQR file",
    "",
    "Output formats are selected from the -o/--output filename extension:",
    "  .nc  : binary NetCDF spring-network file",
    "  .cdl : text NetCDF/CDL spring-network file",
    "  .pdb : Protein Data Bank file",
    "  .pqr : PQR file",
    "",
    "Use -cutoff/--cutoff to create springs between merged structures.",
    "",
    "A cutoff wide enough to hold two fragments together is also wide enough to",
    "weld the base pairs the dynamic hydrogen bond term is meant to form, and such",
    "a spring is mechanically negligible while silently pre-empting the bond. Pass",
    "-hbond/--hbond with the same donor/acceptor table the simulation will use and",
    "those pairs are left alone.",
};

int main(int argc, char ** argv)
{
    // Command-line parsing.
    CommandLineArguments args(std::string(argv[0]), PROGRAM_DESCRIPTION, PROGRAM_VERSION);
    args.parseCommandLine(argc, argv);
    args.printArgumentValues();

    topology::Topology output;

    // Iterates over all input topologies.
    for (size_t topology_id = 0; topology_id < args.pathTopologyList.size(); ++topology_id)
    {
        const std::string & path = args.pathTopologyList[topology_id];
        logging::status("Reading %s", path.c_str());
        auto top = io::readTopology(path);
        top.particles().set_topology_ids(topology_id);
        output = output.merge(top);
    }

    if (args.useUserCutoff())
    {
        // When a donor/acceptor table is given, a pair that could hydrogen bond
        // is left unwelded.
        //
        // WHY THIS MATTERS, measured on a 70S assembled by this very tool: a
        // 7 A join put 1780 springs across the ribosome's own base-edge N/O
        // atoms at the default stiffness of 1.0, against 650 for the all-atom
        // mesh. They held nothing -- removing all of them moved the spring
        // energy from 40.03 to 40.07 kJ/mol -- and yet they took the real base
        // pairs the hydrogen bond term holds from 48 of 49 down to 16 of 49,
        // because the assignment refuses any pair whose atoms are already
        // sprung. Mechanically invisible, functionally decisive.
        //
        // The cutoff is not reconnecting cut bonds here: across every
        // multi-fragment component of that 70S the closest inter-fragment
        // contact is 2.18 A, never a covalent 1.6 A, so no narrower cutoff
        // would do. What it provides is cohesion, and cohesion does not need
        // the base edges.
        topology::Topology::SpringVeto veto = nullptr;
        io::DonorAcceptorRuleReader hbond;
        if (args.useHydrogenBondTable())
        {
            logging::status("Reading donor/acceptor table %s", args.pathHydrogenBond.c_str());
            hbond.setFileName(args.pathHydrogenBond);
            hbond.read();
            // BOTH atoms being hydrogen-bond sites is enough; it does not have to
            // be a donor facing an acceptor.
            //
            // The narrower rule was tried first and it leaves a real residue. On
            // the 70S it removed the 212 donor-to-acceptor springs welding native
            // base pairs, which is what pre-empts a bond, but left 238 that join
            // two donors or two acceptors of the SAME paired bases. Those block
            // nothing -- pre-emption is tested on the donor/acceptor pair at
            // hand -- yet they still tie two bases together mechanically, and two
            // bases should meet only through the hydrogen bond term. Roles also
            // turn over as a structure moves: today's donor-to-donor contact is
            // tomorrow's bond when the partner rotates, so a rule written on
            // today's roles welds a pairing the run would otherwise have formed.
            //
            // Measured before widening it, on that 70S: 3335 of 17501 cutoff
            // springs refused instead of 1657, and the spring graph keeps
            // EXACTLY its 58 connected components with the same sizes and no
            // isolated particle. The cohesion a join provides was never carried
            // by these springs.
            veto = [&hbond](const topology::Particle & a, const topology::Particle & b) {
                const auto * ra = hbond.roleFor(a.properties().residue_name(), a.properties().name(),
                                                a.properties().original_name());
                const auto * rb = hbond.roleFor(b.properties().residue_name(), b.properties().name(),
                                                b.properties().original_name());
                if (ra == nullptr || rb == nullptr)
                    return false;
                const bool a_binds = ra->donorCapacity > 0 || ra->acceptorCapacity > 0;
                const bool b_binds = rb->donorCapacity > 0 || rb->acceptorCapacity > 0;
                return a_binds && b_binds;
            };
        }

        logging::status("Creating springs between structures with a %.2f Angstroms cutoff", args.cutoff);
        const size_t before = output.number_of_springs();
        output.add_springs_between_topologies_from_cutoff(args.cutoff, 1.0, veto);
        logging::info("    springs created: %zu", output.number_of_springs() - before);
        if (!args.useHydrogenBondTable())
            logging::warning("No donor/acceptor table given, so this cutoff may have welded pairs the dynamic "
                             "hydrogen bond term is meant to form. Pass --hbond with the table the simulation "
                             "will use to leave them alone.");
    }

    // Writes output files.
    io::writeTopology(args.pathOutputList, output);
    logging::info("Number of particles: %d", output.number_of_particles());
    logging::info("Number of springs: %d", output.number_of_springs());

    return EXIT_SUCCESS;
}

CommandLineArguments::CommandLineArguments(const std::string & name, const argparse::description_t & description,
                                           const std::string & version)
    : CommandLineArgumentsBase(name, description, version), pathTopologyList(), pathOutputList({"output.nc"}),
      cutoff(-1.0), pathHydrogenBond()
{
    argparse::Argument topology = argparse::Argument()
                                      .name_short("-s")
                                      .name_long("--topology")
                                      .description("input topology file(s); format is selected by extension: .nc, .pdb or .pqr")
                                      .metavar("INPUT_FILE")
                                      .number_of_arguments("+")
                                      .argument_type(argparse::ArgumentType::PATH_INPUT)
                                      .required(true);

    argparse::Argument output = argparse::Argument()
                                    .name_short("-o")
                                    .name_long("--output")
                                    .description("output file name(s); format is selected by extension: .nc, .cdl, .pdb or .pqr")
                                    .metavar("OUTPUT_FILE")
                                    .number_of_arguments("+")
                                    .argument_type(argparse::ArgumentType::PATH_OUTPUT)
                                    .default_value("system.nc");

    argparse::Argument cutoff = argparse::Argument()
                                    .name_short("-cutoff")
                                    .name_long("--cutoff")
                                    .description("cutoff in Angstroms for spring creation (< 0 means no spring)")
                                    .argument_type(argparse::ArgumentType::REAL)
                                    .default_value("-1.0");

    argparse::Argument hbond = argparse::Argument()
                                   .name_short("-hbond")
                                   .name_long("--hbond")
                                   .description("donor/acceptor table; a pair of hydrogen-bond sites is never welded "
                                                "by --cutoff")
                                   .metavar("HBOND_FILE")
                                   .argument_type(argparse::ArgumentType::PATH_INPUT);

    _parser.add_argument(topology);
    _parser.add_argument(output);
    _parser.add_argument(cutoff);
    _parser.add_argument(hbond);
}

void CommandLineArguments::parseCommandLine(int argc, const char * const argv[])
{
    _parser.parse_arguments(argc, argv);
    pathTopologyList = _parser.get_option_values<std::string>("-s");
    pathOutputList = _parser.get_option_values<std::string>("-o");
    cutoff = _parser.get_option_value<float>("--cutoff");
    if (_parser.get_option("--hbond").is_set())
        pathHydrogenBond = _parser.get_option_value<std::string>("--hbond");
}

void CommandLineArguments::printArgumentValues() const
{
    logging::status("Running biospring with arguments:");
    logging::info("    input topologies: ");
    for (const auto & path : pathTopologyList)
        logging::info("        - %s", path.c_str());

    if (pathOutputList.size() > 1)
    {
        logging::info("    output files:");
        for (const auto & path : pathOutputList)
        {
            logging::info("        - %s", path.c_str());
        }
    }
    else
    {
        logging::info("    output file: %s", pathOutputList[0].c_str());
    }
    if (useUserCutoff())
    {
        logging::info("    cutoff: %.2f", cutoff);
    }
    if (useHydrogenBondTable())
    {
        logging::info("    donor/acceptor table: %s", pathHydrogenBond.c_str());
    }
}

} // namespace mergespn
} // namespace biospring
