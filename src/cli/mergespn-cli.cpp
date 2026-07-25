
#include "mergespn-cli.h"
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
        logging::status("Creating springs between structures with a %.2f Angstroms cutoff", args.cutoff);
        output.add_springs_between_topologies_from_cutoff(args.cutoff);
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
      cutoff(-1.0)
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
                                    .description("cutoff for spring creation (< 0 means no spring)")
                                    .argument_type(argparse::ArgumentType::REAL)
                                    .default_value("-1.0");

    _parser.add_argument(topology);
    _parser.add_argument(output);
    _parser.add_argument(cutoff);
}

void CommandLineArguments::parseCommandLine(int argc, const char * const argv[])
{
    _parser.parse_arguments(argc, argv);
    pathTopologyList = _parser.get_option_values<std::string>("-s");
    pathOutputList = _parser.get_option_values<std::string>("-o");
    cutoff = _parser.get_option_value<float>("--cutoff");
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
}

} // namespace mergespn
} // namespace biospring
