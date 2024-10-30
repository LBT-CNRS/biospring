

#include "pdb2spn-cli.h"
#include "IO/ForceFieldReader.h"
#include "IO/io.h"
#include "logging.h"
#include "reduce/Reducer.h"
#include "utils.hpp"

#include <string>
#include <vector>

const std::string PROGRAM_VERSION = "0.1.0";
const biospring::argparse::description_t PROGRAM_DESCRIPTION = {
    "pdb2spn creates a spring network from a topology file.",
};

namespace biospring
{
namespace pdb2spn
{

int main(int argc, char ** argv)
{
    CommandLineArguments args(std::string(argv[0]), PROGRAM_DESCRIPTION, PROGRAM_VERSION);
    args.parseCommandLine(argc, argv);
    args.printArgumentValues();

    // Reads topology and copies it to SpringNetwork.
    auto topology = biospring::io::readTopology(args.pathTopology);

    if (!args.pathGroup.empty())
    {
        std::string extension = utils::path::getExtension(args.pathTopology);
        if (extension == "nc")
        {
            logging::warning("topology initialized from nc file...ignoring group and force field arguments.");
        }
        else
        {
            reduceToCoarseGrain(topology, args);
        }
    }

    // If the PDB contains "CONECT" lines, some springs have already been parsed
    // in the readTopology function. We just set their stiffness here.
    if (topology.number_of_springs() > 0)
    {
        // Update all springs with stiffness
        for (size_t i=0; i<topology.number_of_springs(); i++)
        {
            auto & spring = topology.get_spring(i);
            spring.set_stiffness(args.stiffness);
        }
    }

    if (args.cutoff > 0)
    {
        logging::status("Creating springs with distance cutoff %.2f and stiffness %.2f.", args.cutoff, args.stiffness);
        topology.add_springs_from_cutoff(args.cutoff, args.stiffness);
        logging::info("Created %d spring(s).", topology.number_of_springs());
    }

    // Sets the particle charge to user-defined value.
    if (args.useUserCharge())
    {
        logging::status("Setting particle charge to user defined value: %.f.", args.charge);
        topology.particles().set_charges(args.charge);
    }

    // Sets particle dynamic state to false is user said so.
    if (args.isStatic)
    {
        logging::status("Setting particle mode to static.");
        topology.particles().set_dynamic(false);
    }

    // Writes output files.
    biospring::io::writeTopology(args.pathOutputList, topology);

    return EXIT_SUCCESS;
}

CommandLineArguments::CommandLineArguments(const std::string & name, const argparse::description_t & description,
                                           const std::string & version)
    : CommandLineArgumentsBase(name, description, version), pathTopology(""), pathForceField(""), pathGroup(""),
      pathOutputList(0), cutoff(-1.0), stiffness(1.0), charge(0.0), isStatic(false), ignoreDuplicates(false),
      ignoreMissing(false)
{
    argparse::Argument topology = argparse::Argument()
                                      .name_short("-s")
                                      .name_long("--topology")
                                      .description("input topology.")
                                      .metavar("INPUT_FILE")
                                      .argument_type(argparse::ArgumentType::PATH_INPUT)
                                      .required(true);

    argparse::Argument output = argparse::Argument()
                                    .name_short("-o")
                                    .name_long("--output")
                                    .description("output file name(s).")
                                    .metavar("OUTPUT_FILE")
                                    .number_of_arguments("+")
                                    .argument_type(argparse::ArgumentType::PATH_OUTPUT)
                                    .default_value("system.nc");

    argparse::Argument forcefield = argparse::Argument()
                                        .name_long("--ff")
                                        .description("force field file")
                                        .metavar("INPUT_FILE")
                                        .argument_type(argparse::ArgumentType::PATH_INPUT);

    argparse::Argument grp = argparse::Argument()
                                 .name_long("--grp")
                                 .description("particles definition file")
                                 .metavar("INPUT_FILE")
                                 .argument_type(argparse::ArgumentType::PATH_INPUT);

    argparse::Argument cutoff = argparse::Argument()
                                    .name_long("--cutoff")
                                    .description("cutoff for spring creation (< 0 means no spring)")
                                    .argument_type(argparse::ArgumentType::REAL)
                                    .default_value("-1.0");

    argparse::Argument stiffness = argparse::Argument()
                                       .name_long("--stiffness")
                                       .description("spring stiffness")
                                       .argument_type(argparse::ArgumentType::REAL)
                                       .default_value("1.0");

    argparse::Argument charge = argparse::Argument()
                                    .name_long("--charge")
                                    .description("override force field values for particle charge")
                                    .default_value("0.0")
                                    .argument_type(argparse::ArgumentType::REAL);

    argparse::Argument static_ =
        argparse::StoreTrueArgument("", "--static", "should the particles be freezed during the simulation");

    argparse::Argument ignore_duplicate = argparse::StoreTrueArgument(
        "", "--ignore-duplicate", "ignore duplicate particles when reducing to coarse grain");

    argparse::Argument ignore_missing =
        argparse::StoreTrueArgument("", "--ignore-missing", "ignore missing particles when reducing to coarse grain");

    _parser.add_argument(topology);
    _parser.add_argument(output);
    _parser.add_argument(forcefield);
    _parser.add_argument(grp);
    _parser.add_argument(cutoff);
    _parser.add_argument(stiffness);
    _parser.add_argument(charge);
    _parser.add_argument(static_);
    _parser.add_argument(ignore_duplicate);
    _parser.add_argument(ignore_missing);
}

void CommandLineArguments::parseCommandLine(int argc, const char * const argv[])
{
    _parser.parse_arguments(argc, argv);
    pathTopology = _parser.get_option_value<std::string>("--topology");
    pathOutputList = _parser.get_option_value<std::vector<std::string>>("--output");
    pathForceField = _parser.get_option_value<std::string>("--ff");
    pathGroup = _parser.get_option_value<std::string>("--grp");

    try
    {
        cutoff = _parser.get_option_value<float>("--cutoff");
    }
    catch (const std::invalid_argument & e)
    {
        std::string message = "Invalid argument: --cutoff must be a float (got '" +
                              _parser.get_option_value<std::string>("--cutoff") + "')";
        _parser.die(message);
    }

    try
    {
        charge = _parser.get_option_value<float>("--charge");
    }
    catch (const std::exception & e)
    {
        std::string message = "Invalid argument: --charge must be a float (got '" +
                              _parser.get_option_value<std::string>("--charge") + "')";
        _parser.die(message);
    }

    try
    {
        stiffness = _parser.get_option_value<float>("--stiffness");
    }
    catch (const std::exception & e)
    {
        std::string message = "Invalid argument: --stiffness must be a float (got '" +
                              _parser.get_option_value<std::string>("--stiffness") + "')";
        _parser.die(message);
    }

    isStatic = _parser.get_option("--static").is_set();
    ignoreDuplicates = _parser.get_option("--ignore-duplicate").is_set();
    ignoreMissing = _parser.get_option("--ignore-missing").is_set();

    // Reduce file and force field should be provided together.
    // Dies if not the case.
    if (!pathGroup.empty() && pathForceField.empty())
    {
        _parser.print_help();
        _parser.die("--ff <file> is mandatory when --grp is provided");
    }
    if (!pathForceField.empty() && pathGroup.empty())
    {
        _parser.print_help();
        _parser.die("--grp <file> is mandatory when --ff is provided");
    }
}

void biospring::pdb2spn::CommandLineArguments::printArgumentValues() const
{
    logging::status("Running pdb2spn with arguments:");
    logging::info("    topology: %s", pathTopology.c_str());
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
    if (cutoff < 0)
        logging::info("    cutoff: %.1f (no spring will be created)", cutoff);
    else
        logging::info("    cutoff: %.1f", cutoff);

    logging::info("    spring stiffness: %.2f", stiffness);
    logging::info("    particle charge: %.2f", charge);

    if (isStatic)
        logging::info("    static particles: yes");
    else
        logging::info("    static particles: no");

    if (pathGroup.size())
    {
        logging::info("    coarse grain reduction: true");
        logging::info("        reduce file: %s", pathGroup.c_str());
        logging::info("        force field file: %s", pathForceField.c_str());
        if (ignoreDuplicates)
            logging::info("        ignore duplicate particles: yes");
        else
            logging::info("        ignore duplicate particles: no");
        if (ignoreMissing)
            logging::info("        ignore missing particles: yes");
        else
            logging::info("        ignore missing particles: no");
    }
    else
    {
        logging::info("    coarse grain reduction: false");
    }
}

} // namespace pdb2spn
} // namespace biospring

//
// All-atom to coarse grain reduction.
// Modifies the topology in place.
//
void biospring::pdb2spn::reduceToCoarseGrain(topology::Topology & top, const CommandLineArguments & args)
{
    biospring::reduce::ReductionParameters params;
    params.pathForceField = args.pathForceField;
    params.pathGroup = args.pathGroup;
    params.ignoreDuplicate = args.ignoreDuplicates;
    params.ignoreMissing = args.ignoreMissing;

    logging::status("Reducing input topology to coarse grain");
    biospring::reduce::Reducer reducer(top);
    reducer.reduce(params);
    top = reducer.target_topology();
    logging::info("Number of particles after reduction: %d", top.number_of_particles());
}
