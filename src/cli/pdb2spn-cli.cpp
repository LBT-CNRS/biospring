

#include "pdb2spn-cli.h"
#include "IO/ForceFieldReader.h"
#include "IO/ReduceRuleReader.h"
#include "IO/RigidBodyRuleReader.h"
#include "IO/io.h"
#include "logging.h"
#include "reduce/Reducer.h"
#include "rigidbodygroup/RigidBodyBuilder.h"
#include "utils.hpp"

#include <string>
#include <vector>

const std::string PROGRAM_VERSION = "0.1.0";
const biospring::argparse::description_t PROGRAM_DESCRIPTION = {
    "pdb2spn creates a spring network from a topology file.",
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
    "When the output is a PDB file, -pdbconect/--pdbconect adds CONECT",
    "records at the end of the PDB to visualize the springs of the network.",
    "",
    "-rigidbody/--rigidbody creates springs from a .rbody rigid-body-group",
    "file instead of a distance cutoff: every pairwise spring within each",
    "matching group is created directly (see data/reducerules/*.rbody).",
    "Exclusive with -cutoff/--cutoff -- the two are incompatible strategies",
    "for building the spring network. Uses -stiffness/--stiffness as-is.",
    "If -grp/--grp is also given, it is reused to translate renamed atoms",
    "back to their original name (it must be an all-atom identity mapping,",
    "one atom per rule, like amber.grp -- not a real coarse-grain reduction).",
    "",
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

    if (!args.pathRigidBody.empty())
    {
        logging::status("Creating springs from rigid-body groups %s with stiffness %.2f.", args.pathRigidBody.c_str(),
                         args.stiffness);
        biospring::rigidbodygroup::RigidBodyRuleReader reader(args.pathRigidBody);
        reader.read();

        // If --grp renamed the particles (e.g. amber.grp: CA -> ACA for
        // ALA), re-read that same file as a naming-translation table so
        // atoms can still be found by their original name.
        biospring::reduce::ReduceRuleReader naming_reader;
        biospring::reduce::ReduceRuleReader * naming_reader_ptr = nullptr;
        if (!args.pathGroup.empty())
        {
            naming_reader.setFileName(args.pathGroup);
            naming_reader.read();
            naming_reader_ptr = &naming_reader;
        }

        biospring::rigidbodygroup::RigidBodyBuilder builder(topology, reader.rules(), args.stiffness,
                                                              naming_reader_ptr != nullptr ? &naming_reader_ptr->rules()
                                                                                            : nullptr);
        builder.build();
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
    biospring::io::writeTopology(args.pathOutputList, topology, args.writePdbConect);

    return EXIT_SUCCESS;
}

CommandLineArguments::CommandLineArguments(const std::string & name, const argparse::description_t & description,
                                           const std::string & version)
    : CommandLineArgumentsBase(name, description, version), pathTopology(""), pathForceField(""), pathGroup(""),
      pathRigidBody(""), pathOutputList(0), cutoff(-1.0), stiffness(1.0), charge(0.0),
      isStatic(false), ignoreDuplicates(false), ignoreMissing(false), writePdbConect(false)
{
    argparse::Argument topology = argparse::Argument()
                                      .name_short("-s")
                                      .name_long("--topology")
                                      .description("input topology file; format is selected by extension: .nc, .pdb or .pqr")
                                      .metavar("INPUT_FILE")
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

    argparse::Argument forcefield = argparse::Argument()
                                        .name_short("-ff")
                                        .name_long("--ff")
                                        .description("non-bonded interaction parameter file: per-atom-type charge/"
                                                      "radius/epsilon/mass (steric, electrostatic, hydrophobicity, "
                                                      "IMP), a .ff file")
                                        .metavar("INPUT_FILE")
                                        .argument_type(argparse::ArgumentType::PATH_INPUT);

    argparse::Argument grp = argparse::Argument()
                                 .name_short("-grp")
                                 .name_long("--grp")
                                 .description("particles definition file")
                                 .metavar("INPUT_FILE")
                                 .argument_type(argparse::ArgumentType::PATH_INPUT);

    argparse::Argument cutoff = argparse::Argument()
                                    .name_short("-cutoff")
                                    .name_long("--cutoff")
                                    .description("cutoff in Angstroms for spring creation (< 0 means no spring)")
                                    .argument_type(argparse::ArgumentType::REAL)
                                    .default_value("-1.0");

    argparse::Argument rigidbody = argparse::Argument()
                                       .name_short("-rigidbody")
                                       .name_long("--rigidbody")
                                       .description("create springs from a .rbody rigid-body-group file instead of "
                                                     "a cutoff; exclusive with -cutoff/--cutoff")
                                       .metavar("INPUT_FILE")
                                       .argument_type(argparse::ArgumentType::PATH_INPUT);


    argparse::Argument stiffness = argparse::Argument()
                                       .name_short("-stiffness")
                                       .name_long("--stiffness")
                                       .description("spring stiffness, in kJ.mol-1.A-2")
                                       .argument_type(argparse::ArgumentType::REAL)
                                       .default_value("1.0");

    argparse::Argument charge = argparse::Argument()
                                    .name_short("-charge")
                                    .name_long("--charge")
                                    .description("override force field values for particle charge, in elementary charge units (e)")
                                    .default_value("0.0")
                                    .argument_type(argparse::ArgumentType::REAL);

    argparse::Argument static_ =
        argparse::StoreTrueArgument("-static", "--static", "should the particles be freezed during the simulation");

    argparse::Argument ignore_duplicate = argparse::StoreTrueArgument(
        "-ignore-duplicate", "--ignore-duplicate", "ignore duplicate particles when reducing to coarse grain");

    argparse::Argument ignore_missing =
        argparse::StoreTrueArgument("-ignore-missing", "--ignore-missing", "ignore missing particles when reducing to coarse grain");

    argparse::Argument pdbconect = argparse::StoreTrueArgument(
        "-pdbconect", "--pdbconect",
        "when writing PDB output, add CONECT records for the springs at the end of the PDB file");

    _parser.add_argument(topology);
    _parser.add_argument(output);
    _parser.add_argument(forcefield);
    _parser.add_argument(grp);
    _parser.add_argument(cutoff);
    _parser.add_argument(rigidbody);
    _parser.add_argument(stiffness);
    _parser.add_argument(charge);
    _parser.add_argument(static_);
    _parser.add_argument(ignore_duplicate);
    _parser.add_argument(ignore_missing);
    _parser.add_argument(pdbconect);
}

void CommandLineArguments::parseCommandLine(int argc, const char * const argv[])
{
    _parser.parse_arguments(argc, argv);
    pathTopology = _parser.get_option_value<std::string>("--topology");
    pathOutputList = _parser.get_option_value<std::vector<std::string>>("--output");
    pathForceField = _parser.get_option_value<std::string>("--ff");
    pathGroup = _parser.get_option_value<std::string>("--grp");
    pathRigidBody = _parser.get_option_value<std::string>("--rigidbody");

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
    writePdbConect = _parser.get_option("--pdbconect").is_set();

    // --rigidbody and --cutoff are two incompatible strategies for building
    // the spring network: dies if both are given.
    if (!pathRigidBody.empty() && useUserCutoff() && cutoff > 0)
    {
        _parser.print_help();
        _parser.die("-rigidbody/--rigidbody is exclusive with -cutoff/--cutoff");
    }



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

    if (writePdbConect)
    {
        bool hasPdbOutput = false;
        for (const auto & path : pathOutputList)
        {
            if (biospring::utils::path::getExtension(path) == "pdb")
            {
                hasPdbOutput = true;
                break;
            }
        }

        if (!hasPdbOutput)
        {
            _parser.print_help();
            _parser.die("--pdbconect requires at least one PDB output file in -o/--output");
        }
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
    if (!pathRigidBody.empty())
    {
        logging::info("    rigid-body groups: %s", pathRigidBody.c_str());
    }
    else if (cutoff < 0)
        logging::info("    cutoff: %.1f (no spring will be created)", cutoff);
    else
        logging::info("    cutoff: %.1f", cutoff);

    logging::info("    spring stiffness: %.2f", stiffness);
    logging::info("    particle charge: %.2f", charge);

    if (isStatic)
        logging::info("    static particles: yes");
    else
        logging::info("    static particles: no");

    if (writePdbConect)
        logging::info("    PDB CONECT records: enabled for PDB output file(s)");
    else
        logging::info("    PDB CONECT records: disabled");

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
