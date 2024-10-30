
#include <fstream>
#include <sstream>

#include "IO/io.h"
#include "editspn-cli.h"

namespace biospring
{
namespace editspn
{

const std::string PROGRAM_VERSION = "0.1.0";
const argparse::description_t PROGRAM_DESCRIPTION = {
    "editspn is a tool that aims to edit and create spring networks.",
};

int main(int argc, char ** argv)
{
    CommandLineArguments args(std::string(argv[0]), PROGRAM_DESCRIPTION, PROGRAM_VERSION);
    args.parseCommandLine(argc, argv);
    args.printArgumentValues();

    // Reads the topology.
    auto topology = biospring::io::readTopology(args.pathTopology);

    if (args.useUserCutoff())
    {
        logging::status("Creating springs from user-defined cutoff distance %.2f.", args.cutoff);
        if (args.useUserStiffness())
        {
            logging::info("Using user-defined stiffness %.2f.", args.stiffness);
            topology.springs().set_stiffnesses(args.stiffness);
        }
        topology.add_springs_from_cutoff(args.cutoff);
    }

    // if (not args.pathUpdateSprings.empty())
    // {
    //     logging::status("Updating spring parameters using spring file %s.", args.pathUpdateSprings.c_str());
    //     addSpringsToSpn(args.pathUpdateSprings, spn);
    // }

    // logging::info("Final spring network has %d particles.", spn.getNumberOfParticles());
    // logging::info("Final spring network has %d springs.", spn.getNumberOfSprings());

    // biospring::io::writeTopology(args.pathOutputList, &spn);

    return EXIT_SUCCESS;
}

//
// Reads a spring description text file and store the descriptors in sdesc.
//
void readSpringFile(const std::string & filename, vector<SpringDescriptor> & sdesc)
{
    std::ifstream filein;
    biospring::utils::file::openread(filename, filein);

    istringstream iss;

    for (unsigned lineid = 1; not filein.eof(); ++lineid)
    {
        std::string buffer;
        std::getline(filein, buffer);

        if (not(buffer.empty() or buffer[0] == '#'))
        {
            istringstream iss(buffer);
            std::string id1, id2, sstiff, sequi;
            iss >> id1 >> id2 >> sstiff >> sequi;

            if (sequi.empty())
            {
                logging::die("spring format error at line %d: \"%s\" (expected: <pidentifier1> <pidentifier2> "
                             "<stiffness> <equilibrium>)",
                             lineid, buffer.c_str());
            }

            checkParticleIdFormat(id1, lineid);
            checkParticleIdFormat(id2, lineid);

            SpringDescriptor sd;
            sd.p1 = id1;
            sd.p2 = id2;
            sd.stiffness = atof(sstiff.c_str());
            sd.equilibrium = atof(sequi.c_str());
            sdesc.push_back(sd);
        }
    }
}

//
// Adds springs described in sdesc to the topology.
// If some springs already exist, change their property to those described in
// sdesc.
//
void addSpringsToSpn(const vector<SpringDescriptor> & springDescriptors, biospring::topology::Topology & top)
{
    std::unordered_map<std::string, unsigned> pidmap;

    // First, maps all the particles std::string descriptor and ids
    for (size_t pid = 0; pid < top.number_of_particles(); pid++)
    {
        pidmap[top.get_particle(pid).string_description()] = pid;
    }

    // Then, for each spring descriptor, if the spring does not exist, create it
    // otherwise change its properties.
    for (const SpringDescriptor & descriptor : springDescriptors)
    {
        // Checks particles are actually known to the spring network.
        checkParticleInSpn(descriptor.p1, pidmap);
        checkParticleInSpn(descriptor.p2, pidmap);

        // If so, retrieves them.
        topology::Particle & p1 = top.get_particle(pidmap[descriptor.p1]);
        topology::Particle & p2 = top.get_particle(pidmap[descriptor.p2]);

        // If a spring exists between the two particles, updates its properties.
        if (top.has_spring_between(p1, p2))
        {
            topology::Spring & s = top.get_spring(p1, p2);
            s.set_stiffness(descriptor.stiffness);

            if (descriptor.equilibrium >= 0.0)
                s.set_equilibrium(descriptor.equilibrium);
        }

        // Else, adds a new spring with ad-hoc properties.
        else
        {
            if (descriptor.equilibrium == -1)
                top.add_spring(p1, p2, -1, descriptor.stiffness);

            else
                top.add_spring(p1, p2, descriptor.equilibrium, descriptor.stiffness);
        }
    }
}

//
// Adds springs described in file to the topology.
// If some springs already exist, change their property to those described in
// sdesc.
//
void addSpringsToSpn(const std::string & path, biospring::topology::Topology & top)
{
    std::vector<SpringDescriptor> descriptors;
    readSpringFile(path, descriptors);
    addSpringsToSpn(descriptors, top);
}

//
// Dies if the particle id does not have 4 tokens separated by "::"
//
void checkParticleIdFormat(const std::string & pid, const unsigned & lineid)
{
    vector<std::string> tokens;
    std::string delimiters = "::";

    std::string::size_type lastPos = pid.find_first_not_of(delimiters, 0);
    std::string::size_type pos = pid.find_first_of(delimiters, lastPos);

    while (std::string::npos != pos || std::string::npos != lastPos)
    {
        tokens.push_back(pid.substr(lastPos, pos - lastPos));
        lastPos = pid.find_first_not_of(delimiters, pos);
        pos = pid.find_first_of(delimiters, lastPos);
    }

    if (tokens.size() != 4)
    {
        logging::die(
            "spring format error at line %d: misformatted particle id: '%s' (expected chain::resname::resid::name)",
            lineid, pid.c_str());
    }
}

//
// Dues if the particle id is not in the map containing all particle ids
//
void checkParticleInSpn(const std::string & pid, const unordered_map<std::string, unsigned> & pidmap)
{
    if (pidmap.count(pid) == 0)
    {
        logging::die("spring file error: particle \"%s\" not found in the spring network.", pid.c_str());
    }
}

////////////////////////////////////////////////////////////////////////////////////////
//
// Command-line parsing.
//
////////////////////////////////////////////////////////////////////////////////////////

CommandLineArguments::CommandLineArguments(const std::string & name, const argparse::description_t & description,
                                           const std::string & version)
    : CommandLineArgumentsBase(name, description, version), pathTopology(""), pathOutputList(0), cutoff(-1.0),
      stiffness(0.0)
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

    argparse::Argument springs = argparse::Argument()
                                     .name_short("-u")
                                     .name_long("--update-springs")
                                     .description("add springs defined in file.")
                                     .metavar("SPRING_FILE")
                                     .argument_type(argparse::ArgumentType::PATH_INPUT);

    argparse::Argument cutoff = argparse::Argument()
                                    .name_short("-c")
                                    .name_long("--cutoff")
                                    .description("cutoff for spring creation (< 0 means no spring)")
                                    .argument_type(argparse::ArgumentType::REAL)
                                    .default_value("-1.0");

    argparse::Argument stiffness = argparse::Argument()
                                       .name_short("-d")
                                       .name_long("--stiffness")
                                       .description("when used with -c, sets the spring stiffness")
                                       .argument_type(argparse::ArgumentType::REAL)
                                       .default_value("-1.0");

    _parser.add_argument(topology);
    _parser.add_argument(output);
    _parser.add_argument(springs);
    _parser.add_argument(cutoff);
    _parser.add_argument(stiffness);
}

void CommandLineArguments::parseCommandLine(int argc, const char * const argv[])
{
    _parser.parse_arguments(argc, argv);
    pathTopology = _parser.get_option_value<std::string>("--topology");
    pathOutputList = _parser.get_option_value<std::vector<std::string>>("--output");
    pathUpdateSprings = _parser.get_option_value<std::string>("--update-springs");
    cutoff = _parser.get_option_value<float>("--cutoff");
    stiffness = _parser.get_option_value<float>("--stiffness");
}

void CommandLineArguments::printArgumentValues() const
{
    logging::status("Running editspn with arguments:");
    logging::info("    topology: %s", pathTopology.c_str());
    if (not pathUpdateSprings.empty())
        logging::info("    spring update file: %s", pathUpdateSprings.c_str());
    if (useUserCutoff())
        logging::info("    cutoff: %.1f", cutoff);
    if (useUserStiffness())
        logging::info("    stiffness: %.1f", stiffness);
    if (pathOutputList.size() > 1)
    {
        logging::info("    output files: %s", pathTopology.c_str());
        for (const auto & path : pathOutputList)
        {
            logging::info("        - %s", path.c_str());
        }
    }
    else
    {
        logging::info("    output file: %s", pathOutputList[0].c_str());
    }
}

} // namespace editspn
} // namespace biospring
