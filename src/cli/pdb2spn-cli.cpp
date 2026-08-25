

#include "pdb2spn-cli.h"
#include "IO/ForceFieldReader.h"
#include "IO/ReduceRuleReader.h"
#include "IO/BondedForceFieldReader.h"
#include "IO/RigidBodyRuleReader.h"
#include "IO/io.h"
#include "logging.h"
#include "reduce/Reducer.h"
#include "rigidbodygroup/RigidBodyBuilder.h"
#include "staticbond/StaticBondBuilder.h"
#include "utils.hpp"

#include <optional>
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
    "-bondedinteraction/--bondedinteraction reads a .bi.ff file (real or",
    "ghost-spring parameters translated from a real force field, e.g. AMBER)",
    "on top of -rigidbody/--rigidbody's springs, but applies none of it by",
    "itself: -stretching/--stretching (real 1-2 bonds retuned in place),",
    "-bending/--bending (real 1-3 angles retuned in place, requires",
    "-stretching too), -dihedralbackbone/--dihedralbackbone and",
    "-dihedralsidechain/--dihedralsidechain (dihedral ghost springs added on",
    "top, one flag per proper-dihedral family) each independently opt in;",
    "-dihedral/--dihedral is a shorthand for both dihedral flags together.",
    "Combine only the ones you want -- e.g. -dihedral alone on top of plain",
    "-rigidbody keeps bonds/angles at the uniform -stiffness value while",
    "still adding real dihedral wells.",
    "",
    "-static-hbond/--static-hbond and -static-disulfide/--static-disulfide add",
    "springs for the hydrogen bonds and disulfide bridges ALREADY PRESENT in",
    "the input structure, calibrated to the force constant of the real bond.",
    "Both are additive: they compose with -cutoff, with -rigidbody, or with a",
    "network that has neither, and never replace a spring a pair already has.",
    "A static spring holds a known structure together; it cannot break, and it",
    "cannot form a bond the input did not contain.",
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

    // Reads the .bi.ff file (parsing only, no Spring/particle interaction
    // yet) and reserves enough particle capacity for every ghost particle
    // it might create, BEFORE any Spring exists anywhere in `topology`
    // (cutoff-, CONECT-, or --rigidbody-created -- all of those hold a
    // real Particle& reference, not an index; see
    // Topology::reserve_particles's own comment). Building ghost particles
    // this early would be wrong (their anchors need --rigidbody/--grp
    // resolution first), but reserving capacity only needs to know an
    // upper bound on how many will eventually be added, which the parsed
    // (not yet applied) file already tells us.
    std::optional<biospring::rigidbodygroup::BondedForceFieldReader> bondedReader;
    if (!args.pathBondedInteraction.empty())
    {
        bondedReader.emplace(args.pathBondedInteraction);
        bondedReader->read();
        size_t expected_ghosts = bondedReader->countExpectedGhostParticles(topology);
        topology.reserve_particles(topology.number_of_particles() + expected_ghosts);
        logging::info("Reserved capacity for up to %zu ghost particle(s) before any spring is created.",
                      expected_ghosts);
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

    if (!args.pathBondedInteraction.empty())
    {
        logging::status("Applying bonded interaction parameters from %s (stretching=%s, bending=%s, "
                        "dihedralbackbone=%s, dihedralsidechain=%s, dihedralplanarity=%s).",
                         args.pathBondedInteraction.c_str(), args.stretching ? "on" : "off",
                         args.bending ? "on" : "off", args.dihedralBackbone ? "on" : "off",
                         args.dihedralSidechain ? "on" : "off", args.dihedralPlanarity ? "on" : "off");
        // Already constructed and read() above (before any spring existed,
        // see that call site) -- reused here, not re-read.

        // Same naming-translation table as --rigidbody above (see there):
        // the .bi.ff file uses original atom names even if --grp renamed
        // the particles.
        biospring::reduce::ReduceRuleReader naming_reader;
        biospring::reduce::ReduceRuleReader * naming_reader_ptr = nullptr;
        if (!args.pathGroup.empty())
        {
            naming_reader.setFileName(args.pathGroup);
            naming_reader.read();
            naming_reader_ptr = &naming_reader;
        }

        bondedReader->buildSprings(topology, naming_reader_ptr != nullptr ? &naming_reader_ptr->rules() : nullptr,
                                  args.stretching, args.bending, args.dihedralBackbone, args.dihedralSidechain,
                                  args.dihedralPlanarity);
    }

    // Hydrogen bonds and disulfides as ordinary springs. Deliberately AFTER
    // both network builders and independent of either: these describe bonds
    // the structure already has, whatever strategy built the rest of the
    // network -- or none at all.
    if (!args.pathStaticHydrogenBond.empty())
    {
        logging::status("Adding static hydrogen-bond springs from %s.", args.pathStaticHydrogenBond.c_str());
        const auto table = biospring::staticbond::readDonorAcceptorTable(args.pathStaticHydrogenBond);
        biospring::staticbond::addHydrogenBondSprings(topology, table, biospring::staticbond::HYDROGEN_BOND_CUTOFF,
                                                      biospring::staticbond::HYDROGEN_BOND_STIFFNESS);
    }

    if (args.addStaticDisulfide)
    {
        logging::status("Adding static disulfide springs.");
        biospring::staticbond::addDisulfideSprings(topology, biospring::staticbond::DISULFIDE_CUTOFF,
                                                   biospring::staticbond::DISULFIDE_STIFFNESS);
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
      pathRigidBody(""), pathBondedInteraction(""), pathStaticHydrogenBond(""), pathOutputList(0), cutoff(-1.0),
      stiffness(1.0), charge(0.0), isStatic(false), addStaticDisulfide(false), ignoreDuplicates(false),
      ignoreMissing(false), writePdbConect(false), stretching(false), bending(false), dihedral(false),
      dihedralBackbone(false), dihedralSidechain(false), dihedralPlanarity(false)
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
                                                      "IMP), a .ff or .nbi.ff file (synonym: -nonbondedinteraction/"
                                                      "--nonbondedinteraction)")
                                        .metavar("INPUT_FILE")
                                        .argument_type(argparse::ArgumentType::PATH_INPUT);

    // Synonym for --ff, spelled out for clarity against --bondedinteraction:
    // both write to pathForceField (see parseCommandLine), specifying both
    // at once is a conflict.
    argparse::Argument nonbondedinteraction =
        argparse::Argument()
            .name_short("-nonbondedinteraction")
            .name_long("--nonbondedinteraction")
            .description("synonym for -ff/--ff: a .ff or .nbi.ff file with non-bonded interaction parameters "
                          "(charge/radius/epsilon/mass -- steric, electrostatic, hydrophobicity, IMP)")
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

    argparse::Argument bondedinteraction =
        argparse::Argument()
            .name_short("-bondedinteraction")
            .name_long("--bondedinteraction")
            .description("bonded interaction parameter file: real or virtual/ghost spring parameters (stretching, "
                          "bending, dihedral) translated from a real force field (e.g. AMBER), a .bi.ff file. "
                          "Applies none of it by itself -- see -stretching/-bending/-dihedralbackbone/"
                          "-dihedralsidechain (or -dihedral, a shorthand for both dihedral flags), each an "
                          "independent opt-in; requires -rigidbody/--rigidbody")
            .metavar("INPUT_FILE")
            .argument_type(argparse::ArgumentType::PATH_INPUT);

    argparse::Argument stretching_ = argparse::StoreTrueArgument(
        "-stretching", "--stretching",
        "with -bondedinteraction, retune STRETCH springs (real 1-2 bond lengths/stiffness) in place; has no effect "
        "without -bondedinteraction");

    argparse::Argument bending_ = argparse::StoreTrueArgument(
        "-bending", "--bending",
        "with -bondedinteraction, retune BEND springs (real 1-3 distance-spring equivalent) in place; requires "
        "-stretching too (bending's conversion needs STRETCH's real bond lengths); has no effect without "
        "-bondedinteraction");

    argparse::Argument dihedral_ = argparse::StoreTrueArgument(
        "-dihedral", "--dihedral",
        "with -bondedinteraction, add all dihedral ghost springs on top of -rigidbody's mesh -- shorthand for "
        "-dihedralbackbone, -dihedralsidechain and -dihedralplanarity together; independent of -stretching/-bending -- can be combined "
        "alone for the intermediate model (bonds/angles stay at -rigidbody's uniform value, only dihedral wells "
        "become real); has no effect without -bondedinteraction");

    argparse::Argument dihedralbackbone_ = argparse::StoreTrueArgument(
        "-dihedralbackbone", "--dihedralbackbone",
        "with -bondedinteraction, add backbone (phi/psi/omega) dihedral ghost springs only; has no effect without "
        "-bondedinteraction");

    argparse::Argument dihedralsidechain_ = argparse::StoreTrueArgument(
        "-dihedralsidechain", "--dihedralsidechain",
        "with -bondedinteraction, add side-chain (chi1-4) dihedral ghost springs only; has no effect without "
        "-bondedinteraction");

    argparse::Argument dihedralplanarity_ = argparse::StoreTrueArgument(
        "-dihedralplanarity", "--dihedralplanarity",
        "with -bondedinteraction, add PLANARITY improper ghost springs (aromatic-ring/His hub planarity) only; "
        "has no effect without -bondedinteraction");

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

    argparse::Argument static_hbond =
        argparse::Argument()
            .name_short("-static-hbond")
            .name_long("--static-hbond")
            .description("add a spring for every hydrogen bond already present in the structure, using the "
                         "given .hbond donor/acceptor table; additive to -cutoff and -rigidbody")
            .metavar("INPUT_FILE")
            .argument_type(argparse::ArgumentType::PATH_INPUT);

    argparse::Argument static_disulfide = argparse::StoreTrueArgument(
        "-static-disulfide", "--static-disulfide",
        "add a spring for every disulfide bridge already present in the structure (cysteine sulfurs within "
        "2.5 A); additive to -cutoff and -rigidbody");

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
    _parser.add_argument(nonbondedinteraction);
    _parser.add_argument(grp);
    _parser.add_argument(cutoff);
    _parser.add_argument(rigidbody);
    _parser.add_argument(bondedinteraction);
    _parser.add_argument(stretching_);
    _parser.add_argument(bending_);
    _parser.add_argument(dihedral_);
    _parser.add_argument(dihedralbackbone_);
    _parser.add_argument(dihedralsidechain_);
    _parser.add_argument(dihedralplanarity_);
    _parser.add_argument(stiffness);
    _parser.add_argument(charge);
    _parser.add_argument(static_);
    _parser.add_argument(static_hbond);
    _parser.add_argument(static_disulfide);
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
    {
        // --nonbondedinteraction is a plain synonym for --ff (see argument
        // definitions above): merge it in, dying on a conflict rather than
        // silently picking one.
        const std::string synonym = _parser.get_option_value<std::string>("--nonbondedinteraction");
        if (!synonym.empty())
        {
            if (!pathForceField.empty())
            {
                _parser.print_help();
                _parser.die("-ff/--ff and -nonbondedinteraction/--nonbondedinteraction are synonyms -- specify "
                            "only one");
            }
            pathForceField = synonym;
        }
    }
    pathGroup = _parser.get_option_value<std::string>("--grp");
    pathRigidBody = _parser.get_option_value<std::string>("--rigidbody");
    pathBondedInteraction = _parser.get_option_value<std::string>("--bondedinteraction");

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

    pathStaticHydrogenBond = _parser.get_option_value<std::string>("--static-hbond");
    addStaticDisulfide = _parser.get_option("--static-disulfide").is_set();
    ignoreDuplicates = _parser.get_option("--ignore-duplicate").is_set();
    ignoreMissing = _parser.get_option("--ignore-missing").is_set();
    writePdbConect = _parser.get_option("--pdbconect").is_set();
    stretching = _parser.get_option("--stretching").is_set();
    bending = _parser.get_option("--bending").is_set();
    dihedral = _parser.get_option("--dihedral").is_set();
    dihedralBackbone = _parser.get_option("--dihedralbackbone").is_set();
    dihedralSidechain = _parser.get_option("--dihedralsidechain").is_set();
    dihedralPlanarity = _parser.get_option("--dihedralplanarity").is_set();

    // --rigidbody and --cutoff are two incompatible strategies for building
    // the spring network: dies if both are given.
    if (!pathRigidBody.empty() && useUserCutoff() && cutoff > 0)
    {
        _parser.print_help();
        _parser.die("-rigidbody/--rigidbody is exclusive with -cutoff/--cutoff");
    }

    // -bondedinteraction only retunes/extends -rigidbody's springs in place;
    // it does not build an independent spring network, so it cannot be used
    // without it.
    if (!pathBondedInteraction.empty() && pathRigidBody.empty())
    {
        _parser.print_help();
        _parser.die("-bondedinteraction/--bondedinteraction requires -rigidbody/--rigidbody");
    }

    // -stretching/-bending/-dihedral[backbone|sidechain] only mean something
    // alongside -bondedinteraction, and apply nothing unless at least one is
    // given. -dihedral is a pure convenience alias, resolved here rather
    // than passed down as its own concept: BondedForceFieldReader only ever
    // sees the two resolved per-family booleans.
    dihedralBackbone = dihedralBackbone || dihedral;
    dihedralSidechain = dihedralSidechain || dihedral;
    dihedralPlanarity = dihedralPlanarity || dihedral;

    if ((stretching || bending || dihedralBackbone || dihedralSidechain || dihedralPlanarity) && pathBondedInteraction.empty())
    {
        _parser.print_help();
        _parser.die("-stretching/-bending/-dihedral/-dihedralbackbone/-dihedralsidechain require "
                    "-bondedinteraction/--bondedinteraction");
    }
    if (!pathBondedInteraction.empty() && !stretching && !bending && !dihedralBackbone && !dihedralSidechain && !dihedralPlanarity)
        logging::warning("-bondedinteraction/--bondedinteraction given without -stretching/-bending/-dihedral* -- "
                         "nothing from the .bi.ff file will actually be applied.");

    // Bending's 1-3 conversion needs STRETCH's real bond lengths (see
    // BondedForceFieldReader::buildSprings) -- requesting BEND alone is not
    // a valid configuration.
    if (bending && !stretching)
    {
        _parser.print_help();
        _parser.die("-bending/--bending requires -stretching/--stretching too (bending's conversion needs "
                    "stretching's real bond lengths)");
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
        if (!pathBondedInteraction.empty())
        {
            logging::info("    bonded interaction (stretching/bending retune, dihedral ghost springs): %s",
                          pathBondedInteraction.c_str());
            logging::info("        stretching: %s", stretching ? "enabled" : "disabled");
            logging::info("        bending: %s", bending ? "enabled" : "disabled");
            logging::info("        dihedral backbone: %s", dihedralBackbone ? "enabled" : "disabled");
            logging::info("        dihedral sidechain: %s", dihedralSidechain ? "enabled" : "disabled");
            logging::info("        dihedral planarity: %s", dihedralPlanarity ? "enabled" : "disabled");
        }
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
