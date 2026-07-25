#include "pdb2spn-cli.h"

#include <cstdlib>
#include <initializer_list>
#include <iostream>
#include <string>
#include <vector>

namespace
{

enum class LegacyTool
{
    Pdb2Cdl,
    Pdb2Pqr,
    Pdb2PdbConect,
    Unknown,
};

std::string basename(std::string path)
{
    const std::string separators = "/\\";
    const auto pos = path.find_last_of(separators);
    if (pos != std::string::npos)
        path = path.substr(pos + 1);

#ifdef _WIN32
    const std::string suffix = ".exe";
    if (path.size() >= suffix.size() && path.substr(path.size() - suffix.size()) == suffix)
        path.resize(path.size() - suffix.size());
#endif

    return path;
}

LegacyTool detect_tool(const char * argv0)
{
    const std::string program = basename(argv0 ? argv0 : "");
    if (program == "pdb2cdl")
        return LegacyTool::Pdb2Cdl;
    if (program == "pdb2pqr")
        return LegacyTool::Pdb2Pqr;
    if (program == "pdb2pdbconect" || program == "pdb2pdbconnect")
        return LegacyTool::Pdb2PdbConect;
    return LegacyTool::Unknown;
}

bool is_option(const std::string & value)
{
    return !value.empty() && value[0] == '-';
}

bool is_one_of(const std::string & value, std::initializer_list<const char *> names)
{
    for (const char * name : names)
    {
        if (value == name)
            return true;
    }
    return false;
}

void print_help(LegacyTool tool, const std::string & program)
{
    std::cerr << program << " is a compatibility wrapper around pdb2spn.\n\n";
    switch (tool)
    {
        case LegacyTool::Pdb2Cdl:
            std::cerr << "Usage: " << program << " -pdb INPUT.pdb -cdl OUTPUT.cdl [-grp GROUP.grp -ff FORCEFIELD.ff] [-springcutoff CUTOFF]\n\n"
                      << "Creates a text NetCDF/CDL spring-network file.\n";
            break;
        case LegacyTool::Pdb2Pqr:
            std::cerr << "Usage: " << program << " -pdb INPUT.pdb -pqr OUTPUT.pqr [-grp GROUP.grp -ff FORCEFIELD.ff]\n\n"
                      << "Creates a PQR file.\n";
            break;
        case LegacyTool::Pdb2PdbConect:
            std::cerr << "Usage: " << program << " -pdb INPUT.pdb -pdbconect OUTPUT.pdb [-grp GROUP.grp -ff FORCEFIELD.ff] [-springcutoff CUTOFF]\n\n"
                      << "Creates a PDB file with CONECT records representing the spring network.\n"
                      << "If -springcutoff/--springcutoff is omitted, the historical default cutoff is 9.0 Angstrom.\n";
            break;
        default:
            std::cerr << "Usage: " << program << " [legacy pdb2spn options]\n";
            break;
    }

    std::cerr << "\nAccepted legacy options:\n"
              << "  -pdb, --pdb INPUT          input PDB file\n"
              << "  -grp, --grp GROUP          coarse-grain grouping file\n"
              << "  -ff, --ff FORCEFIELD       force-field file\n"
              << "  -springcutoff, --springcutoff CUTOFF\n"
              << "                              cutoff in Angstroms for spring creation\n"
              << "  -h, -help, --help          show this help message\n\n"
              << "The wrapper translates these options to pdb2spn and reuses the current implementation.\n";
}

[[noreturn]] void die_with_help(LegacyTool tool, const std::string & program, const std::string & message)
{
    print_help(tool, program);
    std::cerr << "\nError: " << message << '\n';
    std::exit(EXIT_FAILURE);
}

std::string require_value(LegacyTool tool, const std::string & program, int argc, char ** argv, int & index)
{
    if (index + 1 >= argc)
        die_with_help(tool, program, std::string("option '") + argv[index] + "' requires an argument");

    const std::string value = argv[++index];
    if (is_option(value))
        die_with_help(tool, program, std::string("option '") + argv[index - 1] + "' requires an argument");

    return value;
}

std::vector<std::string> translate_argv(LegacyTool tool, int argc, char ** argv)
{
    const std::string program = basename(argv[0] ? argv[0] : "legacy-wrapper");
    std::vector<std::string> translated;
    translated.emplace_back(program);

    bool has_input = false;
    bool has_output = false;
    bool has_cutoff = false;

    for (int i = 1; i < argc; ++i)
    {
        const std::string arg = argv[i];

        if (is_one_of(arg, {"-h", "-help", "--help"}))
        {
            print_help(tool, program);
            std::exit(EXIT_SUCCESS);
        }
        if (arg == "--version")
        {
            translated.push_back(arg);
            continue;
        }

        if (is_one_of(arg, {"-pdb", "--pdb", "-s", "--topology", "-topology"}))
        {
            translated.emplace_back("--topology");
            translated.push_back(require_value(tool, program, argc, argv, i));
            has_input = true;
        }
        else if ((tool == LegacyTool::Pdb2Cdl && is_one_of(arg, {"-cdl", "--cdl"})) ||
                 (tool == LegacyTool::Pdb2Pqr && is_one_of(arg, {"-pqr", "--pqr"})) ||
                 (tool == LegacyTool::Pdb2PdbConect && is_one_of(arg, {"-pdbconect", "--pdbconect"})) ||
                 is_one_of(arg, {"-o", "--output", "-output"}))
        {
            translated.emplace_back("--output");
            translated.push_back(require_value(tool, program, argc, argv, i));
            has_output = true;
        }
        else if (is_one_of(arg, {"-grp", "--grp"}))
        {
            translated.emplace_back("--grp");
            translated.push_back(require_value(tool, program, argc, argv, i));
        }
        else if (is_one_of(arg, {"-ff", "--ff"}))
        {
            translated.emplace_back("--ff");
            translated.push_back(require_value(tool, program, argc, argv, i));
        }
        else if (is_one_of(arg, {"-springcutoff", "--springcutoff", "-cutoff", "--cutoff"}))
        {
            translated.emplace_back("--cutoff");
            translated.push_back(require_value(tool, program, argc, argv, i));
            has_cutoff = true;
        }
        else if (is_one_of(arg, {"-stiffness", "--stiffness", "-charge", "--charge"}))
        {
            translated.push_back(arg.rfind("--", 0) == 0 ? arg : "--" + arg.substr(1));
            translated.push_back(require_value(tool, program, argc, argv, i));
        }
        else if (is_one_of(arg, {"-static", "--static", "-ignore-duplicate", "--ignore-duplicate",
                                "-ignore-missing", "--ignore-missing"}))
        {
            translated.push_back(arg.rfind("--", 0) == 0 ? arg : "--" + arg.substr(1));
        }
        else if (is_one_of(arg, {"-enableprobe", "--enableprobe", "-nospringbetweenchain", "--nospringbetweenchain"}))
        {
            die_with_help(tool, program, "legacy option '" + arg + "' is not supported by the current pdb2spn engine");
        }
        else
        {
            die_with_help(tool, program, "unknown option '" + arg + "'");
        }
    }

    if (!has_input)
        die_with_help(tool, program, "missing input PDB file; provide -pdb INPUT.pdb");
    if (!has_output)
    {
        switch (tool)
        {
            case LegacyTool::Pdb2Cdl:
                die_with_help(tool, program, "missing output CDL file; provide -cdl OUTPUT.cdl");
            case LegacyTool::Pdb2Pqr:
                die_with_help(tool, program, "missing output PQR file; provide -pqr OUTPUT.pqr");
            case LegacyTool::Pdb2PdbConect:
                die_with_help(tool, program, "missing output PDB file; provide -pdbconect OUTPUT.pdb");
            default:
                die_with_help(tool, program, "missing output file");
        }
    }

    if (tool == LegacyTool::Pdb2PdbConect)
    {
        if (!has_cutoff)
        {
            translated.emplace_back("--cutoff");
            translated.emplace_back("9.0");
        }
        translated.emplace_back("--pdbconect");
    }

    return translated;
}

} // namespace

int main(int argc, char ** argv)
{
    const LegacyTool tool = detect_tool(argv[0]);
    if (tool == LegacyTool::Unknown)
    {
        std::cerr << "Unknown compatibility wrapper name: " << basename(argv[0] ? argv[0] : "") << '\n';
        return EXIT_FAILURE;
    }

    std::vector<std::string> translated = translate_argv(tool, argc, argv);
    std::vector<char *> translated_argv;
    translated_argv.reserve(translated.size());
    for (std::string & value : translated)
        translated_argv.push_back(value.data());

    return biospring::pdb2spn::main(static_cast<int>(translated_argv.size()), translated_argv.data());
}
