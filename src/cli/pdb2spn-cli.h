
#include "argparse.hpp"
#include "topology.hpp"

#include <string>
#include <vector>

namespace biospring
{

namespace pdb2spn
{

// Handles command-line arguments storage for pdb2spn.
class CommandLineArguments : argparse::CommandLineArgumentsBase
{
  public:
    // I/O paths.
    std::string pathTopology;
    std::string pathForceField;
    std::string pathGroup;
    std::vector<std::string> pathOutputList;

    // User options.
    float cutoff;
    float stiffness;
    float charge;
    bool isStatic;
    bool ignoreDuplicates;
    bool ignoreMissing;

    // Constructor (inherited from argparse::CommandLineArgumentsBase).
    CommandLineArguments(const std::string & name, const argparse::description_t & description,
                         const std::string & version = "");

    // Implements abstract methods from parent class.
    void printArgumentValues() const;
    void parseCommandLine(int argc, const char * const argv[]);

    bool useUserCutoff() const { return _parser.get_option("--cutoff").is_set(); }
    bool useUserStiffness() const { return _parser.get_option("--stiffness").is_set(); }
    bool useUserCharge() const { return _parser.get_option("--charge").is_set(); }
};

void reduceToCoarseGrain(topology::Topology & top, const CommandLineArguments & args);

int main(int argc, char ** argv);

} // namespace pdb2spn

} // namespace biospring
