
#include "argparse.hpp"

#include <string>
#include <vector>

namespace biospring
{

namespace mergespn
{

//
// Handles command-line parsing as well as argument storage for program mergespn.
//
class CommandLineArguments : public argparse::CommandLineArgumentsBase
{
  public:
    // I/O paths.
    std::vector<std::string> pathTopologyList;
    std::vector<std::string> pathOutputList;

    // User options.
    float cutoff;

    // Constructor (inherited from argparse::CommandLineArgumentsBase).
    CommandLineArguments(const std::string & name, const argparse::description_t & description,
                         const std::string & version = "");

    // Implements abstract methods from parent class.
    void printArgumentValues() const;
    void parseCommandLine(int argc, const char * const argv[]);

    bool useUserCutoff() const { return _parser.get_option("--cutoff").is_set(); }
};

int main(int argc, char ** argv);

} // namespace mergespn

} // namespace biospring
