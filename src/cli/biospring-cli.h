
#include "cli/argparse.hpp"
#include "version.h"

#include <string>
#include <set>

namespace biospring
{

namespace biospringcli
{

const std::string PROGRAM_VERSION = "biospring " + biospring::VERSION_STRING;

// Stores MDDriver parameters.
struct MDDriverParameters
{
    bool wait = true;
    unsigned port = 8888;
    unsigned debug = 0;
    std::string logpath = "";
    float forcescale = 1.0;
};

// Stores FreeSASA parameters.
struct FreeSASAParameters
{
    bool sasa_dynamic = false;
    unsigned sasa_sleep = 1000; // in microseconds
    std::string sasa_alg = "lr"; // Default Lee-Richards algorithm.
    unsigned sasa_lr_n = 20;
    unsigned sasa_sr_n = 100;
    double sasa_probe_radius = 1.4;
    std::string sasa_classifier = "default";
    unsigned sasa_n_threads = 2;
};
//
// Handles command-line parsing as well as argument storage for program biospring.
//
class CommandLineArguments : public argparse::CommandLineArgumentsBase
{
  public:
    // I/O paths.
    std::string pathTopology;
    std::string pathConfig;
#ifdef OPENCL_SUPPORT
    bool openclenabled = true;
#else
    bool openclenabled = false;
#endif
    MDDriverParameters mddriverParam;
    FreeSASAParameters freesasaParam;

    // Constructor (inherited from argparse::CommandLineArgumentsBase).
    CommandLineArguments(const std::string & name, const argparse::description_t & description, const std::string & version = "");

    // Implements abstract methods from parent class.
    void printArgumentValues() const;
    void parseCommandLine(int argc, const char * const argv[]);
};

int main(int argc, char ** argv);

} // namespace biospringcli
} // namespace biospring


