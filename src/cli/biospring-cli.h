
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
    // Whether the user ASKED for interactive MD, as opposed to the defaults
    // below simply existing. Set only when --port or --wait is given on the
    // command line.
    //
    // It gates registration of the interactor, and that gate is the point: the
    // interactor used to be registered unconditionally, so every batch run
    // copied all its positions into the transmission buffer once per step for
    // a viewer that was never going to connect. Measured on a 389-particle
    // system, that was about 6 % of the whole step -- and syncSystemStateData's
    // own comment calls the transport "more expensive than the physics it
    // carries". A viewer can now only connect when interactive MD was asked
    // for, which is also the safer behaviour: an unattended run does not open
    // a port.
    bool enabled = false;
    bool wait = true;
    unsigned port = 8888;
    unsigned debug = 0;
    std::string logpath = "";
    float forcescale = 1.0;
    // Send one frame every N steps. 1 keeps the previous behaviour; a client
    // asking for a different rate over IMD_TRATE overrides it at run time.
    unsigned transmissionrate = 1;
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
    // Off unless asked for, or unless IMPALA is on -- the only term that reads
    // the surface FreeSASA computes. See biospring-cli.cpp for why.
    bool sasa_enable = false;
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


