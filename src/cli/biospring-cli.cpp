

#ifdef MDDRIVER_SUPPORT
#include "interactor/mddriver/InteractorMDDriver.h"
#endif

#ifdef FREESASA_SUPPORT
#include "interactor/freesasa/InteractorFreeSASA.h"
#endif

#ifdef OPENMP_SUPPORT
#include "omp.h"
#endif

#ifdef OPENCL_SUPPORT
#include "SpringNetworkOpenCL.h"
#endif

#include "SpringNetwork.h"

#include "IO/NetCDFReader.h"
#include "IO/PDBReader.h"
#include "logging.h"

#include <iostream>
#include <string>
#include <vector>

#include "biospring-cli.h"

#include "configuration/Configuration.hpp"
#include "configuration/SafeConfigurationReader.hpp"

#include "measure.hpp"

namespace biospring
{
namespace biospringcli
{

const argparse::description_t PROGRAM_DESCRIPTION = {
    "biospring spring network engine.",
};

int main(int argc, char ** argv)
{
    // Command-line parsing.
    CommandLineArguments args(std::string(argv[0]), PROGRAM_DESCRIPTION, PROGRAM_VERSION);
    args.parseCommandLine(argc, argv);
    args.printArgumentValues();

    // SpringNetwork initialization: type of SpringNetwork depends on whether
    // OpenCL support is active or not.
    biospring::spn::SpringNetwork * spn = nullptr;

#ifdef OPENCL_SUPPORT
    if (args.openclenabled)
        spn = new SpringNetworkOpenCL();
    else
        spn = new biospring::spn::SpringNetwork();
#else
    spn = new biospring::spn::SpringNetwork();
#endif

        // Interactor initialization.
#if defined(MDDRIVER_SUPPORT)
    biospring::interactor::InteractorMDDriver iMDDriver;
    iMDDriver.setDebug(args.mddriverParam.debug);
    iMDDriver.setPort(args.mddriverParam.port);
    iMDDriver.setLog(args.mddriverParam.logpath.c_str());
    iMDDriver.setWait(args.mddriverParam.wait);
    iMDDriver.setForceScale(args.mddriverParam.forcescale);
    iMDDriver.setSpringNetwork(spn);
    spn->addInteractor(&iMDDriver);
#endif

        // FreeSASA interactor initialization.
#if defined(FREESASA_SUPPORT)
    biospring::interactor::InteractorFreeSASA iFreeSASA;
    iFreeSASA.setDynamic(args.freesasaParam.sasa_dynamic);
    iFreeSASA.setAlg(args.freesasaParam.sasa_alg);
    iFreeSASA.set_lr_n(args.freesasaParam.sasa_lr_n);
    iFreeSASA.set_sr_n(args.freesasaParam.sasa_sr_n);
    iFreeSASA.setProbeRad(args.freesasaParam.sasa_probe_radius);
    iFreeSASA.setRadiiClassifier(args.freesasaParam.sasa_classifier);
    iFreeSASA.setNthreads(args.freesasaParam.sasa_n_threads);
    iFreeSASA.setSpringNetwork(spn);
    spn->addInteractor(&iFreeSASA);
#endif

    // Reads configuration file.
    logging::status("Reading MSP file %s.", args.pathConfig.c_str());
    auto configReader =
        biospring::configuration::SafeConfigurationReader(biospring::configuration::defaultConfiguration());
    configReader.setFileName(args.pathConfig);
    configReader.read();

    logging::status("Using configuration parameters:");
    auto config = configReader.getConfiguration();
    config.print();

    // Reads topology file.
    logging::status("Reading Nc file %s.", args.pathTopology.c_str());
    NetCDFReader netcdfreader;
    netcdfreader.setFileName(args.pathTopology);
    netcdfreader.read();

    netcdfreader.getTopology().to_spring_network(*spn);

    spn->setup(config);

    // if(startposok)
    // {
    // 	PDBReader * startposreader = new PDBReader();
    // 	startposreader->setSpringNetwork(spn);
    // 	logging::status("Reading Starting Positions in %s.", startposname);
    // 	startposreader->setFileName(startposname);
    // 	startposreader->updatePositions();
    // 	delete startposreader;
    // }

    // Launching the simulation.
    logging::status("Running simulation...");
#ifdef OPENMP_SUPPORT
    if (omp_get_thread_num() == 0)
        logging::status("OpenMP running thread %d", omp_get_max_threads());
#endif
    spn->run();

    delete spn;
    return 0;
}

CommandLineArguments::CommandLineArguments(const std::string & name, const argparse::description_t & description,
                                           const std::string & version)
    : CommandLineArgumentsBase(name, description, version), pathTopology(""), pathConfig(""), mddriverParam()
{
    argparse::Argument topology = argparse::Argument()
                                      .name_short("-s")
                                      .name_long("--nc")
                                      .description("input topology (nc format).")
                                      .metavar("NC")
                                      .argument_type(argparse::ArgumentType::PATH_INPUT)
                                      .required(true);

    argparse::Argument config = argparse::Argument()
                                    .name_short("-c")
                                    .name_long("--msp")
                                    .description("input configuration (msp format).")
                                    .metavar("MSP")
                                    .argument_type(argparse::ArgumentType::PATH_INPUT)
                                    .required(true);

    _parser.add_argument(topology);
    _parser.add_argument(config);

    // == MDDriver-specific options ==

#ifdef MDDRIVER_SUPPORT
    argparse::Argument wait =
        argparse::StoreTrueArgument("", "--wait", "Waiting for a connection before starting the simulation.");

    argparse::Argument port = argparse::Argument()
                                  .name_long("--port")
                                  .description("Listening and outcoming port for MDDriver.")
                                  .argument_type(argparse::ArgumentType::INTEGER)
                                  .default_value("3000");

    argparse::Argument debug = argparse::Argument()
                                   .name_long("--debug")
                                   .description("Debug level of the MDDriver simulation (0, 1 or 2).")
                                   .argument_type(argparse::ArgumentType::INTEGER)
                                   .default_value("0");

    argparse::Argument log = argparse::Argument()
                                 .name_long("--log")
                                 .description("Path to MDDriver log file.")
                                 .metavar("LOG")
                                 .default_value("stdout")
                                 .argument_type(argparse::ArgumentType::PATH_OUTPUT);

    _parser.add_argument(port);
    _parser.add_argument(wait);
    _parser.add_argument(debug);
    _parser.add_argument(log);
#endif // MDDRIVER_SUPPORT

#ifdef FREESASA_SUPPORT

    argparse::Argument sasa_dynamic = argparse::Argument()
        .name_long("--sasa-dynamic")
        .argument_type(argparse::ArgumentType::BOOLEAN)
        .default_value("false")
        .description("Update SASA according to sleep duration parameter in sasa calculation thread.");

    argparse::Argument sasa_sleep = argparse::Argument()
        .name_long("--sasa-sleep")
        .argument_type(argparse::ArgumentType::INTEGER)
        .default_value("1000")
        .description("Set sleep duration (in microseconds) for the SASA calculation thread.");

    argparse::Argument sasa_alg = argparse::Argument()
        .name_long("--sasa-alg")
        .argument_type(argparse::ArgumentType::STRING)
        .default_value("lr")
        .description("Choose between Lee-Richards (lr) or Shrake-Rupley (sr) algorithm.");

    argparse::Argument sasa_lr_n = argparse::Argument()
        .name_long("--sasa-lr-n")
        .argument_type(argparse::ArgumentType::INTEGER)
        .default_value("20")
        .description("Number of slices per atom in L&R.");

    argparse::Argument sasa_sr_n = argparse::Argument()
        .name_long("--sasa-sr-n")
        .argument_type(argparse::ArgumentType::INTEGER)
        .default_value("100")
        .description("Number of test points in S&R.");

    argparse::Argument sasa_probe_radius = argparse::Argument()
        .name_long("--sasa-probe-radius")
        .argument_type(argparse::ArgumentType::REAL)
        .default_value("1.4")
        .description("Probe radius (in Ångström).");

    argparse::Argument sasa_classifier = argparse::Argument()
        .name_long("--sasa-classifier")
        .argument_type(argparse::ArgumentType::STRING)
        .default_value("default")
        .description("Struct used to determine radii of atoms (protor, naccess, oons, default, biospring).");

    argparse::Argument sasa_n_threads = argparse::Argument()
        .name_long("--sasa-n-threads")
        .argument_type(argparse::ArgumentType::INTEGER)
        .default_value("1")
        .description("Number of threads to use by FreeSASA.");

    _parser.add_argument(sasa_dynamic);
    _parser.add_argument(sasa_sleep);
    _parser.add_argument(sasa_alg);
    _parser.add_argument(sasa_lr_n);
    _parser.add_argument(sasa_sr_n);
    _parser.add_argument(sasa_probe_radius);
    _parser.add_argument(sasa_classifier);
    _parser.add_argument(sasa_n_threads);

#endif // FREESASA_SUPPORT

    // == OpenCL-specific options ==

#ifdef OPENCL_SUPPORT
    argparse::Argument opencl = argparse::StoreTrueArgument("", "--opencl", "Uses OpenCL.");
    _parser.add_argument(opencl);
#endif // OPENCL_SUPPORT
}

void CommandLineArguments::parseCommandLine(int argc, const char * const argv[])
{
    _parser.parse_arguments(argc, argv);
    pathTopology = _parser.get_option_value<std::string>("--nc");
    pathConfig = _parser.get_option_value<std::string>("--msp");

#ifdef MDDRIVER_SUPPORT
    mddriverParam.wait = _parser.get_option("--wait").is_set();
    mddriverParam.port = _parser.get_option_value<unsigned>("--port");
    mddriverParam.debug = _parser.get_option_value<unsigned>("--debug");
    mddriverParam.logpath = _parser.get_option_value<std::string>("--log");

    if (mddriverParam.debug > 2)
    {
        logging::die("Debug level should be 0, 1 or 2.");
    }
#endif // MDDRIVER_SUPPORT

#ifdef FREESASA_SUPPORT
    freesasaParam.sasa_dynamic = _parser.get_option("--sasa-dynamic").is_set();
    freesasaParam.sasa_sleep = _parser.get_option_value<unsigned>("--sasa-sleep");
    freesasaParam.sasa_alg = _parser.get_option_value<std::string>("--sasa-alg");
    freesasaParam.sasa_lr_n = _parser.get_option_value<unsigned>("--sasa-lr-n");
    freesasaParam.sasa_sr_n = _parser.get_option_value<unsigned>("--sasa-sr-n");
    freesasaParam.sasa_probe_radius = _parser.get_option_value<double>("--sasa-probe-radius");
    freesasaParam.sasa_classifier = _parser.get_option_value<std::string>("--sasa-classifier");
    freesasaParam.sasa_n_threads = _parser.get_option_value<unsigned>("--sasa-n-threads");

    if (freesasaParam.sasa_alg !="lr" && freesasaParam.sasa_alg !="sr")
        logging::die("FreeSASA algorithm option should be lr or sr.");
    if (freesasaParam.sasa_probe_radius < 0.0)
        logging::die("Probe radius should not be a negative number.");
#endif // FREESASA_SUPPORT

#ifdef OPENCL_SUPPORT
    openclenabled = _parser.get_option("--opencl").is_set();
#endif // OPENCL_SUPPORT
}

void CommandLineArguments::printArgumentValues() const
{
    logging::status("Running biospring with arguments:");
    logging::info("    topology: %s", pathTopology.c_str());
    logging::info("    configuration: %s", pathConfig.c_str());

#ifdef MDDRIVER_SUPPORT
    logging::info("    MDDriver parameters:");
    if (mddriverParam.wait)
        logging::info("      wait: ON");
    else
        logging::info("      wait: OFF");
    logging::info("      port: Try to open %d. Check bellow.", mddriverParam.port);
    logging::info("      debug: %d", mddriverParam.debug);
    if (mddriverParam.logpath.empty())
        logging::info("      log file: none (logging to stdout)");
    else
        logging::info("      log file: %s", mddriverParam.logpath.c_str());
#endif // MDDRIVER_SUPPORT

#ifdef FREESASA_SUPPORT
    logging::info("    FreeSASA parameters:");
    if (freesasaParam.sasa_dynamic)
        logging::info("      dynamic: ON");
    else
        logging::info("      dynamic: OFF");
    logging::info("      algorithm: %s", freesasaParam.sasa_alg.c_str());
    if (freesasaParam.sasa_alg =="lr")
        logging::info("      n-slices: %d", freesasaParam.sasa_lr_n);
    else if (freesasaParam.sasa_alg =="sr")
        logging::info("      n-points: %d", freesasaParam.sasa_sr_n);
    logging::info("      probe_radius (Å): %f", freesasaParam.sasa_probe_radius);
    logging::info("      classifier: %s", freesasaParam.sasa_classifier.c_str());
    logging::info("      n_threads: %d", freesasaParam.sasa_n_threads);
#endif // FREESASA_SUPPORT

#ifdef OPENCL_SUPPORT
    if (openclenabled)
        logging::info("    use OpenCL: ON");
    else
        logging::info("    use OpenCL: OFF");
#endif // OPENCL_SUPPORT
}

} // namespace biospringcli
} // namespace biospring