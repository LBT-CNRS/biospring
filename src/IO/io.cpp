
#include "IO/io.h"
#include "IO/NetCDFReader.h"
#include "IO/NetCDFWriter.h"
#include "IO/PDBReader.h"
#include "IO/PDBWriter.h"
#include "IO/PQRReader.h"
#include "IO/PQRWriter.h"
#include "IO/ReaderBase.h"
#include "logging.h"
#include "utils.hpp"

#include <memory>
#include <string>
#include <unordered_set>
#include <vector>

namespace biospring
{
namespace io
{

static void _writeTopology(const std::vector<std::string> & outputfiles, const spn::SpringNetwork * const spn,
                           bool writePdbConect);
static void _writeTopology(const std::string & path, const spn::SpringNetwork * const spn, bool writePdbConect);

//
// Reads topology file.
//
// Format is guessed based on file name extension.
//
topology::Topology readTopology(const std::string & path)
{
    logging::status("Reading input structure %s.", path.c_str());

    std::string extension = biospring::utils::path::getExtension(path);
    std::unique_ptr<TopologyReaderBase> reader;

    if (extension == "pdb")
        reader = std::make_unique<PDBReader>();
    else if (extension == "pqr")
        reader = std::make_unique<PQRReader>();
    else if (extension == "nc")
        reader = std::make_unique<NetCDFReader>();
    else
        logging::die("Unrecognized topology format '%s'", extension.c_str());

    reader->setFileName(path);
    reader->read();

    logging::info("Read %d particles from topology.", reader->getTopology().number_of_particles());
    logging::info("Read %d springs from topology.", reader->getTopology().number_of_springs());

    // Dies if file contained no particle.
    if (reader->getTopology().number_of_particles() == 0)
        logging::die("%s: no particle read from file", path.c_str());

    return reader->getTopology();
}

void writeTopology(const std::string & path, const topology::Topology & topology, bool writePdbConect)
{
    spn::SpringNetwork spn;
    topology.to_spring_network(spn);
    _writeTopology(path, &spn, writePdbConect);
}

void writeTopology(const std::vector<std::string> & outputfiles, const topology::Topology & topology, bool writePdbConect)
{
    spn::SpringNetwork spn;
    topology.to_spring_network(spn);
    _writeTopology(outputfiles, &spn, writePdbConect);
}

//
// Write topology to output file.
// Format is guessed according to the file extension.
//
static void _writeTopology(const std::string & path, const spn::SpringNetwork * const spn, bool writePdbConect)
{
    class TopologyWriter
    {
      public:
        std::string path;
        const spn::SpringNetwork * spn;
        std::string format;

        const std::unordered_set<std::string> allowedFormats{"pdb", "pqr", "cdl", "nc"};

        bool writePdbConect;

        TopologyWriter(const std::string & p, const spn::SpringNetwork * const s, bool conect)
            : path(p), spn(s), writePdbConect(conect)
        {
            format = biospring::utils::path::getExtension(path);

            if (biospring::utils::path::hasNoExtension(path))
            {
                logging::warning("%s: No extension. Writing in NetCDF binary format.", path.c_str());
                format = "nc";
            }
            else if (not isValidFormat())
            {
                logging::warning("%s: Unrecognized format \"%s\". Writing in NetCDF binary format.", path.c_str(),
                                 format.c_str());
                format = "nc";
            }
        }

        bool isValidFormat() const { return allowedFormats.find(format) != allowedFormats.end(); }

        void write() const
        {
            if (format == "pdb")
                writePDB();
            else if (format == "pqr")
                writePQR();
            else if (format == "nc")
                writeNC();
            else if (format == "cdl")
                writeCDL();
        }

        void writePDB() const
        {
            logging::status("Writing spring network to PDB file %s.", path.c_str());
            PDBWriter writer(path, spn);
            writer.setIsConnect(writePdbConect);
            writer.write();
        }

        void writePQR() const
        {
            logging::status("Writing spring network to PQR file %s (springs are ignored).", path.c_str());
            PQRWriter(path, spn).write();
        }

        void writeNC() const
        {
            logging::status("Writing spring network to nc file %s.", path.c_str());
            NetCDFWriter(path, spn).writeBinary();
        }

        void writeCDL() const
        {
            logging::status("Writing spring network to ASCII cdl file %s.", path.c_str());
            NetCDFWriter(path, spn).write();
        }
    };

    TopologyWriter(path, spn, writePdbConect).write();
}

// Write multiple topology files.
// See also: writeTopology.
static void _writeTopology(const std::vector<std::string> & outputfiles, const spn::SpringNetwork * const spn,
                           bool writePdbConect)
{
    for (const std::string & path : outputfiles)
        _writeTopology(path, spn, writePdbConect);
}
} // namespace io
} // namespace biospring