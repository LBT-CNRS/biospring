
#include "IO/PQRReader.h"
#include "logging.h"
#include "utils/string.hpp"

using namespace biospring;

void PQRReader::read()
{
    // Opens file for reading, dies if error occurs.
    safeOpen();

    // Loop over lines.
    std::string buffer;
    std::getline(_instream, buffer);
    do
    {
        if (PDBReader::isAtomLine(buffer))
        {
            topology::Particle p = parseAtomLine(buffer);
            // p.setId(serial++);
            p.properties().set_occupancy(0.0);
            p.properties().set_temperature_factor(0.0);

            if (isInFilter(p.properties().name()))
            {
                int atom_id = p.properties().atom_id();
                if (_extidtoindex.find(atom_id) != _extidtoindex.end())
                {
                    logging::error("%s", buffer.c_str());
                    logging::die("particle with id '%d' already exists", atom_id);
                }
                _extidtoindex[atom_id] = _topology.number_of_particles();
                _topology.add_particle(p);
            }
        }
    } while (std::getline(_instream, buffer));
}

//
// Parses an atom line, creates a particle with according data and adds it to spring network.
//
topology::Particle PQRReader::parseAtomLine(const std::string & line)
{
    // The nominal fixed-width columns run through 69 (radius, cols 62-69),
    // but many real-world PQR files (including some of BioSpring's own
    // historical examples) write charge/radius without padding to the full
    // field width, e.g. "1.8" instead of "1.800  ". Only require the line
    // to reach into the radius field at all (col 62); std::stof on
    // whatever is actually there still parses a short-but-complete value,
    // and throws its own clear error on a genuinely missing/malformed one.
    if (line.size() < 63)
    {
        logging::error("PQR format requires a line that is at least 63 characters long");
        logging::die("line too short: '%s'", line.c_str());
    }
    int id = std::stoi(line.substr(6, 5));
    std::string name = biospring::utils::string::trim(line.substr(12, 4));
    std::string resname = biospring::utils::string::trim(line.substr(17, 3));
    std::string chain = biospring::utils::string::trim(line.substr(21, 1));
    int resid = std::stoi(line.substr(22, 4));
    float x = std::stof(line.substr(30, 7));
    float y = std::stof(line.substr(38, 7));
    float z = std::stof(line.substr(46, 7));
    float charge = std::stof(line.substr(54, 8));
    float radius = std::stof(line.substr(62, 7));

    topology::ParticleProperties properties = topology::ParticleProperties::build()
                                                  .name(name)
                                                  .atom_id(id)
                                                  .residue_name(resname)
                                                  .residue_id(resid)
                                                  .chain_name(chain)
                                                  .position(Vector3f(x, y, z))
                                                  .charge(charge)
                                                  .radius(radius);

    return topology::Particle(properties);
}