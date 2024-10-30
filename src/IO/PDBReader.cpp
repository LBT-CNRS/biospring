
#include "IO/PDBReader.h"
#include "logging.h"
#include "utils/string.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>

using namespace biospring;

//
// Returns true if two particles are the same.
//
// Compares particle name, residue name and residue identifier.
//
static bool sameParticles(const topology::Particle & left, const topology::Particle & right)
{
    return (left.properties().name() == right.properties().name() &&
            left.properties().residue_name() == right.properties().residue_name() &&
            left.properties().residue_id() == right.properties().residue_id());
}

//
// Parses an atom line, creates a particle with according data and adds it to spring network.
//
topology::Particle PDBReader::parseAtomLine(const std::string & line)
{
    if (line.size() < 78)
    {
        logging::error("PDB format requires a line that is at least 78 characters long");
        logging::die("line too short: '%s'", line.c_str());
    }
    int id = std::stoi(line.substr(6, 5));
    std::string name = biospring::utils::string::trim(line.substr(12, 4));
    std::string resname = biospring::utils::string::trim(line.substr(17, 3));
    std::string chain = biospring::utils::string::trim(line.substr(21, 1));
    int resid = std::stoi(line.substr(22, 4));
    float x = std::stof(line.substr(30, 8));
    float y = std::stof(line.substr(38, 8));
    float z = std::stof(line.substr(46, 8));
    float occupancy = std::stof(line.substr(54, 6));
    float tempfactor = std::stof(line.substr(60, 6));
    std::string elementname = biospring::utils::string::trim(line.substr(76, 2));
    float charge = 0.0;
    if (line.size() > 78 and biospring::utils::string::trim(line.substr(78, 2)) != "")
        charge = std::stof(line.substr(78, 2));

    topology::ParticleProperties properties = topology::ParticleProperties::build()
                                                  .name(name)
                                                  .atom_id(id)
                                                  .residue_name(resname)
                                                  .residue_id(resid)
                                                  .chain_name(chain)
                                                  .element_name(elementname)
                                                  .position(Vector3f(x, y, z))
                                                  .occupancy(occupancy)
                                                  .temperature_factor(tempfactor)
                                                  .charge(charge);

    return topology::Particle(properties);
}

std::vector<std::pair<size_t, size_t>> PDBReader::parseConectLine(const std::string & line)
{
    std::vector<std::pair<size_t, size_t>> pairs_indexes;

    size_t serial = std::stoi(line.substr(6, 5));

    // Helper function to try parsing and add to vector if successful
    auto try_add_pair = [&](int start_index) {
        try {
            size_t connected_serial = std::stoi(line.substr(start_index, 5));
            pairs_indexes.emplace_back(serial, connected_serial);
        } catch (const std::exception&) {
            // Catches std::invalid_argument and std::out_of_range
            // Do nothing if parsing fails or substring is out of range
        }
    };

    // Add each connected serial if valid
    try_add_pair(11); // serial1
    try_add_pair(16); // serial2
    try_add_pair(21); // serial3
    try_add_pair(26); // serial4

    return pairs_indexes;
}

void PDBReader::read()
{
    // Opens file for reading, dies if error occurs.
    safeOpen();

    // Loop over lines.
    std::string buffer;
    std::getline(_instream, buffer);
    do
    {
        if (isAtomLine(buffer))
        {
            topology::Particle p = parseAtomLine(buffer);

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
        else if (isConectLine(buffer))
        {
            // https://pdb2pqr.readthedocs.io/en/v3.5.0/_modules/pdb2pqr/pdb.html#CONECT.__init__
            auto pairs_indexes = parseConectLine(buffer);
            for (auto & pair_indexes : pairs_indexes)
            {
                _topology.add_spring(_extidtoindex[pair_indexes.first], _extidtoindex[pair_indexes.second]);
            }
        }
    } while (std::getline(_instream, buffer));
}

//
// Updates a spring network particle positions using data read from a PDB file.
//
void PDBReader::updatePositions()
{
    // Dies if topology contains no particle.
    if (_topology.number_of_particles() == 0)
        logging::die("No valid Topology to update.");

    // Opens files (dies if error occurs).
    safeOpen();

    // Loop over lines
    std::string buffer;
    size_t currentParticleIndex = 0;
    do
    {
        if (PDBReader::isAtomLine(buffer))
        {
            topology::Particle p = parseAtomLine(buffer);
            topology::Particle & fromSPN = _topology.get_particle(currentParticleIndex);
            if (sameParticles(p, fromSPN))
            {
                fromSPN.set_position(p.position());
                currentParticleIndex++;
            }
            else
            {
                logging::warning("Update Positions : %d %s %s %d not found in SPN.", p.properties().atom_id(),
                                 p.properties().name().c_str(), p.properties().residue_name().c_str(),
                                 p.properties().residue_id());
            }
        }
    } while (std::getline(_instream, buffer));

    if (currentParticleIndex != _topology.number_of_particles())
    {
        logging::warning("Positions update file and SpringNetwork have different number of particles (%d vs %d)",
                         currentParticleIndex, _topology.number_of_particles());
    }
}

//
// Returns true if filter is found in the internal filter collection.
//
bool PDBReader::isInFilter(const std::string & filter) const
{
    if (_atomfilter.size() == 0)
        return true;
    return std::count(_atomfilter.begin(), _atomfilter.end(), filter) > 0;
}

int PDBReader::getIdFromExtid(size_t extid) const
{
    if (_extidtoindex.find(extid) == _extidtoindex.end())
    {
        return UNDEFINEDINDEX;
    }
    return _extidtoindex.at(extid);
}
