
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
    // Occupancy (cols 55-60), temperature factor (cols 61-66), element name
    // (cols 77-78) and charge (cols 79-80) are all optional trailing columns
    // of a PDB ATOM/HETATM record; a number of real-world files (including
    // some of BioSpring's own historical examples, e.g. bare xyz-only DNA
    // models with no occupancy/B-factor at all, sometimes not even padded
    // to the nominal 8-character-wide z column) omit some or all of them.
    // Only require the line to reach into the z column (position 47) so
    // std::stof below has at least one digit to parse, and treat each
    // trailing field as absent/defaulted rather than dying when the line
    // doesn't reach that far.
    if (line.size() < 47)
    {
        logging::error("PDB format requires a line that is at least 47 characters long");
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
    float occupancy = line.size() >= 60 && !biospring::utils::string::trim(line.substr(54, 6)).empty()
                           ? std::stof(line.substr(54, 6))
                           : 1.0f;
    float tempfactor = line.size() >= 66 && !biospring::utils::string::trim(line.substr(60, 6)).empty()
                            ? std::stof(line.substr(60, 6))
                            : 0.0f;
    std::string elementname = line.size() >= 78 ? biospring::utils::string::trim(line.substr(76, 2)) : "";
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

    // The connected serials below are each parsed defensively; the record's own
    // serial must be too, otherwise a malformed CONECT line throws
    // std::invalid_argument straight out of the reader.
    size_t serial = 0;
    try
    {
        serial = static_cast<size_t>(std::stoi(line.substr(6, 5)));
    }
    catch (const std::exception &)
    {
        logging::warning("skipping malformed CONECT record: %s", line.c_str());
        return pairs_indexes;
    }

    // Helper function to try parsing and add to vector if successful
    auto try_add_pair = [&](size_t start_index) {
        try {
            size_t connected_serial = static_cast<size_t>(std::stoi(line.substr(start_index, 5)));
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
                if (_extidtoindex.find(atom_id) == _extidtoindex.end())
                {
                    _extidtoindex[atom_id] = _topology.number_of_particles();
                }
                else
                {
                    // Some real-world PDB files (e.g. biological-assembly dumps that
                    // concatenate several symmetry copies without renumbering) reuse
                    // atom serials across otherwise-distinct atoms. CONECT resolution
                    // for such a serial is inherently ambiguous, so only the first
                    // occurrence stays addressable; every atom is still added below.
                    BIOSPRING_WARN_ONCE("PDB file reuses atom serial '%d' for multiple atoms "
                                        "(e.g. unrenumbered symmetry copies in a biological "
                                        "assembly): CONECT records referring to it will only "
                                        "resolve to its first occurrence",
                                        atom_id);
                }
                _topology.add_particle(p);
            }
        }
        else if (isSSBondLine(buffer) || isLinkLine(buffer))
        {
            ResidueAtomRef first, second;
            const bool ok = isSSBondLine(buffer) ? parseSSBondLine(buffer, first, second)
                                                 : parseLinkLine(buffer, first, second);
            if (ok)
                _residue_bonds.emplace_back(first, second);
            else
                logging::warning("skipping malformed %s record: %s", buffer.substr(0, 6).c_str(), buffer.c_str());
        }
        else if (isConectLine(buffer))
        {
            // https://pdb2pqr.readthedocs.io/en/v3.5.0/_modules/pdb2pqr/pdb.html#CONECT.__init__
            auto pairs_indexes = parseConectLine(buffer);
            for (auto & pair_indexes : pairs_indexes)
            {
                // _extidtoindex only holds atoms that passed the filter above.
                // Looking a missing serial up with operator[] would default-insert
                // 0 and silently bond the atom to particle 0 -- a wrong topology
                // that no error reports. Real PDB files routinely CONECT to
                // HETATMs the filter dropped, so skip the pair instead of dying.
                const auto first = _extidtoindex.find(pair_indexes.first);
                const auto second = _extidtoindex.find(pair_indexes.second);
                if (first == _extidtoindex.end() || second == _extidtoindex.end())
                {
                    BIOSPRING_WARN_ONCE("CONECT record refers to atom serials absent from the model "
                                        "(e.g. %zu-%zu): skipping those bonds",
                                        pair_indexes.first, pair_indexes.second);
                    continue;
                }
                // PDB CONECT records are listed symmetrically (a bond between A
                // and B appears once from A's own record and once from B's),
                // so the same spring is routinely seen twice.
                try
                {
                    _topology.add_spring(first->second, second->second);
                }
                catch (const topology::SpringAlreadyExistsException &)
                {
                }
            }
        }
    } while (std::getline(_instream, buffer));

    _resolveResidueBonds();
}

// A fixed-column field, empty when the line is too short to hold it.
static std::string column(const std::string & line, size_t begin, size_t length)
{
    if (line.size() < begin + length)
        return {};
    std::string field = line.substr(begin, length);
    const size_t first = field.find_first_not_of(' ');
    if (first == std::string::npos)
        return {};
    return field.substr(first, field.find_last_not_of(' ') - first + 1);
}

static bool readResidueId(const std::string & line, size_t begin, size_t length, int & out)
{
    const std::string field = column(line, begin, length);
    if (field.empty())
        return false;
    try
    {
        out = std::stoi(field);
    }
    catch (const std::exception &)
    {
        return false;
    }
    return true;
}

// SSBOND: chain 16, seqNum 18-21, then chain 30, seqNum 32-35 (1-based
// columns, PDB v3.3). The atom is always SG on both ends.
bool PDBReader::parseSSBondLine(const std::string & line, ResidueAtomRef & first, ResidueAtomRef & second)
{
    first.chain = column(line, 15, 1);
    second.chain = column(line, 29, 1);
    first.name = "SG";
    second.name = "SG";
    return readResidueId(line, 17, 4, first.residue_id) && readResidueId(line, 31, 4, second.residue_id);
}

// LINK: name 13-16, chain 22, seqNum 23-26, then name 43-46, chain 52,
// seqNum 53-56 (1-based columns, PDB v3.3).
bool PDBReader::parseLinkLine(const std::string & line, ResidueAtomRef & first, ResidueAtomRef & second)
{
    first.name = column(line, 12, 4);
    second.name = column(line, 42, 4);
    first.chain = column(line, 21, 1);
    second.chain = column(line, 51, 1);
    if (first.name.empty() || second.name.empty())
        return false;
    return readResidueId(line, 22, 4, first.residue_id) && readResidueId(line, 52, 4, second.residue_id);
}

void PDBReader::_resolveResidueBonds()
{
    if (_residue_bonds.empty())
        return;

    std::map<std::tuple<std::string, int, std::string>, size_t> by_atom;
    for (size_t i = 0; i < _topology.number_of_particles(); ++i)
    {
        const auto & p = _topology.get_particle(i).properties();
        by_atom.emplace(std::make_tuple(p.chain_name(), p.residue_id(), p.name()), i);
    }

    size_t added = 0, unresolved = 0, duplicate = 0;
    for (const auto & bond : _residue_bonds)
    {
        const auto first = by_atom.find(std::make_tuple(bond.first.chain, bond.first.residue_id, bond.first.name));
        const auto second = by_atom.find(std::make_tuple(bond.second.chain, bond.second.residue_id, bond.second.name));
        if (first == by_atom.end() || second == by_atom.end())
        {
            // Routine on a real file: the record may name a HETATM the filter
            // dropped, or a hydrogen the structure does not carry.
            ++unresolved;
            continue;
        }
        try
        {
            _topology.add_spring(first->second, second->second);
            ++added;
        }
        catch (const topology::SpringAlreadyExistsException &)
        {
            // Routine, and worth its own count: a real PDB often declares the
            // same bridge twice, once by SSBOND and once by CONECT. Folding
            // that into "added = 0" made the log read as a parse failure when
            // the record had in fact been understood -- 1S4Q does exactly this.
            ++duplicate;
        }
    }
    logging::info("Read %zu new bond(s) from SSBOND/LINK records; %zu already declared by a CONECT, %zu named "
                  "an atom absent from the model.",
                  added, duplicate, unresolved);
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
