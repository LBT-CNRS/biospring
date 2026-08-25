#include "staticbond/StaticBondBuilder.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "logging.h"
#include "measure.hpp"

namespace logging = biospring::logging;

namespace biospring
{
namespace staticbond
{

// See the header for where each of these numbers comes from.
const float HYDROGEN_BOND_STIFFNESS = 60.0f;
const float HYDROGEN_BOND_CUTOFF = 3.5f;
const float DISULFIDE_STIFFNESS = 1389.1f;
const float DISULFIDE_CUTOFF = 2.5f;
const int HYDROGEN_BOND_MIN_SEQUENCE_GAP = 2;

DonorAcceptorTable readDonorAcceptorTable(const std::string & path)
{
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("cannot open donor/acceptor table '" + path + "'");

    DonorAcceptorTable table;
    std::string line;
    size_t lineno = 0;
    while (std::getline(input, line))
    {
        ++lineno;
        const size_t comment = line.find('#');
        if (comment != std::string::npos)
            line.erase(comment);

        std::istringstream fields(line);
        std::string resname, name;
        int donor = 0, acceptor = 0;
        if (!(fields >> resname >> name >> donor >> acceptor))
            continue; // blank, comment-only, or too short to be a rule

        // Any further column is an antecedent, which orients a dynamic bond
        // and means nothing to a spring. Read past it without complaint so
        // one table serves both mechanisms.
        table[{resname, name}] = {donor, acceptor};
    }

    if (table.empty())
        logging::warning("StaticBondBuilder: '%s' declares no donor or acceptor -- no hydrogen-bond spring "
                         "can be created from it.",
                         path.c_str());
    return table;
}

size_t addHydrogenBondSprings(topology::Topology & topology, const DonorAcceptorTable & table, float cutoff,
                              float stiffness)
{
    // Collect the two roles once. Scanning every pair of particles would be
    // quadratic in the whole system; donors and acceptors are a small
    // fraction of it, and a particle may be both.
    std::vector<size_t> donors, acceptors;
    for (size_t i = 0; i < topology.number_of_particles(); ++i)
    {
        const auto & properties = topology.get_particle(i).properties();
        const auto rule = table.find({properties.residue_name(), properties.name()});
        if (rule == table.end())
            continue;
        if (rule->second.first > 0)
            donors.push_back(i);
        if (rule->second.second > 0)
            acceptors.push_back(i);
    }

    if (donors.empty() || acceptors.empty())
    {
        logging::warning("StaticBondBuilder: the topology holds %zu donor(s) and %zu acceptor(s) -- no "
                         "hydrogen-bond spring created. A table keyed on names the particles do not carry "
                         "(a .grp renames them) would look exactly like this.",
                         donors.size(), acceptors.size());
        return 0;
    }

    size_t created = 0, skipped_existing = 0, skipped_close = 0;
    for (const size_t d : donors)
    {
        auto & donor = topology.get_particle(d);
        for (const size_t a : acceptors)
        {
            if (d == a)
                continue;
            auto & acceptor = topology.get_particle(a);

            // Neighbours in sequence are held by the backbone itself.
            if (donor.properties().chain_name() == acceptor.properties().chain_name() &&
                std::abs(donor.properties().residue_id() - acceptor.properties().residue_id()) <
                    HYDROGEN_BOND_MIN_SEQUENCE_GAP)
            {
                ++skipped_close;
                continue;
            }

            if (measure::distance(donor, acceptor) > cutoff)
                continue;

            // Already sprung by -cutoff, by -rigidbody, or by a CONECT
            // record: leave that spring alone rather than doubling the pair.
            if (topology.springs().exists(donor.unique_id(), acceptor.unique_id()))
            {
                ++skipped_existing;
                continue;
            }

            topology.add_spring(donor, acceptor, -1.0, stiffness);
            ++created;
        }
    }

    logging::info("Hydrogen-bond springs: %zu created from %zu donor(s) and %zu acceptor(s) within %.2f A "
                  "(stiffness %.1f kJ.mol-1.A-2); %zu pair(s) already sprung, %zu too close in sequence.",
                  created, donors.size(), acceptors.size(), cutoff, stiffness, skipped_existing, skipped_close);
    return created;
}

size_t addDisulfideSprings(topology::Topology & topology, float cutoff, float stiffness)
{
    std::vector<size_t> sulfurs;
    size_t cysteines = 0;
    for (size_t i = 0; i < topology.number_of_particles(); ++i)
    {
        const auto & properties = topology.get_particle(i).properties();
        const std::string & resname = properties.residue_name();
        if (resname.rfind("CY", 0) != 0)
            continue;
        ++cysteines;
        const std::string & name = properties.name();
        if (name.size() >= 2 && name.compare(name.size() - 2, 2, "SG") == 0)
            sulfurs.push_back(i);
    }

    size_t created = 0, recalibrated = 0;
    for (size_t i = 0; i < sulfurs.size(); ++i)
    {
        auto & first = topology.get_particle(sulfurs[i]);
        for (size_t j = i + 1; j < sulfurs.size(); ++j)
        {
            auto & second = topology.get_particle(sulfurs[j]);
            if (measure::distance(first, second) > cutoff)
                continue;

            // A disulfide often arrives as a CONECT record, which pdb2spn has
            // already turned into a spring carrying the GLOBAL -stiffness.
            // That pair is the bridge, so give it the bridge's force constant
            // rather than leaving it at a mesh value it never earned. This is
            // the opposite choice from a hydrogen bond, and for a reason: a
            // pre-existing spring there is a structural neighbour, not the
            // same bond seen twice.
            if (topology.springs().exists(first.unique_id(), second.unique_id()))
            {
                topology.get_spring(first, second).set_stiffness(stiffness);
                ++recalibrated;
                continue;
            }
            topology.add_spring(first, second, -1.0, stiffness);
            ++created;
        }
    }

    // A structure with cysteines but no bridge is perfectly normal -- free
    // thiols exist. A structure with cysteines and no SULFUR is not: it means
    // the name test missed, and saying so is the difference between a model
    // without disulfides and a model that silently lost them.
    if (cysteines > 0 && sulfurs.empty())
        logging::warning("StaticBondBuilder: %zu cysteine particle(s) but not one named SG -- no disulfide "
                         "can be found. Check what the .grp renamed the sulfur to.",
                         cysteines);

    logging::info("Disulfide springs: %zu created and %zu recalibrated (a CONECT record had already sprung "
                  "them) from %zu cysteine sulfur(s) within %.2f A, stiffness %.1f kJ.mol-1.A-2.",
                  created, recalibrated, sulfurs.size(), cutoff, stiffness);
    return created + recalibrated;
}

} // namespace staticbond
} // namespace biospring
