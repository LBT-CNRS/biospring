#include "staticbond/StaticBondBuilder.h"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "logging.h"

namespace logging = biospring::logging;

namespace biospring
{
namespace staticbond
{

// See the header for where each of these numbers comes from.
const float HYDROGEN_BOND_STIFFNESS = 60.0f;
const float DISULFIDE_STIFFNESS = 1389.1f;

namespace
{

bool isCysteineSulfur(const topology::Particle & particle)
{
    const auto & properties = particle.properties();
    const std::string & resname = properties.residue_name();
    const std::string & name = properties.name();
    return resname.rfind("CY", 0) == 0 && name.size() >= 2 && name.compare(name.size() - 2, 2, "SG") == 0;
}

const DonorAcceptorRole * roleOf(const DonorAcceptorTable & table, const topology::Particle & particle)
{
    const auto found = table.find({particle.properties().residue_name(), particle.properties().name()});
    return found == table.end() ? nullptr : &found->second;
}

} // namespace

DonorAcceptorTable readDonorAcceptorTable(const std::string & path)
{
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("cannot open donor/acceptor table '" + path + "'");

    DonorAcceptorTable table;
    std::string line;
    while (std::getline(input, line))
    {
        const size_t comment = line.find('#');
        if (comment != std::string::npos)
            line.erase(comment);

        std::istringstream fields(line);
        std::string resname, name;
        DonorAcceptorRole role;
        if (!(fields >> resname >> name >> role.donorCapacity >> role.acceptorCapacity))
            continue; // blank, comment-only, or too short to be a rule

        table[{resname, name}] = role;
    }

    if (table.empty())
        logging::warning("StaticBondBuilder: '%s' declares no donor or acceptor -- no hydrogen-bond spring "
                         "can be retuned from it.",
                         path.c_str());
    return table;
}

size_t retuneHydrogenBondSprings(topology::Topology & topology, const DonorAcceptorTable & table, float stiffness)
{
    size_t retuned = 0;
    for (size_t i = 0; i < topology.number_of_springs(); ++i)
    {
        auto & spring = topology.get_spring(i);
        const DonorAcceptorRole * first = roleOf(table, spring.first());
        const DonorAcceptorRole * second = roleOf(table, spring.second());
        if (first == nullptr || second == nullptr)
            continue;

        // Either orientation counts: a CONECT record does not say which end
        // donates, and the spring does not care.
        const bool donorAcceptor = (first->donorCapacity > 0 && second->acceptorCapacity > 0) ||
                                   (second->donorCapacity > 0 && first->acceptorCapacity > 0);
        if (!donorAcceptor)
            continue;

        spring.set_stiffness(stiffness);
        ++retuned;
    }

    if (retuned == 0)
        logging::warning("StaticBondBuilder: not one declared bond joins a donor to an acceptor -- no "
                         "hydrogen-bond spring retuned. A structure whose CONECT records do not describe "
                         "its hydrogen bonds, or a table keyed on names the particles do not carry (a .grp "
                         "renames them), both look exactly like this.");
    else
        logging::info("Hydrogen-bond springs: %zu declared bond(s) retuned to %.1f kJ.mol-1.A-2, equilibrium "
                      "left at the declared length.",
                      retuned, stiffness);
    return retuned;
}

size_t retuneDisulfideSprings(topology::Topology & topology, float stiffness)
{
    size_t retuned = 0, cysteines = 0;
    for (size_t i = 0; i < topology.number_of_particles(); ++i)
        if (topology.get_particle(i).properties().residue_name().rfind("CY", 0) == 0)
            ++cysteines;

    for (size_t i = 0; i < topology.number_of_springs(); ++i)
    {
        auto & spring = topology.get_spring(i);
        if (!isCysteineSulfur(spring.first()) || !isCysteineSulfur(spring.second()))
            continue;
        spring.set_stiffness(stiffness);
        ++retuned;
    }

    // A structure with cysteines and no declared bridge is perfectly normal --
    // free thiols exist, and so do PDB files without CONECT records. Saying so
    // is the difference between a model without disulfides and a model that
    // silently lost them.
    if (retuned == 0 && cysteines > 0)
        logging::warning("StaticBondBuilder: %zu cysteine particle(s) but not one declared S-S bond -- no "
                         "disulfide spring retuned. Check that the structure carries CONECT records for its "
                         "bridges, and that the .grp kept the sulfur named SG.",
                         cysteines);
    else
        logging::info("Disulfide springs: %zu declared bridge(s) retuned to %.1f kJ.mol-1.A-2, equilibrium "
                      "left at the declared length.",
                      retuned, stiffness);
    return retuned;
}

} // namespace staticbond
} // namespace biospring
