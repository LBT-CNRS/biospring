#ifndef __PDB_FORMAT_HPP__
#define __PDB_FORMAT_HPP__

#include "Particle.h"
#include "Spring.h"
#include "topology.hpp"

//#include <format>
#include <string>

namespace biospring
{
namespace io
{
namespace pdbfmt
{

std::string model_header(const size_t model_id) { return "MODEL    " + std::to_string(model_id); }

std::string model_footer() { return "ENDMDL"; }

std::string atom_record(const spn::Particle & p)
{
    static const std::string fmt = "%-6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s";

    // Particle name formatting.
    // If the name has less than 4 letters, it is left aligned on column 14, else it's left align on column 13.
    std::string name;
    if (p.getName().size() < 4)
        name = utils::string::format(" %-3s", p.getName().c_str());
    else
        name = utils::string::format("%-4s", p.getName().c_str());

    // Charge formatting.
    std::string charge = "";
    if (p.getElectronCharge() > 0)
        charge = utils::string::format("%+d", p.getElectronCharge());
    else if (p.getElectronCharge() < 0)
        charge = utils::string::format("%d", p.getElectronCharge());

    return utils::string::format(fmt, "ATOM", p.getId() + 1, name.c_str(), " ", p.getResName().c_str(),
                                 p.getChainName().c_str(), p.getResId(), " ", p.getPosition().getX(),
                                 p.getPosition().getY(), p.getPosition().getZ(), p.getOccupancy(), p.getTempFactor(),
                                 p.getElementName().c_str(), charge.c_str());
}

std::string conect_record(const spn::Spring & s)
{
    const std::string fmt = "CONECT%5d%5d";
    return utils::string::format(fmt, s.getParticle1().getId() + 1, s.getParticle2().getId() + 1);
}

// =====================================================================================
//
//     C++20 std::format is not yet supported by GCC 10.2.0, so we have to use biospring::utils::string::format.
//
// =====================================================================================

// std::string atom_record(const topology::Particle & p)
// {
//     const auto & properties = p.properties();
//     return std::format(
//         "{:6s}{:5d} {:4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:2s}{:2s}",
//         "ATOM", properties.atom_id(), properties.name(), " ", properties.residue_name(), properties.chain_name(),
//         properties.residue_id(), " ", p.getX(), p.getY(), p.getZ(), properties.occupancy(),
//         properties.temperature_factor(), properties.element_name(), properties.charge());
// }

// std::string atom_record(const spn::Particle & p)
// {
//     return std::format(
//         "{:6s}{:5d} {:4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:2s}{:2s}",
//         "ATOM", p.getExtid(), p.getName(), " ", p.getResName(), p.getChainName(), p.getResId(), " ", p.getX(),
//         p.getY(), p.getZ(), p.getOccupancy(), p.getTempFactor(), p.getElementName(), p.getElectronCharge());
// }

// std::string conect_record(const topology::Spring & s)
// {
//     return std::format("CONECT{:5d}{:5d}", s.first().properties().atom_id() + 1, s.second().properties().atom_id() +
//     1);
// }

// std::string conect_record(const spn::Spring & s)
// {
//     return std::format("CONECT{:5d}{:5d}", s.getParticle1().getExtid() + 1, s.getParticle2().getExtid() + 1);
// }

std::string end_of_file() { return "END"; }

} // namespace pdbfmt
} // namespace io
} // namespace biospring

#endif // __PDB_FORMAT_HPP__