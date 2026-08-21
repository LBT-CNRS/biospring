#include "IO/PDBWriter.h"

#include "Particle.h"
#include "SpringNetwork.h"
#include "utils/string.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>

static std::string model_header(const size_t model_id)
{
    return biospring::utils::string::format("MODEL    %5zu", model_id);
}

static std::string model_footer() { return biospring::utils::string::format("ENDMDL"); }

static std::string atom_record(const biospring::spn::Particle & p)
{
    static const std::string fmt = "%-6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s";

    // Particle name formatting.
    // If the name has less than 4 letters, it is left aligned on column 14, else it's left align on column 13.
    // Names of 5+ characters (e.g. a renamed particle from a coarse-grain
    // .grp file) are truncated to the standard PDB 4-character atom-name
    // field ("%-4.4s": width 4, precision 4) -- without this, a longer name
    // pushes every following column (resname, chain, resid, coordinates) one
    // character to the right, which PDBReader's fixed-column parsing cannot
    // recover from (it throws trying to parse a shifted, non-numeric resid).
    std::string name;
    if (p.getName().size() < 4)
        name = biospring::utils::string::format(" %-3.3s", p.getName().c_str());
    else
        name = biospring::utils::string::format("%-4.4s", p.getName().c_str());

    // Charge formatting.
    std::string charge = "";
    if (p.getElectronCharge() > 0)
        charge = biospring::utils::string::format("%+d", p.getElectronCharge());
    else if (p.getElectronCharge() < 0)
        charge = biospring::utils::string::format("%d", p.getElectronCharge());

    return biospring::utils::string::format(fmt, "ATOM", p.getId() + 1, name.c_str(), " ", p.getResName().c_str(),
                                            p.getChainName().c_str(), p.getResId(), " ", p.getPosition().getX(),
                                            p.getPosition().getY(), p.getPosition().getZ(), p.getOccupancy(),
                                            p.getTempFactor(), p.getElementName().c_str(), charge.c_str());
}

static std::string conect_record(const biospring::spn::Spring & s)
{
    static const std::string fmt = "CONECT%5d%5d";
    return biospring::utils::string::format(fmt, s.getParticle1().getId() + 1, s.getParticle2().getId() + 1);
}

void PDBWriter::writeModel(size_t modelid)
{
    // Writes MODEL record.
    _ostream << model_header(modelid) << std::endl;

    // Writes ATOM records.
    for (const auto & particle : _spn->getParticles())
    {
        _ostream << atom_record(particle) << std::endl;
    }

    // Writes CONECT records if required.
    if (_isconnect)
    {
        for (const auto & spring : _spn->getSprings())
        {
            _ostream << conect_record(spring) << std::endl;
        }
    }

    // Writes MODEL footer.
    _ostream << model_footer() << std::endl;
}

void PDBWriter::write()
{
    safeOpen();
    writeModel(1);
    _ostream << "END";
    _ostream.close();
}
