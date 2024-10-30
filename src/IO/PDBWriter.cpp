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
    std::string name;
    if (p.getName().size() < 4)
        name = biospring::utils::string::format(" %-3s", p.getName().c_str());
    else
        name = biospring::utils::string::format("%-4s", p.getName().c_str());

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

    // Writes PROBE record if required.
    if (_spn->isProbeEnabled())
    {
        biospring::spn::Particle probe;
        probe.setId(99999);
        probe.setName("PRB");
        probe.setResName("PRB");
        probe.setResId(9999);
        probe.setPosition(Vector3f(0.0, 0.0, 0.0));
        probe.setOccupancy(0.0);
        probe.setTempFactor(0.0);
        probe.setElementName("");
        _ostream << atom_record(probe) << std::endl;
    }

    // Writes CONECT records if required.
    _isconnect = _spn->isSpringEnabled() || _isconnect;
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
