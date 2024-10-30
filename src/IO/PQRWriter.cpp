#include "IO/PQRWriter.h"

#include "Particle.h"
#include "SpringNetwork.h"
#include "utils/string.hpp"

#include <vector>

static std::string atom_record(const biospring::spn::Particle & p)
{
    const std::string pqr_fmt = "%-6s%5i %4s%1s%3s %1s%4i%1s   %8.3f%8.3f%8.3f%8.4f%7.4f            ";

    const float x = p.getPosition().getX();
    const float y = p.getPosition().getY();
    const float z = p.getPosition().getZ();

    return biospring::utils::string::format(pqr_fmt, "ATOM", p.getId() + 1, p.getName().c_str(), " ",
                                            p.getResName().c_str(), " ", p.getResId(), " ", x, y, z, p.getCharge(),
                                            p.getRadius());
}

void PQRWriter::writeModel(size_t modelid)
{
    for (const auto & p : _spn->getParticles())
    {
        _ostream << atom_record(p) << std::endl;
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
        probe.setCharge(0.0);
        probe.setRadius(0.0);
        _ostream << atom_record(probe) << std::endl;
    }
}
