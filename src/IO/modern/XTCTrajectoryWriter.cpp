#include "XTCTrajectoryWriter.hpp"
#include "logging.h"
#include "topology.hpp"
#include "xdrfile_xtc.h"

#include <vector>

namespace biospring
{
namespace io
{
namespace modern
{

void XTCTrajectoryWriter::write_step()
{
    // Copy atoms coordinates into a float array.
    size_t natoms = _topology.getNumberOfParticles();
    std::vector<float> x_xtc(3 * natoms);
    for (size_t i = 0; i < natoms; ++i)
    {
        const auto & particle = _topology.getParticle(i);
        x_xtc[3 * i + 0] = particle.getX() / 10.0;
        x_xtc[3 * i + 1] = particle.getY() / 10.0;
        x_xtc[3 * i + 2] = particle.getZ() / 10.0;
    }

    // The box.
    matrix box = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};

    // Write the step in the xdr file.
    write_xtc(_xdr, natoms, _current_frame, static_cast<float>(_current_frame), box, reinterpret_cast<rvec *>(x_xtc.data()), 1000.0);
    _current_frame++;
}

void XTCTrajectoryWriter::_open_xdr()
{
    if (_path.empty())
        biospring::logging::die("xtctrajectory.enable is set but xtctrajectory.path is empty.");

    if (!_xdr)
    {
        _xdr = xdrfile_open(_path.c_str(), "w");
        if (!_xdr)
            biospring::logging::die("Cannot open %s", _path.c_str());
    }
}

XTCTrajectoryWriter::~XTCTrajectoryWriter()
{
    if (_xdr)
        xdrfile_close(_xdr);
}

} // namespace modern
} // namespace io
} // namespace biospring
