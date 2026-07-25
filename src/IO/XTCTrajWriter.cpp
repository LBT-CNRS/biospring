

#include "XTCTrajWriter.h"

#include "Particle.h"
#include "SpringNetwork.h"
#include "logging.h"

#include "xdrfile.h"
#include "xdrfile_xtc.h"

#include <vector>

void XTCTrajWriter::writeNextStep()
{
    /* copy the atom coordinates into a float array */
    size_t natoms = _spn->getNumberOfParticles();
    std::vector<float> x_xtc(3 * natoms);
    for (size_t i = 0; i < natoms; ++i)
    {
        const biospring::spn::Particle & p = _spn->getParticle(i);
        x_xtc[3 * i + 0] = p.getX() / 10.0;
        x_xtc[3 * i + 1] = p.getY() / 10.0;
        x_xtc[3 * i + 2] = p.getZ() / 10.0;
    }

    /* the box */
    matrix box = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};

    /* write the step in the xdr file */
    write_xtc(_xdr, natoms, _step, static_cast<float>(_step), box, reinterpret_cast<rvec *>(x_xtc.data()), 1000.0);
    _step++;
}

void XTCTrajWriter::safeOpen()
{
    if (not _filename.empty())
    {
        _xdr = xdrfile_open(_filename.c_str(), "w");
        if (not _xdr)
            biospring::logging::die("Cannot open %s", _filename.c_str());
    }
}

void XTCTrajWriter::update() { safeOpen(); }

void XTCTrajWriter::close() { xdrfile_close(_xdr); }
