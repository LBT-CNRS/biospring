
#include <vector>
#include "XTCTrajWriter.h"

#include "Particle.h"
#include "SpringNetwork.h"
#include "logging.h"

#include "xdrfile.h"
#include "xdrfile_xtc.h"

void XTCTrajWriter::writeNextStep()
{
    /* copy the atom coordinates into a float array */
    size_t natoms = _spn->getNumberOfParticles();
    std::vector<rvec> x_xtc(natoms);
    for (size_t i = 0; i < natoms; ++i)
    {
        const biospring::spn::Particle & p = _spn->getParticle(i);
        x_xtc[i][0] = p.getX() / 10.0;
        x_xtc[i][1] = p.getY() / 10.0;
        x_xtc[i][2] = p.getZ() / 10.0;
    }

    /* the box */
    matrix box = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};

    /* write the step in the xdr file */
    write_xtc(_xdr, natoms, _step, (float)_step, box, x_xtc, 1000.0);
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
