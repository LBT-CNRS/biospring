#include "CSVTrajectoryWriter.hpp"
#include "SpringNetwork.h"

namespace biospring
{
namespace io
{
namespace modern
{

void CSVTrajectoryWriter::write_header()
{
    _ostream << "# Step"
             << "\t"
             << "FrameRate(Hz)"
             << "\t"
             << "Kinetic energy (kJ.mol-1)";
    if (_topology.isSpringEnabled())
        _ostream << "\t"
                 << "Spring energy (kJ.mol-1)";
    if (_topology.isStericEnabled())
        _ostream << "\t"
                 << "Steric energy (kJ.mol-1)";
    if (_topology.isElectrostaticEnabled())
        _ostream << "\t"
                 << "Electrostatic energy (kJ.mol-1)";
    if (_topology.isIMPEnabled())
        _ostream << "\t"
                 << "IMP energy (kJ.mol-1)";
    if (_topology.isInsertionVectorEnabled())
        _ostream << "\t"
                 << "Insertion Angle (degrees) \tInsertion Depth (A)";
    _ostream << std::endl;
}

void CSVTrajectoryWriter::write_step()
{
    _ostream << _topology.getNbIterations() << "\t" << _topology.getFrameRate() << "\t" << _topology.getKineticEnergy();
    if (_topology.isSpringEnabled())
        _ostream << "\t" << _topology.getSpringEnergy();
    if (_topology.isStericEnabled())
        _ostream << "\t" << _topology.getStericEnergy();
    if (_topology.isElectrostaticEnabled())
        _ostream << "\t" << _topology.getElectrostaticEnergy();
    if (_topology.isIMPEnabled())
        _ostream << "\t" << _topology.getIMPEnergy();
    if (_topology.isInsertionVectorEnabled())
    {
        const InsertionVector & iv = _topology.getInsertionVector();
        _ostream << "\t" << iv.getAngle() << "\t" << iv.getInsertionDepth();
    }

    _ostream << std::endl;
}

} // namespace modern
} // namespace io
} // namespace biospring
