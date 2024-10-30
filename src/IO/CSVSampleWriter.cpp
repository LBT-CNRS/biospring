
#include "CSVSampleWriter.h"
#include "InsertionVector.h"
#include "SpringNetwork.h"

#include <iostream>

void CSVSampleWriter::write()
{
    const InsertionVector & iv = _spn->getInsertionVector();
    if (_firstline)
    {
        _firstline = false;
        _ostream << "# Step"
                 << "\t"
                 << "FrameRate(Hz)"
                 << "\t"
                 << "Kinetic energy (kJ.mol-1)";
        if (_spn->isSpringEnabled())
            _ostream << "\t"
                     << "Spring energy (kJ.mol-1)";
        if (_spn->isStericEnabled())
            _ostream << "\t"
                     << "Steric energy (kJ.mol-1)";
        if (_spn->isElectrostaticEnabled())
            _ostream << "\t"
                     << "Electrostatic energy (kJ.mol-1)";
        if (_spn->isIMPEnabled())
            _ostream << "\t"
                     << "IMP energy (kJ.mol-1)";
        if (_spn->isInsertionVectorEnabled())
            _ostream << "\t"
                     << "Insertion Angle (degrees) \tPenetration (A)";


        _ostream << std::endl;
    }

    _ostream << _spn->getNbIterations() << "\t" << _spn->getFrameRate() << "\t"
             << _spn->getKineticEnergy();
    if (_spn->isSpringEnabled())
        _ostream << "\t" << _spn->getSpringEnergy();
    if (_spn->isStericEnabled())
        _ostream << "\t" << _spn->getStericEnergy();
    if (_spn->isElectrostaticEnabled())
        _ostream << "\t" << _spn->getElectrostaticEnergy();
    if (_spn->isIMPEnabled())
        _ostream << "\t" << _spn->getIMPEnergy();
    if (_spn->isInsertionVectorEnabled())
        _ostream << "\t" << iv.getAngle() << "\t" << iv.getInsertionDepth();

    _ostream << endl;
}
