#include "ParticleProperty.h"

namespace biospring
{
namespace spn
{

void ParticleProperty::setCharge(float charge)
{
    if (charge != 0.0)
        _ischarged = true;
    else
        _ischarged = false;
    _charge = charge;
    _electroncharge = (int)_charge;
}

void ParticleProperty::setElectronCharge(int electroncharge)
{
    if (electroncharge != 0)
        _ischarged = true;
    else
        _ischarged = false;
    _electroncharge = electroncharge;
    _charge = (float)electroncharge;
}

void ParticleProperty::setHydrophobicity(float hydrophobicity)
{
    if (hydrophobicity != 0.0)
        _ishydrophobic = true;
    else
        _ishydrophobic = false;
    _hydrophobicity = hydrophobicity;
}

} // namespace spn
} // namespace biospring
