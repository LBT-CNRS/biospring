#ifndef __PARTICLEPROPERTY_H__
#define __PARTICLEPROPERTY_H__

namespace biospring
{
namespace spn
{

// Per-particle physical properties read from the .ff force field file.
// Units (matching the simulation's internal unit system, see
// forcefield/constants.hpp): mass in Dalton (Da), charge in elementary
// charge (e), radius in Angstrom (A), epsilon in kJ.mol-1 (Lennard-Jones
// well depth), solvent accessibility surface in A^2, and two quantities
// that are NOT the same and whose names invite confusing them:
//   - hydrophobicity, the PAIRWISE term's per-particle amplitude, in
//     sqrt(kJ.mol-1.A-1): the law multiplies two of them together and the
//     product is a force amplitude (forcefield/shared/hydrophobic_shared.h).
//     It is written to the .nc as "hydrophobicity".
//   - transfer energy by accessible surface, IMPALA's, in kJ.mol-1.A-2 --
//     per unit of surface, since imp.hpp multiplies it by one
//     (forcefield/energy/imp.hpp). It is written to the .nc as
//     "hydrophobicityscale", which is the trap: the name with "scale" in it
//     is IMPALA's, not the pairwise term's.
class ParticleProperty
{
  public:
    ParticleProperty()
        : _mass(1.0), _charge(0.0), _electroncharge(0), _radius(1.0), _epsilon(0.0), _tempfactor(0.0), _occupancy(0.0),
          _hydrophobicity(0.0), _solventaccessibilitysurface(0.0), _transferenergybyaccessiblesurface(0.0),
          _ischarged(false), _ishydrophobic(false), _burying(1.0)
    {
    }

    void setCharge(float charge);
    float getCharge() const { return _charge; }

    void setElectronCharge(int charge);
    int getElectronCharge() const { return _electroncharge; }

    bool isCharged() const { return _ischarged; }

    void setRadius(float radius) { _radius = radius; }
    float getRadius() const { return _radius; }

    void setMass(float mass) { _mass = mass; }
    float getMass() const { return _mass; }

    float getTempFactor() const { return _tempfactor; }
    void setTempFactor(float temp) { _tempfactor = temp; }

    float getOccupancy() const { return _occupancy; }
    void setOccupancy(float occupancy) { _occupancy = occupancy; }

    float getEpsilon() const { return _epsilon; }
    void setEpsilon(float epsilon) { _epsilon = epsilon; }

    float getHydrophobicity() const { return _hydrophobicity; }
    void setHydrophobicity(float hydrophobicity);
    bool isHydrophobic() const { return _ishydrophobic; }

    float getSolventAccessibilitySurface() const { return _solventaccessibilitysurface; }
    void setSolventAccessibilitySurface(float sas) { _solventaccessibilitysurface = sas; }

    float getTransferEnergyByAccessibleSurface() const { return _transferenergybyaccessiblesurface; }
    void setTransferEnergyByAccessibleSurface(float ener) {_transferenergybyaccessiblesurface = ener;};

    float getBurying() const { return _burying; }
    void setBurying(float burying) { _burying = burying; }

  protected:
  private:
    float _mass;
    float _charge;
    int _electroncharge;
    float _radius;
    float _epsilon;
    float _tempfactor;
    float _occupancy;
    float _hydrophobicity;
    float _solventaccessibilitysurface;
    float _transferenergybyaccessiblesurface;
    bool _ischarged;
    bool _ishydrophobic;
    float _burying;
};

} // namespace spn
} // namespace biospring

#endif // __PARTICLEPROPERTY_H__
