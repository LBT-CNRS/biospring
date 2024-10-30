#ifndef __PARTICLEPROPERTY_H__
#define __PARTICLEPROPERTY_H__

namespace biospring
{
namespace spn
{

class ParticleProperty
{
  public:
    ParticleProperty()
        : _mass(1.0), _charge(0.0), _electroncharge(0), _radius(1.0), _epsilon(0.0), _tempfactor(0.0), _occupancy(0.0),
          _solventaccessibilitysurface(0.0), _transferenergybyaccessiblesurface(0.0), _ischarged(false), _burying(1.0)
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
