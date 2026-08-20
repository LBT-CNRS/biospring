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
// well depth), hydrophobicity in kJ.mol-1 (see NetCDFWriter's
// "hydrophobicityscale" units), solvent accessibility surface in A^2,
// transfer energy by accessible surface in kJ.mol-1.A-2 (see IMPALA,
// forcefield/energy/imp.hpp).
class ParticleProperty
{
  public:
    ParticleProperty()
        : _mass(1.0), _charge(0.0), _electroncharge(0), _radius(1.0), _epsilon(0.0), _tempfactor(0.0), _occupancy(0.0),
          _hydrophobicity(0.0), _solventaccessibilitysurface(0.0), _transferenergybyaccessiblesurface(0.0),
          _ischarged(false), _ishydrophobic(false), _burying(1.0), _donorcapacity(0), _acceptorcapacity(0), _antecedentindex(-1)
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

    // Hydrogen-bond donor/acceptor CAPACITY: how many bonds this atom can
    // hold at once in each role, not a yes/no. Chemistry sets it -- a
    // donatable hydrogen each for a donor (an amino nitrogen has two), a
    // lone pair each for an acceptor (a carbonyl oxygen has two) -- and it
    // is read from the .hbond table, whose columns are counts. 0/1 remains
    // valid and means exactly what it always did, so a table written before
    // capacities existed keeps its behaviour. A particle can be both (a
    // Ser/Thr/Tyr hydroxyl donates one and accepts two).
    unsigned donorCapacity() const { return _donorcapacity; }
    void setDonorCapacity(unsigned n) { _donorcapacity = n; }

    unsigned acceptorCapacity() const { return _acceptorcapacity; }
    void setAcceptorCapacity(unsigned n) { _acceptorcapacity = n; }

    bool isDonor() const { return _donorcapacity > 0; }
    bool isAcceptor() const { return _acceptorcapacity > 0; }

    // Index of the heavy atom this donor/acceptor hangs off, or -1 when the
    // .hbond table names none. It is what gives a hydrogen bond a direction
    // without an explicit hydrogen: the antecedent->self vector stands in
    // for where the H (or the lone pair) points, and the angle between it
    // and self->partner weights the Morse well (see
    // forcefield::hydrogen_bond_angular_factor). Measured on a B-DNA duplex,
    // that weight is ~0.25 on a real Watson-Crick bond and 0.003 on a
    // stacked same-strand pair that the distance criterion alone accepts.
    int antecedentIndex() const { return _antecedentindex; }
    void setAntecedentIndex(int index) { _antecedentindex = index; }

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
    unsigned _donorcapacity;
    unsigned _acceptorcapacity;
    int _antecedentindex;
};

} // namespace spn
} // namespace biospring

#endif // __PARTICLEPROPERTY_H__
