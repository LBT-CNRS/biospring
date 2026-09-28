#ifndef __PARTICLEPROPERTY_H__
#define __PARTICLEPROPERTY_H__

#include <cmath>

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
          _ischarged(false), _ishydrophobic(false), _burying(1.0), _donorcapacity(0), _acceptorcapacity(0), _antecedentindex(-1), _antecedentindex2(-1), _hblobemode(0), _hblobecos(1.0f), _hblobesin(0.0f)
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

    // A SECOND antecedent, or -1. One is not enough for the commonest donor
    // there is. A backbone amide nitrogen is planar with two heavy
    // neighbours, CA and the previous residue's C, and its hydrogen points
    // opposite their bisector -- 58 degrees away from the CA->N direction a
    // single antecedent gives. Measured on ubiquitin's alpha helix, that
    // costs a factor 3.2 on cos^2(theta): weight 0.28 where the real N-H
    // direction gives 0.91, on bonds whose geometry is ideal.
    //
    // With both set, the direction is -(u1 + u2) normalised, which is exact
    // for any planar sp2 centre -- the protein amide, and equally a guanine
    // N1 or a thymine N3 sitting between two ring carbons. It reproduces the
    // true N-H direction to 0.7 degrees on ubiquitin (18.1 against 17.4),
    // recovering 99 % of the correct weight. This is the geometric form of
    // the rule DSSP has used since Kabsch & Sander 1983.
    int antecedentIndex2() const { return _antecedentindex2; }
    void setAntecedentIndex2(int index) { _antecedentindex2 = index; }

    // Does the second antecedent define the PLANE rather than a second bond?
    //
    // The bisector above is exact for a site with TWO heavy neighbours and ONE
    // hydrogen or lone pair. A site with ONE heavy neighbour and TWO of them is
    // the other common case, and there the bisector is the worst possible
    // answer: an exocyclic amine (guanine N2, adenine N6, cytosine N4) or a
    // carbonyl oxygen (guanine O6, cytosine O2, thymine O4) is planar sp2, so
    // its two hydrogens -- or its two lone pairs -- sit in the base plane at
    // +/- 62 degrees of the C->self axis, and "away from the antecedent" aims
    // exactly BETWEEN them. Measured over 56 amine hydrogens of a B-DNA duplex
    // the angle is 62.0, 62.1 and 62.9 degrees for DG, DA and DC, so the best
    // weight such a site can reach is cos^2(62) = 0.22, and 0.05 is what the
    // running duplex reports.
    //
    // Three of the five Watson-Crick bond types have such a site at BOTH ends,
    // so every base pair kept one bond and lost the rest: 13 of 24 complementary
    // pairs had two bonds at the right distance, 0 of 24 had two at a usable
    // angle. With the two lobes the three broken types go from 0.06-0.08 to
    // 0.89-0.94, which is what the two working types already had.
    //
    // When this is set, antecedentIndex2 is not a bonded neighbour: it is any
    // other atom of the same planar group, and it serves only to fix which
    // plane the lobes lie in. The .hbond table marks it with a '~' prefix.
    // How the second antecedent is read, and at what angle the lobes sit.
    //
    //   0  NONE       antecedent2 is a second bond; the direction is the bisector
    //   1  IN PLANE   the two lobes lie in the plane of the three atoms, which is
    //                 where an sp2 centre's two hydrogens or two lone pairs are:
    //                 an exocyclic amine (62 deg), a carbonyl oxygen
    //   2  OUT OF PLANE  the lobes leave that plane, along its normal. This is the
    //                 Burgi-Dunitz approach to a carbonyl CARBON, at 105 deg of the
    //                 C=O axis and perpendicular to the sp2 plane -- a nucleophile
    //                 attacks a carbonyl from one of its two FACES, and there are
    //                 exactly two such trajectories, so it is the same two-lobe
    //                 problem with a different angle and a different plane.
    //
    // The cosine and sine are stored rather than the angle, because the direction is
    // rebuilt once per candidate pair and a transcendental there would be paid
    // millions of times a step.
    static const int HBOND_LOBES_NONE = 0;
    static const int HBOND_LOBES_IN_PLANE = 1;
    static const int HBOND_LOBES_OUT_OF_PLANE = 2;

    int lobeMode() const { return _hblobemode; }
    float lobeCos() const { return _hblobecos; }
    float lobeSin() const { return _hblobesin; }
    bool hasLobes() const { return _hblobemode != HBOND_LOBES_NONE; }
    void setLobes(int mode, float degrees)
    {
        _hblobemode = mode;
        const float r = degrees * 3.14159265358979323846f / 180.0f;
        _hblobecos = std::cos(r);
        _hblobesin = std::sin(r);
    }

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
    int _antecedentindex2;
    int _hblobemode;
    float _hblobecos;
    float _hblobesin;
};

} // namespace spn
} // namespace biospring

#endif // __PARTICLEPROPERTY_H__
