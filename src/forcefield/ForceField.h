#ifndef _FORCEFIELD_H_
#define _FORCEFIELD_H_

#include "ParticleProperty.h"
#include <string>

#include <map>
#include <unordered_map>

#include "Vector3f.h"
#include "energy/spring.hpp"
#include "energy/electrostatic.hpp"
#include "energy/hydrophobic.hpp"

namespace biospring
{
namespace forcefield
{

typedef std::map<std::string, biospring::spn::ParticleProperty> propertiesmap;

class ForceField
{
  public:
    ForceField()
        : _stericscale(1.0), _springscale(1.0), _impscale(1.0), _impuppermebraneoffset(0.0), _implowermembraneoffset(0.0),
          _uppermembtubecurv(0.0), _lowermembtubecurv(0.0),
          _forcefieldscale(1.0), _coulombscale(1.0), _hydrophobicityscale(1.0),
          _hydrophobicitydecaylength(10.0), _dielectric(1.0), _distancedependentdielectric(false),
          _hydrogenbondscale(1.0), _hydrogenbondwelldepth(16.7), _hydrogenbondequilibrium(2.9),
          _hydrogenbondwidth(1.34)
    {
    }
    virtual ~ForceField() {}

    // Assignement operator.
    ForceField & operator=(const ForceField & other)
    {
        _stericscale = other._stericscale;
        _springscale = other._springscale;
        _impscale = other._impscale;
        _forcefieldscale = other._forcefieldscale;
        _coulombscale = other._coulombscale;
        _hydrophobicityscale = other._hydrophobicityscale;
        _hydrophobicitydecaylength = other._hydrophobicitydecaylength;
        _dielectric = other._dielectric;
        _distancedependentdielectric = other._distancedependentdielectric;
        _hydrogenbondscale = other._hydrogenbondscale;
        _hydrogenbondwelldepth = other._hydrogenbondwelldepth;
        _hydrogenbondequilibrium = other._hydrogenbondequilibrium;
        _hydrogenbondwidth = other._hydrogenbondwidth;
        _propertiesfromname = other._propertiesfromname;
        return *this;
    }

    void addPropertiesFromName(const std::string name, const biospring::spn::ParticleProperty & pp);
    biospring::spn::ParticleProperty getPropertiesFromName(const std::string & name) const;
    void print();

    size_t getNumberOfProperties() const { return _propertiesfromname.size(); }

    bool hasProperty(const std::string & name) const
    {
        return _propertiesfromname.find(name) != _propertiesfromname.end();
    }

    float computeElectrostaticFieldEnergy(float potential, float charge) const;

    // NOT virtual, and defined here, for the same reason as the spring pair
    // below: no subclass in the tree overrides either, and the four that exist
    // override the STERIC pair instead. Every charged pair of every step was
    // paying an indirect call, and an uninlinable one, to reach a multiply.
    float computeElectrostaticEnergy(float charge1, float charge2, float distance) const
    {
        const float eps = _distancedependentdielectric ? _dielectric * distance : _dielectric;
        return _coulombscale * electrostatic_energy(charge1, charge2, distance, eps);
    }
    float computeElectrostaticForceModule(float charge1, float charge2, float distance) const
    {
        if (!_distancedependentdielectric)
            return _coulombscale * electrostatic_force_module(charge1, charge2, distance, _dielectric);
        // epsilon(r) = e0*r makes the energy go as 1/r^2, so -dE/dr is TWICE what
        // substituting e0*r into the 1/r^2 force expression gives. See
        // ForceField::isDistanceDependentDielectric.
        return 2.0f * _coulombscale *
               electrostatic_force_module(charge1, charge2, distance, _dielectric * distance);
    }

    virtual float computeStericEnergy(float radius1, float radius2, float epsilon1, float epsilon2,
                                      float distance) const;
    virtual float computeStericForceModule(float radius1, float radius2, float epsilon1, float epsilon2,
                                           float distance) const;

    // NOT virtual, and defined here rather than in ForceField.cpp, unlike their
    // steric neighbours above.
    //
    // Those are genuinely polymorphic -- four subclasses override
    // computeStericEnergy/ForceModule to switch Lennard-Jones variants. These
    // two never were: no subclass in the tree overrides either, so every spring
    // force in every step paid an indirect call, and an uninlinable one, to
    // reach a multiply and a subtraction. A sampling profile of RecA put the
    // pair at about 4 % of the process's own work, on top of what it cost by
    // blocking inlining at the call site.
    //
    // Making them virtual again would be a deliberate act: a subclass that
    // needs a different spring law can have it, but it should be added with
    // the override, not left standing as an option nobody took.
    float computeSpringEnergy(float distance, float stiffness, float equilibrium) const
    {
        return _springscale * spring_energy(distance, stiffness, equilibrium);
    }

    float computeSpringForceModule(float distance, float stiffness, float equilibrium) const
    {
        return _springscale * spring_force_module(distance, stiffness, equilibrium);
    }

    float computeIMPEnergy(float x, float y, float z, float surface, float transfer) const;
    Vector3f computeIMPForceVector(float x, float y, float z, float surface, float transfer) const;

    // Same again: nothing overrides these, and they sit on the hydrophobic
    // term's per-pair path.
    float computeHydrophobicityEnergy(float hydrophobicity1, float hydrophobicity2, float distance) const
    {
        return _hydrophobicityscale *
               hydrophobic_energy(hydrophobicity1, hydrophobicity2, distance, _hydrophobicitydecaylength);
    }
    float computeHydrophobicityForceModule(float hydrophobicity1, float hydrophobicity2, float distance) const
    {
        return _hydrophobicityscale *
               hydrophobic_force_module(hydrophobicity1, hydrophobicity2, distance, _hydrophobicitydecaylength);
    }

    // Morse potential between a donor and an acceptor heavy atom (see
    // energy/hydrogenbond.hpp for the derivation of the default De/re/a).
    virtual float computeHydrogenBondEnergy(float distance) const;
    virtual float computeHydrogenBondForceModule(float distance) const;

    // ================================================================================
    // Getters and setters
    //
    // All *Scale members below are dimensionless multipliers applied
    // uniformly to the corresponding energy/force term (see msp options
    // spring.scale, steric.gridscale, coulomb.scale, imp.scale,
    // hydrophobicity.scale in doc/MSP_Options.md); they do not change units,
    // only magnitude. _dielectric is likewise dimensionless (relative
    // permittivity, applied on top of vacuum permittivity in
    // electrostatic_energy/electrostatic_force_module).

    float getCoulombScale() const { return _coulombscale; }
    void setCoulombScale(float coulombscale) { _coulombscale = coulombscale; }

    float getForceFieldScale() const { return _forcefieldscale; }
    void setForceFieldScale(float forcefieldscale) { _forcefieldscale = forcefieldscale; }

    float getStericScale() const { return _stericscale; }
    void setStericScale(float stericscale) { _stericscale = stericscale; }

    float getSpringScale() const { return _springscale; }
    void setSpringScale(float springscale) { _springscale = springscale; }

    float getIMPScale() const { return _impscale; }
    void setIMPScale(float impscale) { _impscale = impscale; }

    float getImpDoubleMembraneUpperMembOffset() const { return _impuppermebraneoffset; }
    void setImpDoubleMembraneUpperMembOffset(float offset) { _impuppermebraneoffset = offset; }

    float getImpDoubleMembraneLowerMembOffset() const { return _implowermembraneoffset; }
    void setImpDoubleMembraneLowerMembOffset(float offset) { _implowermembraneoffset = offset; }

    float getImpDoubleMembraneUpperMembTubeCurv() const { return _uppermembtubecurv; }
    void setImpDoubleMembraneUpperMembTubeCurv(float curv) { _uppermembtubecurv = curv; }

    float getImpDoubleMembraneLowerMembTubeCurv() const { return _lowermembtubecurv; }
    void setImpDoubleMembraneLowerMembTubeCurv(float curv) { _lowermembtubecurv = curv; }

    float getHydrophobicityScale() const { return _hydrophobicityscale; }
    void setHydrophobicityScale(float hydrophobicityscale) { _hydrophobicityscale = hydrophobicityscale; }

    // The distance over which the pairwise hydrophobic attraction decays, in A
    // (msp: hydrophobicity.decaylength). A solvent property the modeller picks,
    // like _dielectric below, not a universal constant.
    float getHydrophobicityDecayLength() const { return _hydrophobicitydecaylength; }
    void setHydrophobicityDecayLength(float decaylength) { _hydrophobicitydecaylength = decaylength; }

    float getDielectric() const { return _dielectric; }
    void setDielectric(float dielectric) { _dielectric = dielectric; }

    // Distance-dependent dielectric: epsilon(r) = dielectric * r, with r in
    // Angstrom (Warshel & Levitt 1976). A constant dielectric of 78 screens a
    // pair 3 A apart exactly as hard as one 30 A apart, which is wrong in the
    // direction that matters here: no water molecule fits between two atoms
    // at contact, so the short-range interaction is barely screened in
    // reality. Measured on a native alpha-helix hydrogen bond of ubiquitin,
    // the amide/carbonyl attraction is -13.61 kJ/mol unscreened and -0.17 at
    // 78 -- annihilated.
    //
    // The FORCE is not the same substitution. With epsilon(r) = e0*r the
    // energy goes as 1/r^2, so the force goes as 2/r^3: putting e0*r into the
    // force formula alone would give half the right answer, and nothing but a
    // gradient check would notice.
    bool isDistanceDependentDielectric() const { return _distancedependentdielectric; }
    void setDistanceDependentDielectric(bool on) { _distancedependentdielectric = on; }

    // Hydrogen-bond Morse potential parameters. Defaults derived from
    // literature (see energy/hydrogenbond.hpp); settable mainly for testing.
    float getHydrogenBondScale() const { return _hydrogenbondscale; }
    void setHydrogenBondScale(float scale) { _hydrogenbondscale = scale; }

    float getHydrogenBondWellDepth() const { return _hydrogenbondwelldepth; }
    void setHydrogenBondWellDepth(float welldepth) { _hydrogenbondwelldepth = welldepth; }

    float getHydrogenBondEquilibrium() const { return _hydrogenbondequilibrium; }
    void setHydrogenBondEquilibrium(float equilibrium) { _hydrogenbondequilibrium = equilibrium; }

    float getHydrogenBondWidth() const { return _hydrogenbondwidth; }
    void setHydrogenBondWidth(float width) { _hydrogenbondwidth = width; }

  protected:
    float _stericscale;
    float _springscale;
    float _impscale;
    float _impuppermebraneoffset;
    float _implowermembraneoffset;
    float _uppermembtubecurv;
    float _lowermembtubecurv;
    float _forcefieldscale;
    float _coulombscale;
    float _hydrophobicityscale;
    float _hydrophobicitydecaylength;
    float _dielectric;
    bool _distancedependentdielectric;
    float _hydrogenbondscale;
    float _hydrogenbondwelldepth;
    float _hydrogenbondequilibrium;
    float _hydrogenbondwidth;

  private:
    propertiesmap _propertiesfromname;
};

} // namespace forcefield
} // namespace biospring

#endif
