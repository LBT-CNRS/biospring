#ifndef _FORCEFIELD_H_
#define _FORCEFIELD_H_

#include "ParticleProperty.h"
#include <string>

#include <map>
#include <unordered_map>

#include "Vector3f.h"
#include "energy/spring.hpp"

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
          _forcefieldscale(1.0), _coulombscale(1.0), _hydrophobicityscale(1.0), _hydrophobicitydecaylength(10.0), _dielectric(1.0)
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

    virtual float computeElectrostaticFieldEnergy(float potential, float charge) const;

    virtual float computeElectrostaticEnergy(float charge1, float charge2, float distance) const;
    virtual float computeElectrostaticForceModule(float charge1, float charge2, float distance) const;

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

    virtual float computeIMPEnergy(float x, float y, float z, float surface, float transfer) const;
    virtual Vector3f computeIMPForceVector(float x, float y, float z, float surface, float transfer) const;

    virtual float computeHydrophobicityEnergy(float hydrophobicity1, float hydrophobicity2, float distance) const;
    virtual float computeHydrophobicityForceModule(float hydrophobicity1, float hydrophobicity2, float distance) const;

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

  private:
    propertiesmap _propertiesfromname;
};

} // namespace forcefield
} // namespace biospring

#endif
