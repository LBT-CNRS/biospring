#ifndef _FORCEFIELD_H_
#define _FORCEFIELD_H_

#include "ParticleProperty.h"
#include <string>

#include <map>
#include <unordered_map>

#include "Vector3f.h"

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
          _forcefieldscale(1.0), _coulombscale(1.0), _hydrophobicityscale(1.0), _dielectric(1.0)
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

    virtual float computeSpringEnergy(float distance, float stiffness, float equilibrium) const;
    virtual float computeSpringForceModule(float distance, float stiffness, float equilibrium) const;

    virtual float computeIMPEnergy(float x, float y, float z, float surface, float transfer) const;
    virtual Vector3f computeIMPForceVector(float x, float y, float z, float surface, float transfer) const;

    virtual float computeHydrophobicityEnergy(float hydrophobicity1, float hydrophobicity2, float distance) const;
    virtual float computeHydrophobicityForceModule(float hydrophobicity1, float hydrophobicity2, float distance) const;

    // ================================================================================
    // Getters and setters

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
    float _dielectric;

  private:
    propertiesmap _propertiesfromname;
};

} // namespace forcefield
} // namespace biospring

#endif
