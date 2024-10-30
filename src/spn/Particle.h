#ifndef __PARTICLE_H__
#define __PARTICLE_H__

#include "forcefield/ForceField.h"

#include "ParticleProperty.h"
#include "Vector3f.h"

namespace biospring
{
namespace spn
{

class Spring;
class SpringNetwork;

class Particle : public ParticleProperty
{
  public:
    Particle()
        : ParticleProperty(), _force(Vector3f()), _position(Vector3f()), _velocity(Vector3f()), _isstatic(false),
          _name(""), _resid(0), _resname(""), _chainname(""), _internalstructid(0), _rigidbodyid(0),
          _isrigid(false), _springneighbors(), _extid(-1), _springnetwork(NULL)
    {
    }

    void setSpringNetwork(SpringNetwork * springnetwork) { _springnetwork = springnetwork; }
    SpringNetwork * getSpringNetwork() const { return _springnetwork; }

    void setId(int id) { _id = id; }
    int getId() const { return _id; }

    void setExtid(unsigned id) { _extid = id; }
    unsigned getExtid() const { return _extid; }

    bool isDynamic() const { return !isStatic(); }
    void setDynamic(bool isdynamic) { _isstatic = !isdynamic; }

    bool isStatic() const { return _isstatic; }
    void setStatic(bool isstatic) { _isstatic = isstatic; }

    void setElectrostaticEnergy(float energy) { _electrostaticenergy = energy; }
    float getElectrostaticEnergy() const { return _electrostaticenergy; }

    void setStericEnergy(float energy) { _stericenergy = energy; }
    float getStericEnergy() const { return _stericenergy; }

    void setKineticEnergy(float energy) { _kineticenergy = energy; }
    float getKineticEnergy() const { return _kineticenergy; }

    void setPreviousPosition(const Vector3f & v) { _previousposition = v; }
    Vector3f getPreviousPosition() const { return _previousposition; }

    void setIMPEnergy(float energy) { _impenergy = energy; }
    float getIMPEnergy() const { return _impenergy; };

    void setHydrophobicityEnergy(float energy) { _hydrophobicityenergy = energy; }
    float getHydrophobicityEnergy() const { return _hydrophobicityenergy; }

    void setPosition(const Vector3f & v);
    Vector3f getPosition() const { return _position; }

    void setForce(const Vector3f & v) { _force = v; }
    Vector3f getForce() const { return _force; }
    void addForce(const Vector3f & v) { _force = _force + v; }
    void setPreviousForce() { _previousForce = _force;}
    Vector3f getPreviousForce() const { return _previousForce; }

    void setVelocity(const Vector3f & v) { _velocity = v; }
    Vector3f getVelocity() const { return _velocity; }
    void addVelocity(const Vector3f& v) { _velocity = _velocity + v; }

    float distance(const Particle & p) const;
    static float distance(const Particle & p1, const Particle & p2);

    float getX() const { return _position.getX(); }
    float getY() const { return _position.getY(); }
    float getZ() const { return _position.getZ(); }

    void setX(float x) { _position.setX(x); }
    void setY(float y) { _position.setY(y); }
    void setZ(float z) { _position.setZ(z); }

    unsigned getResId() const { return _resid; }
    void setResId(unsigned resid) { _resid = resid; }

    const std::string & getResName() const { return _resname; }
    void setResName(const std::string & resname) { _resname = resname; }

    const std::string & getName() const { return _name; }
    void setName(const std::string & name) { _name = name; }

    const std::string & getCGName() const { return _cgname; }
    void setCGName(const std::string & cgname) { _cgname = cgname; }

    const std::string & getChainName() const { return _chainname; }
    void setChainName(const std::string & chainname) { _chainname = chainname; }

    const std::string & getElementName() const { return _elementname; }
    void setElementName(const std::string & elementname) { _elementname = elementname; }

    const std::unordered_map<unsigned, Spring *> getSpringNeighbors() const { return _springneighbors; }
    void removeSpringNeighbor(unsigned index) { _springneighbors.erase(index); }
    void addToSpringNeighbors(unsigned index, Spring * spring);

    unsigned getNumberOfSprings() const { return _springneighbors.size(); }

    bool isInSpringNeighbors(unsigned index) const { return _springneighbors.find(index) != _springneighbors.end(); }

    void addElectrostaticFieldForce();
    void addDensityFieldForce();
    void addElectrostaticForceNoGrid(float cutoff);
    void addElectrostaticForce();
    void addStericForce();
    void addIMPForce();
    void addHydrophobicityForce();
    void addElectrostaticProbeForce(Particle & probe);
    void addStericProbeForce(Particle & probe);

    void resetForce();

    void applyViscosity(float viscosity);

    void computeStericNeighbors();
    void computeElectrostaticNeighbors();
    void computeHydrophobicNeighbors();

    void updateFromForceField(const biospring::forcefield::ForceField & ff);
    void IntegrateVelocityVerlet(float timestep);
    void IntegrateEuler(float timestep);

    inline void clearSpringNeighbors() { _springneighbors.clear(); }

    inline unsigned getInternalStructId() const { return _internalstructid; }
    inline void setInternalStructId(const unsigned structid) { _internalstructid = structid; }

    inline unsigned getRigidBodyId() const { return _rigidbodyid; }
    inline void setRigidBodyId(const unsigned rbid) { _rigidbodyid = rbid; }

    inline bool isRigid() const { return _isrigid; }
    inline void setRigid(bool isrigid) { _isrigid = isrigid; }

    std::string tostr() const;

  private:
    void _integrateForce(float timestep);
    void _integrateVelocity(float timestep);

    Vector3f _force;
    Vector3f _previousForce;
    Vector3f _position;
    Vector3f _previousposition;
    Vector3f _velocity;
    bool _isstatic;

    std::string _name;   // name in pdb
    std::string _cgname; // coarse grain name in pdb
    unsigned _resid;
    std::string _resname;
    std::string _chainname;
    std::string _elementname;
    unsigned _internalstructid;
    unsigned _rigidbodyid;
    bool _isrigid;

    std::unordered_map<unsigned, Spring *> _springneighbors;

    float _kineticenergy;
    float _electrostaticenergy;
    float _stericenergy;
    float _impenergy;
    float _hydrophobicityenergy;

    unsigned _extid;
    int _id;
    SpringNetwork * _springnetwork;
};

} // namespace spn
} // namespace biospring

#endif // __PARTICLE_H__
