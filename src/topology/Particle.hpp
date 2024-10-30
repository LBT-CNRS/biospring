#ifndef __PARTICLE_HPP__
#define __PARTICLE_HPP__

#include "ParticleProperties.hpp"
#include <sstream>

namespace biospring
{
namespace topology
{

using pid_t = size_t; // Particle ID type.

class Particle
{
  protected:
    static pid_t uid_; // Unique ID for each particle.

  private:
    static pid_t new_uid() { return Particle::uid_++; }

  public:
    ParticleProperties _properties;
    pid_t _uid;

    Particle() : _properties(ParticleProperties()), _uid(new_uid()) {}
    Particle(const ParticleProperties & properties) : _properties(properties), _uid(new_uid()) {}

    // Copy constructor.
    // Used when a particle is copied to a new memory location, i.e. containers.
    // The unique id is not changed.
    // To create a copy of a `Particle` outside this context, use `Particle::copy()`.
    Particle(const Particle & other) : _properties(other._properties), _uid(other._uid) {}

    // Create a copy of this particle (assigns a new unique id).
    Particle copy() const
    {
        Particle p = Particle(*this);
        p._uid = new_uid();
        return p;
    }

    // Assignment operator.
    Particle & operator=(const Particle & other)
    {
        if (this != &other)
            _properties = other._properties;
        return *this;
    }

    // Getters.
    pid_t uid() const { return unique_id(); }
    pid_t unique_id() const { return _uid; }

    const ParticleProperties & properties() const { return _properties; }
    ParticleProperties & properties() { return _properties; }

    // Generates a string descriptor for this particle.
    // The format is: "<chain>::<residue_name>::<residue_id>::<atom_name>".
    std::string string_description() const
    {
        std::ostringstream oss;
        oss << _properties.chain_name() << "::" << _properties.residue_name() << "::" << _properties.residue_id()
            << "::" << _properties.name();
        return oss.str();
    }

    // =================================================================================
    // Shortcuts to ParticleProperties::position.
    // =================================================================================

    auto position() const { return _properties.position(); }

    template <typename T> void set_position(const T & pos) { _properties.set_position(pos); }

    // getX(), getY(), getZ() functions to make Particle compliant with `Locatable`.
    double getX() const { return _properties.position().getX(); }
    double getY() const { return _properties.position().getY(); }
    double getZ() const { return _properties.position().getZ(); }
};

} // namespace topology
} // namespace biospring

#endif // __PARTICLE_HPP__