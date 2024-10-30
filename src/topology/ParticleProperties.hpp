#ifndef __PARTICLEPROPERTIES_HPP__
#define __PARTICLEPROPERTIES_HPP__

#include "Vector3f.h"
#include <cmath>
namespace biospring
{
namespace topology
{

class IMPPropertiesBuilder;
class ParticlePropertiesBuilder;

// Float equality.
static bool float_eq(float a, float b) { return std::abs(a - b) < 1e-6; }

// Holds the Impala properties of a particle.
// IMPPropertiesBuilder is builder-pattern compliant.
class IMPProperties
{
  protected:
    float _transfert_energy_by_accessible_surface = 0.0;
    float _solvent_accessible_surface = 0.0;

  public:
    friend class IMPPropertiesBuilder;
    static IMPPropertiesBuilder build();

    float transfert_energy_by_accessible_surface() const { return _transfert_energy_by_accessible_surface; }
    float solvent_accessible_surface() const { return _solvent_accessible_surface; }

    void set_transfert_energy_by_accessible_surface(float transfert_energy_by_accessible_surface)
    {
        _transfert_energy_by_accessible_surface = transfert_energy_by_accessible_surface;
    }
    void set_solvent_accessible_surface(float solvent_accessible_surface)
    {
        _solvent_accessible_surface = solvent_accessible_surface;
    }

    // Equality operator.
    bool operator==(const IMPProperties & other) const
    {
        return float_eq(_transfert_energy_by_accessible_surface, other._transfert_energy_by_accessible_surface) &&
               float_eq(_solvent_accessible_surface, other._solvent_accessible_surface);
    }
};

// IMPPropertiesBuilder.
class IMPPropertiesBuilder
{
  protected:
    IMPProperties _imp;

  public:
    IMPPropertiesBuilder & transfert_energy_by_accessible_surface(float transfert_energy_by_accessible_surface)
    {
        _imp.set_transfert_energy_by_accessible_surface(transfert_energy_by_accessible_surface);
        return *this;
    }
    IMPPropertiesBuilder & solvent_accessible_surface(float solvent_accessible_surface)
    {
        _imp.set_solvent_accessible_surface(solvent_accessible_surface);
        return *this;
    }

    operator IMPProperties() const { return _imp; }
};

// Holds the properties of a particle.
class ParticleProperties
{
  protected:
    float _mass = 1.0;
    float _charge = 0.0;
    float _radius = 1.0;
    float _epsilon = 0.0;
    float _temperature_factor = 0.0;
    float _occupancy = 0.0;
    float _hydrophobicity = 0.0;
    float _burying = 1.0;
    IMPProperties _imp;

    Vector3f _position = Vector3f();
    Vector3f _velocity = Vector3f();
    Vector3f _force = Vector3f();
    Vector3f _previous_position = Vector3f();

    std::string _name;
    std::string _residue_name;
    std::string _chain_name;
    std::string _element_name;

    int _residue_id = 0;
    int _atom_id = 0; // stores the atom id read from the topology files.

    bool _is_static = false;

    size_t _topology_id = 0; // stores the topology id of the particle (needed to merge topologies).

  public:
    friend class ParticlePropertiesBuilder;
    static ParticlePropertiesBuilder build();

    // ==================================================================================
    // Getters.
    // ==================================================================================
    float mass() const { return _mass; }
    float charge() const { return _charge; }
    float radius() const { return _radius; }
    float epsilon() const { return _epsilon; }
    float temperature_factor() const { return _temperature_factor; }
    float occupancy() const { return _occupancy; }
    float hydrophobicity() const { return _hydrophobicity; }
    float burying() const { return _burying; }

    IMPProperties imp() const { return _imp; }

    Vector3f position() const { return _position; }
    Vector3f velocity() const { return _velocity; }
    Vector3f force() const { return _force; }
    Vector3f previous_position() const { return _previous_position; }

    std::string name() const { return _name; }
    std::string residue_name() const { return _residue_name; }
    std::string chain_name() const { return _chain_name; }
    std::string element_name() const { return _element_name; }

    int residue_id() const { return _residue_id; }
    int atom_id() const { return _atom_id; }

    bool is_charged() const { return std::abs(_charge) > 1e-6; }
    bool is_hydrophobic() const { return std::abs(_hydrophobicity) > 1e-6; }
    bool is_static() const { return _is_static; }
    bool is_dynamic() const { return !_is_static; }

    size_t topology_id() const { return _topology_id; }

    // ==================================================================================
    // Setters.
    // ==================================================================================

    // clang-format off
    void set_mass(float mass) { _mass = mass; }
    void set_charge(float charge) { _charge = charge; }
    void set_radius(float radius) { _radius = radius; }
    void set_epsilon(float epsilon) { _epsilon = epsilon; }
    void set_temperature_factor(float temperature_factor) { _temperature_factor = temperature_factor; }
    void set_occupancy(float occupancy) { _occupancy = occupancy; }
    void set_hydrophobicity(float hydrophobicity) { _hydrophobicity = hydrophobicity; }
    void set_burying(float burying) { _burying = burying; }
    void set_imp(IMPProperties imp) { _imp = imp; }

    void set_position(const Vector3f & position) { _previous_position = _position; _position = position; }
    void set_velocity(const Vector3f & velocity) { _velocity = velocity; }
    void set_force(const Vector3f & force) { _force = force; }

    void set_name(const std::string & name) { _name = name; }
    void set_residue_name(const std::string & residue_name) { _residue_name = residue_name; }
    void set_chain_name(const std::string & chain_name) { _chain_name = chain_name; }

    void set_residue_id(int residue_id) { _residue_id = residue_id; }
    void set_atom_id(int atom_id) { _atom_id = atom_id; }

    void set_static(bool is_static = true) { _is_static = is_static; }
    void set_dynamic(bool dynamic = true) { _is_static = !dynamic; }

    void set_topology_id(size_t topology_id) { _topology_id = topology_id; }
    // clang-format on

    // ==================================================================================

    // Equality operator.
    bool operator==(const ParticleProperties & other) const
    {
        return float_eq(_mass, other._mass) && float_eq(_charge, other._charge) && float_eq(_radius, other._radius) &&
               float_eq(_epsilon, other._epsilon) && float_eq(_temperature_factor, other._temperature_factor) &&
               float_eq(_occupancy, other._occupancy) && float_eq(_hydrophobicity, other._hydrophobicity) &&
               float_eq(_burying, other._burying) && _imp == other._imp && _position == other._position &&
               _velocity == other._velocity && _force == other._force &&
               _previous_position == other._previous_position && _name == other._name &&
               _residue_name == other._residue_name && _chain_name == other._chain_name &&
               _element_name == other._element_name && _residue_id == other._residue_id && _atom_id == other._atom_id &&
               _is_static == other._is_static && _topology_id == other._topology_id;
    }
};

// ParticlePropertiesBuilder.
class ParticlePropertiesBuilder
{
  protected:
    ParticleProperties _properties;

  public:
    // clang-format off
    ParticlePropertiesBuilder & mass(float mass) { _properties._mass = mass; return *this;}
    ParticlePropertiesBuilder & charge(float charge) { _properties._charge = charge; return *this;}
    ParticlePropertiesBuilder & radius(float radius) { _properties._radius = radius; return *this;}
    ParticlePropertiesBuilder & epsilon(float epsilon) { _properties._epsilon = epsilon; return *this;}
    ParticlePropertiesBuilder & temperature_factor(float temperature_factor) { _properties._temperature_factor = temperature_factor; return *this;}
    ParticlePropertiesBuilder & occupancy(float occupancy) { _properties._occupancy = occupancy; return *this;}
    ParticlePropertiesBuilder & hydrophobicity(float hydrophobicity) { _properties._hydrophobicity = hydrophobicity; return *this;}
    ParticlePropertiesBuilder & burying(float burying) { _properties._burying = burying; return *this;}
    ParticlePropertiesBuilder & imp(IMPProperties imp) { _properties._imp = imp; return *this;}

    ParticlePropertiesBuilder & position(const Vector3f & position) { _properties._position = position; return *this;}
    ParticlePropertiesBuilder & velocity(const Vector3f & velocity) { _properties._velocity = velocity; return *this;}
    ParticlePropertiesBuilder & force(const Vector3f & force) { _properties._force = force; return *this;}

    ParticlePropertiesBuilder & name(const std::string & name) { _properties._name = name; return *this;}
    ParticlePropertiesBuilder & residue_name(const std::string & residue_name) { _properties._residue_name = residue_name; return *this;}
    ParticlePropertiesBuilder & chain_name(const std::string & chain_name) { _properties._chain_name = chain_name; return *this;}
    ParticlePropertiesBuilder & element_name(const std::string & element_name) { _properties._element_name = element_name; return *this;}

    ParticlePropertiesBuilder & residue_id(int residue_id) { _properties._residue_id = residue_id; return *this;}
    ParticlePropertiesBuilder & atom_id(int atom_id) { _properties._atom_id = atom_id; return *this;}

    ParticlePropertiesBuilder & dynamic(bool dynamic = true) { _properties._is_static = !dynamic; return *this;}

    ParticlePropertiesBuilder & topology_id(size_t topology_id) { _properties._topology_id = topology_id; return *this;}

    operator ParticleProperties() const { return _properties; }
    // clang-format on
};

} // namespace topology
} // namespace biospring

#endif // __PARTICLEPROPERTIES_HPP__