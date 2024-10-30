#ifndef __PARTICLECOLLECTION_HPP__
#define __PARTICLECOLLECTION_HPP__

// Defines the `ParticleCollection` class.
//
// A `ParticleCollection` is a vector-like container for `Particle` objects.
// It provides essential `std::vector` methods such as `push_back`, `empty`, as well as
// several helper methods for accessing and manipulating the particles.

#include "Particle.hpp"
#include "concepts.hpp"

#include <cassert>
#include <limits>
#include <type_traits>
#include <unordered_map>

namespace biospring
{
namespace topology
{

class ParticleCollection
{
  public:
    using value_type = Particle;
    using size_type = std::size_t;
    using reference = value_type &;
    using const_reference = const value_type &;
    using pointer = value_type *;
    using iterator = std::vector<value_type>::iterator;
    using const_iterator = std::vector<value_type>::const_iterator;

  protected:
    // The particles in the collection.
    std::vector<Particle> _data;

    // Maps particle unique ids to their position in the collection.
    std::unordered_map<pid_t, size_t> _by_uid;

  public:
    // =============================================================================
    // Returned by `find` if the particle is not found.
    static const size_t PARTICLE_NOT_FOUND;

    // =============================================================================
    // Constructors.
    // =============================================================================

    // Default constructor.
    ParticleCollection() = default;

    // Constructs a collection from a collection of particles.
    template <typename container> ParticleCollection(const container & particles) { push_back(particles); }

    // Construct from initializer lists.
    ParticleCollection(const std::initializer_list<Particle> & particles) { push_back(particles); }

    // Assignment operator.
    // Copies particles from `other` to `this`.
    ParticleCollection & operator=(const ParticleCollection & other)
    {
        _data.clear();
        _by_uid.clear();
        for (const Particle & particle : other._data)
            push_back(particle);
        return *this;
    }

    // =============================================================================
    // Adds particles.
    // =============================================================================

    // Adds a (copy of a) particle to the end of the collection.
    void push_back(const Particle & particle)
    {
        _data.push_back(particle.copy());
        _by_uid[_data.back().unique_id()] = _data.size() - 1;
    }

    // Adds a severak (copies of) particle to the end of the collection.
    template <typename container> void push_back(const container & particles)
    {
        for (const Particle & particle : particles)
            push_back(particle);
    }

    // Template specialization for initializer lists.
    void push_back(const std::initializer_list<Particle> & particles)
    {
        for (const Particle & particle : particles)
            push_back(particle);
    }

    // The += operator is an alias for the `push_back` methods.
    ParticleCollection & operator+=(const Particle & particle)
    {
        push_back(particle);
        return *this;
    }

    ParticleCollection & operator+=(const std::initializer_list<Particle> & particles)
    {
        push_back(particles);
        return *this;
    }

    template <typename container> ParticleCollection & operator+=(const container & particles)
    {
        push_back(particles);
        return *this;
    }

    // The + operator returns a copy of the collection with the given `Particle` or
    // `ParticleCollection` appended.
    ParticleCollection operator+(const Particle & particle) const
    {
        ParticleCollection copy = *this;
        copy += particle;
        return copy;
    }

    ParticleCollection operator+(const ParticleCollection & container) const
    {
        ParticleCollection copy = *this;
        copy += container;
        return copy;
    }

    // =============================================================================
    // Remove particles.
    // =============================================================================

    // Removes the particle at the given position.
    void remove_particle(size_type pos)
    {
        _data.erase(_data.begin() + pos);
        _by_uid.clear();
        for (size_type i = 0; i < _data.size(); ++i)
            _by_uid[_data[i].unique_id()] = i;
    }

    // Removes the particle with the given unique id.
    void remove_particle_with_uid(pid_t uid) { remove_particle(_by_uid.at(uid)); }

    // Removes all particles from the collection.
    void clear()
    {
        _data.clear();
        _by_uid.clear();
    }

    // =============================================================================
    // Access to particles.
    // =============================================================================

    // Returns a reference to the particle at the given position.
    reference operator[](size_type pos) { return _data[pos]; }
    const_reference operator[](size_type pos) const { return _data[pos]; }
    reference at(size_type pos) { return _data.at(pos); }
    const_reference at(size_type pos) const { return _data.at(pos); }

    // Returns a reference to the particle with the given unique id.
    reference at_uid(pid_t uid) { return _data[_by_uid.at(uid)]; }
    const_reference at_uid(pid_t uid) const { return _data[_by_uid.at(uid)]; }

    // Returns an iterator to the first particle in the collection.
    iterator begin() { return _data.begin(); }
    const_iterator begin() const { return _data.begin(); }

    // Returns an iterator to the last particle in the collection.
    iterator end() { return _data.end(); }
    const_iterator end() const { return _data.end(); }

    auto & data() { return _data; }
    const auto & data() const { return _data; }

    auto & by_uid() { return _by_uid; }
    const auto & by_uid() const { return _by_uid; }

    // =============================================================================
    // Lookup.
    // =============================================================================

    // Returns the position of the given particle in the collection.
    size_type find(const Particle & particle) const { return find(particle.unique_id()); }

    // Returns the position of the particle with the given unique id in the collection.
    // Returns std::numeric_limits<size_type>::max() if the particle is not found.
    size_type find(pid_t uid) const
    {
        auto it = _by_uid.find(uid);
        return it == _by_uid.end() ? PARTICLE_NOT_FOUND : it->second;
    }

    // =============================================================================
    // Capacity.
    // =============================================================================

    // Returns the number of particles in the collection.
    size_type size() const { return _data.size(); }

    // Returns whether the collection is empty.
    bool empty() const { return _data.empty(); }

    // =============================================================================
    // Lookup.
    // =============================================================================

    // Returns whether the collection contains the given particle.
    bool contains(const Particle & particle) const { return contains(particle.unique_id()); }

    // Returns whether the collection contains a particle with the given unique id.
    bool contains(pid_t uid) const { return _by_uid.contains(uid); }

    // =============================================================================
    // Getters and setters for ParticleProperties.
    // =============================================================================

    // Mass.
    std::vector<float> masses() const { return get_property(&ParticleProperties::mass); }
    void set_masses(float mass) { set_property(mass, &ParticleProperties::set_mass); }
    void set_masses(const std::vector<float> & masses) { set_property(masses, &ParticleProperties::set_mass); }

    // Charge.
    std::vector<float> charges() const { return get_property(&ParticleProperties::charge); }
    void set_charges(float charge) { set_property(charge, &ParticleProperties::set_charge); }
    void set_charges(const std::vector<float> & charges) { set_property(charges, &ParticleProperties::set_charge); }

    // Radius.
    std::vector<float> radii() const { return get_property(&ParticleProperties::radius); }
    void set_radii(float radius) { set_property(radius, &ParticleProperties::set_radius); }
    void set_radii(const std::vector<float> & radii) { set_property(radii, &ParticleProperties::set_radius); }

    // Epsilon.
    std::vector<float> epsilons() const { return get_property(&ParticleProperties::epsilon); }
    void set_epsilons(float epsilon) { set_property(epsilon, &ParticleProperties::set_epsilon); }
    void set_epsilons(const std::vector<float> & epsilons) { set_property(epsilons, &ParticleProperties::set_epsilon); }

    // Temperature factor.
    std::vector<float> temperature_factors() const { return get_property(&ParticleProperties::temperature_factor); }
    void set_temperature_factors(float temperature_factor)
    {
        set_property(temperature_factor, &ParticleProperties::set_temperature_factor);
    }
    template <template <typename...> class Container>
    void set_temperature_factors(const std::vector<float> & temperature_factors)
    {
        set_property(temperature_factors, &ParticleProperties::set_temperature_factor);
    }

    // Occupancy.
    std::vector<float> occupancies() const { return get_property(&ParticleProperties::occupancy); }
    void set_occupancies(float occupancy) { set_property(occupancy, &ParticleProperties::set_occupancy); }
    void set_occupancies(const std::vector<float> & occupancies)
    {
        set_property(occupancies, &ParticleProperties::set_occupancy);
    }

    // Hydrophobicity.
    std::vector<float> hydrophobicities() const { return get_property(&ParticleProperties::hydrophobicity); }
    void set_hydrophobicities(float hydrophobicity)
    {
        set_property(hydrophobicity, &ParticleProperties::set_hydrophobicity);
    }
    void set_hydrophobicities(const std::vector<float> & hydrophobicities)
    {
        set_property(hydrophobicities, &ParticleProperties::set_hydrophobicity);
    }

    // Burying.
    std::vector<float> buryings() const { return get_property(&ParticleProperties::burying); }
    void set_buryings(float burying) { set_property(burying, &ParticleProperties::set_burying); }
    void set_buryings(const std::vector<float> & buryings) { set_property(buryings, &ParticleProperties::set_burying); }

    // IMP properties.
    std::vector<IMPProperties> imps() const { return get_property(&ParticleProperties::imp); }
    void set_imps(IMPProperties imp) { set_property(imp, &ParticleProperties::set_imp); }
    void set_imps(const std::vector<IMPProperties> & imps) { set_property(imps, &ParticleProperties::set_imp); }

    // Position.
    std::vector<Vector3f> positions() const { return get_property(&ParticleProperties::position); }
    void set_positions(const Vector3f & position) { set_property(position, &ParticleProperties::set_position); }
    void set_positions(const std::vector<Vector3f> & positions)
    {
        set_property(positions, &ParticleProperties::set_position);
    }

    // Velocity
    std::vector<Vector3f> velocities() const { return get_property(&ParticleProperties::velocity); }
    void set_velocities(const Vector3f & velocity) { set_property(velocity, &ParticleProperties::set_velocity); }
    void set_velocities(const std::vector<Vector3f> & velocities)
    {
        set_property(velocities, &ParticleProperties::set_velocity);
    }

    // Force.
    std::vector<Vector3f> forces() const { return get_property(&ParticleProperties::force); }
    void set_forces(const Vector3f & force) { set_property(force, &ParticleProperties::set_force); }
    void set_forces(const std::vector<Vector3f> & forces) { set_property(forces, &ParticleProperties::set_force); }

    // Name.
    std::vector<std::string> names() const { return get_property(&ParticleProperties::name); }
    void set_names(const std::string & name) { set_property(name, &ParticleProperties::set_name); }
    void set_names(const std::vector<std::string> & names) { set_property(names, &ParticleProperties::set_name); }

    // Residue name.
    std::vector<std::string> residue_names() const { return get_property(&ParticleProperties::residue_name); }
    void set_residue_names(const std::string & residue_name)
    {
        set_property(residue_name, &ParticleProperties::set_residue_name);
    }
    void set_residue_names(const std::vector<std::string> & residue_names)
    {
        set_property(residue_names, &ParticleProperties::set_residue_name);
    }

    // Chain name.
    std::vector<std::string> chain_names() const { return get_property(&ParticleProperties::chain_name); }
    void set_chain_names(const std::string & chain_name)
    {
        set_property(chain_name, &ParticleProperties::set_chain_name);
    }
    void set_chain_names(const std::vector<std::string> & chain_names)
    {
        set_property(chain_names, &ParticleProperties::set_chain_name);
    }

    // Residue id.
    std::vector<int> residue_ids() const { return get_property(&ParticleProperties::residue_id); }
    void set_residue_ids(int residue_id) { set_property(residue_id, &ParticleProperties::set_residue_id); }
    void set_residue_ids(const std::vector<int> & residue_ids)
    {
        set_property(residue_ids, &ParticleProperties::set_residue_id);
    }

    // Atom id.
    std::vector<int> atom_ids() const { return get_property(&ParticleProperties::atom_id); }
    void set_atom_ids(int atom_id) { set_property(atom_id, &ParticleProperties::set_atom_id); }
    void set_atom_ids(const std::vector<int> & atom_ids) { set_property(atom_ids, &ParticleProperties::set_atom_id); }

    // Dynamic.
    std::vector<bool> dynamic() const { return get_property(&ParticleProperties::is_dynamic); }
    void set_dynamic(bool is_dynamic) { set_property(is_dynamic, &ParticleProperties::set_dynamic); }
    void set_dynamic(const std::vector<bool> & is_dynamic)
    {
        set_property(is_dynamic, &ParticleProperties::set_dynamic);
    }

    // Topology id.
    std::vector<size_t> topology_ids() const { return get_property(&ParticleProperties::topology_id); }
    void set_topology_ids(size_t topology_id) { set_property(topology_id, &ParticleProperties::set_topology_id); }
    void set_topology_ids(const std::vector<size_t> & topology_ids)
    {
        set_property(topology_ids, &ParticleProperties::set_topology_id);
    }

  private:
    // ParticleProperties getter and setter template functions.
    template <typename T> std::vector<T> get_property(T (ParticleProperties::*getProperty)() const) const
    {
        std::vector<T> values;
        values.reserve(_data.size());

        for (const Particle & particle : _data)
            values.push_back((particle.properties().*getProperty)());

        return values;
    }

    template <typename T, typename Setter>
    void set_property(const T & value, Setter setter, std::type_identity_t<const T> * = nullptr)
    {
        for (Particle & particle : _data)
            (particle.properties().*setter)(value);
    }

    template <typename T, typename Setter>
    void set_property(const std::vector<T> & values, Setter setter, std::type_identity_t<const T> * = nullptr)
    {
        assert(_data.size() == values.size());
        for (size_t i = 0; i < _data.size(); ++i)
            (_data[i].properties().*setter)(values[i]);
    }
};

// Checks that `ParticleCollection` is a `biospring::concepts::LocatableContainer`.
static_assert(concepts::LocatableContainer<ParticleCollection>);

} // namespace topology
} // namespace biospring

#endif // __PARTICLECOLLECTION_HPP__