#ifndef __SPRINGCOLLECTION_HPP__
#define __SPRINGCOLLECTION_HPP__

#include "ParticleCollection.hpp"
#include "Spring.hpp"

#include <unordered_map>
#include <utility> // std::hash
#include <vector>

namespace biospring
{
namespace topology
{

// =============================================================================
//
// Exceptions.
//
// =============================================================================

// Exception raised when trying to add a spring between particles that are not in the
// collection.
class InvalidParticleException : public std::exception
{
  protected:
    std::string _msg;

  public:
    InvalidParticleException(const std::string & msg) : _msg(msg) {}

    const char * what() const noexcept override { return _msg.c_str(); }
};

// Exception raised when trying to add a spring that already exists.
class SpringAlreadyExistsException : public std::exception
{
  protected:
    std::string _msg;

  public:
    SpringAlreadyExistsException(const std::string & msg) : _msg(msg) {}

    const char * what() const noexcept override { return _msg.c_str(); }
};

// Exception raised when trying to add a spring between a particle and itself.
class SelfSpringException : public std::exception
{
  protected:
    std::string _msg;

  public:
    SelfSpringException(const std::string & msg) : _msg(msg) {}

    const char * what() const noexcept override { return _msg.c_str(); }
};

// Specialization of `std::hash` for spring identifiers.
// This allows us to use `sid_t` as keys in `std::unordered_map`.
struct sid_hash
{
    std::size_t operator()(const std::pair<pid_t, pid_t> & uid) const
    {
        return std::hash<pid_t>()(uid.first) ^ std::hash<pid_t>()(uid.second);
    }
};

// =============================================================================
//
// SpringCollection class.
//
// A `SpringCollection` is a container for `Spring` objects.
//
// A `SpringCollection` can only exists within the context of a partner `ParticleCollection`.
// This is because the `SpringCollection` is a container for `Spring` objects, which are
// defined by pairs of `Particle` objects.
//
// =============================================================================

class SpringCollection
{
  public:
    using value_type = Spring;
    using size_type = std::size_t;
    using reference = value_type &;
    using const_reference = const value_type &;
    using pointer = value_type *;
    using iterator = std::vector<value_type>::iterator;
    using const_iterator = std::vector<value_type>::const_iterator;

  protected:
    std::vector<Spring> _data;
    std::unordered_map<sid_t, size_t, sid_hash> _by_uid;
    ParticleCollection & _particles;

  public:
    SpringCollection(ParticleCollection & particles) : _particles(particles) {}

    // =============================================================================
    // Access to springs.
    // =============================================================================

    // Returns a reference to the spring at the given position.
    reference operator[](size_type pos) { return _data[pos]; }
    const_reference operator[](size_type pos) const { return _data[pos]; }

    // Returns a reference to the spring with the given id.
    reference at_uid(sid_t uid) { return _data[_by_uid.at(uid)]; }
    const_reference at_uid(sid_t uid) const { return _data[_by_uid.at(uid)]; }

    // Returns an iterator to the first particle in the collection.
    iterator begin() { return _data.begin(); }
    const_iterator begin() const { return _data.begin(); }

    // Returns an iterator to the last particle in the collection.
    iterator end() { return _data.end(); }
    const_iterator end() const { return _data.end(); }

    // Returns a reference to the underlying vector of springs.
    auto & data() { return _data; }
    const auto & data() const { return _data; }

    // Returns a reference to the underlying map of springs.
    auto & by_uid() { return _by_uid; }
    const auto & by_uid() const { return _by_uid; }

    // Return a reference to the underlying particle collection.
    auto & particles() { return _particles; }
    const auto & particles() const { return _particles; }

    // =============================================================================
    // Modifiers.
    // =============================================================================

    // Adds a spring to the collection.
    // The spring must be initialized with two particles that are already in the
    // collection.
    auto & add_spring(Particle & p1, Particle & p2, double equilibrium = -1.0, double stiffness = 1.0)
    {
        // Checks that the particles are in the collection.
        if (!_particles.contains(p1) || !_particles.contains(p2))
            throw InvalidParticleException("Cannot add a spring between particles that are not in the collection.");

        // Checks that the particles are not the same.
        if (p1.unique_id() == p2.unique_id())
            throw SelfSpringException("Cannot add a spring between a particle and itself.");

        // Checks that the spring does not already exist.
        if (exists(p1, p2))
            throw SpringAlreadyExistsException("Cannot add a spring that already exists.");

        // Adds the spring to the vector of springs.
        _data.push_back(Spring(p1, p2, equilibrium, stiffness));

        // Adds the spring to the map of springs.
        _by_uid[_data.back().uid()] = _data.size() - 1;

        return _data.back();
    }

    // + operator.
    // Returns a new SpringCollection that contains the springs of both collections.
    SpringCollection operator+(const SpringCollection & other)
    {
        SpringCollection result(_particles);

        result += *this;
        result += other;

        return result;
    }

    SpringCollection & operator+=(const SpringCollection & other)
    {
        // Checks that the collections are compatible.
        if (&_particles != &other._particles)
            throw std::invalid_argument(
                "Cannot add two SpringCollection objects that are not made from the same particle collections.");

        // Iterates over the springs in the other collection.
        for (const auto & [uid, pos] : other.by_uid())
        {
            const Spring & source = other[pos];

            // Gets the particles from the collection.
            auto & p1 = _particles.at_uid(source.first().unique_id());
            auto & p2 = _particles.at_uid(source.second().unique_id());

            add_spring(p1, p2, source.equilibrium(), source.stiffness());
        }

        return *this;
    }

    void clear()
    {
        _data.clear();
        _by_uid.clear();
    }

    // =============================================================================
    // Getters / Setters.
    // =============================================================================
    void set_stiffnesses(double stiffness)
    {
        for (auto & spring : _data)
            spring.set_stiffness(stiffness);
    }

    // =============================================================================
    // Capacity.
    // =============================================================================

    // Returns the number of springs in the collection.
    size_t size() const { return _data.size(); }

    // Returns true if the collection is empty.
    bool empty() const { return _data.empty(); }

    // =============================================================================
    // Lookup.
    // =============================================================================

    // Returns true if the spring exists.
    bool exists(sid_t uid) const { return _by_uid.contains(uid) || _by_uid.contains({uid.second, uid.first}); }

    bool exists(pid_t uid1, pid_t uid2) const
    {
        return _by_uid.contains({uid1, uid2}) || _by_uid.contains({uid2, uid1});
    }

    bool exists(const Particle & p1, const Particle & p2) const { return exists(p1.unique_id(), p2.unique_id()); }
};

} // namespace topology
} // namespace biospring

#endif // __SPRINGCOLLECTION_HPP__