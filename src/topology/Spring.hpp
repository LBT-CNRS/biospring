#ifndef __TOPOLOGY_SPRING_HPP__
#define __TOPOLOGY_SPRING_HPP__

#include "Particle.hpp"
#include "measure.hpp"

#include <utility> // std::pair

namespace biospring
{
namespace topology
{

// Spring identifier is made of the two unique ids of its particles.
using sid_t = std::pair<pid_t, pid_t>;

// Spring class.
//
// A `Spring` is a connection between two `Particle`s.
// It has a rest length (`equilibrium`) and a spring constant (`stiffness`).
//
// The two particles are reference to as `first` and `second`.
// They are stored as references to the particles, so they must exist, and
// they must not be deleted while the spring exists.
// It also ensures that the two particles cannot be set to nullptr.

class Spring
{
  protected:
    // The two particles that the spring connects.
    Particle & _first;
    Particle & _second;

    // The equilibrium length of the spring.
    double _equilibrium;

    // The spring constant of the spring.
    double _stiffness;

    // The spring identifier is made of the two unique ids of its particles.
    sid_t _uid;


  public:
    // Initializes a spring with the given particles, equilibrium length, and
    // spring constant.
    Spring(Particle & first, Particle & second, double equilibrium = -1.0, double stiffness = 1.0):
        _first(first), _second(second), _equilibrium(equilibrium), _stiffness(stiffness),
        _uid(generate_uid(first, second))
    {
        // Checks that the particles are not the same.
        if (first.unique_id() == second.unique_id())
            throw std::invalid_argument("Cannot add a spring between a particle and itself.");

        if (equilibrium < 0.0)
            _equilibrium = measure::distance(first, second);
    }

    // ================================================================================
    // Getters and setters
    // ================================================================================

    // Returns the first particle.
    const Particle & first() const { return _first; }
    Particle & first() { return _first; }

    // Returns the second particle.
    const Particle & second() const { return _second; }
    Particle & second() { return _second; }

    // Gets/Sets the equilibrium length.
    double equilibrium() const { return _equilibrium; }
    void set_equilibrium(double equilibrium) { _equilibrium = equilibrium; }

    // Gets/Sets the spring constant.
    double stiffness() const { return _stiffness; }
    void set_stiffness(double stiffness) { _stiffness = stiffness; }


    // ================================================================================
    // Methods.
    // ================================================================================

    // Computes and returns the length of the spring.
    double length() const { return measure::distance(_first, _second); }

    // Returns the spring unique identifier which is the pair of the unique
    // identifiers of the two particles.
    sid_t uid() const { return unique_id(); }
    sid_t unique_id() const { return _uid; }


    // =============================================================================
    // Identifier generation.
    // =============================================================================

    // Returns a unique id for a spring between the two given particles.
    static sid_t generate_uid(const Particle & p1, const Particle & p2)
    {
        return {p1.unique_id(), p2.unique_id()};
    }

};

} // namespace topology
} // namespace biospring


#endif // __TOPOLOGY_SPRING_HPP__