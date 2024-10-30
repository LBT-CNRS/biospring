#ifndef __BOX_HPP__
#define __BOX_HPP__

#include <array>
#include <ostream>

#include "concepts.hpp"

namespace biospring
{

// Defines a box in 3D space.
class Box
{
  protected:
    // The coordinates of the box boundaries.
    std::array<double, 6> _boundaries;
    std::array<double, 3> _length;

  public:
    // Empty constructor.
    Box() {}

    // Initializes the box from dimensions.
    Box(const std::array<double, 6> & box) { set_boundaries(box); }

    // Initializes the box from dimensions.
    Box(const std::array<double, 3> & origin, const std::array<double, 3> & length)
    {
        double max_x = origin[0] + length[0];
        double max_y = origin[1] + length[1];
        double max_z = origin[2] + length[2];
        set_boundaries({origin[0], origin[1], origin[2], max_x, max_y, max_z});
    }

    // Initializes the box from a container of `Locatable`.
    template <concepts::LocatableContainer ContainerType> Box(const ContainerType & particles)
    {
        // Finds the minimum and maximum coordinates of the system.
        double min_x = concepts::locatable::getX(particles[0]);
        double min_y = concepts::locatable::getY(particles[0]);
        double min_z = concepts::locatable::getZ(particles[0]);

        double max_x = concepts::locatable::getX(particles[0]);
        double max_y = concepts::locatable::getY(particles[0]);
        double max_z = concepts::locatable::getZ(particles[0]);

        for (size_t i = 1; i < particles.size(); i++)
        {
            double x = concepts::locatable::getX(particles[i]);
            double y = concepts::locatable::getY(particles[i]);
            double z = concepts::locatable::getZ(particles[i]);

            if (x < min_x)
                min_x = x;
            if (y < min_y)
                min_y = y;
            if (z < min_z)
                min_z = z;
            if (x > max_x)
                max_x = x;
            if (y > max_y)
                max_y = y;
            if (z > max_z)
                max_z = z;
        }

        set_boundaries({min_x, min_y, min_z, max_x, max_y, max_z});
    }

    // Sets the boundaries of the box.
    // Calculates the length of the box in each dimension.
    void set_boundaries(const std::array<double, 6> & boundaries)
    {
        _boundaries = boundaries;
        _length[0] = _boundaries[3] - _boundaries[0];
        _length[1] = _boundaries[4] - _boundaries[1];
        _length[2] = _boundaries[5] - _boundaries[2];
    }

    // Returns the length of the box in each dimension.
    const std::array<double, 3> & length() const { return _length; }
    double length_x() const { return _length[0]; }
    double length_y() const { return _length[1]; }
    double length_z() const { return _length[2]; }

    // Returns the origin of the box.
    const std::array<double, 3> origin() const { return {origin_x(), origin_y(), origin_z()}; }
    double origin_x() const { return _boundaries[0]; }
    double origin_y() const { return _boundaries[1]; }
    double origin_z() const { return _boundaries[2]; }

    // Returns the minimum and maximum coordinates of the box.
    const std::array<double, 6> & boundaries() const { return _boundaries; }
    std::array<double, 3> min() const { return {min_x(), min_y(), min_z()}; }
    std::array<double, 3> max() const { return {max_x(), max_y(), max_z()}; }

    double min_x() const { return _boundaries[0]; }
    double min_y() const { return _boundaries[1]; }
    double min_z() const { return _boundaries[2]; }
    double max_x() const { return _boundaries[3]; }
    double max_y() const { return _boundaries[4]; }
    double max_z() const { return _boundaries[5]; }

    // Self-printing.
    friend std::ostream & operator<<(std::ostream & os, const Box & box)
    {
        os << "Box(" << box.min_x() << ", " << box.min_y() << ", " << box.min_z() << ", " << box.max_x() << ", "
           << box.max_y() << ", " << box.max_z() << ")";
        return os;
    }
};

} // namespace biospring

#endif // __BOX_HPP__