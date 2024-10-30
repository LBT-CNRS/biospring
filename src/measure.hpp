#ifndef __MEASURE_HPP__
#define __MEASURE_HPP__

#include <array>
#include <cmath>
#include <vector>

#include "Particle.h"
#include "Vector3f.h"

#include "box.hpp"
#include "concepts.hpp"

class Vector3f;
namespace biospring
{
namespace measure
{

// ======================================================================================
//
// Function prototypes.
//
// ======================================================================================

bool almost_equal(double a, double b, double epsilon = 1e-6);
template <typename Container> bool almost_equal(const Container & lhs, const Container & rhs, double epsilon = 1e-6);

template <concepts::Locatable T> double norm(const T & obj);
template <concepts::Locatable T> double distance(const T & lhs, const T & rhs);
template <concepts::Locatable T1, concepts::Locatable T2> double distance(const T1 & lhs, const T2 & rhs);
template <concepts::Locatable T> double distance(const T & lhs, const std::array<double, 3> & rhs);
double distance(const std::array<double, 3> & lhs, const std::array<double, 3> & rhs);

template <concepts::LocatableContainer ContainerType>
std::array<double, 3> closest_to(const ContainerType & container, const std::array<double, 3> & point);
template <concepts::LocatableContainer ContainerType>
std::array<double, 3> closest_to_origin(const ContainerType & container);
template <concepts::LocatableContainer ContainerType>
std::array<double, 3> closest_to_centroid(const ContainerType & container);

template <concepts::LocatableContainer ContainerType> double rmsd(const ContainerType & lhs, const ContainerType & rhs);
template <concepts::LocatableContainer ContainerType> std::array<double, 3> centroid(const ContainerType & container);
template <concepts::LocatableContainer ContainerType> double radius(const ContainerType & container);
template <concepts::LocatableContainer ContainerType> Box box(const ContainerType & particles);

// ======================================================================================
//
// `closest_to` functions.
//
// ======================================================================================

// Returns the coordinates of the point that is the closest to a given point.
template <concepts::LocatableContainer ContainerType>
std::array<double, 3> closest_to(const ContainerType & container, const std::array<double, 3> & point)
{
    if (container.empty())
        return {0, 0, 0};

    double min_distance = std::numeric_limits<double>::max();
    std::array<double, 3> closest_point{0.0, 0.0, 0.0};

    for (const auto & element : container)
    {
        double d = distance(element, point);
        if (d < min_distance)
        {
            min_distance = d;
            closest_point[0] = concepts::locatable::getX(element);
            closest_point[1] = concepts::locatable::getY(element);
            closest_point[2] = concepts::locatable::getZ(element);
        }
    }

    return closest_point;
}

// Returns the coordinates of the point that is the closest to the origin.
template <concepts::LocatableContainer ContainerType>
std::array<double, 3> closest_to_origin(const ContainerType & container)
{
    return closest_to(container, {0.0, 0.0, 0.0});
}

// Returns the coordinates of the point that is the closest to the centroid of the particles.
template <concepts::LocatableContainer ContainerType>
std::array<double, 3> closest_to_centroid(const ContainerType & container)
{
    return closest_to(container, centroid(container));
}

// ======================================================================================
//
// `rmsd`, `centroid` and `radius` functions.
//
// ======================================================================================

// Returns the RMSD between two containers of `Locatable` objects.
template <concepts::LocatableContainer ContainerType> double rmsd(const ContainerType & lhs, const ContainerType & rhs)
{
    float N = lhs.size();

    if (rhs.size() != N)
        throw std::invalid_argument("cannot calculate RMSD from containers with different sizes");

    if (N == 0)
        return 0.0;

    float sum_distances = 0.0;

    for (size_t i = 0; i < N; i++)
    {
        float d = distance(lhs[i], rhs[i]);
        sum_distances += d * d;
    }

    return std::sqrt(sum_distances / float(N));
}

// Returns the centroid of a container of `Locatable` objects.
template <concepts::LocatableContainer ContainerType> std::array<double, 3> centroid(const ContainerType & container)
{
    if (container.empty())
        return {0, 0, 0};

    std::array<double, 3> cumsum{0.0, 0.0, 0.0};

    for (const auto & element : container)
    {
        cumsum[0] += concepts::locatable::getX(element);
        cumsum[1] += concepts::locatable::getY(element);
        cumsum[2] += concepts::locatable::getZ(element);
    }

    double f = 1.0 / container.size();
    cumsum[0] *= f;
    cumsum[1] *= f;
    cumsum[2] *= f;
    return cumsum;
}

// Returns the radius of the smallest sphere that can contain all the particles in the container.
template <concepts::LocatableContainer ContainerType> double radius(const ContainerType & container)
{
    double max_distance = 0.0;
    std::array<double, 3> center = centroid(container);
    for (const auto & element : container)
    {
        double d = distance(element, center);
        if (d > max_distance)
            max_distance = d;
    }
    return max_distance;
}

// ======================================================================================
//
// `Locatable` functions.
//
// ======================================================================================

// Returns the distance between two points in 3D space.
template <concepts::Locatable T> double distance(const T & lhs, const T & rhs)
{
    const double dx = concepts::locatable::getX(lhs) - concepts::locatable::getX(rhs);
    const double dy = concepts::locatable::getY(lhs) - concepts::locatable::getY(rhs);
    const double dz = concepts::locatable::getZ(lhs) - concepts::locatable::getZ(rhs);
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

template <concepts::Locatable T1, concepts::Locatable T2> double distance(const T1 & lhs, const T2 & rhs)
{
    const double dx = concepts::locatable::getX(lhs) - concepts::locatable::getX(rhs);
    const double dy = concepts::locatable::getY(lhs) - concepts::locatable::getY(rhs);
    const double dz = concepts::locatable::getZ(lhs) - concepts::locatable::getZ(rhs);
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

template <concepts::Locatable T> double distance(const T & lhs, const std::array<double, 3> & rhs)
{
    const double dx = concepts::locatable::getX(lhs) - rhs[0];
    const double dy = concepts::locatable::getY(lhs) - rhs[1];
    const double dz = concepts::locatable::getZ(lhs) - rhs[2];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

// Returns the norm of a vector.
template <concepts::Locatable T> double norm(const T & obj)
{
    const double x = concepts::locatable::getX(obj);
    const double y = concepts::locatable::getY(obj);
    const double z = concepts::locatable::getZ(obj);
    if (std::abs(x) < 1e-40 && std::abs(y) < 1e-40 && std::abs(z) < 1e-40)
        return 0.0;
    return std::sqrt(x * x + y * y + z * z);
}

// ======================================================================================
// measure::box.
//
// Returns the bounding box of a container of `Locatable` objects.
template <concepts::LocatableContainer ContainerType> Box box(const ContainerType & particles)
{
    return Box(particles);
}

} // namespace measure
} // namespace biospring

#endif // __MEASURE_HPP__