#include "measure.hpp"

namespace biospring
{
namespace measure
{

double distance(const std::array<double, 3> & lhs, const std::array<double, 3> & rhs)
{
    const double dx = lhs[0] - rhs[0];
    const double dy = lhs[1] - rhs[1];
    const double dz = lhs[2] - rhs[2];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

// ======================================================================================
//
// `almost_equal` functions.
//
// ======================================================================================

bool almost_equal(double a, double b, double epsilon) { return std::abs(a - b) < epsilon; }

template <typename Container> bool almost_equal(const Container & lhs, const Container & rhs, double epsilon)
{
    if (lhs.size() != rhs.size())
        return false;

    for (size_t i = 0; i < lhs.size(); i++)
    {
        if (!almost_equal(lhs[i], rhs[i], epsilon))
            return false;
    }

    return true;
}

} // namespace measure
} // namespace biospring