#include "GridCoordinatesSystem.hpp"
#include "measure.hpp"

// Equality operator: template specialization for double.
template <> bool biospring::grid::_coordinates<double>::operator==(const _coordinates<double> & other) const
{
    return measure::almost_equal(x, other.x) && measure::almost_equal(y, other.y) && measure::almost_equal(z, other.z);
}