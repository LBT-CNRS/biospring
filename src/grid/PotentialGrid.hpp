#ifndef __POTENTIAL_GRID_HPP__
#define __POTENTIAL_GRID_HPP__

#include "DenseGrid.hpp"
#include "Vector3f.h"

namespace biospring
{
namespace grid
{

struct PotentialCell
{
    float scalar;
    Vector3f vector;
};

class PotentialGrid : public DenseGrid<PotentialCell>
{
  public:
    static const double GRADIENT_SCALE;

    // Computes the gradient of each cell of the grid.
    void compute_gradient();
};

} // namespace grid
} // namespace biospring

#endif // __POTENTIAL_GRID_HPP__