#include "PotentialGrid.hpp"
#include "forcefield/constants.hpp"
#include "logging.h"

namespace biospring
{
namespace grid
{

// Initialize gradient scale.
const double PotentialGrid::GRADIENT_SCALE = biospring::forcefield::BOLTZMANJPERK *
                                             biospring::forcefield::METER_TO_ANGSTROM *
                                             biospring::forcefield::NEWTON_TO_DALTON_ANGSTROM_PER_FEMTOSECOND_2;

static float compute_gradient_(float current, float previous, float next, float cell_size)
{
    return ((current - previous) + (next - current)) / (cell_size * 2.0);
}

void PotentialGrid::compute_gradient()
{
    if (size() == 0)
        biospring::logging::die("cannot compute gradient of grid which size is 0");

    float scale = biospring::forcefield::BOLTZMANJPERK * biospring::forcefield::METER_TO_ANGSTROM *
                  biospring::forcefield::NEWTON_TO_DALTON_ANGSTROM_PER_FEMTOSECOND_2;

    double dx = cell_size()[0], dy = cell_size()[1], dz = cell_size()[2];

    for (int i = 0; i < shape()[0]; ++i)
    {
        for (int j = 0; j < shape()[1]; ++j)
        {
            for (int k = 0; k < shape()[2]; ++k)
            {
                Vector3f gradient = Vector3f(0.0f, 0.0f, 0.0f);

                // Computes the gradient along the x axis.
                if (i > 0 && i < shape()[0] - 1)
                    gradient[0] = compute_gradient_(_data[i][j][k].scalar, _data[i - 1][j][k].scalar,
                                                    _data[i + 1][j][k].scalar, dx);

                // Computes the gradient along the y axis.
                if (j > 0 && j < shape()[1] - 1)
                    gradient[1] = compute_gradient_(_data[i][j][k].scalar, _data[i][j - 1][k].scalar,
                                                    _data[i][j + 1][k].scalar, dy);

                // Computes the gradient along the z axis.
                if (k > 0 && k < shape()[2] - 1)
                    gradient[2] = compute_gradient_(_data[i][j][k].scalar, _data[i][j][k - 1].scalar,
                                                    _data[i][j][k + 1].scalar, dz);

                // Stores the gradient in the cell.
                _data[i][j][k].vector = gradient * -scale;
            }
        }
    }
}

} // namespace grid
} // namespace biospring