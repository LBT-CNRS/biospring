#include "Vector3f.h"
#include "measure.hpp"

float Vector3f::distance(const Vector3f & v1, const Vector3f & v2)
{
    return biospring::measure::distance(v1, v2);
}

