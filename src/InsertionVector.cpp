#include "InsertionVector.h"
#include "SpringNetwork.h"
#include "logging.h"
#include "measure.hpp"

#include <cmath>


biospring::spn::Particle & InsertionVector::getParticle(unsigned i) 
{
    if (i > 1)
        throw std::out_of_range("i should be 0 or 1");
    if (i == 0)
        return _p1;
    return _p2;
}

void InsertionVector::computeVector()
{
    _vector = _p2.getPosition() - _p1.getPosition();
}

void InsertionVector::computeAngle()
{
    _angle = std::acos(-_vector.getZ() / _vector.norm());
    _angle *= 180.0 / M_PI;
    _angle -= 90.0;
}

void InsertionVector::computeRollAngle()
{
    // FOR NOW ROLL ANGLE ONLY CORRECT FOR PLANE MEMBRANE
    // if (!_sp->getForceField()->isMeshMembraneVerticesArrayEmpty())
    // {
    //     // Need to be implemented...
    // } else {
    // }

    // Center of Mass
    Vector3f bary = _sp.getCentroid();
    // Position of reference particle 1
    Vector3f p1Pos = _p1.getPosition();
    // Position of reference particle 2
    Vector3f p2Pos = _p2.getPosition();
    // Vector betwenn p1 and p2 : insertion vector
    Vector3f p2p1 = p2Pos - p1Pos;
    // Projection of barycentre on insertion vector
    Vector3f _projOnIVPos = p1Pos + p2p1*(((bary - p1Pos) * p2p1)/pow(p2p1.norm(), 2.0));

    // First vector to define the roll rotation circle
    Vector3f n1 = (_projOnIVPos - bary);
    n1.normalize();
    // Second vector to define the roll rotation circle
    Vector3f n3 = (n1 ^ p2p1);
    n3.normalize();
    // Radian angle to find the position on the roll rotation circle that has the minimum z component.
    float am = std::atan2(n3.getZ(), n1.getZ()) - M_PI;
    // Radius of the roll rotation circle
    float r = Vector3f::distance(bary, _projOnIVPos);
    // Position on the roll rotation circle that has the minimum z component.
    Vector3f _refRollRotPos = _projOnIVPos + n1*r*std::cos(am) + n3*r*std::sin(am);
    // Reference vector to compute angle
    Vector3f n2 = (_projOnIVPos - _refRollRotPos);
    n2.normalize();
    // Compute angle in radian between n1 and n2 : between 0 and PI
    float a_rad = std::acos((n1.dot(n2))/(n1.norm()*n2.norm()));
    // Determine the z component of n1n2 normal vector to get the side of bary on circle
    // to allow calculate roll angle between 0 and 360°.
    float s = (n1 ^ n2).dot(p2Pos - _projOnIVPos);
    // if s > 0.0 : 0 to 180 else: 360 - (180 to 0)
    float a_deg_1 = a_rad * 180/M_PI;
    float a_deg_final = a_deg_1;
    if (s > 0.0)
    {
        a_deg_final = 360 - a_deg_1;
    }

    _rollAngle = a_deg_final;
}


float InsertionVector::getInsertionDepth() const
{
    return _sp.getCentroid()[2];
}
