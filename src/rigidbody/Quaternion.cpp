#include "Quaternion.h"
#include "logging.h"

namespace logging = biospring::logging;

namespace biospring
{
namespace rigidbody
{

void Quaternion::normalize()
{
    float magnitude = sqrt(_x*_x+_y*_y+_z*_z+_w*_w);
    if (magnitude > 0.005)
    {
        _x /= magnitude;
        _y /= magnitude;
        _z /= magnitude;
        _w /= magnitude;
    }
}

Quaternion Quaternion::inverse() const
{
    float absoluteValue = norm();
    absoluteValue *= absoluteValue;
    absoluteValue = 1 / absoluteValue;

    Quaternion conjugateValue = conjugate();

    // logging::info("absoluteValue: %f", absoluteValue);
    // logging::info("conjugateValue-> x: %f; y: %f; z: %f; w: %f", conjugateValue.getX(), conjugateValue.getY(), conjugateValue.getZ(), conjugateValue.getW());

    float scalar = conjugateValue.getW() * absoluteValue;
    Vector3f imaginary = conjugateValue.getV() * absoluteValue;

    return Quaternion(imaginary, scalar);
}

Quaternion Quaternion::rotateVectorAboutAxisAndAngle(Vector3f axis, float angle) const
{
    //logging::info("axis before -> x: %f; y: %f; z: %f", axis.getX(), axis.getY(), axis.getZ());
    Quaternion qrot = Quaternion();
    qrot.fromAxisAngle(axis.getX(), axis.getY(), axis.getZ(), angle);
    Quaternion qrot_inv = qrot.inverse();
    return qrot * Quaternion(_x, _y, _z, 0.) * qrot_inv;
}


void Quaternion::fromAxisAngle(float x, float y, float z, float angle)
{
    float sin_a = std::sin(angle / 2.0);
    float cos_a = std::cos(angle / 2.0);

    _x = x * sin_a;
    _y = y * sin_a;
    _z = z * sin_a;
    _w = cos_a;
    normalize();
}

void Quaternion::toAxisAngle(Vector3f * axis, float * angle) const
{
    *angle = std::acos(_w) * 2.0;
    *axis = Vector3f(_x, _y, _z);
    axis->normalize();
}

void Quaternion::fromEuler(float ax, float ay, float az)
{
    Vector3f vx = Vector3f(1.0, 0.0, 0.0);
    Vector3f vy = Vector3f(0.0, 1.0, 0.0);
    Vector3f vz = Vector3f(0.0, 0.0, 1.0);

    Quaternion qx;
    qx.fromAxisAngle(vx, ax);
    Quaternion qy;
    qy.fromAxisAngle(vy, ay);
    Quaternion qz;
    qz.fromAxisAngle(vz, az);

    *this = qx * qy * qz;
}

Quaternion Quaternion::operator*(const Quaternion & v) const
{
    Vector3f v2(v._x, v._y, v._z);
    Vector3f v1(_x, _y, _z);

    Vector3f axis1 = v2 * _w;
    Vector3f axis2 = v1 * v._w;
    Vector3f axis3 = v1 ^ v2;
    Vector3f axis = axis1 + axis2 + axis3;

    return Quaternion(axis.getX(), axis.getY(), axis.getZ(), _w * v._w - v1.dot(v2));
}

Quaternion Quaternion::operator+(const Quaternion & v) const
{
    return Quaternion(_x + v.getX(), _y + v.getY(), _z + v.getZ(), _w + v.getW());
}

Quaternion& Quaternion::operator+=(const Quaternion & v)
{
    _x += v.getX();
    _y += v.getY();
    _z += v.getZ();
    _w += v.getW();
    return *this;
}

Quaternion Quaternion::operator*(float d) const
{
    return Quaternion(_x*d, _y*d, _z*d, _w*d);
}

void Quaternion::rotationUsing2Vectors(Vector3f & v1, Vector3f & v2)
{
    Vector3f v3 = v1 ^ v2;
    _x = v3.getX();
    _y = v3.getY();
    _z = v3.getZ();
    float norm1 = v1.norm();
    float norm2 = v2.norm();
    _w = norm1 * norm1 * norm2 * norm2 + v1.dot(v2);
    normalize();
}

Matrix Quaternion::quaternionToRotationMatrix() const
{
    /*
    Covert a quaternion into a full three-dimensional rotation matrix.

    Input
    :param Q: A 4 element array representing the quaternion (q0,q1,q2,q3) 

    Output
    :return: A 3x3 element matrix representing the full 3f rotation matrix. 
                This rotation matrix converts a point in the local reference 
                frame to a point in the global reference frame.
    */
    Matrix m = Matrix(3,3);

    // // First row of the rotation matrix
    // m(0,0) = 2 * (_x * _x + _y * _y) - 1;
    // m(0,1) = 2 * (_y * _z - _x * _w);
    // m(0,2) = 2 * (_y * _w + _x * _z);
        
    // // Second row of the rotation matrix
    // m(1,0) = 2 * (_y * _z + _x * _w);
    // m(1,1) = 2 * (_x * _x + _z * _z) - 1;
    // m(1,2) = 2 * (_z * _w - _x * _y);
        
    // // Third row of the rotation matrix
    // m(2,0) = 2 * (_y * _w - _x * _z);
    // m(2,1) = 2 * (_z * _w + _x * _y);
    // m(2,2) = 2 * (_x * _x + _w * _w) - 1;

    // First row of the rotation matrix
    m(0,0) = 1 - 2 * (_y * _y - _z * _z);
    m(0,1) = 2 * (_x * _y - _w * _z);
    m(0,2) = 2 * (_x * _z + _w * _y);
        
    // Second row of the rotation matrix
    m(1,0) = 2 * (_x * _y + _w * _z);
    m(1,1) = 1 - 2 * (_x * _x - _z * _z);
    m(1,2) = 2 * (_y * _z - _w * _x);
        
    // Third row of the rotation matrix
    m(2,0) = 2 * (_x * _z - _w * _y);
    m(2,1) = 2 * (_y * _z + _w * _x);
    m(2,2) = 1 - 2 * (_x * _x - _y * _y);

    return m;
}

Matrix Quaternion::quaternionToMatrix() const
{
    Matrix m = Matrix(1,4);
    m(0,0) = _x;
    m(0,1) = _y;
    m(0,2) = _z;
    m(0,3) = _w;

    return m;
}

    

} // namespace rigidbody    
} // namespace biospring