#include "Quaternion.h"

void Quaternion::normalize()
{
    float magnitude = norm();
    if (magnitude > 0.005)
    {
        _x /= magnitude;
        _y /= magnitude;
        _z /= magnitude;
        _w /= magnitude;
    }
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

    return Quaternion(axis.getX(), axis.getY(), axis.getZ(), _w * v._w - v1 * v2);
}

void Quaternion::rotationUsing2Vectors(Vector3f & v1, Vector3f & v2)
{
    Vector3f v3 = v1 ^ v2;
    _x = v3.getX();
    _y = v3.getY();
    _z = v3.getZ();
    float norm1 = v1.norm();
    float norm2 = v2.norm();
    _w = norm1 * norm1 * norm2 * norm2 + v1 * v2;
    normalize();
}
