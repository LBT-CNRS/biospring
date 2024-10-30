#ifndef __QUATERNION_H__
#define __QUATERNION_H__

#include <Vector3f.h>
#include <cmath>
#include <iostream>

class Quaternion
{
  public:
    Quaternion(float x, float y, float z, float w) : _x(x), _y(y), _z(z), _w(w) { normalize(); }

    Quaternion() : Quaternion(0.0, 0.0, 0.0, 0.0) {}

    Quaternion(const Quaternion & v) : _x(v._x), _y(v._y), _z(v._z), _w(v._w) { normalize(); }

    void setX(float x) { _x = x; }
    float getX() const { return _x; }

    void setY(float y) { _y = y; }
    float getY() const { return _y; }

    void setZ(float z) { _z = z; }
    float getZ() const { return _z; }

    void setW(float w) { _w = w; }
    float getW() const { return _w; }

    float norm() const { return std::sqrt(_w * _w + _x * _x + _y * _y + _z * _z); }

    void normalize();

    void fromAxisAngle(float x, float y, float z, float angle);
    void fromAxisAngle(Vector3f & axis, float angle) { fromAxisAngle(axis.getX(), axis.getY(), axis.getZ(), angle); }

    void toAxisAngle(Vector3f * axis, float * angle) const;

    void fromEuler(float ax, float ay, float az);

    Quaternion operator*(const Quaternion & v) const;
    void rotationUsing2Vectors(Vector3f & v1, Vector3f & v2);

  private:
    float _x;
    float _y;
    float _z;
    float _w;
};

#endif // __QUATERNION_H__
