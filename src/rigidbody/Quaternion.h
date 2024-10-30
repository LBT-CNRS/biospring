#ifndef __QUATERNION_H__
#define __QUATERNION_H__

#include <Vector3f.h>
#include <cmath>
#include <iostream>
#include "matrix.h"

namespace biospring
{
namespace rigidbody
{

class Quaternion
{
  public:
    // Quaternion(float x, float y, float z, float w) : _x(x), _y(y), _z(z), _w(w) { normalize(); }
    Quaternion(float x, float y, float z, float w) : _x(x), _y(y), _z(z), _w(w) {  }

    Quaternion() : Quaternion(0.0, 0.0, 0.0, 0.0) {}

    // Quaternion(const Quaternion & v) : _x(v._x), _y(v._y), _z(v._z), _w(v._w) { normalize(); }
    Quaternion(const Quaternion & v) : _x(v._x), _y(v._y), _z(v._z), _w(v._w) {  }

    // Quaternion(Vector3f v, float w) : _x(v.getX()), _y(v.getY()), _z(v.getZ()), _w(w) { normalize(); }
    Quaternion(Vector3f v, float w) : _x(v.getX()), _y(v.getY()), _z(v.getZ()), _w(w) {  }

    void setX(float x) { _x = x; }
    float getX() const { return _x; }

    void setY(float y) { _y = y; }
    float getY() const { return _y; }

    void setZ(float z) { _z = z; }
    float getZ() const { return _z; }

    void setW(float w) { _w = w; }
    float getW() const { return _w; }

    Vector3f getV() const { return Vector3f(_x, _y, _z); }

    /**Return the components of this quaternion to string format
		@return the components in string format
		*/
		inline std::string to_string() 
    {
    return "x: "+std::to_string(_x)+" y: "+std::to_string(_y)+" z: "+std::to_string(_z)+" w: "+std::to_string(_w);
    }

    inline int sign() const { return (_w >= 0) ? 1 : -1; }


    float norm() const { return std::sqrt(_w * _w + _x * _x + _y * _y + _z * _z); }

    Quaternion conjugate() const { return Quaternion(-_x, -_y, -_z, _w); }

    Quaternion inverse() const;

    Quaternion rotateVectorAboutAxisAndAngle(Vector3f axis, float angle) const;

    void normalize();

    void fromAxisAngle(float x, float y, float z, float angle);
    void fromAxisAngle(Vector3f & axis, float angle) { fromAxisAngle(axis.getX(), axis.getY(), axis.getZ(), angle); }

    void toAxisAngle(Vector3f * axis, float * angle) const;

    void fromEuler(float ax, float ay, float az);

    Quaternion operator*(const Quaternion & v) const;
    Quaternion operator+(const Quaternion & v) const;
    Quaternion& operator+=(const Quaternion & v);
    Quaternion operator*(float d) const;
    void rotationUsing2Vectors(Vector3f & v1, Vector3f & v2);

    Matrix quaternionToRotationMatrix() const;
    Matrix quaternionToMatrix() const;

  private:
    float _x;
    float _y;
    float _z;
    float _w;
};

} // namespace rigidbody
} // namespace biospring

#endif // __QUATERNION_H__
