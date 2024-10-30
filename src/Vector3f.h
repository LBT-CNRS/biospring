#ifndef _VECTOR3f_H_
#define _VECTOR3f_H_

#include <array>
#include <cmath>
#include <ostream>
#include <string>

// Three dimensions vector with single floating point precision, used for position, velocity, or force.
//
// Author: Nicolas Ferey

class Vector3f
{
  public:
    Vector3f() : _x(0.0), _y(0.0), _z(0.0) {}
    Vector3f(float x, float y, float z) : _x(x), _y(y), _z(z) {}
    Vector3f(std::array<double, 3> x) : _x(x[0]), _y(x[1]), _z(x[2]) {}

    // ==================================================================================
    // Gets/Sets coordinates.
    float getX() const { return _x; }
    float getY() const { return _y; }
    float getZ() const { return _z; }

    void setX(const float x) { _x = x; }
    void setY(const float y) { _y = y; }
    void setZ(const float z) { _z = z; }

    // Subscript operator.
    float & operator[](int i)
    {
        if (i == 0)
            return _x;
        else if (i == 1)
            return _y;
        else if (i == 2)
            return _z;
        throw std::out_of_range("Vector3f subscript out of range");
    }

    // ==================================================================================

    // Normalizes the vector.
    void normalize()
    {
        float norm = this->norm();
        if (norm > 1.0E-40)
        {
            _x /= norm;
            _y /= norm;
            _z /= norm;
        }
        else
        {
            _x = 0.0;
            _y = 0.0;
            _z = 0.0;
        }
    }

    // Returns the module of this vector.
    float norm() const;

    // Returns the opposite vector.
    Vector3f operator-() const { return Vector3f(-_x, -_y, -_z); }

    // Returns the difference between this vector and v.
    Vector3f operator-(const Vector3f & v) const { return Vector3f(_x - v._x, _y - v._y, _z - v._z); }

    // Returns the addition between this vector and v.
    Vector3f operator+(const Vector3f & v) const { return Vector3f(_x + v._x, _y + v._y, _z + v._z); }

    // Returns the vector scaled with a factor.
    Vector3f operator*(float scale) const { return Vector3f(_x * scale, _y * scale, _z * scale); }

    // Returns the term by term product between this vector and v.
    Vector3f operator*(const Vector3f & v) const { return Vector3f(_x * v._x, _y * v._y, _z * v._z); }

    // Returns the cross product between this vector and v.
    Vector3f operator^(const Vector3f & v) const
    {
        return Vector3f(_y * v._z - _z * v._y, _z * v._x - _x * v._z, _x * v._y - _y * v._x);
    }

    // Returns a new vector down scaled by scale.
    Vector3f operator/(float scale) const { return Vector3f(_x / scale, _y / scale, _z / scale); }

    // Adds v to this vector.
    void operator+=(const Vector3f & v)
    {
        _x += v._x;
        _y += v._y;
        _z += v._z;
    }

    // Subtracts v to this vector.
    void operator-=(const Vector3f & v)
    {
        _x -= v._x;
        _y -= v._y;
        _z -= v._z;
    }

    // Scales this vector with a factor.
    void operator*=(float scale)
    {
        _x *= scale;
        _y *= scale;
        _z *= scale;
    }

    // Term-by-term product between this vector and v.
    void operator*=(const Vector3f & v)
    {
        _x *= v._x;
        _y *= v._y;
        _z *= v._z;
    }

    // Scales this vector with a factor.
    void operator/=(float scale)
    {
        _x /= scale;
        _y /= scale;
        _z /= scale;
    }

    // Sets the components of this vector to 0.0
    void reset()
    {
        _x = 0;
        _y = 0;
        _z = 0;
    }

    // Vector equality.
    bool operator==(const Vector3f & v) const
    {
        return (std::abs(_x - v._x) < 1e-6) && (std::abs(_y - v._y) < 1e-6) && (std::abs(_z - v._z) < 1e-6);
    }

    // Returns the dot product between this vector and v.
    float dot(const Vector3f & v) const { return _x * v._x + _y * v._y + _z * v._z; }

    // Returns the distance between this vector and v
    float distance(const Vector3f & v) const { return distance(*this, v); }

    // Returns the distance between two vectors
    static float distance(const Vector3f & v1, const Vector3f & v2);

    friend std::ostream & operator<<(std::ostream & os, const Vector3f & v)
    {
        return os << "[" << v.getX() << ", " << v.getY() << ", " << v.getZ() << "]";
    }

    std::array<double, 3> to_array() const { return {_x, _y, _z}; }

    void to_array(float arr[3]) const {
        arr[0] = _x;
        arr[1] = _y;
        arr[2] = _z;
    }

    /**Return the components of this vector to string format
	@return the components in string format
	*/
	inline std::string to_string()
	{
		return "x: "+std::to_string(_x)+" y: "+std::to_string(_y)+" z: "+std::to_string(_z);
	}

  private:
    float _x;
    float _y;
    float _z;
};

#endif
