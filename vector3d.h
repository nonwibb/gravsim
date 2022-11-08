// This is 3d vector class
// VECTOR3D is a vector in cartesian coordinate system

#ifndef _VECTOR3D_H
#define _VECTOR3D_H

#include <cmath>

class VECTOR3D 
{
public:
    double x, y, z;
    
    // make a 3d vector, this is called the constructor function
    // note: it has the same name as the class
    VECTOR3D(double initial_x = 0.0, double initial_y = 0.0, double initial_z = 0.0)
    {
      x = initial_x;
      y = initial_y;
      z = initial_z;
    }
    
    // make any member functions you need
    double Magnitude()
    {
      double r;
      r = sqrt(x*x + y*y + z*z);
      return r;
    }

    // add two 3d vectors
    VECTOR3D operator+(const VECTOR3D& v)  // v is the second vector you want to add to first
    {
      return (VECTOR3D(x + v.x, y + v.y, z + v.z));
    }
    
    // subtract two 3d vectors
    VECTOR3D operator-(const VECTOR3D& vec)
    {
      return VECTOR3D(x - vec.x, y - vec.y, z - vec.z);
    }

    // multiply a 3d vector by a scalar
    // note: scalar has to follow vector so vector*scalar when you call, NOT scalar*vec.
    VECTOR3D operator*(double scalar)
    {
        return (VECTOR3D(x*scalar, y*scalar, z*scalar));
    }

    // divide a 3d vector by a scalar
    // note: scalar has to follow a 3d vector so 3dvector/scalar when you call, NOT scalar/3dvector (this does not make any sense, even mathematically).
    VECTOR3D operator/(double scalar)
    {
        return (VECTOR3D(x/scalar, y/scalar, z/scalar));
    }
};

#endif