// This is the particle class
// It lists features of the particle as members and also lists how its position, velocity, etc are updated

#ifndef _PARTICLE_H
#define _PARTICLE_H

#include <iostream>
#include "vector3d.h"

class PARTICLE 
{
  public:

  // members
  int ty;           // type of the particle
  double diameter;	// diameter of the particle
  double m; 		// mass of the particle
  VECTOR3D position;	// position vector of the particle
  VECTOR3D velocity;	// velocity vector of the particle
  VECTOR3D force;	// force vector on the particle
  double pe;		// potential energy
  double ke;	    // kinetic energy
  
  // member functions
  
  // make a particle
  PARTICLE(int initial_type = 1, double initial_diameter = 0, double initial_mass = 0, VECTOR3D initial_position = VECTOR3D(0,0,0), VECTOR3D initial_velocity = VECTOR3D(0,0,0))
  {
    ty = initial_type;
    diameter = initial_diameter;
    m = initial_mass;
    position = initial_position;
    velocity = initial_velocity;
  }
  
  // the next two functions are central to the velocity-Verlet algorithm:

  // update position of the particle
  void update_position(double dt, double L)		// dt is the time-step, L is the box length
  {
    position = ( position + (velocity * dt) );	// position updated to a full time-step
      if (position.x > L/2.0)
          position.x = position.x - L;
      if (position.x < -L/2.0)
          position.x = position.x + L;
      if (position.y > L/2)
          position.y = position.y - L;
      if (position.y < -L/2.0)
          position.y = position.y + L;
      if (position.z > L/2.0)
          velocity.z = -velocity.z;
      if (position.z < -L/2.0)
          velocity.z = -velocity.z;
  }
  
  // update velocity of the particle
  void update_velocity(double half_dt)
  {
    velocity = ( velocity + ( (force) * ( half_dt / m ) ) );	// notice the half time-step
  }
  
  // calculate kinetic energy of a particle
  void kinetic_energy()				
  {
    ke = 0.5 * m * velocity.Magnitude() * velocity.Magnitude();
  }
};

#endif