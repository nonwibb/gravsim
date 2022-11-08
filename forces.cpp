// This file contains the routine that computes the force on the particle exerted by all other particles

using namespace std;

#include <vector>
#include "particle.h"
#include<iostream>

void compute_forces(vector<PARTICLE>& particle, double lj_epsilon, double r_g, double r_c, double L)
{
    double rij;           // magnitude of distance between particle i and j
    VECTOR3D distance_ij; // distance vector between particle i and j
    VECTOR3D Fij;         // Force vector ON i BY J

    for (unsigned int i = 0; i < particle.size(); i++)
    {
        Fij = VECTOR3D(0,0,0);
        for (unsigned int j = 0; j < particle.size(); j++)
        {
            if (j == i) continue;
            distance_ij = particle[i].position-particle[j].position;
            
            if (distance_ij.x > L/2)
                distance_ij.x = distance_ij.x - L;
            if (distance_ij.x < -L/2)
                distance_ij.x = distance_ij.x + L;
            if (distance_ij.y > L/2)
                distance_ij.y = distance_ij.y - L;
            if (distance_ij.y < -L/2)
                distance_ij.y = distance_ij.y + L;

            rij = distance_ij.Magnitude();
            if (rij <= r_c) {
                Fij = Fij + (distance_ij / rij) * (48 * (1 / rij) * ((1 / pow(rij, 12)) - 0.5 * (1 / pow(rij, 6))));
            }
            else
            {
                Fij = Fij + VECTOR3D(0, 0, 0);
            }
        }
        particle[i].force =  Fij - VECTOR3D(0,0,r_g);
    }
}