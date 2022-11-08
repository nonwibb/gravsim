// This function has the potential energy evaluations

using namespace std;
#include <vector>
#include "particle.h"
#include<iostream>

double compute_potential_energies(vector<PARTICLE>& particle, double lj_epsilon, double r_g, double r_c, double L)
{
    // what you need to compute energy
    double rij;           // magnitude of distance between particle i and j
    VECTOR3D distance_ij; // distance vector between particle i and j
    double energy_shift = 4 * (1 / pow(r_c, 12) - 1 / pow(r_c, 6));

    // energy is computed as pair-wise sums
    for (unsigned int i = 0; i < particle.size(); i++)
    {
        particle[i].pe = 0.0;
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
                particle[i].pe = particle[i].pe + 4 * (1 / pow(rij, 12) - 1 / pow(rij, 6)) - energy_shift;
            }
            else
                particle[i].pe = particle[i].pe + 0;
            particle[i].pe = particle[i].pe - r_g;
        }
    }

    double total_pe = 0;
    for (unsigned int i = 0; i < particle.size(); i++) {
        total_pe += particle[i].pe;
    }
    total_pe = 0.5 * total_pe;    // factor of 0.5 to account for double counting

    return total_pe;
}