// This file contains functions that are useful to:
// initialize many particle systems
// monitor and analyze the dynamics of many particle systems

#include "functions.h"
#include <random>
#include <chrono>

// make movie
void make_movie(int num, vector<PARTICLE>& particle, ofstream& outdump, double L)
{
  outdump << "ITEM: TIMESTEP" << endl;
  outdump << num << endl;
  outdump << "ITEM: NUMBER OF ATOMS" << endl;
  outdump << particle.size() << endl;
  outdump << "ITEM: BOX BOUNDS" << endl;
  outdump << -L/2 << "\t" << L/2 << endl;
  outdump << -L/2 << "\t" << L/2 << endl;
  outdump << -L/2 << "\t" << L/2 << endl;
  outdump << "ITEM: ATOMS index type x y z v f" << endl;
  for (unsigned int i = 0; i < particle.size(); i++)
  {
    outdump << i+1 << "   " << particle[i].ty << "   " << particle[i].position.x << "   " << particle[i].position.y << "   " << particle[i].position.z << "   " << particle[i].velocity.Magnitude() << "   " << particle[i].force.Magnitude() << endl;
  }
}

// display progress
void ProgressBar (double progress)
{
    int val = (int) (progress * 100);
    int lpad = (int) (progress * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf ("\r%3d%% |%.*s%*s|", val, lpad, PBSTR, rpad, "");
    fflush (stdout);
}

// initialize velocities of particles to start simulation
void initialize_particle_velocities(vector<PARTICLE>& particle, double initial_temperature)
{
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator (seed);
    std::normal_distribution<double> distribution (0.0,initial_temperature);
    for (unsigned int i = 0; i < particle.size(); i++) {
        double vx = distribution(generator);
        double vy = distribution(generator);
        double vz = distribution(generator);
        particle[i].velocity = VECTOR3D(vx, vy, vz);
    }

    // average velocity should be 0; as there is no net flow of the system in any particular direction; we do this next
    VECTOR3D average_velocity_vector = VECTOR3D(0,0,0);
    for (unsigned int i = 0; i < particle.size(); i++)
        average_velocity_vector = average_velocity_vector + particle[i].velocity;
    average_velocity_vector = average_velocity_vector*(1.0/particle.size());

    // subtract this computed average_velocity_vector from the velocity of each particle to ensure that the total average after this operation is 0
    for (unsigned int i = 0; i < particle.size(); i++)
        particle[i].velocity = particle[i].velocity - average_velocity_vector;
}

void initialize_particle_positions(vector<PARTICLE>& particle, double density, int N, double sigma, double m, double L)
{
    // initialize the positions of the particles on a lattice (not randomly)

    unsigned int N_axis = int(ceil(pow((double(N)),1.0/3.0)));

    for (unsigned int i = 0; i < N_axis; i++)
    {
        for (unsigned int j = 0; j < N_axis; j++)
        {
            for (unsigned int k = 0; k < N_axis; k++)
            {
                if (int(particle.size()) < N)	// stop if atoms created are equal to the requested number
                {
                    double x = (-L/2 + 1.0/2) + i;
                    double y = (-L/2 + 1.0/2) + j;
                    double z = (-L/2 + 1.0/2) + k;
                    VECTOR3D position = VECTOR3D(x,y,z);
                    PARTICLE fresh_particle = PARTICLE(1,sigma,m,position,VECTOR3D(0,0,0));
                    particle.push_back(fresh_particle);
                }
            }
        }
    }
}