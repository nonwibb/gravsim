#ifndef _FUNCTIONS_H
#define _FUNCTIONS_H

#include<iostream>
#include<iomanip>
#include<fstream>
#include<vector>
#include<cmath>
#include "vector3d.h"
#include "particle.h"

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

using namespace std;

// make movie
void make_movie(int, vector<PARTICLE>&, ofstream&, double);

// initialize particle velocities
void initialize_particle_velocities(vector<PARTICLE>&, double);

// initialize particle positions
void initialize_particle_positions(vector<PARTICLE>&, double, int, double, double, double);

// display progress of the simulation
void ProgressBar(double);

#endif