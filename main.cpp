// This is main.
// This is MD simulation of many particles interacting with pair potentials

using namespace std;

#include "functions.h"

void compute_forces(vector<PARTICLE>&, double, double, double, double);
double compute_potential_energies(vector<PARTICLE>&, double, double, double, double);

int main(int argc, char* argv[]) 
{
  // we begin with defining the Boltzmann constant
  const double kB = 1.38e-23;	// Joules per Kelvin

  cout << "\n---Simulating fluid (Argon)---\n";

  // key descriptors of the many particle system
  double diameter = 3.405e-10;  // what is the size of each particle? (in meters)
  double mass = 6.634e-26;  // what is the mass of each particle? (in kg)
  double ljenergy = 120 * kB; // what is the strength of typical particle-particle collisions? (in Joules)
  double gacceleration = -9.81; // (m/s/s)

  // reduced units
  double unitlength = diameter; // unit of length is the simulation (diameter of the particle)
  double unitmass = mass; // unit of mass in the simulation (mass of the particle)
  double unitenergy = ljenergy; // unit of energy in the simulation (characteristic pair potential strength)
  double unitgravity = mass*gacceleration;
  double unittemperature = ljenergy/kB;

  double unittime = sqrt(unitmass * unitlength * unitlength / unitenergy); // Unit of time

  cout << "unit of length is " << unitlength << " meters" << endl;
  cout << "unit of mass is " << unitmass << " kilograms" << endl;
  cout << "unit of energy is " << unitenergy << " Joules" << endl;
  cout << "unit of time is " << unittime << " seconds" << endl;
  cout << "unit of temperature is " << unittemperature << " Kelvin" << endl;
  cout << "non-reduced gravity on particle is " << unitgravity << " Newtons" << endl;
  cout << "\n";


  double reduced_gravity = unitgravity / ljenergy * mass;
  double reduced_diameter = diameter / unitlength;
  double reduced_mass = mass / unitmass;
  double reduced_ljenergy = ljenergy / unitenergy;
  cout << "reduced gravity on particle is " << reduced_gravity << endl;
  int Nparticles = 216;
  cout << "initializing 216 particles on a lattice (simple cubic). density will be changed by changing box size" << endl;

  double density;
  cout << "enter density (in reduced units; suggested 0.8442) " << endl;
  cin >> density;


  double initial_temperature = 0.728;   // in reduced units

  double distance_cutoff = 2.5; // in units of the diameter (chosen as unitlength)

  double boxL = pow(Nparticles/density,1.0/3.0); // box length
  cout << "cubic box edge length (in reduced units) calculated using " << Nparticles << " particles and " << density << " density: " << boxL << endl;

  if (boxL <= 2*distance_cutoff)
  {
      cout << "the box edge length is not large enough to make PBC work" << endl;
      cout << "not proceeding further. increase number of particles to resolve the issue. choose a cubic number." << endl;
      return 0;
  }

  // Different parts of the system
  vector<PARTICLE> particle;		// all particles in the system

  initialize_particle_positions(particle, density, Nparticles, reduced_diameter, reduced_mass, boxL);

  initialize_particle_velocities(particle, initial_temperature);

  // output to screen
  cout << "particle diameter in meters: " << diameter << " and in reduced units: " << reduced_diameter << endl;
  cout << "particle mass in kg: " << mass << " and in reduced units: " << reduced_mass << endl;
  cout << "pair potential energy strength in Joules: " << ljenergy << " and in reduced units: " << reduced_ljenergy << endl;
  cout << "Number of particles inside the box: " << particle.size() << endl;
  cout << "density of the fluid is " << particle.size()/(boxL*boxL*boxL) << endl;
  cout << "volume packing fraction of the fluid is " << (M_PI/6.0)*particle.size()/(boxL*boxL*boxL);
  cout << "\n";

  // initial energies and forces computation
  double totalke = 0.0;
  for (unsigned int i = 0; i < particle.size(); i++)
  {
      particle[i].kinetic_energy();
      totalke += particle[i].ke;
  }

  double totalpe = compute_potential_energies(particle, reduced_ljenergy, reduced_gravity, distance_cutoff, boxL);
  compute_forces(particle, reduced_ljenergy, reduced_gravity, distance_cutoff, boxL);

  // create files for storing movie and energies

  char file_movie[200], file_energy[200];

  // file movie
  sprintf(file_movie, "movie_rho%f_N%d.out", density, Nparticles);
  ofstream list_propagation(file_movie, ios::out); // create a file to store and visualize 3D data
  make_movie(0,particle,list_propagation, boxL);

  // file ke, pe, and total energies of all particles
  sprintf(file_energy, "energy_rho%f_N%d.out", density, Nparticles);
  ofstream output_energy(file_energy, ios::out);
  output_energy << 0 << "  " << totalke/particle.size() << "  " << totalpe/particle.size() << "  " << (totalke+totalpe)/particle.size() << endl;

  // print energies
  cout << "initial kinetic energy per particle (in reduced units): " << totalke/particle.size() << endl;
  cout << "initial potential energy per particle (in reduced units): " << totalpe/particle.size() << endl;
  cout << "initial total system energy per particle (in reduced units): " << (totalke+totalpe)/particle.size() << endl;
  cout << "\n";

  double totaltime = 50;
  int steps = 50000;		// number of time discretizations (slices)
  double delta_t = totaltime/steps;	// choose steps carefully
  int movie_step = 100;  // movie will be made every movie_step
  int energycalc_step = 100; // energies will be computed every energycalc_step

  cout << "Fluid simulation time (in seconds): " << totaltime*unittime << " and in reduced units: " << totaltime << endl;

  // code to setup computing averages from simulations
  int hit_eqm = 1000; // this is your choice of where you think the system hit equilibrium
  double average_pe = 0.0;
  double average_ke = 0.0;
  int samples = 0;

  // Molecular Dynamics
  cout << "progress..." << endl;
  for (int num = 1; num <= steps; num++)
  {
      // velocity-Verlet
      for (unsigned int i = 0; i < particle.size(); i++) {
          particle[i].update_velocity(delta_t/2);  // update velocity half timestep
      }

      for (unsigned int i = 0; i < particle.size(); i++) {
          particle[i].update_position(delta_t, boxL);  // update position full timestep
      }

      compute_forces(particle, reduced_ljenergy, reduced_gravity, distance_cutoff, boxL);  // expensive step

      for (unsigned int i = 0; i < particle.size(); i++) {
          particle[i].update_velocity(delta_t/2);  // update velocity half timestep
      }

      // calling a movie function to get a movie of the simulation every movie_step
      if (num%movie_step == 0)
          make_movie(num,particle,list_propagation,boxL);

      // calculating energies every energycalc_step
      if (num%energycalc_step == 0) {
          totalke = 0.0;
          for (unsigned int i = 0; i < particle.size(); i++) {
              particle[i].kinetic_energy();
              totalke += particle[i].ke;
          }
          totalpe = compute_potential_energies(particle, reduced_ljenergy, reduced_gravity, distance_cutoff, boxL);

          // outputting the energy (per particle) to make sure simulation can be trusted
          output_energy << num << "  " << totalke/particle.size() << "  " << totalpe/particle.size() << "  " << (totalke+totalpe)/particle.size() << endl;

          // if equilibrium (steady-state) is reached, we measure the average equilibrium pe and ke by utilizing the generated samples (totalpe, total ke).
          // other equilibrium properties can be measured similarly
          if (num > hit_eqm)
          {
              average_pe = average_pe + totalpe;
              average_ke = average_ke + totalke;
              samples++;
          }
      }

      // monitoring progress of the simulation
      double progress = ((num)/(double)steps);
      ProgressBar(progress);
  }

  cout << endl;
  cout << "equilibrium properties: ..." << endl;
  cout << "average equilibrium pe per particle is " << (average_pe/samples/particle.size()) << endl;
  cout << "average equilibrium ke per particle is " << (average_ke/samples/particle.size()) << endl;

  double equilibrium_temperature = (2.0/3.0)*(average_ke/samples)/(particle.size());
  cout << "average equilibrium temperature is " << equilibrium_temperature << endl;

  return 0;
} 
// End of main