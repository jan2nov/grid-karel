//**********************************************//
// Code developed by Jan & Karel in Oxford 2013 //
//**********************************************//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>

using namespace std;


// some startup constant
//
const int N_shell = 20;
const int N_grids = 4;
const int N_time_steps = 5;
const float Ratio_mass_core_ambient = .05;
const float Ratio_InEnergy_core_ambient = 50.0;
const float Pi = 3.14159265359;
const float Gamma = 5./3.;
const float CFLfactor = 0.1;
const int T_max = 50;
const float dfactor = 2.25;
const float Art_viscosity = 2.0;
// ********************************************

// false or true debug mode
//
const bool DEBUG_MODE = true;
//const bool DEBUG_MODE = false;
// ********************************************

// just a function to call a variable (array) and print it
//
void debug_print(const double *variable,
				 const char *name,
				 const int count){	
	for (int i = 0; i < count; i++)	{
		printf("%i. Your number #%s is: %lf\n", i, name, variable[i]);
	}
	printf("\n");
}
// *****************************************************

// function to set all initial variable
//
void initial_condition(double *f_rho, 
					   double *f_dmass12, 
					   double *f_energy, 
					   double *f_pressure,
					   double *f_velocity,
					   double *f_mass, 
					   double *f_dmass, 
					   double *f_radius,
					   double *f_volume, 
					   double *f_radius12,
					   double *f_A,
					   double *f_A12, 
					   double *f_Ak12,
					   double *f_viscosity){
	
	// initial condition for the core 
	double rho_core = (pow((double)(N_shell / N_grids),3.0) - 1.0) * Ratio_mass_core_ambient;
	double delta_mass12_core = 4.0/3.0 * Pi * pow((double)N_grids/(double)N_shell,3.0) * rho_core/(double)N_grids;
	
	// the loop is going inside the core ...
	for (int i = 0; i < N_grids; i++){
		f_rho[i]      = rho_core;
		f_dmass12[i]  = delta_mass12_core;
		f_energy[i]   = Ratio_InEnergy_core_ambient;
        f_pressure[i] = (Gamma - 1.) * f_rho[i] * f_energy[i];
        f_velocity[i] = 0.0;
	}

	// initial condition for the ambient part
	double rho_amb = 1.0;
	double delta_mass12_ambient = 4.0/3.0 * Pi * (1.0 - pow((double)N_grids / (double)N_shell, 3.0)) / (double)(N_shell - N_grids) * rho_amb;

	// loop begin at the end of the core ...
	for (int i = N_grids; i < N_shell; i++){
		f_rho[i]      = rho_amb;
		f_dmass12[i]  = delta_mass12_ambient;
		f_energy[i]   = 1.0;
        f_pressure[i] = (Gamma - 1.0) * f_rho[i] * f_energy[i];
        f_velocity[i] = 0.0;
	}
	
	// f_dmass no need to have number on position 0
	// not in cartesian coordinates, 
	f_mass[0] = f_dmass12[0];
	f_dmass[0] = 0.00000001;
	f_radius[0] = 0.0;
	f_volume[0] = 0.0;
	f_radius12[0] = 0.00000000001;
	f_A[0] = 0.000001;
	f_A12[0] = 0.0;
	f_viscosity[0] = 0.0;
	for (int i = 1; i < N_shell; i++){
		f_mass[i]   = f_mass[i - 1] + f_dmass12[i];
		f_dmass[i]  = 0.5 * (f_dmass12[i] + f_dmass12[i - 1]);
		f_volume[i] = f_volume[i - 1] + f_dmass12[i - 1] / f_rho[i - 1];
		f_radius[i] = pow(f_volume[i] / (4.0/3.0 * Pi), (1.0/3.0));
		f_A[i]      = 4 * Pi * pow(f_radius[i],2);
		f_A12[i]    = 0.0;
		f_viscosity[i] = 0.0;
	}

	for (int i = 0; i < N_shell - 1; i++){
		f_radius12[i] = f_radius[i + 1] - f_radius[i];
	}

	
	// for checking the variables
	if (DEBUG_MODE){
		debug_print(f_rho, "rho", N_shell);
		debug_print(f_dmass12, "mass12", N_shell);
		debug_print(f_energy, "energy", N_shell);
		debug_print(f_pressure, "pres", N_shell);
		debug_print(f_velocity, "vel", N_shell);
		debug_print(f_mass, "f_mass", N_shell);
		debug_print(f_dmass, "f_dmass", N_shell);
		debug_print(f_radius, "f_radius", N_shell);
		debug_print(f_volume, "f_volume", N_shell);
		debug_print(f_radius12, "f_radius12", N_shell - 1);
		debug_print(f_A, "f_A", N_shell);
	}	
}
// *******************************************


// setting the time scale for each step you need to evaluate
// need to take care of sound in medium ie sound velocity and
// no time in respect to changing velocity and radius and volume 
// of the star.
//
void time_scale_step(double *f_radius12,
					 double *f_velocity,
					 double *f_energy,
					 double *f_volume,
					 double *f_A12, 
					 double f_time, 
					 double *f_dtime){

	double delta_tc = 1.0e30;
	
	for (int i = 0; i < N_shell - 1; i++)
	 delta_tc = min(delta_tc, f_radius12[i] / (abs(f_velocity[i]) + sqrt( (Gamma - 1) * f_energy[i])));
	
	delta_tc = delta_tc * CFLfactor;
	if (f_time + delta_tc > T_max ) delta_tc = T_max - f_time;
	
	// diffusion limit...speed of sound
	double delta_td = 1.0e-30;
	for (int i=0; i< N_shell - 1; i++)
	 delta_td = max(delta_td, abs(f_A12[i + 1] * f_velocity[i + 1] - f_A12[i] * f_velocity[i]) / (f_volume[i + 1] - f_volume[i]));

	delta_td = 0.5 / delta_td / dfactor;

	//check what is minimum if speed of shell or sound
	delta_tc = min(delta_tc,delta_td);

	f_dtime[0] = 0.5 * (f_dtime[1] + delta_tc);
    f_dtime[1] = delta_tc;
	if (DEBUG_MODE) printf("Time: %.8lf Time12: %.10lf \n",f_dtime[0], f_dtime[1]);
}
// **************************************************************


//updating function of all variables
//
void update_step(double *f_velocity,
				 double *f_A,
				 double *f_A12,
				 double *f_Ak12,
				 double *f_pressure,
				 double *f_rho,
				 double *f_energy,
				 double *f_dtime,
				 double *f_dmass,
				 double *f_dmass12,
				 double *f_volume,
				 double *f_viscosity, 
				 double *f_radius,
				 double *f_radius12){

	double old_radius[N_shell];
	double old_energy[N_shell];
	
	//in no cartesian set the v(0) to zero
	f_velocity[0] = 0.0;
	// safe the radius before update; swap process?
	// velocity update
	for(int i = 1; i < N_shell; i++){
		f_velocity[i] = f_velocity[i] - f_A[i] * f_dtime[0] / f_dmass[i] * (f_pressure[i] - f_pressure[i - 1]) 
						- 0.5 * ( f_viscosity[i] * (3.0 * f_Ak12[i] - f_A[i]) 
								- f_viscosity[i - 1] * (3.0 * f_Ak12[i - 1] - f_A[i])) * f_dtime[0] / f_dmass[i] ;
		//printf("%i %lf = %.8lf - %.8lf * %.8lf - %.8lf * %.8lf / %.8lf \n", i, f_velocity[i], f_velocity[i], f_A[i], f_pressure[i], f_pressure[i-1], f_dtime[0], f_dmass[i]);
	}

	//radius update, volume, rho, a, a12
	for (int i = 0; i < N_shell; i++){
		old_radius[i] = f_radius[i];
		f_radius[i]   = old_radius[i] + f_velocity[i] * f_dtime[1];
		f_A12[i]      = 4.0 * Pi * pow(0.5 * ( f_radius[i] + old_radius[i]), 2 );
        f_A[i]        = 4.0 * Pi * pow(f_radius[i],2); 
        f_volume[i]   = 4.0/3.0 * Pi * pow(f_radius[i], 3);
		//printf("A: %lf vol: %lf\n", f_A12[i], f_volume[i]);
	}
	for (int i = 0; i < N_shell - 1; i++){
		f_radius12[i] = f_radius[i + 1] - f_radius[i];
		f_Ak12[i]     = 0.5 * (f_A12[i + 1] + f_A12[i]);
		f_radius12[i] = f_radius[i + 1] - f_radius[i];
		f_rho[i]      = f_dmass12[i] / (f_volume[i + 1] - f_volume[i]);
		f_viscosity[i]= - Art_viscosity * Art_viscosity * f_rho[i] * fabs(f_velocity[i + 1] - f_velocity[i]) 
						* ( f_velocity[i + 1] * (1.0 - f_A12[i + 1] / 3.0 / f_Ak12[i]) 
						  - f_velocity[i] * (1.0 - f_A12[i] / 3.0 / f_Ak12[i]));
		if (f_velocity[i + 1] > f_velocity[i]) f_viscosity[i] = 0.0;
	}
	// boundary of the density?
	f_rho[N_shell] = f_rho[N_shell - 1];
	
	// internal energies and pressure
	for (int i = 0; i < N_shell - 1; i++){
		old_energy[i] = f_energy[i] - f_pressure[i] * (f_A12[i + 1] * f_velocity[i + 1] - f_A12[i] * f_velocity[i]) * f_dtime[1] / f_dmass12[i];
		f_pressure[i] = 0.5 * (f_pressure[i] + (Gamma - 1.0) * f_rho[i] * old_energy[i]);
	}

	
}
//***************************************************************

// main function of hydrocode to evaluate a star collapse
//
int main(int argc, char *argv[]){
	
	FILE *output;

	double *rho, *pressure, *energy, *velocity, *viscosity;
	double *mass, *dmass, *dmass12;
	double *radius, *radius12, *A, *A12, *Ak12, *volume;
	double time = 0.0, *dtime;

	//allocation of variable memory for initial condition and later use
	 rho      = (double*)malloc(sizeof(double)*(N_shell)); 
	 dmass12  = (double*)malloc(sizeof(double)*(N_shell));
	 dmass    = (double*)malloc(sizeof(double)*(N_shell));
	 pressure = (double*)malloc(sizeof(double)*(N_shell)); 
	 energy   = (double*)malloc(sizeof(double)*(N_shell)); 
	 velocity = (double*)malloc(sizeof(double)*(N_shell)); 
	 mass     = (double*)malloc(sizeof(double)*(N_shell));
	 radius   = (double*)malloc(sizeof(double)*(N_shell));
	 volume   = (double*)malloc(sizeof(double)*(N_shell));
	 viscosity= (double*)malloc(sizeof(double)*(N_shell));
	 radius12 = (double*)malloc(sizeof(double)*(N_shell));
	 A        = (double*)malloc(sizeof(double)*(N_shell));
	 A12	  = (double*)malloc(sizeof(double)*(N_shell));
	 Ak12	  = (double*)malloc(sizeof(double)*(N_shell));
	 dtime    = (double*)malloc(sizeof(double)*(2));

	// setup the initial data, like number of grid in the ejecta (core) and number of shells
	initial_condition(rho, dmass12, energy, pressure, velocity, mass, dmass, radius, volume, radius12, A, A12, Ak12, viscosity);

	// set the time pointer to zero's
	dtime[0] = 0.0;
	dtime[1] = 0.0;

	//begin the loop over time
	for (int i = 0; i < 1; i++){

	// launch the time scale algorithm
	time_scale_step(radius12, velocity, energy, volume, A12, time, dtime);
	
	// update the values to next step
	debug_print(pressure, "Press", N_shell - 1);
	update_step(velocity, A, A12, Ak12, pressure, rho, energy, dtime, dmass, dmass12, volume, viscosity, radius, radius12);
	printf("*****************************\n");
	debug_print(pressure, "Press", N_shell - 1);
	}
	/*
	output = fopen("output.dat","w");

	if (DEBUG_MODE){
		for (int i = 0; i< N_shell; i++){
			fprintf(output, "%lf %lf %lf %lf %lf %lf\n", radius[i], rho[i], energy[i], volume[i], dmass12[i], dmass[i]);
		}
	}*/
	
//	system("PAUSE");
	
	// clean up process; be a good kid!

	//fclose(output);
	
	free(rho);
	free(energy);
	free(dmass12);
	free(velocity);
	free(pressure);
	free(mass);
	free(dmass);
	free(radius);
	free(volume);
	free(radius12);
	free(A);
	free(A12);
	free(dtime);
	free(viscosity);
	free(Ak12);

	return 0;
}