/* ----------------------------------------------------------------------
    Radiation Penetration Depth for Particle Beds
    
    This is a Monte Carlo code to find the Radiation Distribution Factors 
    or View Factors in groups of uniform-sized spheres. 
    
    This code and more information can be found at:    
    https://github.com/ef-johnson/Radiation-Penetration-Depth-in-Particle-Beds
    
    Written by Evan Johnson
    
    This code is provided as-is, as a tool to help other researchers, but
    without any guarantees. It is provided under the GNU General Public
    License, v3.0. 
     
------------------------------------------------------------------------- */



#include <math.h>
#include "mpi.h"
#include <iostream>
#include <cmath>
#include <random>
#include <cstdio>
#include <fstream>
#include <string>
#include <iomanip>
#include <time.h>
#include <sstream>
#include <cstring>


using namespace std;


//------------------------------- MAIN -----------------------------------------------------

int main(int argc, char *argv[]) {

	// Start MPI
	MPI_Status status;
	MPI::Init(argc, argv);
	int n_procs = MPI::COMM_WORLD.Get_size();
	int rank = MPI::COMM_WORLD.Get_rank();

	// Frequently Changed Parameters

	unsigned long int n_photons;
	char char_n_photons[4]; // char strings need one more than the number of digits! More is ok though
	double vol_frac; 		
	int n_particles; 
	double abs;
	int LTS_ts;
	double radius;
	double xy_max;
	double z_plane;
	double inc_angle;
	bool specular;

if (rank == 0){



cout << endl << "Welcome to the Monte-Carlo code for finding photon extinction in a particle bed." << endl << endl;
cout << "Enter information - n_particles and vol_frac must match the pos file!" << endl << endl ;
cout << "Enter number of photons in scientific notation (ex: 1e5): ";
cin >> char_n_photons;
cout << "Enter nominal volume fraction (must be in format 0.XX): ";
cin >> vol_frac;
cout << "Enter number of particles: ";
cin >> n_particles;
cout << "Enter absorptivity of the particles (format 0.XX): ";
cin >> abs;
cout << "Enter LIGGGHTS timestep number (format XXXXXX): ";
cin >> LTS_ts;
cout << "Enter radius of the particles (meters): ";
cin >> radius;
cout << "Enter the xy max to emit photons from. Square is from x=+/- this value and y=+/- this value.  (meters): ";
cin >> xy_max;
cout << "Enter the z plane to emit photons from (meters):";
cin >> z_plane;
cout << "Enter the incidence angle (angle between the Z=0 plane the incident rays) in degrees: ";
cin >> inc_angle;
} // end if rank = 0


// Hard coding the answers to the initial questions in here, to make it easier to troubleshoot while developing!
/*
// For specular validation, using SC packed bed:
vol_frac = 0.52;
n_particles = 17500;
abs = 0.86;
LTS_ts = 1000;
radius = 0.01;
xy_max = 0.25;
z_plane = 0.14+1.1*radius;
inc_angle = 45;
specular = false; // true is for specular, false is for diffuse reflections
*/


// Use this to run multiple simulations to test different absorptivities and incidence angles:
//int n_abs = 4;
//double abs_all[4] = {0.2, 0.4, 0.6, 0.8}; 

//int n_inc_angles = 8;
//double inc_angles_all[8] = {10, 20, 30, 40, 50, 60, 70, 80};



// Big loop to iterate over all the absorptivity values
//for (int abs_i = 0; abs_i<n_abs; abs_i++){
//abs = abs_all[abs_i];

// Big loop to iterate over the whole code for different angles.
//for (int angle_i = 0; angle_i<n_inc_angles; angle_i++){

//inc_angle = inc_angles_all[angle_i];




// Send info to all processors
// Broadcast: Address of the array, number of elements to pass, data type, root processor, communication group
MPI::COMM_WORLD.Bcast(&char_n_photons, 4, MPI_CHAR, 0);
MPI::COMM_WORLD.Bcast(&vol_frac, 1, MPI_DOUBLE, 0);
MPI::COMM_WORLD.Bcast(&n_particles, 1, MPI_INT, 0);
MPI::COMM_WORLD.Bcast(&abs, 1, MPI_DOUBLE, 0);
MPI::COMM_WORLD.Bcast(&LTS_ts, 1, MPI_INT, 0);
MPI::COMM_WORLD.Bcast(&radius, 1, MPI_DOUBLE, 0);
MPI::COMM_WORLD.Bcast(&xy_max, 1, MPI_DOUBLE, 0);
MPI::COMM_WORLD.Bcast(&z_plane, 1, MPI_DOUBLE, 0);
MPI::COMM_WORLD.Bcast(&inc_angle, 1, MPI_DOUBLE, 0);



string str_n_photons(char_n_photons); // convert to string for later use in the file name

int base = char_n_photons[0] - 48; // converting from char to int requires subtracting 48 since it's from ASCII table!!
int exp  = char_n_photons[2] - 48;

n_photons = base * pow(10,exp);



if (rank == 0){
printf("\nCheck that the values were read in correctly and file names are right: \n\n");
}

printf(" --------------- Processor: %d -------------------\n", rank);
printf("n_photons: %s, vol_frac: %1.2f, n_particles: %d, abs: %f \n",char_n_photons, vol_frac, n_particles, abs);




	// print out the zPos_ctr?
	bool print_zPos_ctr = true;

	

	double cutoff_radii = 1000; // radii. lowering the cutoff distance may save time but could lead to error. 
	double cutoff_dist = cutoff_radii * radius;
		

	//View Factor Output Files
	ostringstream stream_abs;
	stream_abs << fixed << setprecision(2) << abs;
	string string_abs = stream_abs.str();
	
	ostringstream stream_vol_frac;
	stream_vol_frac << fixed << setprecision(2) << vol_frac;
	string string_vol_frac = stream_vol_frac.str();


	ostringstream stream_LTS_ts;
	stream_LTS_ts << LTS_ts;
	string string_LTS_ts = stream_LTS_ts.str();

	ostringstream stream_n_particles;
	stream_n_particles << n_particles;
	string string_n_particles = stream_n_particles.str();

	ostringstream stream_rank;
	stream_rank << rank;
	string string_rank = stream_rank.str();
	
	ostringstream stream_inc_angle;
	stream_inc_angle << inc_angle;
	string string_inc_angle = stream_inc_angle.str();

	// Position file (input)
	string pos_filename = "pos/pos_" + string_n_particles + "_" + string_vol_frac + "_" + string_LTS_ts + ".txt";
	ifstream file(pos_filename);



	// For the temporary files which write out results from each processor
	string cc_tmp = "tmp_data/zPos_ctr_" + string_vol_frac + "_" + string_n_particles + "_" + str_n_photons + "_" + string_abs + "_" + string_LTS_ts + "_" + string_inc_angle + "_proc"  +  string_rank + ".txt";
		
	string abs_loc_tmp = "tmp_data/abs_loc_all_" + string_vol_frac + "_" + string_n_particles + "_" + str_n_photons + "_" + string_abs + "_"  + string_LTS_ts + "_" + string_inc_angle + "_proc" + string_rank + ".txt";

	// For final file after combining data from all processors
	string cc_filename =    "outputs/zPos_ctr_" + string_vol_frac + "_" + string_n_particles + "_" + str_n_photons + "_" + string_abs  + "_" + string_LTS_ts + "_" + string_inc_angle + ".txt";


	// For radiation penetration: print out the location of the absorption
	string abs_loc_filename = "outputs/abs_loc_all_" + string_vol_frac + "_" + string_n_particles + "_" + str_n_photons + "_" + string_abs + "_" + string_LTS_ts + "_" + string_inc_angle + ".txt";

	
	cout << "File names:"  << endl;
	cout << pos_filename << endl;
	cout << cc_tmp << endl;
	cout << cc_filename << endl;
	cout << abs_loc_filename << endl << endl;



	// CC vs View Factor output
	ofstream out_file2;	
	if (print_zPos_ctr == true){ 	
	out_file2.open (cc_tmp);}
	
	// Absorption Location output
	ofstream out_file_abs_loc_tmp;
	out_file_abs_loc_tmp.open (abs_loc_tmp);
	
	// Emission Location output
	// Good for visualizing how the photon reflects around the particle bed
	//ofstream out_file_emiss_loc;
	//out_file_emiss_loc.open ("emiss_loc.txt");

		



	constexpr int n_read_col = 3;




	// Random number generator
	random_device rd;   // non-deterministic generator
	mt19937 gen(rd());  // rd() is a random seed. Replace rd() with a number to get a consistent seed
	uniform_real_distribution<> dist(0,1); // distribute results between 0 and 1 inclusive
	


	double pos[n_particles][n_read_col]; // x,y,z position and temps of each particle		

	double cc_dist[n_particles]; // center-center distance from home particle to all other particles
	double pi = atan(1)*4; // find pi

	int particles_simulated = 0;
	int particles_simulated_all[n_procs];

	clock_t t_before;
	clock_t t_after;

	// Declare variables for later
	double emiss_loc_rel[3];
	double emiss_loc[3];
	double last_emiss_loc[3]; // needed for specular case
	
	int current_emit_id; // needed for specular case
	double alpha;
	double beta;

	// Geometry info:
	// V1: vector from emission location to center of particle being checked
	// V1_mag: distance from emission location to center of particle being checked
	// V2: Vector from emission location to the destination of the photon
	// V2_mag: Magnitude of V2 doesn't matter - but it's necessary to have a value for angle calculations
	// zeta: angle formed between V1 and V2
	// zeta_max: maximum angle between V1 and V2 to hit the sphere



if (rank == 0){


	// ----------------------- Read file with 3 columns: x, y, z -------------------------------

	int counter = 0;
	if(file.is_open()){
			for (int row = 0; row<n_particles; row++){
				for(int col = 0; col < n_read_col; col++){
					
					file >> setprecision(8) >> pos[row][col];
										
					counter++;
				}
			}
	}
	
	else {printf("\n\nCould not open pos/pos_xxxxxx file. Exiting.\n\n");
		exit (EXIT_FAILURE);
	}


			

} // end if rank = 0


// Broadcast: Address of the array, number of elements to pass, data type, root processor, communication group
MPI::COMM_WORLD.Bcast(&pos, n_particles*3, MPI_DOUBLE, 0);


// *********************START HOME PARTICLE LOOP ************************  //




		// ------ Create neighbor list - just one neighbor list now for photon extinction  -------------------
		

		double NL_xPos[n_particles+1];
		double NL_yPos[n_particles+1];
		double NL_zPos[n_particles+1];
		double NL_hits[n_particles+1]; // Make this a double so it can be divided easily later
		

		
		for (int nei_id=0; nei_id<n_particles; nei_id++){

			double p_check[3];
			p_check[0] = pos[nei_id][0];
			p_check[1] = pos[nei_id][1];
			p_check[2] = pos[nei_id][2];
			
			



			// Build the neighbor list with the id, cc_dist, XYZ coordinates, and space for the view factor calcs which come later
			// Changing the neighbor_list into:
			
			NL_xPos[nei_id] = p_check[0]; // x-pos
			NL_yPos[nei_id] = p_check[1]; // y-pos
			NL_zPos[nei_id] = p_check[2]; // z-pos
			NL_hits[nei_id] = 0.0; // total hits


		} // end of building neighbor list
		

			NL_xPos[n_particles] = 0; // x-pos
			NL_yPos[n_particles] = 0; // y-pos
			NL_zPos[n_particles] = 0; // z-pos
			NL_hits[n_particles] = 0.0; // total hits
			



		//----------------- Monte Carlo for this home particle --------------------------


		// For this home particle, send photons, then check with each neighbor. Use Unsigned long int since i and n_photons are very large.
		for (unsigned long int i = 0; i<n_photons; ++i){

			// This loops through as many reflections as needed until the photon absorbs.
			// The first time through, refl is set to true so that the original emission happens, and the emission location is set by alpha and beta
			// On successive runs through (until refl is set to false, representing an absorption) the emission location is calculated
			// based on the incident location on the sphere surface.
			bool refl = true;
			bool first_run = true;
			bool initial_emiss = true;
			bool overlapping = false; // set to true if the photon is released inside a neighboring sphere

			int refl_count = 0;

			// EMISSION/REFLECTION LOOP
			while (refl == true){
			
	
				// On the first run, the emission location is a random place on the surface
				if (first_run == true) {

					// Choose angles for calculation of emission location and direction
					// Beta is angle down from vertical axis - choose a random angle between 0 and pi, where
					// Alpha is azimuth angle from x axis - choose a random angle between 0 and 2*pi
					
					//beta = acos(2*dist(gen) - 1);
					//alpha = 2*pi*dist(gen);
					


					// Find emission location relative to the emitting particle center, using alpha and beta

					
					// Choose the random x, y location of emission of the photon
					emiss_loc[0] = (dist(gen)*2-1) * xy_max; // x location is between -xy_max and xy_max
					emiss_loc[1] = (dist(gen)*2-1) * xy_max; // y location is between -xy_max and xy_max
					emiss_loc[2] = z_plane; 
					
			
					// Track the id of the particle where photon is emitted from. It starts here with the home particle.
					current_emit_id = -999; // set this for the original emission

					// The first emission location has been set. Make first_run false now so next time
					// the emission location will be determined by the incident point on the sphere hit by the photon
					first_run = false;
				}
				// If this is NOT the first run, the incident location, alpha, and beta are already known from the last iteration.
				// They are set at the bottom of this while(refl==TRUE) loop
				
				//else {
					// find incident location from the V2, home_coords, the id and xyz of the particle that was hit, etc.
					// find the alpha and beta from the sphere that photon is leaving from
					// This is done whether or not it's the first run, so taking out this else statement

				//}


				// Choose angles for direction of emission

				// Initially, choose angles in temporary axes, where axis is mounted on surface of sphere at emission location
				// with z in line with vector from center of sphere to emission location, and x parallel to the original x-axis
				// Eta is angle down from vertical axis - choose a random angle between 0 and pi/2
				// Gamma is azimuth angle, measured from x-axis - choose a random angle between 0 and 2*pi
				double eta;
				double gamma;
				
				
				// For initial emission directions:
				
				// Set photon emission angles. On the initial emission of this photon (for extinction code) manually set the angles. 
				// If not, emission is a re-emission due to a reflection, so choose the angles by specular or diffuse angles.
				if (initial_emiss == true){
					
					// alpha = 0 means all rays are aligned to y-normal planes (from the emission loc it leaves in the xz plane)			
					alpha = 0; 
					
					// beta is angle from +z axis to the ray. So since incident angle (inc_angle) input by user is the angle between xy plane
					// and the rays, beta is pi plus inc_angle. 
					beta = pi/2 + inc_angle*pi/180; // inc_angle is converted to radians. Beta measured from Z+ axis down to the ray
					
					// eta and gamma for releasing first ray should be exactly aligned with beta and alpha
					eta = 0; // 
					gamma = 0; // 
	

					
					
					initial_emiss = false; // now that initial emission is ready, set it to false for the next run through the emission/reflection loop

	
	
				}
				
				// For reflections:
				else{	
				
					
					// Diffuse reflection angles
					if (specular == false){
					eta = asin( sqrt(dist(gen)) );
					gamma = 2*pi*dist(gen);
					
					}
					
					// Specular reflection angles
					else{
					
					// Convert the incident vector from the global XYZ coord sys to the local x'y'z' which is at the reflecting particle's surface (and tilted at alpha and beta)
					// This is the other way from the original Euler Angles implemented, where I go from x'y'z' and convert to the XYZ coord system
					// Then, use this incoming ray to find the incident eta and incident gamma values
					// This requires another use of Euler Angles:
					
					// Calculate vector inc_ray in the XYZ (overall) coord sys. The incident ray comes from last_emiss_loc to emiss_loc, but I actually want the vector
					// going from the origin of the x'y'z' coord sys outward towards last_emiss_loc. This way, I can calculate the azimuth and elevation angles. 
					double inc_ray_XYZ[3];
					
					inc_ray_XYZ[0] = last_emiss_loc[0] - emiss_loc[0];
					inc_ray_XYZ[1] = last_emiss_loc[1] - emiss_loc[1];
					inc_ray_XYZ[2] = last_emiss_loc[2] - emiss_loc[2];
					


					// Use Euler angles to rotate from the temporary x-,y-,z-prime coord sys back to a coord sys aligned
					// with the overall XYZ axis where Z is up, X is towards you, and Y is to the right.
					// This is the "y-convention" euler angles, where the rotation is around the z axis first, then
					// around y axis, then around z again. In this case, the rotations are:
					// 1) 0 radians around Z
					// 2) -beta around y’ axis, which brings z’ into alignment with Z
					// 3) -alpha around z” (aka Z) so x’ aligns with X and y’ aligns with Y

					// For reflection case, I have the vector (inc_ray_XYZ) in the overall XYZ axes to start with, and I want them in the x'y'z'.

					// Alpha and Beta should be known from the previous run, where the relfection was determined
					
					//Euler transform matrix- going from XYZ to x'y'z' system. 
					// This is the ORIGINAL ZYZ equation, not the transposed matrix.
					
					beta = -beta;
					alpha = -alpha;
					double a11= cos(alpha)*cos(beta);
					double a21= sin(alpha);
					double a31=-cos(alpha)*sin(beta);
					
					double a12=-sin(alpha)*cos(beta);
					double a22= cos(alpha);
					double a32= sin(alpha)*sin(beta);
					
					double a13= sin(beta);
					double a23= 0;
					double a33= cos(beta);
					beta =  -beta;
					alpha = -alpha;

				
					//Multiply transform matrix by the original inc_ray_XYZ vector to get the vector in the new x'y'z' coord system 
					double inc_ray_prime[3];
					inc_ray_prime[0] = a11*inc_ray_XYZ[0] + a12*inc_ray_XYZ[1] + a13*inc_ray_XYZ[2];
					inc_ray_prime[1] = a21*inc_ray_XYZ[0] + a22*inc_ray_XYZ[1] + a23*inc_ray_XYZ[2];
					inc_ray_prime[2] = a31*inc_ray_XYZ[0] + a32*inc_ray_XYZ[1] + a33*inc_ray_XYZ[2];					
				
					// Find the eta and gamma angles, in the x'y'z' coord sys.
							 
					// Find alpha and beta from the hit-particle's center by converting to spherical coordinates
					// Use the x,y,z values of a point defined by the inc_ray_prime vector (from the x'y'z' origin towards the location where photon was emitted)
					
					double inc_ray_mag = sqrt(pow(inc_ray_prime[0],2) + pow(inc_ray_prime[1],2) + pow(inc_ray_prime[2],2));

					
					double x_prime = inc_ray_prime[0];
					double y_prime = inc_ray_prime[1]; 
					double z_prime = inc_ray_prime[2]; 
					


					// Calculate gamma. atan gives an angle in the range -pi/2 to +pi/2, so adjust as necessary
					double inc_gamma = atan2(y_prime, x_prime); // accounts correctly for the diffrerent quadrants, as opposed to atan() which is ambiguous
					double inc_eta  = acos( z_prime / inc_ray_mag );
					
					eta = inc_eta; // eta is measured down from z axis, so incoming and outgoing specular ray have same eta value
					gamma = inc_gamma + 3.1415926535; // gamma (azimuth, measured from temporary x' axis) of the outgoing ray is 180 degrees from the incident gamma angle 
								

					}
					
					//printf("in ELSE: eta and gamma: %f, %f \n", eta, gamma);
				}
				
				
				// Print out the emission location, either from the initial emission or from a reflection
				// Good for visualizing how a photon is reflecting around a particle bed
				//out_file_emiss_loc << std::fixed << std::setprecision(8) <<  emiss_loc[0] << " " << emiss_loc[1] << " " << emiss_loc[2] << endl; 

		
				// FIND EMISSION DIRECTION in XYZ axes
		
				// Calculate vector V2 in the temporary axes (denoted by "prime")
				double V2_mag = 1;
				double V2_prime[]= {V2_mag*sin(eta)*cos(gamma), V2_mag*sin(eta)*sin(gamma), V2_mag*cos(eta)};

//				printf("V2_prime, refl ray in prime coord sys: %f, %f, %f \n", V2_prime[0], V2_prime[1], V2_prime[2]);
				
				
				// Use Euler angles to rotate from the temporary x-,y-,z-prime coord sys back to a coord sys aligned
				// with the overall XYZ axis where Z is up, X is towards you, and Y is to the right.
				// This is the "y-convention" euler angles, where the rotation is around the z axis first, then
				// around y axis, then around z again. In this case, the first rotation around z-axis is 0, and
				// the rotation around the y is negative alpha, and the rotation around z is negative beta.
				// The negative are needed to convert from the current x'y'z' back to the original by the angles
				// alpha and beta chosen for the emission location.

				//Euler transform matrix
				
				// This is actually the TRANSPOSE matrix, since we're going from x'y'z' to XYZ
				
				beta = -beta;
				alpha = -alpha;
				double a11= cos(alpha)*cos(beta);
				double a12= sin(alpha);
				double a13=-cos(alpha)*sin(beta);
				double a21=-sin(alpha)*cos(beta);
				double a22= cos(alpha);
				double a23= sin(alpha)*sin(beta);
				double a31= sin(beta);
				double a32= 0;
				double a33= cos(beta);
				beta = -beta;
				alpha = -alpha;
/*				
				printf("a11 a12 a13: %f, %f, %f \n", a11, a12, a13);
				printf("a21 a22 a23: %f, %f, %f \n", a21, a22, a23);
				printf("a31 a32 a33: %f, %f, %f \n", a31, a32, a33);
*/				
				
				// Multiply transform matrix by x',y',z' to get vector in new coord sys aligned with XYZ
				double V2[3];
				V2[0] = a11*V2_prime[0] + a12*V2_prime[1] + a13*V2_prime[2];
				V2[1] = a21*V2_prime[0] + a22*V2_prime[1] + a23*V2_prime[2];
				V2[2] = a31*V2_prime[0] + a32*V2_prime[1] + a33*V2_prime[2];

//				printf("V2, in XYZ coord sys: %f, %f, %f \n", V2[0], V2[1], V2[2]);




				//----------------- Check neighbors for hits and absorption --------------------------

				// Re-Initialize the ID and distance for the closest neighbor to emission location, for this photon
				double smallest_V1_mag = 1e9;
				int closest_id = -99;

				// Check to see which neighbors the photon would hit
				// Iterate from zero until the number of neighbors has been reached
				for (int nei_id = 0; nei_id < n_particles; nei_id++){

					// Don't check the particle that is currently emitting the photon, or else it would always be closest!
					if (nei_id != current_emit_id){

						// Get the xyz position of the particle we're checking
						double p_check[] = {NL_xPos[nei_id],NL_yPos[nei_id],NL_zPos[nei_id]};

						// Calculate vector V1 from emission location to p_check center
						double V1[3] = {p_check[0]-emiss_loc[0], p_check[1]-emiss_loc[1], p_check[2]-emiss_loc[2]};
						double V1_mag = sqrt(pow(V1[0],2) + pow(V1[1],2) + pow(V1[2],2));
				
			
						
						// Special case: if the photon emits from inside the overlapping spheres, consider it an
						// automatic absorption. In reality, it would hit, but where would the incident location be
						// calculated? It's not physically realistic. It would likely reflect a couple times at
						// most until absorbing.

						//if (V1_mag < radius) { // if V1_mag (the dist from emiss_loc to center of particle we're checking) is less than the radius
						if (V1_mag < radius*0.999999) { // instead of just V1_mag<radius. Before, when I used a SC lattice, the values were basically equal sometimes, even when reflecting... so now V1_mag has to be a very tiny amount less than the radius to actually count as an overlap. 
						// then the photon must be emitting from inside the neighboring sphere.
							overlapping = true;
							closest_id = nei_id; // directly set this as the closest neighbor that hits, and setting overlapping=true prevents 
							// having to do further calculations for all the other neighbors due to the (overlapping==false) condition in the 
							// IF statement below. 
							
							// It should only happen on the first emission, not on a reflection
							
							//printf("Problem in overlapping particle. V1_mag: %1.10f, particle radius: %1.10f \n", V1_mag, radius);
						}
						
						

						// Only check for hits if current neighbor is closer than the current closest hit,
						// AND only check if neighbor is closer than the cutoff distance.
						// AND only check if the photon is not being released inside of a neighboring sphere.
						if ((V1_mag < smallest_V1_mag) && (V1_mag < cutoff_dist) && overlapping == false){

							// Calculate angle between V1 and V2, using dot product
							double zeta = acos( (V1[0]*V2[0] + V1[1]*V2[1] + V1[2]*V2[2])/(V1_mag * V2_mag));

							// Calculate maximum angle to qualify as a hit
							double zeta_max = asin(radius/V1_mag);

						
							// The photon will hit the neighbor particle if the angle is less than the max angle
							if (zeta < zeta_max){
								// If the V1_mag (emission location to neighbor_center) for this neighbor is less than the previous closest value,
								// save the value and id for the new closest hit

								smallest_V1_mag = V1_mag;
								closest_id = nei_id;
								
							
								//}
							} // end if zeta < zeta_max
						
						} // end if current neighbor is closer than the current closest hit
						
					} // end exclude the current emitting particle from the neighbor list
					
				} // end neighbor list for loop, which checks to see which neighbors the photon would hit




				//----------------- Check WALL for hits and absorption --------------------------
				

				
				// For the current home particle, if the current photon didn't hit any neighbors, check if it hits the wall or not
				if (closest_id == -99){
				
				// In this extinction code, if it doesn't hit any neighbors, just set refl = false so it exits the loop, and the particle must be absorbed by the wall
			
					// Wall is assumed to be under all of the particles (at Z=0)
					// Check for wall hit: any photon not hitting another particle AND going in the negative z-direction must hit the wall
						if (V2[2] < 0) { // photon hits wall at Z=0
						
						 
							refl = false; // photon terminated
							NL_hits[n_particles]++; // add one to the count of photons absorbed by the wall, which is the 
							
							
							// Calculate the incident location where the photon crosses the Z=0 plane.
							// With the parametric equation x2i + y2j + z2k = (x1+s*x_vec)i + (y1+s*y_vec)j + (z1+s*z_vec)k
							// where x2 is the point we are looking for where the vector (x_vec i, y_vec j, z_vec k) passes through the z=0 plane
							// Here, x_vec, y_vec, z_vec is the V2[] vector, which has magnitude 10 set with V2_mag at the beginning
							// and (x1, y1, z1) is the emission location where the vector leaves from, and s is the parametric variable.
							// For the special case of the plane Z=0, s can be solved with s=z1/z_vec		
							double inc_loc[3];						
							double s = -emiss_loc[2]/V2[2];
							inc_loc[0] = emiss_loc[0] + s*V2[0];
							inc_loc[1] = emiss_loc[1] + s*V2[1];
							inc_loc[2] = 0;
							
							
							out_file_abs_loc_tmp << std::fixed << std::setprecision(8) <<  inc_loc[0] << " " << inc_loc[1] << " " << inc_loc[2] << endl;
														
							
 
						
					}	
					
					
					
					
						// If photon is not crossing the wall at Z=0, then the photon is done and didn't absorb into any particles or the
						// wall - it must have gone out one of the other sides of the domain.
						else{
						refl = false; // photon terminated
						
						// Need to find the point where we can call the photon terminated, the lost location. From the emission point, follow V2 a distance of V2_mag before recording the point.
						double lost_loc[3];
						lost_loc[0] = emiss_loc[0] + V2[0];
						lost_loc[1] = emiss_loc[1] + V2[1];
						lost_loc[2] = emiss_loc[2] + V2[2];
						
						
						}				
						
	
					} // end if the photon hit NO particles
				





				//----------------- PARTICLE absorption or reflection --------------------------
				
				// If the photon did hit some neighbor(s)
				else  {
				// printf("Particle hit a neighbor!\n");
					// Special case for photon released from inside the neighboring sphere -> an automatic absorption
					if (overlapping == true){
						refl = false; // photon terminated, absorbed into overlapping neighbor
						NL_hits[closest_id]++; // add one to the count of photons absorbed by that neighbor
						printf("Photon release inside sphere: direct absorption\n");	
										
					}
					
					// Photon not released inside an overlapping particle (almost always the case)
					else{ 
						
						
						// Calculate the incident location. Use this either to find the alpha/beta values in the case of a reflection, or to print out the absorption location in the case of absorption.

							// Get the center coordinates of the sphere that it actually hits
							double dest_ctr_coords[] = {pos[closest_id][0], pos[closest_id][1], pos[closest_id][2] };
							

							// Prep for quadratic:
							double a = V2[0]*V2[0] + V2[1]*V2[1] + V2[2]*V2[2];

							double dx = emiss_loc[0]-dest_ctr_coords[0];
							double dy = emiss_loc[1]-dest_ctr_coords[1];
							double dz = emiss_loc[2]-dest_ctr_coords[2];

							double b = 2*(dx*V2[0] + dy*V2[1] + dz*V2[2]);
							double c = -(radius*radius) + dx*dx + dy*dy + dz*dz;

							// Quadratic to find t:
							double t1 = (-b + sqrt( b*b - 4*a*c))/(2*a);
							double t2 = (-b - sqrt( b*b - 4*a*c))/(2*a);

							// Find the xyz locations where it hits
							double x1 = emiss_loc[0] + t1*V2[0];
							double y1 = emiss_loc[1] + t1*V2[1];
							double z1 = emiss_loc[2] + t1*V2[2];

							double x2 = emiss_loc[0] + t2*V2[0];
							double y2 = emiss_loc[1] + t2*V2[1];
							double z2 = emiss_loc[2] + t2*V2[2];

							double mag1 = sqrt ( (x1-emiss_loc[0])*(x1-emiss_loc[0]) + (y1-emiss_loc[1])*(y1-emiss_loc[1]) + (z1-emiss_loc[2])*(z1-emiss_loc[2]) );
							double mag2 = sqrt ( (x2-emiss_loc[0])*(x2-emiss_loc[0]) + (y2-emiss_loc[1])*(y2-emiss_loc[1]) + (z2-emiss_loc[2])*(z2-emiss_loc[2]) );

							
							// Vector will pass thru sphere twice, one going in, and one coming out
							// Take whichever point is closest to the emission location to be the incident location
							double inc_loc[3];
							if ( mag1 < mag2){
								inc_loc[0] = x1;
								inc_loc[1] = y1; // first point is closer,
								inc_loc[2] = z1;
							}
							else{
								inc_loc[0] = x2;
								inc_loc[1] = y2; // second point is closer,
								inc_loc[2] = z2;
							}
												
						
						
						//  Check to see if photon is absorbed
						double rand_num = dist(gen);
						
						if (rand_num<abs){
							refl = false; // photon terminated, absorbed
							
							
							NL_hits[closest_id]++; // add one to the count of photons absorbed by that neighbor	
								
							// If absorption is into a particle, write out the incident location:				
							out_file_abs_loc_tmp << std::fixed << std::setprecision(8) <<  inc_loc[0] << " " << inc_loc[1] << " " << inc_loc[2] << endl; // print xyz position to abs_loc_tmp file 
														
						}

						else{ // If photon is reflected

							refl = true;
							current_emit_id = closest_id; // Set this for the next run through the reflection loop


							// Find Alpha and Beta angles for next run through. 
							// The x_tmp, y_tmp, z_tmp coord sys is aligned with the overall, original XYZ system, but center is at the relfecting particle's center. 
							
							// Because the photon is reflecting from the point where it hit the incident sphere's surface,
							// we need to find the alpha and beta angles of that point for the next run through the reflection loop
							// instead of finding them randomly like we did for the original emitting sphere. 
							// Find alpha and beta from the hit-particle's center by converting to spherical coordinates
							// Use the x,y,z values from the center of the hit-particle
							double x_tmp = inc_loc[0] - dest_ctr_coords[0];
							double y_tmp = inc_loc[1] - dest_ctr_coords[1];
							double z_tmp = inc_loc[2] - dest_ctr_coords[2];


							// Calculate alpha. atan gives an angle in the range -pi/2 to +pi/2, so adjust as necessary
							alpha = atan2( y_tmp, x_tmp );
						
							beta  = acos( z_tmp / radius );
							
							
							// ADDING FOR SPECULAR CASE
							// Save the emission location from this photon so the incidence angle can be computed if the photon is reflected
							last_emiss_loc[0] = emiss_loc[0];
							last_emiss_loc[1] = emiss_loc[1];
							last_emiss_loc[2] = emiss_loc[2];
							
							
							// Set the emission location for the next run as the incident location in this run
							emiss_loc[0] = inc_loc[0];
							emiss_loc[1] = inc_loc[1];
							emiss_loc[2] = inc_loc[2];
							
							

						} // end check if reflected
						
					} // end if not overlapping 

				} // else if this photon hit at least one neighbor

				refl_count++;

			} // end while refl = true loop (aka the reflection loop)

		} // end photon for loop



		// Find total photons absorbed, total missing photons, and total view factor for this particular home particle

		double total_abs = 0;
		for (int i=0;i<n_particles+1;i++){
			 //NL_RDF[i] = NL_hits[i] / n_photons; // Calculate the View Factor for each neighbor
			 total_abs = total_abs + NL_hits[i];
		}



		// Warn if there is an escaped photon
		//if (n_photons != total_abs){
		//	printf("Lost a photon!! %lu photons sent, but only %f absorbed!\n", n_photons,total_abs);
		//}
		//double total_vf = total_hits / n_photons;


		// Print these center-center distances and view factors to the output text file
		if (print_zPos_ctr == true){ 
			for (int row=0;row<n_particles + 1;row++){ // Adding 1 back in for the wall here
			
				// Only print the AVF if it is above zero just to save file space. It can reduce file size by about 200 times!
				//if (neighbor_list[row][6] > 0){ - no that's no good! you need those zeros - if not you only print the non-zero terms wich skews the curve upward slightly (maybe not appreciably though)
				//if (neighbor_list[row][1] < cutoff_dist){

//					a << std::fixed << std::setprecision(0) <<  home_num; // "from" particle id
//					out_file2 <<  " ";
					if (row < n_particles){
					out_file2 << std::fixed << std::setprecision(0) <<  row ; // "to" particle id
					}
					else{
						out_file2 << std::fixed << std::setprecision(0) <<  -1 ; // the wall
					}
					
					out_file2 <<  " ";
					out_file2 << std::fixed << std::setprecision(8) <<  NL_zPos[row]; // z position 
					out_file2 <<  " "; 					
					out_file2 << std::fixed << std::setprecision(0) <<  NL_hits[row]; // Number of hits
					out_file2 <<  endl;
					//} // end only print RDF if it's less than the cutoff distance
				//} // end print RDF if it's above zero
			}
		}
		


		

		// Outputs to console



		if (rank == 0){
			printf("Finished simulating. Sent %lu photons per processor, %1.0f absorbed on proc0.\n", n_photons, total_abs);
		
			
		} // end if rank = 0



MPI_Barrier(MPI::COMM_WORLD);

// MPI_Gather(&particles_simulated, 1, MPI_INT, &particles_simulated_all, MPI_INT, 0, MPI::COMM_WORLD);

if (rank == 0){

	// Open the files to print the entire output
	cout << endl << "Iterating through the temporary files: " << endl;

	// CC vs View Factor output
	ofstream out_file4;
	if (print_zPos_ctr == true){ 
		out_file4.open (cc_filename);
	}
	
		// Absorption Location output
		ofstream out_file_abs_loc;
		out_file_abs_loc.open (abs_loc_filename);
	

	

	for (int rnk = 0; rnk <n_procs;rnk++){

		ostringstream stream_rnk;
		stream_rnk << rnk;
		string string_rnk = stream_rnk.str();
		string inFile_cc = "tmp_data/zPos_ctr_" + string_vol_frac + "_" + string_n_particles + "_" + str_n_photons + "_" + string_abs + "_" + string_LTS_ts +"_" + string_inc_angle + "_proc" + string_rnk + ".txt";
		
		
		// For Absorption Location files
		string inFile_abs_loc = "tmp_data/abs_loc_all_" + string_vol_frac + "_" + string_n_particles + "_" + str_n_photons + "_" + string_abs + "_" + string_LTS_ts +"_" + string_inc_angle + "_proc" + string_rnk + ".txt";
		

		
		string line;

		// Print the zPos ctr output file
		if (print_zPos_ctr == true){
			ifstream cc_filestream(inFile_cc);
			cout << inFile_cc << endl;
			
				// Check to make sure the files to read are open
			    if (!cc_filestream.is_open() )
			    {
				cout << "Path Wrong (in print_zPos_ctr)!!!! Check the file name and location." << endl;
				exit(EXIT_FAILURE);
			    }

			// Iterates through each line of the file, putting all data into "line" string. Then insert that line
			// into the overall output file.			
			while (getline(cc_filestream,line))
			{
			    istringstream iss(line);
			    string lineStream;
			    out_file4 << line << endl;
			}
		} // end if print_zPos_ctr = true
		
		
		
		
		// Print the Absorption Locations

			ifstream abs_loc_filestream(inFile_abs_loc);
			cout << inFile_abs_loc << endl;
			
			// Check to make sure the files to read are open
			    if (!abs_loc_filestream.is_open() )
			    {
					cout << "Path Wrong (in Absorption Location print-out)!!!! Check the file name and location." << endl;
					exit(EXIT_FAILURE);
			    }

			// Iterates through each line of the file, putting all data into "line" string. Then insert that line
			// into the overall output file.			
			while (getline(abs_loc_filestream,line))
			{
			    istringstream iss(line);
			    string lineStream;
			    out_file_abs_loc << line << endl;
			}

		} // end if rank = 0


out_file4.close();



printf("\nSuccessfully finished code execution.\n");
printf("Absorptivity: %1.2f \n", abs);
printf("Photons per particle: %lu \n",n_photons);
printf("Volume Fraction: %1.2f \n\n", vol_frac); 



} // end if rank = 0


//} // Big loop to iterate over the whole code for different angles.
//} // Big loop to iterate over the different absorptivities.

MPI_Barrier(MPI::COMM_WORLD);



MPI::Finalize();

return 0;


}
