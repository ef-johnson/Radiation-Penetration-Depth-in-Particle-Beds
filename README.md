# Radiation-Penetration-Depth-in-Particle-Beds
Monte Carlo Ray Tracing for beam radiation striking a particle bed.

This code is a modified version of Ray-Tracing-Many-Spheres (https://github.com/ef-johnson/Ray-Tracing-Many-Spheres). Ray Tracing Many Spheres finds the RDF or View Factors from one sphere to other spheres and walls. In contrast, "Radiation Penetration Depth in Particle Beds" is for the case of beam radiation (parallel rays) striking a particle bed. Since most of the code and parameters are the same, see Ray Tracing Many Spheres for more background detail. A short summany is given here:

The code is written in parallel C++, developed Ubuntu 16.04. Ubuntu or a similar Linux based operating system is recommended because of the other pieces of software which are useful for working with particles, using as LIGGGHTS for Discrete Element Method (DEM) modeling to generate the particle domains, and ParaView for visualization.

The particle bed must be generated in some other program. It must be defined by a space separated text file with the xyz positions of the particle centers, with one line for each set of particle coordinates. For an ordered packing (such as simple cubic, body centered cubic, etc.) these could be computed and saved as a text file with a program such as Matlab, Octave, etc. For generating realistic particle beds, LIGGGHTS is recommended for DEM. The "position" file should be named: pos_[number of particles]_[solid fraction]_[DEM time step number].txt, and placed in the pos directory. The number of particles must match the actual number of particles in the simulation. The solid fraction and DEM time step number are not actually used in the simulation, but these parameters will be passed to the name of the output file. This is done so you can keep track of many different simulations based on their solid fraction and DEM time step number. If these are unknown, they could be assigned as any numbers which help you keep track of the different simulations and their input/output files.  

Recommended Prerequisites: 
1) LIGGGHTS for doing DEM simulations, to generte particle beds 
2) OpenMPI must be installed to run parallel simulations. This is a prerequisite for LIGGGHTS as well, so most users will already have it installed from that. 
3) ParaView, recommended for visualizing particles


To run a simulation:
1) Make sure the pos file is in the pos folder
2) The source code in the main folder is compiled to make the executable with the command:  mpicxx -std=c++11 Rad_Depth_1.0.cpp -o Rad_Depth_1.0
3) The executable is run using MPI with the command: mpirun -np 10 Rad_Depth_1.0 (where the number of processors is 10 in this case)
4) 



The input parameters will be requrested by the program:

The user specifies the independent 
