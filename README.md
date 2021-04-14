# Radiation-Penetration-Depth-in-Particle-Beds
Monte Carlo Ray Tracing for beam radiation striking a particle bed.
It is written in C++ and does parallel computations with MPI.

The xyz center positions of the particles in the particle bed are supplied by the user. The code simulates emitted photons (rays) from a plane above the bed, with a specified incidence angle. Photons are followed through possibly multiple reflections off particles, until the photon is absorbed in a particle. The output files contain the locations where photons are absorbed. 

Please see the attached PDF for more information.

This code is a modified version of Ray-Tracing-Many-Spheres (https://github.com/ef-johnson/Ray-Tracing-Many-Spheres). Ray Tracing Many Spheres finds the RDF or View Factors from one sphere to other spheres and walls. In contrast, "Radiation Penetration Depth in Particle Beds" is for the case of beam radiation (parallel rays) striking a particle bed. Since most of the code and parameters are the same, see Ray Tracing Many Spheres for more background detail. 

Particle beds of up to ~100,000 particles are acceptable (memory runs out if trying to use too many particles.) 

The code is written in parallel C++, developed Ubuntu 16.04.

If you use the concepts, equations or code from this project, please cite one of the following articles:

[1] E.F. Johnson, İ. Tarı, D. Baker, Radiative heat transfer in the discrete element method using distance based approximations, Powder Technol. 380 (2021) 164–182. doi:10.1016/j.powtec.2020.11.050.

[2] E. Johnson, İ. Tarı, D. Baker, A Monte Carlo method to solve for radiative effective thermal conductivity for particle beds of various solid fractions and emissivities, J. Quant. Spectrosc. Radiat. Transf. 250 (2020). doi:https://doi.org/10.1016/j.jqsrt.2020.107014.


This code was developed for my research, but I am sharing it with others in the spirit of open source collaboration. It has been validated in several peer-reviewed articles, but it is provided without any guarantees. 


Evan
