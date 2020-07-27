#ifndef DFTYPES_HPP
#define DFTYPES_HPP
#include "./include_eigen.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include <math.h>

using namespace std;

/*
particle() is a datatype used to define the particle that will diffuse through the medium or potential conditions imposed on the problem
*/
struct particle{
	Eigen::Vector3i position;//position is defined as lattice indices
	double mass;//mass of the particle
	double stickiness;//the probability that a particle will stick to this particle upon collision
	Eigen::Vector3d velocity;//the translational velocity of this particle
	double charge;//he charge of this particle
	double sigma; double epsilon;//sigma and epsilon are the lennard-jones parameters
	double radius;//radius of the particle
	double penetration_depth;//this metric defines the distance to which the particle is considered bound, if 1.0 then the full radius is used
	int spin;//definition of a spin
	//constructor and destructor zone below:
	particle();
	particle(Eigen::Vector3i indices, double gradius, double gm, double gstick, Eigen::Vector3d gvel, double gq, double gsigma, double geps, double gpendepth, int gspin);
	particle(const particle & o);
	~particle(){}
	particle operator=(const particle & o);
	//member functions
	//1) Assuming a field provides a force we use newtons laws to transform the velocity and position
	Eigen::Vector3d find_acceleration_electrostatic(Eigen::Vector3d fieldvalue, Eigen::Vector3d r);//the 3D vector field at the point of interest is represented by the vector input and the new position 
	//2) Interaction of this particle with another particle
	Eigen::Vector3d interact_to(const particle & o, Eigen::Vector3d dielectric_const, double dx, double dy, double dz);
};

//defining the default constructor
particle::particle(): position(0,0,0), radius(1.0), mass(1.0), stickiness(1.0), velocity(0.0,0.0,0.0), charge(0.0), sigma(1.0), epsilon(0.0), penetration_depth(1.0), spin(0){}
//defining the parametric constructor
particle::particle(Eigen::Vector3i indices, double gradius, double gm, double gstick, Eigen::Vector3d gvel, double gq, double gsigma, double geps, double gpendepth, int gspin): position(0,0,0), radius(gradius), mass(gm), stickiness(gstick), velocity(0.0,0.0,0.0), charge(gq), sigma(gsigma), epsilon(geps), penetration_depth(gpendepth){
	position = indices; velocity = gvel; spin = gspin;
}
//defining the copy constructor
particle::particle(const particle &o): position(0,0,0), radius(1.0), mass(1.0), stickiness(1.0), velocity(0.0,0.0,0.0), charge(0.0), sigma(1.0), epsilon(0.0){
	position = o.position; radius = o.radius; mass = o.mass; stickiness = o.stickiness; velocity = o.velocity; charge = o.charge; sigma = o.sigma; epsilon = o.epsilon; spin = o.spin;
}
//defining the equality(=) operator
particle particle::operator=(const particle & o){
	position = o.position; radius = o.radius; mass = o.mass; stickiness = o.stickiness; velocity = o.velocity; charge = o.charge; sigma = o.sigma; epsilon = o.epsilon; spin = o.spin;
	return *this;
}
//1) Member function to find the acceleration
Eigen::Vector3d particle::find_acceleration_electrostatic(Eigen::Vector3d fieldvalue, Eigen::Vector3d r){
	double Fx,Fy,Fz;
	Fx = ((fieldvalue(0))*charge)/((r(0)));
	Fy = ((fieldvalue(1))*charge)/((r(1)));
	Fz = ((fieldvalue(2))*charge)/((r(2)));
	Eigen::Vector3d acceleration((Fx/mass), (Fy/mass), (Fz/mass));
	return acceleration;
}

//2) Interaction to another particle, having a 3d dielectric constant allows us to account for anisotropies if they arise eventually
Eigen::Vector3d particle::interact_to(const particle & o, Eigen::Vector3d dielectric_const, double dx, double dy, double dz){
	double Eqx, Eqy, Eqz, Eljx, Eljy, Eljz;
	Eigen::Vector3i dr; dr = (position - (o.position));
	Eigen::Vector3d r; r(0) = (dx*((double)dr(0))); r(1) = (dy*((double)dr(1))); r(2) = (dz*((double)dr(2)));
	Eqx = (((o.charge)*charge)/(4.0*3.1415962*(dielectric_const(0))*(r(0))));
	Eqy = (((o.charge)*charge)/(4.0*3.1415962*(dielectric_const(1))*(r(1))));
	Eqz = (((o.charge)*charge)/(4.0*3.1415962*(dielectric_const(2))*(r(2))));
	double sigij,epsij;
	sigij = (sigma+(o.sigma))/2.0; epsij = pow(((o.epsilon)*epsilon),0.5);
	Eljx = 4.0*epsij*( (pow( (sigij/(r(0))),12.0)) - (pow( (sigij/(r(0))),6.0)) );
	Eljy = 4.0*epsij*( (pow( (sigij/(r(1))),12.0)) - (pow( (sigij/(r(1))),6.0)) );
	Eljz = 4.0*epsij*( (pow( (sigij/(r(2))),12.0)) - (pow( (sigij/(r(2))),6.0)) );
	Eigen::Vector3d Enet( (Eqx+Eljx), (Eqy+Eljy), (Eqz+Eljz) );
	return Enet;
}

/*
we define the center_of_mass() function to find the center of mass of a set of particles
*/
Eigen::Vector3d center_of_mass(vector<particle> gparticles, double dx, double dy, double dz){
	vector<particle>::iterator it; it = gparticles.begin();
	Eigen::Vector3d mass_weighted_sum(0.0,0.0,0.0); double mass_sum = 0.0;
	while(it != gparticles.end()){
		particle newparticle(*it);
		double mass_curr = newparticle.mass;
		Eigen::Vector3i ri_curr; ri_curr = newparticle.position;
		Eigen::Vector3d r_curr; r_curr(0) = (dx*((double)ri_curr(0))); r_curr(1) = (dy*((double)ri_curr(1))); r_curr(2) = (dz*((double)ri_curr(2)));
		mass_weighted_sum = (mass_weighted_sum + (mass_curr*(r_curr)));
		mass_sum = (mass_sum + mass_curr);
		++it;
	}
	Eigen::Vector3d Rcom; Rcom = ((1.0/mass_sum)*mass_weighted_sum);
	return Rcom;
}
/*
We define the radius_of_gyration() function to find the radius of gyration given a set of particles and a center of mass, that ideally corresponds to the set
*/
double radius_of_gyration(vector<particle> gparticles, Eigen::Vector3d Rcom, double dx, double dy, double dz){
	vector<particle>::iterator it; it = gparticles.begin();
	double N = 0.0; double rdiffsum = 0.0;
	while(it != gparticles.end()){
		particle newparticle(*it);
		Eigen::Vector3i ri_curr; ri_curr = newparticle.position;
		Eigen::Vector3d r_curr; r_curr(0) = (dx*((double)ri_curr(0))); r_curr(1) = (dy*((double)ri_curr(1))); r_curr(2) = (dz*((double)ri_curr(2)));
		Eigen::Vector3d dr; dr = r_curr - Rcom; double drnorm = dr.norm();
		rdiffsum = (rdiffsum + (drnorm*drnorm));
		N = (N + 1.0);
		++it;
	}
	double Rgyr = sqrt((rdiffsum/N));
	return Rgyr;
}
/*
We define a fuction that calculates the fractal dimension using the radius of gyration relation N ~ Rgyr^Dh
*/
/*
cluster() is a datatype that contains the index positions of all particles that become frozen upon collision and also decides if the collision is achieved and whether a particle will stick
*/
struct cluster{
	vector<particle> particles;
	double bc_in; double bc_on;//the boundary conditions in and on the cluster indices
	Eigen::Vector3d r_com;
	//constructor and destructor below:
	cluster();
	cluster(vector<particle> gparticles, double gbcin, double gbcon, double dx, double dy, double dz);
	cluster(const cluster &o);
	~cluster(){}
	cluster operator=(const cluster &o);
	//member functions of the cluster:
	bool is_colliding(const particle & op, double bondlength);
	double r_gyr(double dx, double dy, double dz);
	double fractal_dimension();
	int find_num_neighbors(const particle & op, double bondlength);
};
//we define the default and parametric constructors for cluster()
cluster::cluster(): bc_in(0.0), bc_on(0.0), r_com(0.0,0.0,0.0){}
cluster::cluster(vector<particle> gparticles, double gbcin, double gbcon, double dx, double dy, double dz): bc_in(gbcin), bc_on(gbcon), r_com(0.0,0.0,0.0){
	particles = gparticles;
	r_com = center_of_mass(gparticles, dx, dy, dz);
}
//we will now define the copy constructor
cluster::cluster(const cluster & o): bc_in(0.0), bc_on(0.0), r_com(0.0,0.0,0.0){
	particles = o.particles;
	bc_in = o.bc_in;
	bc_on = o.bc_on;
	r_com = o.r_com;
}
//then we define the equality(=) operator
cluster cluster::operator=(const cluster & o){
	particles = o.particles;
	bc_in = o.bc_in;
	bc_on = o.bc_on;
	r_com = o.r_com;
	return *this;
}
//1) is_colliding() will test if a particle is close to any of the frozen cluster particles allowing for the new particle to be inserted into the cluster
bool cluster::is_colliding(const particle & op, double bondlength){
	vector<particle>::iterator it; it = particles.begin();
	bool retval = false;
	while(it != particles.end()){
		particle newp; newp = *it;
		double cblx = ((double)(op.position(0))) - ((double)(newp.position(0)));
		double cbly = ((double)(op.position(1))) - ((double)(newp.position(1)));
		double cblz = ((double)(op.position(2))) - ((double)(newp.position(2)));
		Eigen::Vector3d dr(cblx,cbly,cblz);
		//the particle must be less than or equal to bond length distance away from the cluster
		if( islessequal((abs((dr.norm()))),bondlength) != false ){retval = true; break;}
		else{++it;}
	}
	return retval;
}

//2)r_gyr() finds out the radius of gyration
double cluster::r_gyr(double dx, double dy, double dz){
	double Rgyr = radius_of_gyration(particles, r_com, dx, dy, dz);
	return Rgyr;
}

//3)fractal_dimension() finds out the fractal dimension of the cluster, which may have fractality or multifractality

//4)We define a random integer generator
int GenerateRandomInteger(int imin, int imax){
	int range_from  = imin;
	int range_to    = imax;
	std::random_device rand_dev;
	std::mt19937 generator(rand_dev());
	std::uniform_int_distribution<int> distr(range_from, range_to);
	int cointoss = distr(generator);
	return cointoss;
}
//5)We define a random double generator
double GenerateRandomDouble(double imin, double imax){
	double range_from  = imin;
	double range_to    = imax;
	std::random_device rand_dev;
	std::mt19937 generator(rand_dev());
	std::uniform_real_distribution<double> distr(range_from, range_to);
	double cointoss = distr(generator);
	return cointoss;
}
//6)A function to rectify position indices
//a function to rectify the index position
Eigen::Vector3i rectify_position(Eigen::Vector3i ri, int Nx, int Ny, int Nz){
	Eigen::Vector3i r; r = ri;
	Eigen::Vector3i rmax(Nx,Ny,Nz);
	int i = 0;
	while(i < 3){
		int j = i;
		i = i + 1;
		if( (r(j)) == (rmax(j)) ){ r(j) = 0; }
		else if( (r(j)) == -1 ){ r(j) = ((rmax(j)) - 1); }
	}
	return r;
}

//7) A function to find neighbors given a collision exists
int cluster::find_num_neighbors(const particle & op, double bondlength){
	vector<particle>::iterator it; it = particles.begin();
	particle collisionp;
	while(it != particles.end()){
		particle newp; newp = *it;
		double cblx = ((double)(op.position(0))) - ((double)(newp.position(0)));
		double cbly = ((double)(op.position(1))) - ((double)(newp.position(1)));
		double cblz = ((double)(op.position(2))) - ((double)(newp.position(2)));
		Eigen::Vector3d dr(cblx,cbly,cblz);
		//the particle must be less than or equal to bond length distance away from the cluster
		if( islessequal((abs((dr.norm()))),bondlength) != false ){collisionp = newp; break;}
		else{++it;}
	}
	//we now have the colliding particle, we want to find out it's neighbors
	it = particles.begin(); int N = 0;
	while(it != particles.end()){
		particle newp; newp = *it;
		double cblx = ((double)(op.position(0))) - ((double)(newp.position(0)));
		double cbly = ((double)(op.position(1))) - ((double)(newp.position(1)));
		double cblz = ((double)(op.position(2))) - ((double)(newp.position(2)));
		Eigen::Vector3d dr(cblx,cbly,cblz);
		//the particle must be less than or equal to 1 lattice distance away from the cluster
		if( islessequal((abs((dr.norm()))),bondlength) != false ){N = N + 1; ++it; }
		else{++it;}
	}
	N = N - 1;//the reason being that the above loop will always count the particle itself so we have to subtract that
	return N;
}

#endif
