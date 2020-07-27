#include "./CustomLaplace.hpp"
#include <chrono>

using namespace std;
using namespace arma;


//we define a function that initiates the position along one of the borders
Eigen::Vector3i initiate_position(int Nx, int Ny, int Nz){
	Eigen::Vector3i retval; Eigen::Vector3i maxvals(Nx,Ny,Nz);
	Eigen::Vector3d ctvals;
	ctvals(0) = GenerateRandomDouble(0.0, 1.0); ctvals(1) = GenerateRandomDouble(0.0, 1.0); ctvals(2) = GenerateRandomDouble(0.0, 1.0);
	int i = 0;
	while(i < 3){
		if(ctvals(i) < 0.5){ retval(i) = 0; }
		else{ retval(i) = GenerateRandomInteger(0,(maxvals(i)-1)); }
		i = i +1;
	}
	//if neither is zero then one of them has to be randomly set to zero
	//otherwise the initial position is not on the edges or boundaries of the box
	if( (retval(0) != 0) && (retval(1) != 0) && (retval(2) != 0) ){
		double dd = GenerateRandomDouble(0,3);
		if( dd < 1.0 ){ retval(0) = 0; return retval; }
		else if( (islessequal(1.0,dd)!=false) && (dd < 2.0) ){ retval(1) = 0; return retval; }
		else{ retval(2) = 0; return retval; }
	}
	else{ return retval; }
	
}
//a function to gnerate a random displacement in 3D
Eigen::Vector3i generate_displacement_3d(int dx, int dy, int dz){
	Eigen::Vector3i retval;
	retval(0) = GenerateRandomInteger(0,dx); retval(1) = GenerateRandomInteger(0,dy); retval(2) = GenerateRandomInteger(0,dz);
	//if all three zero then one of them has to be randomly set to 1
	if( (retval(0) == 0) && (retval(1) == 0) && (retval(2) == 0) ){
		double dd = GenerateRandomDouble(0,3);
		if( dd < 1.0 ){ retval(0) = dx; return retval; }
		else if( (islessequal(1.0,dd)!=false) && (dd < 2.0) ){ retval(1) = dy; return retval; }
		else{ retval(2) = dz; return retval; }
	}
	return retval;
}
//the commonly used displacement function does not have diagonal moves
Eigen::Vector3i generate_displacement(int dx, int dy, int dz){
	Eigen::Vector3i retval;
	retval(0) = 0; retval(1) = 0; retval(2) = 0;//initiate all the moves at zero
	//change any one of the indices to 1 initiating a 1 dimensional motion
	double dd = GenerateRandomDouble(0,3);
	if( dd < 1.0 ){ retval(0) = dx; return retval; }
	else if( (islessequal(1.0,dd)!=false) && (dd < 2.0) ){ retval(1) = dy; return retval; }
	else{ retval(2) = dz; return retval; }
}

//4) We will define functions to create clusters of some specified symmetry
//4.a) Spherical Geometry
vector<particle> create_spherical_cluster(double radius, particle &clusterparticle, double dx, double dy, double dz, int Nx, int Ny, int Nz){
	vector<particle> returnedparticles; cluster tempcluster;
	particle newparticleB; newparticleB.position(0) = (Nx/2) - 1; newparticleB.position(1) = (Ny/2) - 1; newparticleB.position(2) = (Nz/2) - 1;
	tempcluster.particles.push_back(newparticleB);
	int i = 0;
	while( i < Nx ){
		int j = 0;
		while( j < Ny ){
			int k = 0;
			while( k < Nz ){
				double x = (dx*((double)i)); double y = (dy*((double)j)); double z = (dz*((double)k));
				double x0 = (dx*((double)((Nx/2) - 1))); double y0 = (dy*((double)((Ny/2) - 1))); double z0 = (dz*((double)((Nz/2) - 1)));
				double dr = (((x-x0)*(x-x0)) + ((y-y0)*(y-y0)) + ((z-z0)*(z-z0)));
				double crit = (abs(dr - radius));
				if( islessequal(crit,1.0) ){
					particle newparticle; newparticle = clusterparticle;
					newparticle.position(0) = i; newparticle.position(1) = j; newparticle.position(2) = k;
					Eigen::Vector3i gapdist(1,1,1);//one bond length apart
					//check if the position is in cluster colliding
					if( (i != ((Nx/2) - 1)) && (j != ((Ny/2) - 1)) && (k != ((Nz/2) - 1)) ){
						returnedparticles.push_back(newparticle);
						k = k + 1;
					}
					else{ k = k + 1; }
					//else nothing is done
				}
				else{ k = k + 1; }
			}
			j = j + 1;
		}
		i = i + 1;
	}
	return returnedparticles;
}
//4.b) Ellipse Geometry
vector<particle> create_ellipsoidal_cluster(double a, double b, double c, particle &clusterparticle, double dx, double dy, double dz, int Nx, int Ny, int Nz){
	vector<particle> returnedparticles; cluster tempcluster;
	particle newparticleB; newparticleB.position(0) = (Nx/2) - 1; newparticleB.position(1) = (Ny/2) - 1; newparticleB.position(2) = (Nz/2) - 1;
	tempcluster.particles.push_back(newparticleB);
	int i = 0;
	while( i < Nx ){
		int j = 0;
		while( j < Ny ){
			int k = 0;
			while( k < Nz ){
				double x = (dx*((double)i)); double y = (dy*((double)j)); double z = (dz*((double)k));
				double x0 = (dx*((double)((Nx/2) - 1))); double y0 = (dy*((double)((Ny/2) - 1))); double z0 = (dz*((double)((Nz/2) - 1)));
				double dr_x = ((x-x0)*(x-x0)); double dr_y = ((y-y0)*(y-y0)); double dr_z = ((z-z0)*(z-z0));
				double crit = ((dr_x/(a*a)) + (dr_y/(b*b)) + (dr_z/(c*c)));
				//the criterea is x^2/a^2 + y^2/b^2 + z^2/c^2 <= 1
				if( islessequal(crit,1.0) ){
					particle newparticle; newparticle = clusterparticle;
					newparticle.position(0) = i; newparticle.position(1) = j; newparticle.position(2) = k;
					Eigen::Vector3i gapdist(1,1,1);//one bond length apart
					//check if the position is in cluster colliding
					if( (i != ((Nx/2) - 1)) && (j != ((Ny/2) - 1)) && (k != ((Nz/2) - 1)) ){
						returnedparticles.push_back(newparticle);
						k = k + 1;
					}
					else{ k = k + 1; }
					//else nothing is done
				}
				else{ k = k + 1; }
			}
			j = j + 1;
		}
		i = i + 1;
	}
	return returnedparticles;
}
//4. c) Cuboidal Geometry
vector<particle> create_cuboidal_cluster(double lx, double ly, double lz, particle &clusterparticle, double dx, double dy, double dz, int Nx, int Ny, int Nz){
	vector<particle> returnedparticles; cluster tempcluster;
	particle newparticleB; newparticleB.position(0) = (Nx/2) - 1; newparticleB.position(1) = (Ny/2) - 1; newparticleB.position(2) = (Nz/2) - 1;
	tempcluster.particles.push_back(newparticleB);
	int i = 0;
	while( i < Nx ){
		int j = 0;
		while( j < Ny ){
			int k = 0;
			while( k < Nz ){
				double x = (dx*((double)i)); double y = (dy*((double)j)); double z = (dz*((double)k));
				double x0 = (dx*((double)((Nx/2) - 1))); double y0 = (dy*((double)((Ny/2) - 1))); double z0 = (dz*((double)((Nz/2) - 1)));
				double dr_x = ((x-x0)*(x-x0)); double dr_y = ((y-y0)*(y-y0)); double dr_z = ((z-z0)*(z-z0));
				//the criterea is x^2/a^2 + y^2/b^2 + z^2/c^2 <= 1
				if( (islessequal(dr_x,(lx/2.0)) != false) && (islessequal(dr_y,(ly/2.0)) != false) && (islessequal(dr_z,(lz/2.0)) != false) ){
					particle newparticle; newparticle = clusterparticle;
					newparticle.position(0) = i; newparticle.position(1) = j; newparticle.position(2) = k;
					Eigen::Vector3i gapdist(1,1,1);//one bond length apart
					//check if the position is in cluster colliding
					if( (i != ((Nx/2) - 1)) && (j != ((Ny/2) - 1)) && (k != ((Nz/2) - 1)) ){
						returnedparticles.push_back(newparticle);
						k = k + 1;
					}
					else{ k = k + 1; }
					//else nothing is done
				}
				else{ k = k + 1; }
			}
			j = j + 1;
		}
		i = i + 1;
	}
	return returnedparticles;
}
//A print function for the cluster
void print_cluster(cluster &o){
	vector<particle>::iterator it; it = o.particles.begin();
	string filename = "./trajectory_frames/Nucleation_Cluster."+to_string(o.particles.size())+".xyz";
	ofstream ofile; ofile.open(filename, ofstream::out|ofstream::app );
	ofile<<o.particles.size()<<endl;
	ofile<<endl;
	while(it != o.particles.end()){
		particle particle_current; particle_current = *it;
		ofile<<"H"<<'\t'<<particle_current.position(0)<<'\t'<<particle_current.position(1)<<'\t'<<particle_current.position(2)<<endl;
		++it;
	}
	ofile.close();
}
//A print function for the electric field
void print_field(cube &E_final, int N_particles, int Nx, int Ny, int Nz){
	string filename = "./field_frames/Nucleation_field."+to_string(N_particles)+".xyz";
	ofstream ofileb; ofileb.open(filename, ofstream::out|ofstream::app );
	//ofile<<"i"<<'\t'<<"j"<<'\t'<<"k"<<'\t'<<"F_arb(i,j,k)"<<endl;
	int a = 0;
	while(a < Nx){
		int b = 0;
		while(b < Ny){
			int c = 0;
			while(c < Nz){
				ofileb<<a<<'\t'<<b<<'\t'<<c<<'\t'<<E_final(a,b,c)<<endl;
				c = c+1; 
			}
			b = b + 1;
		}
		a = a + 1;
	}
	ofileb.close();
}

/*
This is successfully tested for a point particle seed, it seems to be working fine
*/
int main(){
	//now we initialize a cluster with phi = 0.0 in and on the boundary
	//we set up a new cluster
	cluster newcluster; newcluster.bc_in = 0.0; newcluster.bc_on = 0.0;
	double eps = 0.001;
	double bc_on = 1.0;//we define phi = 1.0 at infinity
	int Nx = 50;
	int Ny = 50;
	int Nz = 50;
	double Q_cluster = 1.0;
	double Q_insert = -1.0;
	Eigen::Vector3d epsconst(1.0,1.0,1.0);//the three dimensional dielectric constants along x, y and z
	int Nparticle_max = 300;//the total number of particles in the cluster
	int xstep = 1; int ystep = 1; int zstep = 1;//the steps in all three directions
	double bondlength = 1.0;//this is the perception to which the cluster can hold another particle
	double dx = 1.0; double dy = 1.0; double dz = 1.0;//these are the lattice spacings of the lattice
	//shape parameters
	int shape = 1;//0 = point, 1 = sphere, 2 = ellipse, 3 = cuboid
	double r = 10.0;//radius of a sphere
	double aa = 5.0; double bb = 5.0; double cc = 5.0;//a,b,c parameters of an ellipsoid
	double lx = 5.0; double ly = 5.0; double lz = 5.0;//x,y,z lengths of a cuboid
	vector<particle> seedparticles;
	particle seed_particle; seed_particle.position(0) = (Nx/2) - 1; seed_particle.position(1) = (Ny/2) - 1; seed_particle.position(2) = (Nz/2) - 1; seed_particle.charge = Q_cluster;//create the cluster particle
	if( shape == 0 ){
		//do not generate a nucleation cluster
		seedparticles.push_back(seed_particle);
		cout<<"Point cluster has been created!"<<endl;
	}
	else if( shape == 1 ){
		//generate a spherical nucleation cluster
		seedparticles = create_spherical_cluster(r, seed_particle, dx, dy, dz, Nx, Ny, Nz);
		cout<<"Spherical Cluster of radius "<<r<<" lattice units has been created!"<<endl;
	}
	else if( shape == 2 ){
		//generate an ellipsoidal nucleation cluster
		seedparticles = create_ellipsoidal_cluster(aa, bb, cc, seed_particle, dx, dy, dz, Nx, Ny, Nz);
		cout<<"Ellipsoidal cluster with a = "<<aa<<" lattice units, b = "<<bb<<" lattice units, c = "<<cc<<" lattice units, has been created!"<<endl;
	}
	else if( shape == 3 ){
		//generate a cuboidal nucleation cluster
		seedparticles = create_cuboidal_cluster(lx, ly, lz, seed_particle, dx, dy, dz, Nx, Ny, Nz);
		cout<<"Cuboidal cluster with lx = "<<lx<<" lattice units, ly = "<<ly<<" lattice units, lz = "<<lz<<" lattice units, has been created!"<<endl;
	}
	else{
		//do not generate a nucleation cluster
		seedparticles.push_back(seed_particle);
		cout<<"Point cluster has been created!"<<endl;
	}
	//insert the seed particles and calculate the initial quantities
	vector<particle>::iterator partit; partit = seedparticles.begin();
	while(partit != seedparticles.end()){
		particle currparticle; currparticle = *partit;
		newcluster.particles.push_back(currparticle);
		++partit;
	}
	newcluster.r_com = center_of_mass(newcluster.particles, dx, dy, dz);//this will be updated upon insertion
	double R_gyration = newcluster.r_gyr(dx, dy, dz);//this will be updated upon insertion
	cout<<"New cluster created"<<endl;
	//initialize a field object
	auto start = std::chrono::high_resolution_clock::now();//start measuring time from here
	cube F_arb = laplace_solver_3d(newcluster, eps, bc_on, newcluster.bc_on, Nx, Ny, Nz);//initially all points are at 0.0
	auto finish = std::chrono::high_resolution_clock::now();//end measuring time here
	//F_arb.update_field(newcluster, 1.0);//the new cluster is inserted and boundary is maintained at phi = 1.0
	//once updated we print the laplacian field:
	ofstream ofile; ofile.open("test_laplace.dat");
	ofile<<"i"<<'\t'<<"j"<<'\t'<<"k"<<'\t'<<"F_arb(i,j,k)"<<endl;
	int a = 0;
	while(a < Nx){
		int b = 0;
		while(b < Ny){
			int c = 0;
			while(c < Nz){
				ofile<<a<<'\t'<<b<<'\t'<<c<<'\t'<<F_arb(a,b,c)<<endl;
				c = c+1; 
			}
			b = b + 1;
		}
		a = a + 1;
	}
	ofile.close();
	
	std::chrono::duration<double> elapsed = finish - start;//the time recorded
	cout << "Elapsed time: " << elapsed.count() << " s\n";//report the time
	return 0;
}
