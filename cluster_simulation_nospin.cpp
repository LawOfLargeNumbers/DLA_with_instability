#include "./Electrostatic.hpp"

using namespace std;
using namespace arma;


//we define a function that initiates the position along one of the borders
Eigen::Vector3i initiate_position_box(int Nx, int Ny, int Nz){
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
//We define a better function using border planes
Eigen::Vector3i initiate_position(int Nx, int Ny, int Nz){
	vector<Eigen::Vector3i> starterpositions;
	//first the xy planes
	int i = 0;
	while( i < Nx ){
		 int j = 0;
		 while( j < Ny ){
		 	Eigen::Vector3i va(i,j,0); Eigen::Vector3i vb(i,j,(Nz - 1));
		 	starterpositions.push_back(va); starterpositions.push_back(vb);
		 	j = j + 1;
		 }
		 i = i + 1;
	}
	//second xz planes
	i = 0;
	while( i < Nx ){
		 int k = 0;
		 while( k < Nz ){
		 	Eigen::Vector3i va(i,0,k); Eigen::Vector3i vb(i,(Ny - 1),k);
		 	starterpositions.push_back(va); starterpositions.push_back(vb);
		 	k = k + 1;
		 }
		 i = i + 1;
	}
	//third yz planes
	int j = 0;
	while( j < Ny ){
		 int k = 0;
		 while( k < Nz ){
		 	Eigen::Vector3i va(0,j,k); Eigen::Vector3i vb(0,j,(Nz - 1));
		 	starterpositions.push_back(va); starterpositions.push_back(vb);
		 	k = k + 1;
		 }
		 j = j + 1;
	}
	//next we generate a random integer
	int randindx = GenerateRandomInteger(0,(starterpositions.size()));
	Eigen::Vector3i retvec; retvec = starterpositions[randindx];
	return retvec;
}
//a function to gnerate a random displacement in 3D
Eigen::Vector3i generate_displacement_3d(int dx, int dy, int dz){
	Eigen::Vector3i retval;
	retval(0) = GenerateRandomInteger((-1.0*dx),dx); retval(1) = GenerateRandomInteger((-1.0*dy),dy); retval(2) = GenerateRandomInteger((-1.0*dz),dz);
	//UPDATE: The particle may choose to not move at all and thats fine really...
	//if all three zero then one of them has to be randomly set to 1
	/*
	if( (retval(0) == 0) && (retval(1) == 0) && (retval(2) == 0) ){
		double dd = GenerateRandomDouble(0,3);
		if( dd < 0.5 ){ retval(0) = dx; return retval; }
		else if( islessequal(0.5,dd)!=false && (dd < 1.0) ){ retval(0) = (-1.0)*dx; return retval; }
		else if( (islessequal(1.0,dd)!=false) && (dd < 2.0) ){ retval(1) = dy; return retval; }
		else{ retval(2) = dz; return retval; }
	}
	*/
	return retval;
}
//the commonly used displacement function does not have diagonal moves
Eigen::Vector3i generate_displacement(int dx, int dy, int dz){
	Eigen::Vector3i retval;
	retval(0) = 0; retval(1) = 0; retval(2) = 0;
	//generate a stencil motion walker with a +/- motion in each direction
	double dd = GenerateRandomDouble(0,3);
	if( dd < 1.0 ){ retval(0) = GenerateRandomInteger((-1.0*dx),dx); return retval; }
	else if( (islessequal(1.0,dd)!=false) && (dd < 2.0) ){ retval(1) = GenerateRandomInteger((-1.0*dy),dy); return retval; }
	else{ retval(2) = GenerateRandomInteger((-1.0*dz),dz); return retval; }
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
//the main functions:
int main(){
	//we set up a new cluster
	cluster newcluster; newcluster.bc_in = 0.0; newcluster.bc_on = 0.0;
	double eps = 0.00001; double t = 0.5;//base of instability probability
	double bc_on = 1.0;//we define phi = 1.0 at infinity
	int Nx = 50;
	int Ny = 50;
	int Nz = 50;
	double Q_cluster = 1.0;
	double Q_insert = -1.0;
	Eigen::Vector3d epsconst(1.0,1.0,1.0);//the three dimensional dielectric constants along x, y and z
	int Nparticle_max = 3000;//the total number of particles in the cluster
	int xstep = 1; int ystep = 1; int zstep = 1;//the steps in all three directions
	double bondlength = 1.0;//this is the perception to which the cluster can hold another particle
	double dx = 1.0; double dy = 1.0; double dz = 1.0;//these are the lattice spacings of the lattice
	//shape parameters
	int shape = 1;//0 = point, 1 = sphere, 2 = ellipse, 3 = cuboid
	double r = 5.0;//radius of a sphere
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
	cout<<"The number of seed particles in the system is "<<seedparticles.size()<<endl;
	cout<<"parameters set up!"<<endl;
	//insert the seed particles and calculate the initial quantities
	vector<particle>::iterator partit; partit = seedparticles.begin();
	while(partit != seedparticles.end()){
		particle currparticle; currparticle = *partit;
		newcluster.particles.push_back(currparticle);
		++partit;
	}
	newcluster.r_com = center_of_mass(newcluster.particles, dx, dy, dz);//this will be updated upon insertion
	double R_gyration = newcluster.r_gyr(dx, dy, dz);//this will be updated upon insertion
	//create a new mobile particle
	particle mobile_particle; mobile_particle.charge = Q_insert; mobile_particle.position = initiate_position(Nx,Ny,Nz);//this initiates the mobile particle
	//a temporary particle allows for the manipulation of the mobile particle without touching it physically:
	particle temporary_particle; temporary_particle = mobile_particle;//the mobile particle is the one that will move around
	cout<<"System set up!"<<endl;
	cube E_final = randu<cube>(Nx,Ny,Nz);
	//set up a vector for storing the step vs time difference
	vector<tuple<int,int,double>> timeseries;
	//Initiate the field based on the seed cluster
	cube E_tensor = laplace_solver_3d(newcluster, eps, bc_on, newcluster.bc_on, Nx, Ny, Nz);
	//now we can initate a growth loop
	int Nparticle_current = newcluster.particles.size(); int mainstep = 0; int mainstep_old = 0;
	ofstream logfile; logfile.open("cluster_nucleus.log");
	while(Nparticle_current < (Nparticle_max + 1)){
		logfile<<"+++++++++++++++++++++++++++++++++++"<<endl;
		logfile<<"Net Step Number "<<mainstep<<endl;
		logfile<<"+++++++++++++++++++++++++++++++++++"<<endl;
		mainstep = mainstep + 1;
		//we generate a random displacement for the particle:
		Eigen::Vector3i disp; disp = generate_displacement(xstep, ystep, zstep);
		//calculate the Electric field tensor
		//cube E_tensor = laplace_solver_3d(newcluster, eps, bc_on, newcluster.bc_on, Nx, Ny, Nz);
		//E_tensor = ElectricFieldTensor(newcluster, dx, dy, dz, Nx, Ny, Nz, epsconst);//this is the electric field due to the current state of the cluster
		//Print out some quantities to a log file:
		logfile<<"Electric field tensor meanvalue "<<tensor_mean_value(E_tensor, eps, bc_on, newcluster.bc_on, Nx, Ny, Nz)<<endl;//print the mean value of the electric field
		logfile<<"Center of mass of the cluster is at "<<newcluster.r_com<<endl;//print the center of mass
		logfile<<"N_particle = "<<Nparticle_current<<endl;//print the number of particles
		logfile<<"R_gyration = "<<R_gyration<<endl;//print the radius of gyration
		//now using this tensor we calculate the probability of the translation
		Eigen::Vector3i newpos; newpos = ( (temporary_particle.position) + disp);
		//cout<<"displacement = "<<disp<<" old position = "<<temporary_particle.position<<" new position = "<<newpos<<endl;
		newpos = rectify_position(newpos, Nx, Ny, Nz);//the position is rectified based on the reappearance geometry of the periodic box
		//cout<<"rectified position = "<<newpos<<endl;
		//also calculate the electric field potential at some point
		//double Energy_E_new, Energy_E_old;
		//Energy_E_new = ElectricFieldAtPoint(newpos(0), newpos(1), newpos(2), newcluster, dx, dy, dz, epsconst);
		//Energy_E_old = ElectricFieldAtPoint(temporary_particle.position(0), temporary_particle.position(1), temporary_particle.position(2), newcluster, dx, dy, dz, epsconst);
		//double dE = Energy_E_new - Energy_E_old;
		//double P_q = exp((-1.0*dE));
		double P_phi = exp((-1.0*((E_tensor(newpos(0),newpos(1),newpos(2))) - (E_tensor(temporary_particle.position(0),temporary_particle.position(1), temporary_particle.position(2))))));
		double P = ((1.0/6.0)*P_phi);//we have to multiply the field-probability with the generation probability
		//double P = P_electrostatic(E_tensor, (temporary_particle.charge), (temporary_particle.position), newpos, Nx, Ny, Nz);
		//toss a coin and check if we can accept:
		double tossacoin = GenerateRandomDouble(0.0, 1.0);
		if( islessequal(tossacoin,P) ){
			temporary_particle.position = newpos;//the teporary particle is displaced
			//now we check if the particle is in a position to be inserted into the cluster
			if( (newcluster.is_colliding(temporary_particle, bondlength)) != false ){
				int mainstep_new = mainstep;
				//we need to find out how many neighbors there are
				double N_neighbors = (double)(newcluster.find_num_neighbors(temporary_particle, 1.155));
				double s = pow(t,(13.0-N_neighbors));
				//now toss the coin again to see if we will accept the move
				tossacoin = GenerateRandomDouble(0.0, 1.0);
				//if the criterea is met then accept and insert the particle
				if( islessequal(tossacoin,s) ){
					newcluster.particles.push_back(temporary_particle);//the new particle is pushed back
					newcluster.r_com = center_of_mass(newcluster.particles, dx, dy, dz);//the new center of mass is calculated
					R_gyration = newcluster.r_gyr(dx, dy, dz);
					//now we reset the position of the temporary particle
					temporary_particle.position = initiate_position(Nx,Ny,Nz);//a new initial position is assigned
					//the particle count is updated
					E_final = E_tensor;//put in the new field
					Nparticle_current = Nparticle_current + 1;
					int dT_insert = mainstep_new - mainstep_old;
					logfile<<"+++++INSERTION++++ Inserted Particle number "<<Nparticle_current<<endl;
					logfile<<"dT_insert = "<<dT_insert<<endl;
					mainstep_old = mainstep_new;//maintain a record of this step so that we can calculate the time of insertion again
					//create a timeseries object
					tuple<int,int,double> timepair; timepair = make_tuple(Nparticle_current, dT_insert, R_gyration); timeseries.push_back(timepair);
					//update the field object
					E_tensor = laplace_solver_3d(newcluster, eps, bc_on, newcluster.bc_on, Nx, Ny, Nz);//this reduces computational load since the field will remain constant unless the cluster is updated
					//finally we will print the new cluster and the new field
					print_cluster(newcluster);
					print_field(E_tensor,(newcluster.particles.size()), Nx, Ny, Nz);
				}
				else{
					Nparticle_current = Nparticle_current + 0;//we continue the iteration
				}
			}
			else{
				Nparticle_current = Nparticle_current + 0;//we continue the iteration
			}
		}
		else{
			Nparticle_current = Nparticle_current + 0;//the temporary particle is not displaced
		}
		logfile<<"++++++++++++++++++++++++++++++++++++"<<endl;
		logfile<<endl;
		cout<<"|Step number "<<mainstep<<" | Probability of acceptance is = "<<P<<", Cointoss is = "<<tossacoin<<", Number of particles in system = "<<Nparticle_current<<endl;
	}
	cout<<endl;
	cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
	cout<<"Code Finished!"<<endl;
	logfile.close();
	//now we write out the particles
	//print out the final electric field
	
	//print out the time series
	ofstream ofilec; ofilec.open("Nucleation_insertiontime.dat");
	//ofile<<"i"<<'\t'<<"j"<<'\t'<<"k"<<'\t'<<"F_arb(i,j,k)"<<endl;
	vector<tuple<int,int,double>>::iterator titer; titer = timeseries.begin();
	while(titer != timeseries.end()){
		tuple<int,int,double> tcurr; tcurr = *titer;
		ofilec<<get<0>(tcurr)<<'\t'<<get<1>(tcurr)<<'\t'<<get<2>(tcurr)<<endl;
		++titer;
	}
	ofilec.close();
	cout<<endl;
	cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
	cout<<"FILES HAVE BEEN PRINTED!"<<endl;
	return 0;
}
