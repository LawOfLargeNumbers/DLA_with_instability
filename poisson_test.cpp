#include "./CustomLaplace.hpp"
#include <chrono>

using namespace std;
using namespace arma;

/*
This is successfully tested for a point particle seed, it seems to be working fine
*/
int main(){
	//now we initialize a cluster with phi = 0.0 in and on the boundary
	cluster newcluster; newcluster.bc_in = 0.0; newcluster.bc_on = 0.0;
	double eps = 0.0001;
	double bc_on = 1.0;//we define phi = 1.0 at infinity
	int Nx = 50;
	int Ny = 50;
	int Nz = 50;
	double electric_rho = 0.5;//the charge density in the space
	double epsconst = 1.0;//the permittivity of space
	cout<<"parameters set up!"<<endl;
	//initialize a lattice object
	//Eigen::Tensor<int,3> lattice(1000,1000,1000); lattice.setZero();//a 3D 1000 cross dimensional lattice with 1000000000 points
	//initialize the seed particle at the middle of the lattice
	particle seed_particle; seed_particle.position(0) = (Nx/2) - 1; seed_particle.position(1) = (Ny/2) - 1; seed_particle.position(2) = (Nz/2) - 1;
	cout<<"Seed particle created"<<endl;
	//put the seed particle in the lattice and cluster
	//lattice(499,499,499) = 1;
	newcluster.particles.push_back(seed_particle);
	cout<<"New cluster created"<<endl;
	//initialize an external field object
	cube F_external = randu<cube>(Nx,Ny,Nz);
	int m = 0;
	while(m < Nx){
		int mm = 1;
		while(mm < Ny){
			int mmm = 1;
			while(mmm < Nz){
				F_external(m,mm,mmm) = (electric_rho/epsconst);//f = rho/epsilon definition is being applied
				mmm = mmm + 1;
			}
			mm = mm + 1;
		}
		m = m + 1;
	}
	cout<<"External electric field with field density "<<electric_rho<<", and permittivity "<<epsconst<<", has been created!"<<endl;
	//initialize a field object
	auto start = std::chrono::high_resolution_clock::now();//start measuring time from here
	cube F_arb = poisson_solver_3d(newcluster, F_external, eps, bc_on, newcluster.bc_on, Nx, Ny, Nz);//initially all points are at 0.0
	auto finish = std::chrono::high_resolution_clock::now();//end measuring time here
	//F_arb.update_field(newcluster, 1.0);//the new cluster is inserted and boundary is maintained at phi = 1.0
	//once updated we print the laplacian field:
	ofstream ofile; ofile.open("test_poisson.dat");
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
