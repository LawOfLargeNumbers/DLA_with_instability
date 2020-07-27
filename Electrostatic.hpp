#ifndef DFELECTROSTATIC_HPP
#define DFELECTROSTATIC_HPP
#include "./CustomLaplace.hpp"

using namespace std;
using namespace arma;

//We first define a function that iterates through a tensor lattice object and returns the net electric field at some point in a field tensor
double ElectricFieldAtPoint(int i, int j, int k, cluster &c, double dx, double dy, double dz, Eigen::Vector3d epsconst){
	Eigen::Vector3d r_pos; r_pos(0) = (dx*((double)i)); r_pos(1) = (dy*((double)j)); r_pos(2) = (dz*((double)k));
	vector<particle>::iterator it; it = c.particles.begin();
	double Eret = 0.0;
	while(it != c.particles.end()){
		particle p_curr; p_curr = *it;
		Eigen::Vector3d r_particle; r_particle(0) = (dx*((double)(p_curr.position(0)))); r_particle(1) = (dy*((double)(p_curr.position(1)))); r_particle(2) = (dz*((double)(p_curr.position(2))));
		Eigen::Vector3d dr; dr = r_pos - r_particle;
		double Ex, Ey, Ez, Enet;
		Ex = ((p_curr.charge)/(4.0*3.1415962*(epsconst(0))*(dr(0)))); Ey = ((p_curr.charge)/(4.0*3.1415962*(epsconst(1))*(dr(1)))); Ez = ((p_curr.charge)/(4.0*3.1415962*(epsconst(2))*(dr(2))));
		Enet = sqrt( ( (Ex*Ex)+(Ey*Ey)+(Ez*Ez) ) );
		Eret = Eret + Enet;
		++it;
	}
	return Eret;
}
//This is the method to get the analytical solution of the electric field at any given point in a lattice
cube ElectricFieldTensor(cluster &c, double dx, double dy, double dz, int Nx, int Ny, int Nz, Eigen::Vector3d epsconst){
	cube F_q = randu<cube>(Nx,Ny,Nz);//in the beginning the entire field tensor is filled with a random tensor
	//we need to rectify the field along the particles
	int ii = 0;
	while(ii < Nx){
		int jj = 0;
		while(jj < Ny){
			int kk = 0;
			while(kk < Nz){
				if( indices_are_in_cluster(ii,jj,kk, c.particles) ){ F_q(ii,jj,kk) = 0.0; kk = kk + 1; }//we assign the field to be 0 at the surface of the crystal
				else{
					F_q(ii,jj,kk) = 1.0;
					kk = kk + 1;
				}
			}
			jj = jj + 1;
		}
		ii = ii + 1;
	}
	//once the values are rectified along the cluster and on the box boundaries we can iterate to generate the field at each point
	int i = 0;
	while(i < Nx){
		int j = 0;
		while(j < Ny){
			int k = 0;
			while(k < Nz){
				if( indices_are_in_cluster(i,j,k, c.particles) ){ k = k + 1; }//we will not find out the value for the points on the cluster
				else{
					double Ec = ElectricFieldAtPoint(i, j, k, c, dx, dy, dz, epsconst);
					F_q(i,j,k) = Ec;
					k = k + 1;
				}
			}
			j = j + 1;
		}
		i = i + 1;
	}
	//once this calculation is done we will return the tensor
	return F_q;
}

//A function to calculate the probability of a translation
double P_electrostatic(cube fieldtensor, double charge_curr, Eigen::Vector3i ru_o, Eigen::Vector3i ru_new, int Nx, int Ny, int Nz){
	//We need to rectify each of the positions
	Eigen::Vector3i ineg, ipos, jneg, jpos, kneg, kpos, r_o, r_new;
	r_o = rectify_position(ru_o, Nx, Ny, Nz); r_new = rectify_position(ru_new, Nx, Ny, Nz);
	int i,j,k;
	i = r_o(0); j = r_o(1); k = r_o(2);
	Eigen::Vector3i indx_ineg((i-1),j,k); Eigen::Vector3i indx_ipos((i+1),j,k); Eigen::Vector3i indx_jneg(i,(j-1),k);
	Eigen::Vector3i indx_jpos(i,(j+1),k); Eigen::Vector3i indx_kneg(i,j,(k-1)); Eigen::Vector3i indx_kpos(i,j,(k+1));
	ineg = rectify_position(indx_ineg, Nx, Ny, Nz); ipos = rectify_position(indx_ipos, Nx, Ny, Nz);
	jneg = rectify_position(indx_jneg, Nx, Ny, Nz); jpos = rectify_position(indx_jpos, Nx, Ny, Nz);
	kneg = rectify_position(indx_kneg, Nx, Ny, Nz); ipos = rectify_position(indx_kpos, Nx, Ny, Nz);
	//we need to find all the field values around the previous point
	double E_ineg, E_ipos, E_jneg, E_jpos, E_kneg, E_kpos, E_ijk, E_newpos;
	E_ineg = fieldtensor(ineg(0),ineg(1),ineg(2)); E_ipos = fieldtensor(ipos(0),ipos(1),ipos(2));
	E_jneg = fieldtensor(jneg(0),jneg(1),jneg(2)); E_jpos = fieldtensor(jpos(0),jpos(1),jpos(2));
	E_kneg = fieldtensor(kneg(0),kneg(1),kneg(2)); E_kpos = fieldtensor(kpos(0),kpos(1),kpos(2));
	E_ijk = fieldtensor(i,j,k); E_newpos = fieldtensor((r_new(0)), (r_new(1)), (r_new(2)));
	//normalize the difference
	double dE_ineg, dE_ipos, dE_jneg, dE_jpos, dE_kneg, dE_kpos;
	dE_ineg = (abs((E_ineg - E_ijk))); dE_ipos = (abs((E_ipos - E_ijk)));
	dE_jneg = (abs((E_jneg - E_ijk))); dE_jpos = (abs((E_jpos - E_ijk)));
	dE_kneg = (abs((E_kneg - E_ijk))); dE_kpos = (abs((E_kpos - E_ijk)));
	//we then find the absolute sum of the fields multiplied by the absolute value of charge to give us the potential
	double dV_net = ((abs(charge_curr))*(dE_ineg + dE_ipos + dE_jneg + dE_jpos + dE_kneg + dE_kpos));
	double exponent_value = ((charge_curr*(E_newpos - E_ijk))/dV_net);
	//when raised to the negative of the calculated exponent we can obtain a measure of the probability, based on a gaussian local kernel
	double Pq = (exp(((-1.0)*exponent_value)));
	return Pq;
}
#endif
