#include "./datatypes.hpp"
#include <armadillo>

using namespace std;
using namespace arma;

//first we define a method to test if a given index is in a cluster
bool indices_are_in_cluster(int i, int j, int k, vector<particle> gparticles){
	vector<particle>::iterator it; it = gparticles.begin();
	bool retval = false;
	while(it != gparticles.end()){
		particle currparticle; currparticle = *it;
		Eigen::Vector3i currindx; currindx = currparticle.position;
		if( (i == currindx(0)) && (j == currindx(1)) && (k == currindx(2)) ){
			retval = true;
			break;
		}
		else{
			retval = false;
			++it;
		}
	}
	return retval;
}
//figure out the mean value a tensor
double tensor_mean_value(cube latticecube, double eps, double bc_on, double cluster_bc_on, int Nx, int Ny, int Nz){
	double valsum = 0.0; double N = 0.0;
	int l = 0;
	while( l < Nx ){
		int m = 0;
		while( m < Ny ){
			int n = 0;
			while( n < Nz ){
				double valcurr = latticecube(l,m,n);
				if( (abs((valcurr-bc_on)) < 0.00000000001) || (abs((valcurr-cluster_bc_on)) < 0.00000000001) ){ n = n + 1; }
				else{ valsum = valsum + valcurr; N = N + 1.0; n = n + 1; }//only vales neither at the boundaries nor on the cluster, are counted
			}
			m = m + 1;
		}
		l = l + 1;
	}
	double valmean = (valsum/N) ;
	return valmean;
}
//figure out max value in a tensor
double tensor_max_value(cube latticecube, double eps, double bc_on, double cluster_bc_on, int Nx, int Ny, int Nz){
	double valsum = 0.0; double N = 0.0;
	int l = 0; double maxval = latticecube(0,0,0);
	while( l < Nx ){
		int m = 0;
		while( m < Ny ){
			int n = 0;
			while( n < Nz ){
				double valcurr = latticecube(l,m,n);
				if( (abs((valcurr-bc_on)) < 0.00000000001) || (abs((valcurr-cluster_bc_on)) < 0.00000000001) ){ n = n + 1; }
				else{
					if(valcurr > maxval){ maxval = valcurr; n = n + 1; }//assign a new maximum value
					else{ maxval = maxval; n = n + 1; }//no new maximum value
				}
			}
			m = m + 1;
		}
		l = l + 1;
	}
	double valmean = (valsum/N) ;
	return valmean;
}
//we will initalize a tensor and then compute the laplacian, we have a border boundary condition and a cluster boundary condition
cube laplace_solver_3d(cluster &provided_cluster, double eps, double bc_on, double cluster_bc_on, int Nx, int Ny, int Nz){
	cube F_old = randu<cube>(Nx,Ny,Nz);//we initialize a uniform distribution random cube
	cube F_new = randu<cube>(Nx,Ny,Nz);//we initialize another random tensor for the new field
	cube F_placeholder = randu<cube>(Nx,Ny,Nz);
	//cout<<"set up of laplacian is done"<<endl;
	//assign the boundaries their values
	int i = 0;
	while(i < Nx){
		int j = 0;
		while(j < Ny){
			int k = 0;
			while(k < Nz){
				if( i == 0 || j == 0 || k == 0 || i == (Nx-1) || j == (Ny-1) || k == (Nz-1)  ){
					//impose the condition for box boundary
					F_old(i,j,k) = bc_on; F_new(i,j,k) = bc_on; F_placeholder(i,j,k) = bc_on;
					k = k + 1;
				}
				else if( indices_are_in_cluster(i,j,k, provided_cluster.particles) ){
					//impose the boundary condition on the cluster
					F_old(i,j,k) = cluster_bc_on; F_new(i,j,k) = cluster_bc_on; F_placeholder(i,j,k) = cluster_bc_on;
					k = k + 1;
				}
				else{
					k = k + 1;
				}
			}
			j = j + 1;
		}
		i = i + 1;
	}
	//cout<<"boundary conditions are imposed"<<endl;
	//now that the cluster and box boundaries have been set, we can safely iterate till convergence.
	cube dF_tensor = F_new - F_old;
	double valmean = tensor_mean_value(dF_tensor, eps, bc_on, cluster_bc_on, Nx, Ny, Nz);
	//cout<<"pre calculation is done for the mean of the tensor"<<endl;
	int N_laplace_iter = 0;
	//we keep updating the dynamic probability sum
	double p_sum = 0.0; double p_sum_final = 0.0;
	while( (abs(valmean)) > eps ){
		double valmean_old = valmean;//assign the old mean value to another variable to allow us to calculate the second derivative metric once the new value is calculated
		int ii = 0;
		while(ii < Nx){
			int jj = 0;
			while(jj < Ny){
				int kk = 0;
				while(kk < Nz){
					if( ((abs(F_old(ii,jj,kk)-bc_on)) < eps) || ((abs(F_old(ii,jj,kk)-cluster_bc_on)) < eps) ){ kk = kk + 1; }
					else{
						double Vineg, Vjneg, Vkneg,  Vipos, Vjpos, Vkpos, Vijk;
						Vineg = F_old((ii-1),jj,kk); Vjneg = F_old(ii,(jj-1),kk); Vkneg = F_old(ii,jj,(kk-1));
						Vipos = F_old((ii+1),jj,kk); Vjpos = F_old(ii,(jj+1),kk); Vkpos = F_old(ii,jj,(kk+1));
						Vijk = F_old(ii,jj,kk);
						double Vijk_new = (1.0/6.0)*(Vineg + Vipos + Vjneg + Vjpos + Vkneg + Vkpos);//Laplace equation FDM equation
						F_new(ii,jj,kk) = Vijk_new;//the new value is updated in the new tensor field
						p_sum = p_sum + Vijk_new;//this is the probabilistic sum
						kk = kk + 1;
					}
				}
				jj = jj + 1;
			}
			ii = ii + 1;
		}
		//find the difference_type
		cube dF_tensor_new = F_new - F_old;//we calculate the difference between two elements in two tensors
		valmean = tensor_mean_value(dF_tensor_new, eps, bc_on, cluster_bc_on, Nx, Ny, Nz);//recalculate the mean value for loop criterea
		//valmean = tensor_max_value(dF_tensor_new, (0.01*eps), bc_on, cluster_bc_on, Nx, Ny, Nz);
		F_placeholder = F_old;
		F_old = F_new;//we assign the new tensor to the old one, so that the new one can be updated
		//p_sum_final = p_sum;
		N_laplace_iter = N_laplace_iter + 1;
		//cout<<"Laplacian Iteration = "<<N_laplace_iter<<", The tensor average difference is = "<<valmean<<endl;
		//if the difference in mean values is negative then we accept the move, if it is positive we break and the previous value is used
		if( (N_laplace_iter > 1) && (( (abs(valmean)) - (abs(valmean_old)) ) > 0.0) ){ F_new = F_placeholder; valmean = (0.1*eps); break; }//we impose two conditions to stop the full loop, only after one step though!
	}
	//cout<<"Laplace calculated!"<<endl;
	//if the threshhold is reached then the loop will break and we will be left with the laplace equation solution
	//we then normalize the full equation:
	//F_new = (F_new/p_sum_final);
	//then we return the field
	return F_new;//this field is the tensor we want.
}

//we will initalize a tensor and then compute the laplacian, we have a border boundary condition and a cluster boundary condition
cube poisson_solver_3d(cluster &provided_cluster, cube &rhs_cluster,double eps, double bc_on, double cluster_bc_on, int Nx, int Ny, int Nz){
	cube F_old = randu<cube>(Nx,Ny,Nz);//we initialize a uniform distribution random cube
	cube F_new = randu<cube>(Nx,Ny,Nz);//we initialize another random tensor for the new field
	cube F_placeholder = randu<cube>(Nx,Ny,Nz);
	//cout<<"set up of Poisson Equation is done"<<endl;
	//assign the boundaries their values
	int i = 0;
	while(i < Nx){
		int j = 0;
		while(j < Ny){
			int k = 0;
			while(k < Nz){
				if( i == 0 || j == 0 || k == 0 || i == (Nx-1) || j == (Ny-1) || k == (Nz-1)  ){
					//impose the condition for box boundary
					F_old(i,j,k) = bc_on; F_new(i,j,k) = bc_on; F_placeholder(i,j,k) = bc_on;
					k = k + 1;
				}
				else if( indices_are_in_cluster(i,j,k, provided_cluster.particles) ){
					//impose the boundary condition on the cluster
					F_old(i,j,k) = cluster_bc_on; F_new(i,j,k) = cluster_bc_on; F_placeholder(i,j,k) = cluster_bc_on;
					k = k + 1;
				}
				else{
					k = k + 1;
				}
			}
			j = j + 1;
		}
		i = i + 1;
	}
	//cout<<"boundary conditions are imposed"<<endl;
	//now that the cluster and box boundaries have been set, we can safely iterate till convergence.
	cube dF_tensor = F_new - F_old;
	double valmean = tensor_mean_value(dF_tensor, eps, bc_on, cluster_bc_on, Nx, Ny, Nz);
	//cout<<"pre calculation is done for the mean of the tensor"<<endl;
	int N_laplace_iter = 0;
	//we keep updating the dynamic probability sum
	double p_sum = 0.0; double p_sum_final = 0.0;
	while( (abs(valmean)) > eps ){
		double valmean_old = valmean;//assign the old mean value to another variable to allow us to calculate the second derivative metric once the new value is calculated
		int ii = 0;
		while(ii < Nx){
			int jj = 0;
			while(jj < Ny){
				int kk = 0;
				while(kk < Nz){
					if( ((abs(F_old(ii,jj,kk)-bc_on)) < eps) || ((abs(F_old(ii,jj,kk)-cluster_bc_on)) < eps) ){ kk = kk + 1; }
					else{
						double Vineg, Vjneg, Vkneg,  Vipos, Vjpos, Vkpos, Vijk;
						Vineg = F_old((ii-1),jj,kk); Vjneg = F_old(ii,(jj-1),kk); Vkneg = F_old(ii,jj,(kk-1));
						Vipos = F_old((ii+1),jj,kk); Vjpos = F_old(ii,(jj+1),kk); Vkpos = F_old(ii,jj,(kk+1));
						Vijk = F_old(ii,jj,kk);
						double Vijk_new = (1.0/6.0)*(Vineg + Vipos + Vjneg + Vjpos + Vkneg + Vkpos - (rhs_cluster(ii,jj,kk)) );//Poisson equation FDM equation
						F_new(ii,jj,kk) = Vijk_new;//the new value is updated in the new tensor field
						p_sum = p_sum + Vijk_new;//this is the probabilistic sum
						kk = kk + 1;
					}
				}
				jj = jj + 1;
			}
			ii = ii + 1;
		}
		//find the difference_type
		cube dF_tensor_new = F_new - F_old;//we calculate the difference between two elements in two tensors
		valmean = tensor_mean_value(dF_tensor_new, eps, bc_on, cluster_bc_on, Nx, Ny, Nz);//recalculate the mean value for loop criterea
		F_placeholder = F_old;
		F_old = F_new;//we assign the new tensor to the old one, so that the new one can be updated
		//p_sum_final = p_sum;
		N_laplace_iter = N_laplace_iter + 1;
		//cout<<"Poisson Iteration = "<<N_laplace_iter<<", The tensor average difference is = "<<valmean<<endl;
		//if the difference in mean values is negative then we accept the move, if it is positive we break and the previous value is used
		if( (N_laplace_iter > 1) && (( (abs(valmean)) - (abs(valmean_old)) ) > 0.0) ){ F_new = F_placeholder; valmean = (0.1*eps); break; }//we impose two conditions to stop the full loop, only after one step though!
	}
	//cout<<"Poisson calculated!"<<endl;
	//if the threshhold is reached then the loop will break and we will be left with the laplace equation solution
	//we then normalize the full equation:
	//F_new = (F_new/p_sum_final);
	//then we return the field
	return F_new;//this field is the tensor we want.
}
