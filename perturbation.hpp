#include "./Electrostatic.hpp"

using namespace std;
using namespace arma;


vector<Eigen::Vector3i> find_open_positions(vector<particle> gparticles, int i, int dx, int dy, int dz){
	particle chosenp; chosenp = gparticles[i];
	vector<Eigen::Vector3i> genposes,finalposes;
	//first we generate the positions possible in an octahedral stencil around the particle of interest
	int di = 1;
	while(di <= dx){
		Eigen::Vector3i rt((chosenp.position(0)+di),(chosenp.position(1)),(chosenp.position(2))); genposes.push_back(rt);
		Eigen::Vector3i ru((chosenp.position(0)-di),(chosenp.position(1)),(chosenp.position(2))); genposes.push_back(ru);
		di = di + 1;
	}
	int dj = 1;
	while(dj <= dy){
		Eigen::Vector3i rt((chosenp.position(0)),(chosenp.position(1)+dj),(chosenp.position(2))); genposes.push_back(rt);
		Eigen::Vector3i ru((chosenp.position(0)),(chosenp.position(1)-dj),(chosenp.position(2))); genposes.push_back(ru);
		dj = dj + 1;
	}
	int dk = 1;
	while(dk <= dz){
		Eigen::Vector3i rt((chosenp.position(0)),(chosenp.position(1)),(chosenp.position(2)+dk)); genposes.push_back(rt);
		Eigen::Vector3i ru((chosenp.position(0)),(chosenp.position(1)),(chosenp.position(2)-dk)); genposes.push_back(ru);
		dk = dk + 1;
	}
	//then we test which ones are free
	vector<Eigen::Vector3i>::iterator it; it = genposes.begin();
	while(it != genposes.end()){
		Eigen::Vector3i posc; posc = *it;
		bool rispresent; rispresent = indices_are_in_cluster(posc(0), posc(1), posc(2), gparticles);
		if(rispresent){ ++it; }
		else{ finalposes.push_back(posc); ++it; }
	}
	//return this final set of positions
	return finalposes;
	
}

cluster perturb(cluster gclust, int nmaxstep, int dx, int dy, int dz, double P_perturb=0.0){
	cluster retclust; retclust = gclust;
	double toss = GenerateRandomDouble(0.0, 1.0);
	//we attempt a perturbation in the cluster at random probabilities
	if( islessequal(toss,P_perturb) ){
		int maxindx = retclust.particles.size();
		int stepno = 0;
		while(stepno < nmaxstep){
			int moveindx = GenerateRandomInteger(0, (maxindx-1));
			particle movep; movep = retclust.particles[moveindx];
			//figure out the free positions of motions
			vector<Eigen::Vector3i> genposes; genposes = find_open_positions(gclust.particles, moveindx);
			if(genposes.size() > 0){
				int rindx = GenerateRandomInteger(0, ((genposes.size())-1));
				Eigen::Vector3i r_random; r_random = genposes[rindx];
				movep.position = r_random;//modify particle position
				retclust.particles[moveindx] = movep;//modify the particle in the cluster
				stepno = stepno + 1;//increment step number
			}
			else{stepno = stepno + 0;}
		}
		cout<<"Cluster Perturbation was performed"<<endl;
		return retclust;
	}
	else{ return retclust; }
}