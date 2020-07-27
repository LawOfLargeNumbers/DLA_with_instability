#ifndef DFLATTICE_HPP
#define DFLATTICE_HPP

#include "../Eigen/Eigen/Core"
#include "../Eigen/Eigen/Dense"

using namespace std;

/**
 * We want to figure out the electric field in an around a lattice point in a matrix
 * Uniform charge density on a sphere/shell described over a lattice "m". This allows us to make up electric fields in a box
*/
Eigen::MatrixXd generate_electric_field(Eigen::MatrixXd m, double dx, double dy, int ii, int jj){
    int i,j;
    //if ii or jj exced the index of the matrix we will assign it to the beginning
    if( ii == (m.rows()) ){ i = 0; j = jj; }
    else if( jj == (m.cols()) ){ j = 0; i = ii }
    else{ i = ii; j = jj; }
    //next we find the value of the electric field somewhere along the indices
}

#endif