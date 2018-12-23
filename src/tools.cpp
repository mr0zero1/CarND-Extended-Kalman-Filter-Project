#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0, 
          0, 
          0, 
          0;
  // check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size
	// ... your code here

	// check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size
  if ( (estimations.size()!=ground_truth.size() ) || (estimations.size() == 0) ) {
    cout << "*E, Error, invalid estimation or ground_truth data" << endl;
    return rmse;
  }
  
  //accumulate squared residuals
	for(unsigned int i=0; i < estimations.size(); ++i){
		VectorXd residual = estimations[i] - ground_truth[i];
		//coefficient-wise multiplication
		residual = residual.array()*residual.array();
		rmse += residual;
	}
	//calculate the mean
	rmse = rmse/estimations.size();

	//calculate the squared root
	rmse = rmse.array().sqrt();

	//return the result
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  	MatrixXd Hj(3,4);
  
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//TODO: YOUR CODE HERE 
    //pre-compute a set of terms to avoid repeated calculation
    float c1 = px*px + py*py;
    float c2 = sqrt(c1);
    float c3 = (c1*c2);
	//check division by zero
	if (fabs(c1) < 0.0001) {
	    cout << "Error: CalculateJacobian() div by zero " << endl;
	    return Hj;
	}
	//compute the Jacobian matrix
	Hj << (px/c2), (py/c2), 0, 0,
		  -(py/c1), (px/c1), 0, 0,
		  py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

	return Hj;

}

VectorXd Tools::Cartesian2Polar(const VectorXd& x_state)
{
    VectorXd polor_coord(3);
    polor_coord << 0
                 , 0
                 , 0
                 ;
  
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);
  
    float c1 = sqrt(px*px + py*py);
    if (fabs(c1) < 0.0001) {
	    //cout << "Error: CalculatePolar() ro , div by zero " << endl;
	    c1 = 0.0001;
	}
  
    // if ( (py == 0.0 ) && (px == 0.0) ) {
	//    cout << "Error: CalculatePolar() ro_dot , div by zero " << endl;
	//    return polor_coord;
    //	}
  
	float ro = c1;
    float theta = atan2(py, px);
    float ro_dot = (px*vx + py*vy)/c1;
  
  	polor_coord << ro
                 , theta
                 , ro_dot
                 ;
  
  	return polor_coord;
}

