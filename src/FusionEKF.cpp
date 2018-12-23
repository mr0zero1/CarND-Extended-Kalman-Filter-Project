#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  
 
  
  

  //measurement covariance matrix - laser
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
   R_radar_ = MatrixXd(3, 3);
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;
  
  // H - measurement matrix
  H_laser_ = MatrixXd(2, 4);
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;
  
  // 
  H_radar_ = MatrixXd(3, 4);
  H_radar_ << 1, 0, 0, 0,
              0, 1, 0, 0,
              0, 0, 1, 0;
  
  // F - state transition matrix
  ekf_.F_= MatrixXd(4, 4);
  ekf_.F_ <<  1, 0, 1, 0,
              0, 1, 0, 1,
  			  0, 0, 1, 0,
  			  0, 0, 0, 1;

  // P - state covariance matrix
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0,    0,    0,
             0, 1,    0,    0,
             0, 0, 1000,    0,
             0, 0,    0, 1000;
    
  // Q - process covariance Matrix : 
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << 1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 1, 0,
             0, 0, 0, 1;
   
  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  noise_ax = 9.0;
  noise_ay = 9.0;
  
  //ekf_.Init(x, P, F, H, R, Q);

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    
    
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float ro    = measurement_pack.raw_measurements_[0];
      float theta = measurement_pack.raw_measurements_[1];  
      float x = ro * cos( theta );
      float y = ro * sin( theta );
      ekf_.x_   << x,
      			   y,
                   0,
                   0;
      
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      float x = measurement_pack.raw_measurements_[0];
      float y = measurement_pack.raw_measurements_[1];
      ekf_.x_   << x,
      			   y,
                   0,
                   0;
    }
    
    // done initializing, no need to predict or update
    is_initialized_ = true;
    previous_timestamp_ = measurement_pack.timestamp_;
    return;
  }

  //compute the time elapsed between the current and previous measurements
  	float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
	previous_timestamp_ = measurement_pack.timestamp_;
	


  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
	float dt_2 = dt * dt;
	float dt_3 = dt * dt_2;
	float dt_4 = dt * dt_3;
	
    // TODO: YOUR CODE HERE
	//1. Modify the F matrix so that the time is integrated
	ekf_.F_(0,2) = dt;
	ekf_.F_(1,3) = dt;
	//2. Set the process covariance matrix Q
	ekf_.Q_ = MatrixXd(4,4);
	ekf_.Q_ << dt_4/4*noise_ax,    0,               dt_3/2*noise_ax,      0,
               0,                  dt_4/4*noise_ay, 0,                    dt_3/2*noise_ay,
               dt_3/2*noise_ax,    0,               dt_2*noise_ax,        0,
               0,                  dt_3/2*noise_ay, 0,                    dt_2*noise_ay;

    ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    H_radar_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = H_radar_;
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  //cout << "x_ = " << ekf_.x_ << endl;
  //cout << "P_ = " << ekf_.P_ << endl;
}
