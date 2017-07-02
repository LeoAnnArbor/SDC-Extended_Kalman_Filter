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

   /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
    
    // initializing matrices
  R_laser_ = MatrixXd(2, 2);
    
  R_laser_ << 0.0225, 0,
  0, 0.0225;
    
  R_radar_ = MatrixXd(3, 3);
    
  R_radar_ << 0.09, 0, 0,
  0, 0.0009, 0,
  0, 0, 0.09;
    
  H_laser_ = MatrixXd(2, 4);
    
  H_laser_ << 1, 0, 0, 0,
  0, 1, 0, 0;
    
  Hj_ = MatrixXd(3, 4);
    

  noise_ax = 9;
  noise_ay = 9;
  
  //state covariance matrix P
  MatrixXd P_ = MatrixXd(4, 4);
    
  P_ << 1, 0, 0, 0,
  0, 1, 0, 0,
  0, 0, 1000, 0,
  0, 0, 0, 1000;
    
  MatrixXd F_ = MatrixXd(4, 4);
    
  MatrixXd Q_ = MatrixXd(4, 4);
    
  VectorXd x_ = VectorXd(4);
  x_ << 1, 1, 1, 1;
    
  ekf_.Init(x_, P_, F_, H_laser_, R_laser_, Q_);
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
//    ekf_.x_ = VectorXd(4);
//    ekf_.x_ << 1, 1, 1, 1;

    double px = 0;
    double py = 0;


    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      double rho = measurement_pack.raw_measurements_[0];
      double phi = measurement_pack.raw_measurements_[1];
        
      px = rho * cos(phi);
      py = rho * sin(phi);
        
       // If initial values are zero they will set to an initial guess
       // and the uncertainty will be increased.
       // Initial zeros would cause the algorithm to fail when using only Radar data.
      
     // if(fabs(px) < 0.0001)
     // {
     //   px = 1;
     //   ekf_.P_(0,0) = 1000;
     // }
      
     // if(fabs(py) < 0.0001)
     // {
     //   py = 1;
     //   ekf_.P_(1,1) = 1000;
     // }

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      px = measurement_pack.raw_measurements_[0];
      py = measurement_pack.raw_measurements_[1];
        
    }
    
    ekf_.x_  << px, py, 0, 0;

    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

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

  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
    
  ekf_.F_ << 1, 0, dt, 0,
  0, 1, 0, dt,
  0, 0, 1, 0,
  0, 0, 0, 1;
    
  float dt2 = dt * dt;
  float dt3 = dt2 * dt;
  float dt4 = dt3 * dt;
    
  ekf_.Q_ << noise_ax*dt4/4, 0, noise_ax*dt3/2, 0,
  0, noise_ay*dt4/4, 0, noise_ay*dt3/2,
  noise_ax*dt3/2, 0, noise_ax*dt2, 0,
  0, noise_ay*dt3/2, 0, noise_ay*dt2;
    
    
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
      ekf_.R_ = R_radar_;
      ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
      ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
  // Laser updates
      ekf_.R_ = R_laser_;
      ekf_.H_ = H_laser_;
      ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
