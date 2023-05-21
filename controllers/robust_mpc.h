/*
    Unconstrained Robust MPC implementation.
*/

#include "mpc.h"
#include "mpc.cpp"


class robust_MPC: public MPC
{
    /*
        This class implements an unconstrained Robust MPC
    */

public:

    // Constructor
    robust_MPC(Matrix &A, Matrix &B, Matrix &C, double Qk, double Rk,
                          double noise_var_est, double P0, double kld_thresh, double lambda, double eps);

protected:

    // Evaluate KL divergence
    float eval_kld(Matrix Pt, float param_t, float kld_thresh);
    // Bijection algorithm implementation
    float bijection_algo(Matrix Pt);

    // Compute state prediction in Kalman filter
    Matrix state_prediction(Matrix yt);
    // Compute state update in Kalman filter
    void state_update(Matrix utt, Matrix yt);

    // State estimate covariance matrix
    Matrix Vt        {SS_X_LEN, SS_X_LEN};

    // Create identity matrix
    Matrix identity_mat {SS_X_LEN, SS_X_LEN};

    // KL Divergence threshold
    float kld_thresh;

    // Upper bound threshold for bijection algorithm
    float lambda;

    // Bijection algorithm threshold
    float eps;
};
