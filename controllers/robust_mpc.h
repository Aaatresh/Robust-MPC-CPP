/*
    Unconstrained Robust MPC implementation.
*/

#include "model_defs.h"
#include "../matrix.h"
#include "../utils.h"

class robust_MPC
{
    /*
        This class implements an unconstrained Robust MPC
    */

public:

    // Constructor
    robust_MPC(Matrix &A, Matrix &B, Matrix &C, double Qk, double Rk,
                          double noise_var_est, double P0, double kld_thresh, double lambda, double eps);

    // Function to initialize state space model
    void init_model(Matrix &A, Matrix &B, Matrix &C, double Qk, double Rk);

    // Function to calculate the optimal control
    bool get_utt(Matrix , Matrix , Matrix *u);

    // Function to perform a complete optimization step operation of the MPC at any given time step
    void step(Matrix rt, Matrix *x, Matrix *u, Matrix yt);

    // State space matrices
    Matrix A        {SS_X_LEN, SS_X_LEN};
    Matrix B        {SS_X_LEN, SS_U_LEN};
    Matrix C        {SS_Z_LEN, SS_X_LEN};

protected:

    // Evaluate KL divergence
    float eval_kld(Matrix Pt, float param_t, float kld_thresh);
    // Bijection algorithm implementation
    float bijection_algo(Matrix Pt);

    // Compute state prediction in Kalman filter
    Matrix state_prediction(Matrix yt);
    // Compute state update in Kalman filter
    void state_update(Matrix utt, Matrix yt);

    // Utility matrices to compute closed form solution for unconstrained standard MPC
    Matrix psi     {(HP*SS_Z_LEN), SS_X_LEN};
    Matrix omega   {(HP*SS_Z_LEN), SS_U_LEN};
    Matrix gamma   {(HP*SS_Z_LEN), (HU*SS_U_LEN)};

    // Temporary control matrix
    Matrix U_       {(HU*SS_U_LEN), 1};

    // Weighting matrices for prediction and control in loss function
    Matrix Q        {(HP*SS_Z_LEN), (HP*SS_Z_LEN)};
    Matrix R        {(HU*SS_U_LEN), (HU*SS_U_LEN)};

    // Defining utility variables for kalman filtering inside MPC framework
    Matrix Lt        {SS_X_LEN, 1};
    Matrix Kg        {SS_X_LEN, 1};
    Matrix xtt_1     {SS_X_LEN, 1};

    // Process and measurement covariance matrices
    Matrix G1        {SS_X_LEN, SS_X_LEN};
    Matrix D1        {1, SS_X_LEN};

    // State estimate covariance matrices
    Matrix Pt        {SS_X_LEN, SS_X_LEN};
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
