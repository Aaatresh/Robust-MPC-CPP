/*
    Unconstrained standard MPC implementation.
*/

#include "model_defs.h"
#include "../matrix/matrix.h"
#include "../../utils/utils.h"

class controller
{
    public:

        void init_model();

        void step();
};

class MPC: public controller
{
    /*
        This class implements an unconstrained standard MPC
    */

public:
    // Constructor
    MPC(Matrix &A, Matrix &B, Matrix &C, double Qk, double Rk, double noise_var_est, double P0);
    // Function to initialize state space model
    void init_model(Matrix &A, Matrix &B, Matrix &C, double Qk, double Rk);

    // Function to calculate the optimal control
    bool get_utt(Matrix rt, Matrix x, Matrix *u);
    // Function to perform a complete optimization step operation of the MPC at any given time step
    void step(Matrix rt, Matrix *x, Matrix *u, Matrix yt);

protected:

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

    // State space matrices
    Matrix A        {SS_X_LEN, SS_X_LEN};
    Matrix B        {SS_X_LEN, SS_U_LEN};
    Matrix C        {SS_Z_LEN, SS_X_LEN};

    // Weighting matrices for prediction and control in loss function
    Matrix Q        {(HP*SS_Z_LEN), (HP*SS_Z_LEN)};
    Matrix R        {(HU*SS_U_LEN), (HU*SS_U_LEN)};

    // Utility matrices of Kalman filter
    Matrix Lt       {SS_X_LEN, 1};
    Matrix Kg       {SS_X_LEN, 1};
    Matrix xtt_1    {SS_X_LEN, 1};

    // Covariance matrices in Kalman filter
    Matrix G1       {SS_X_LEN, SS_X_LEN};
    Matrix D1       {1, SS_X_LEN};

    // Covariance matrix of state estimate
    Matrix Pt       {SS_X_LEN, SS_X_LEN};

};
