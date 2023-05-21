// Robust MPC

/*
        This program simulates and tests a robust MPC controller on a servo-mechanical system.
            - The robust MPC controller uses a robust Kalman filter to perform state estimation in the MPC, with the
              knowledge that the linear model approximation does not exactly describe the model. This controller has
              been implemented in files robust_mpc.h and robust_mpc.cpp.
            - The servo mechanical system that has been considered here is inherently non-linear due to friction. The
              non-linear and linear models of this system are considered for simulation and testing. These models are
              implemented in servo_config.h.
*/

// Include necessary files that implement matrices, system model and the robust MPC controller.
#include "../controllers/model_defs.h"             // Contains functions associated with configuration
#include "../matrix/matrix.h"                 // Contains code related to matrix operations
#include "../matrix/matrix.cpp"
#include "../controllers/robust_mpc.h"             // Implements the robust MPC
#include "../controllers/robust_mpc.cpp"
#include <fstream>
#include <random>
#include "../servo_config.h"           // Contains functions that implement the linear and non-linear system servo-mechanism models

#define KLD_THRESH 0.1              // KL Divergence threshold (Radius of N-Dim ball around nominal (linear) model)

//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Function prototypes
void sample_gauss_pdf_iid(Matrix &m, float mean, float var, std::mt19937 gen);


int main() {

    // Initialize random number generator for sampling
    std::random_device rd{};
    std::mt19937 gen{rd()};

    // Create A, B, and C matrices for state space model
    Matrix A(SS_X_LEN, SS_X_LEN);
    Matrix B(SS_X_LEN, SS_U_LEN);
    Matrix C(SS_Z_LEN, SS_X_LEN);

    Matrix rt((HP * SS_Z_LEN), 1);          // Matrix to hold set point trajectory
    Matrix xtt(SS_X_LEN, 1);                // State vector
    Matrix utt(SS_U_LEN, 1);                // Input vector
    Matrix yt(SS_Z_LEN, 1);                 // Output vector

    // Weighting factors for MPC
    double Qk = 5000.0;
    double Rk = 0.01;
    double P0 = 1.0;
    double noise_var_est = 0.3E-2L;
    float eps = 1E-7;
    float lam = 2.0;

    // Create set point vector for MPC to track
    rt.vSetHomogen(1.57);

    // Initialize system with state space representative matrices
    initialize_linear_system(A, B, C);

    // Create an instance of robust_MPC type and initialize it with
    robust_MPC controller(A, B, C, Qk, Rk, noise_var_est, P0, KLD_THRESH, lam, eps);
//    robust_MPC controller(Qk, Rk, noise_var_est, P0, KLD_THRESH, lam, eps);

    // Create a file object to write input and output into a data file
    fstream f("data/data.txt", ios::out);
    f << "Time(s), Input Voltage (V), Output Load Angle (rad)\n";

    // Initialize state vector with samples drawn from a Gaussian distribution with certain mean and covariance
    float gau_mean = 0.0;
    float gau_var = 1.0;
    sample_gauss_pdf_iid(xtt, gau_mean, gau_var, gen);


    // Simulate model forward and determine optimal input signal at each time step.
    int e = 0;
    for(float t = 0; t < 20; t += 0.1){

        // Step controller through this time step given set point trajectory, current state estimate, control input,
        // and output from the system
        controller.step(rt, &xtt, &utt, yt);

        // Apply utt to the system - This is a simulation step. When applied to a real world system, this would be a
        // noisy measurement that would be received.
        yt = (C * (xtt + (SAMP_TIME * non_lin_dyn(xtt, utt, e + 1))));
        e++;

        // write time stamped input and output into file for visualization in post processing step
        f << t << ", " << utt[0][0] << ", " << yt[0][0] << "\n";

        cout << "output: " << yt[0][0] << "\n";
        cout << "---------------\n";
    }

    f.close();

    return 0;
}

// Function declarations
void sample_gauss_pdf_iid(Matrix &m, float mean, float var, std::mt19937 gen)
{
    /*
            Function to initialize entries of a matrix using IID (independent and identically distributed)
            samples drawn from a Gaussian distribution of fixed mean and variance.

            Arguments:
                m - A matrix whose elements are to random initialized.
                mean - mean of univariate gaussian distribution
                var - variance of univariate gaussian distribution
                gen - Mersenne Twister pseudo-random generator
    */

    // Create a Gaussian distribution sampler
    std::normal_distribution<> d{mean, var};

    // Iterate through each element of the matrix and sample from sampler for each element
    for(int row = 0; row < m.i32getRow(); row++)
    {
        for(int col=0; col < m.i32getColumn(); col++)
        {
            m[row][col] = d(gen);
        }
    }

}

