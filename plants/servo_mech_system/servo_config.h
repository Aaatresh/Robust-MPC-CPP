/*
    This file contains data and information associated with the servo mechanical system parameters. Both the linear and
    non-linear model can be found here.
*/


void initialize_linear_system(Matrix& A, Matrix& B, Matrix& C)
{
    // Define A
    A[0][0] = 7.63672682E-01L;      A[0][1] = 8.72694126E-02L;      A[0][2] =   1.18163659e-02;      A[0][3] = 3.18363234e-04;
    A[1][0] =  -4.42813522  ;      A[1][1] = 6.76403269E-01L  ;      A[1][2] = 2.21406761E-01L  ;      A[1][3] = 8.56906092E-03L;
    A[2][0] = 4.44371208E-01L ;      A[2][1] =   1.59181617E-02L ;      A[2][2] =  9.77781440E-01L  ;      A[2][3] = 6.20505175E-02L;
    A[3][0] =  7.12857003 ;      A[3][1] =  4.28453046E-01L;      A[3][2] =   -3.56428501E-01L  ;      A[3][3] = 3.44866161E-01L;

    B[0][0] =  -1.73122579E-04L;
    B[1][0] = 3.18363234E-04L;
    B[2][0] =  8.65612893E-06L;
    B[3][0] =  6.20505175E-02L;

    C[0][0] =  1.0000;      C[0][1] =  0.0000;      C[0][2] =  0.0000;      C[0][3] =  0.0000;
}


Matrix non_lin_dyn(Matrix x, Matrix u, int k)
{
    /*
        Function to simulate non-linear dynamics of the servo-mechancial system.

        Arguments:
            - x: Current state (x[k])
            - u: Control signal (u[k])
            - k: Time step index as we are working with a discrete time system

        Returns:
            - Next state (d/dt(x[k]))
    */

    // Model parameters
    float L = 0.8;   // armature coil inductance
    // J_m = 0.5;  // motor inertia
    float J_m = 0.5;  // motor inertia
    float beta_m = 0.1;   // motor viscous friction coefficient
    float R = 20;   // resistance of armature
    float Kt = 10;  // motor constant
    float rho = 20;  // gear ratio
    float k_theta = 1280.2;  // torsional rigidity
    float J_l = 25;  // nominal load inertia
    // J_l = 10;  // nominal load inertia
    float beta_l = 25;  // load viscous friction coefficient

    // Non-linear coefficients
    float alpha_l0 = 0.5;
    float alpha_l1 = 10.0;
    float alpha_l2 = 0.5;

    float alpha_m0 = 0.1;
    float alpha_m1 = 2;
    float alpha_m2 = 0.5;

    // Create matrix to hold d/dt(x)
    Matrix xdot(x.i32getRow(), x.i32getColumn());
    xdot = x;

    xdot[0][0] = x[1][0];
    xdot[2][0] = x[3][0];

    // Compute non-linearities
    float Tfl = (alpha_l0 + (alpha_l1 * exp(-alpha_l2 * abs(x[1][0])))) * sign(x[1][0]);
    float Ts = (k_theta / rho) * ((x[2][0] / rho) - x[0][0]);
    xdot[1][0] = (1 / J_l) * ((rho * Ts) - (beta_l * x[1][0]) - Tfl);

    float Tfm = (alpha_m0 + (alpha_m1 * exp(-alpha_m2 * abs(x[3][0])))) * sign(x[3][0]);
    float Im = ((u[0][0] - (Kt * x[3][0])) / R) * (1 - exp(-R * k * SAMP_TIME / L));
    float Tm = Kt * Im;
    xdot[3][0] = (1 / J_m) * (Tm - Ts - (beta_m * x[3][0]) - Tfm);

    return xdot;
}