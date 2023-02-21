/*
    This program implements the Robust unconstrained MPC. Algorithm information can be found in the following
    references:
        - [1] Predictive Control: With Constraints, by Jan Marian Maciejowski.
        - [2] On the coupling of model predictive control and robust Kalman filtering, by Alberto Zenere and Mattia Zorzi.
*/


robust_MPC::robust_MPC(Matrix &A, Matrix &B, Matrix &C, double Qk, double Rk,
                         double noise_var_est, double P0, double kld_thresh, double lambda, double eps)
{
    /*
        Constructor for the Robust MPC class

        Arguments:
            - A, B, C: State space matrices
            - Qk, Rk: Weighting matrices for the prediction and control loss
            - noise_var_est: Estimate of the noise variance in the linear model
            - P0: Prior covariance estimate
    */

    // Initialize model
    init_model(A, B, C, Qk, Rk);

    // Initialize covariance matrices
    G1.vSetDiag(noise_var_est); G1[0][0] = 0.0;
    D1[0][0] = noise_var_est;

    // Initialize prior on state estimate
    Pt.vSetDiag(P0);

    // Create identity matrix
    identity_mat.vSetHomogen(1.0);

    // Initialize KL Divergence threshold, which is the radius of the ball around the nominal model
    robust_MPC::kld_thresh = kld_thresh;

    // Lambda value for the bijection algorithm
    robust_MPC::lambda = lambda;

    robust_MPC::eps = eps;
}

void robust_MPC::init_model(Matrix &A, Matrix &B, Matrix &C, double Qk, double Rk)
{
     /*
        Function to initialize state space model.

        Arguments:
            - A, B, C: State space matrices
            - Qk, Rk: Weighting matrices for the prediction and control loss
    */


    robust_MPC::A = A;
    robust_MPC::B = B;
    robust_MPC::C = C;
    Q.vSetDiag(Qk);
    R.vSetDiag(Rk);


    // Following the terminology of: Predictive Control: With Constraints, by Jan Marian Maciejowski
    Matrix Apow(SS_X_LEN, SS_X_LEN);
    Apow = A;
    for (int i = 0; i < HP; i++) {
        psi = psi.InsertSubMatrix((C * Apow), i * SS_Z_LEN, 0);
        Apow = Apow * A;
    }

    Matrix tempSigma(SS_X_LEN, SS_U_LEN);
    Apow.vSetIdentity();
    tempSigma = B;
    for (int i = 0; i < HP; i++) {
        omega = omega.InsertSubMatrix((C * Apow * B), i * SS_Z_LEN, 0);
        Apow = Apow * A;
    }

    for (int i = 0; i < HU; i++) {
        gamma = gamma.InsertSubMatrix(omega, i * SS_Z_LEN, i * SS_U_LEN, (HP * SS_Z_LEN) - (i * SS_Z_LEN), SS_U_LEN);
    }
}


bool robust_MPC::get_utt(Matrix SP, Matrix x, Matrix *u)
{
    /*
        Determine the optimal control by following the closed form solution of the unconstrained MPC problem.
        Reference: On the coupling of model predictive control and robust Kalman filtering.
    */

    // Calculate utility matrices to calculate optimal control for unconstrained MPC
    Matrix G((HU*SS_U_LEN), 1);
    Matrix H((HU*SS_U_LEN), (HU*SS_U_LEN));

    G = 2.0 * (gamma.Transpose()) * Q * (SP - (psi * x));
    H = ((gamma.Transpose()) * Q * gamma) + R;

    Matrix H_inv = H.Invers();
    if (!H_inv.bMatrixIsValid()) {
        U_.vSetToZero();
        return false;
    } else {
        U_ = (H_inv) * G * 0.5;
    }

    // Clip input to permissible voltage limit
    for (int i = 0; i < SS_U_LEN; i++) {

        if((U_[i][0] < VOLT_LIMIT) && (U_[i][0] > -VOLT_LIMIT))
            (*u)[i][0] = U_[i][0];
        else
        {
            (*u)[i][0] = sign(U_[i][0]) * VOLT_LIMIT;
        }
    }

    return true;
}



float robust_MPC::eval_kld(Matrix Pt, float param_t, float kld_thresh)
{
    /*
        Function to calculate the KL divergence

        Arguments:
            - Pt: Covariance matrix of state estimate
            - param_t: inverse of lagrange multiplier
            - kld_thresh: KL Divergence threshold

        Returns:
            - KL Divergence
    */

    float term1 = log(determinant(identity_mat - (param_t * Pt)));
    float term2 = trace((identity_mat - (param_t * Pt)).Invers() - identity_mat);

    return term1 + term2 - kld_thresh;
}

float robust_MPC::bijection_algo(Matrix Pt)
{
    /*
        Bijection algorithm to compute the roots of a function.

        Arguments:
            - Pt: Covariance matrix of state estimate

        Returns:
            root of function
    */


    // Define param1 and param2 as two bounds for the function to find a root
//    float eps = 1.0E-7;
    float param1 = eps, param2;

    if(((1 / lambda) > 2 * eps) and (lambda > 1E-5L))
        param2 = (1 / lambda) - eps;
    else
        param2 = eps;

    float param_t = param1;

    // Iterate until param1 and param2 are very close
    while(abs(param1 - param2) > eps)
    {
        param_t = (param1 + param2) / 2;
        float kld = eval_kld(Pt, param_t, kld_thresh);

        if(kld < 0)
            param1 = param_t;
        else
            param2 = param_t;
    }

    return param_t;
}


Matrix robust_MPC::state_prediction(Matrix yt)
{
     /*
            Computing the standard Kalman filter state prediction step by finding the mean prediction.

            Arguments:
                - yt: Output of the system at time step t

            Returns:
                - Mean state prediction
        */

        float param_t = bijection_algo(Pt);
        Vt = (Pt.Invers() - (param_t * identity_mat)).Invers();

        // Finding Lt
        Lt = (Vt * C.Transpose()) * ((C * Vt * C.Transpose()) + (D1 * D1.Transpose())).Invers();

        // Prediction step to find xtt
        return (xtt_1 + (Lt * (yt - (C * xtt_1))));
}

void robust_MPC::state_update(Matrix utt, Matrix yt)
{
        /*
        Standard Kalman filter state update step to update state estimate covariance matrix and estimated mean

        Arguments:
            - utt: Control signal at time step t
            - yt: Output at time step t
        */

        // Kalman gain
        Kg = A * Lt;

        // Update Pt
        Pt = (A * Vt * A.Transpose()) - (Kg * ((C * Vt * C.Transpose()) + (D1 * D1.Transpose())) * Kg.Transpose()) + (G1 * G1.Transpose());

        // Prediction step
        xtt_1 = (A * xtt_1) + (Kg * (yt - (C * xtt_1))) + (B * utt);
}


void robust_MPC::step(Matrix rt, Matrix *xtt, Matrix *utt, Matrix yt)
{
    /*
        Function to step through a single computation of the optimal control in the standard unconstrained Kalman filter, once.

        Arguments:
            - rt: Set point trajectory at time step t
            - xtt: State vector at time step t
            - utt: Control at time step t
            - yt: Output at time step t
    */

    // Compute state prediction using the robust Kalman filter
    *xtt = state_prediction(yt);

    // Determine control signal using the solution to the unconstrained robust MPC
    get_utt(rt, *xtt, utt);

    // Compute state update using the robust Kalman filter
    state_update(*utt, yt);
}
