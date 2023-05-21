/*
    This file is a configuration file for the simulation and testing of the robust MPC and standard MPC.
*/
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

/* State Space dimension */
#define SS_X_LEN    (4)
#define SS_Z_LEN    (1)
#define SS_U_LEN    (1)

// Sampling time of digital controller
#define SAMP_TIME 0.1

#define MATRIX_USE_BOUND_CHECKING

#define MATRIX_MAXIMUM_SIZE 28

#define VOLT_LIMIT 220.0

/* MPC horizon */
#define HP           (7)
#define HU           (3)

/* Set this define to choose math precision of the system */
#define PRECISION_SINGLE    1
#define PRECISION_DOUBLE    2
#define FPU_PRECISION       (PRECISION_DOUBLE)

#if (FPU_PRECISION == PRECISION_SINGLE)
    #define float_prec          float
    #define float_prec_ZERO     (1e-7)
    #define float_prec_ZERO_ECO (1e-5)      /* 'Economical' zero, for noisy calculation where 'somewhat zero' is good enough */
#elif (FPU_PRECISION == PRECISION_DOUBLE)
    #define float_prec          double
    #define float_prec_ZERO     (1e-30)
    #define float_prec_ZERO_ECO (1e-30)      /* 'Economical' zero, for noisy calculation where 'somewhat zero' is good enough */
#else
    #error("FPU_PRECISION has not been defined!");
#endif

