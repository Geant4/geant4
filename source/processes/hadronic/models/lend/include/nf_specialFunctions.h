/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#ifndef specialFunctions_h_included
#define specialFunctions_h_included

#include <math.h>
#include <float.h>
#include <nf_utilities.h>

#ifdef WIN32
#define isfinite _finite
#define M_PI 3.141592653589793238463
/*#define INFINITY (DBL_MAX+DBL_MAX)*/
#endif

#if defined __cplusplus
    extern "C" {
    namespace GIDI {
#endif

double nf_polevl( double x, double coef[], int N );
double nf_p1evl( double x, double coef[], int N );
double nf_exponentialIntegral( int n, double x, nfu_status *status );
double nf_gammaFunction( double x, nfu_status *status );
double nf_logGammaFunction( double x, nfu_status *status );
double nf_incompleteGammaFunction( double a, double x, nfu_status *status );
double nf_incompleteGammaFunctionComplementary( double a, double x, nfu_status *status );

double  nf_amc_log_factorial( int );
double  nf_amc_factorial( int );
double  nf_amc_wigner_3j( int, int, int, int, int, int );
double  nf_amc_wigner_6j( int, int, int, int, int, int );
double  nf_amc_wigner_9j( int, int, int, int, int, int, int, int, int );
double  nf_amc_racah( int, int, int, int, int, int );
double  nf_amc_clebsh_gordan( int, int, int, int, int );
double  nf_amc_z_coefficient( int, int, int, int, int, int );
double  nf_amc_zbar_coefficient( int, int, int, int, int, int );
double  nf_amc_reduced_matrix_element( int, int, int, int, int, int, int );

#if defined __cplusplus
    }
    }
#endif

#endif          /* End of ptwXY_h_included. */
