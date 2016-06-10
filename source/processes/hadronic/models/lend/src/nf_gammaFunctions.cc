/*                            gamma.c
 *
 *    Gamma function
 *
 * DESCRIPTION:
 *
 * Returns gamma function of the argument.  The result is
 * correctly signed, and the sign (+1 or -1) is also
 * returned in a global (extern) variable named sgngam.
 * This variable is also filled in by the logarithmic gamma
 * function lgam().
 *
 * Arguments |x| <= 34 are reduced by recurrence and the function
 * approximated by a rational function of degree 6/7 in the
 * interval (2,3).  Large arguments are handled by Stirling's
 * formula. Large negative arguments are made positive using
 * a reflection formula.  
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE    -170,-33      20000       2.3e-15     3.3e-16
 *    IEEE     -33,  33     20000       9.4e-16     2.2e-16
 *    IEEE      33, 171.6   20000       2.3e-15     3.2e-16
 *
 * Error for arguments outside the test range will be larger
 * owing to error amplification by the exponential function.
 *
 */
/*                            lgam()
 *
 *    Natural logarithm of gamma function
 *
 *
 * DESCRIPTION:
 *
 * Returns the base e (2.718...) logarithm of the absolute
 * value of the gamma function of the argument.
 * The sign (+1 or -1) of the gamma function is returned in a
 * global (extern) variable named sgngam.
 *
 * For arguments greater than 13, the logarithm of the gamma
 * function is approximated by the logarithmic version of
 * Stirling's formula using a polynomial approximation of
 * degree 4. Arguments between -33 and +33 are reduced by
 * recurrence to the interval [2,3] of a rational approximation.
 * The cosecant reflection formula is employed for arguments
 * less than -33.
 *
 * Arguments greater than MAXLGM return DBL_MAX and an error
 * message.  MAXLGM = 2.556348e305 for IEEE arithmetic.
 *
 *
 * ACCURACY:
 *
 * arithmetic      domain        # trials     peak         rms
 *    IEEE    0, 3                 28000     5.4e-16     1.1e-16
 *    IEEE    2.718, 2.556e305     40000     3.5e-16     8.3e-17
 * The error criterion was relative when the function magnitude
 * was greater than one but absolute when it was less than one.
 *
 * The following test used the relative error criterion, though
 * at certain points the relative error could be much higher than
 * indicated.
 *    IEEE    -200, -4             10000     4.8e-16     1.3e-16
 *
 */
/*                            gamma.c    */
/*    gamma function    */

/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1989, 1992, 2000 by Stephen L. Moshier
*/

#include "nf_specialFunctions.h"

#if defined __cplusplus
#include "G4Exp.hh"
#include "G4Log.hh"
#include "G4Pow.hh"
namespace GIDI {
using namespace GIDI;
using namespace std;
#endif

static double P[] = { 1.60119522476751861407E-4, 1.19135147006586384913E-3, 1.04213797561761569935E-2, 4.76367800457137231464E-2,
                      2.07448227648435975150E-1, 4.94214826801497100753E-1, 9.99999999999999996796E-1 };
static double Q[] = { -2.31581873324120129819E-5, 5.39605580493303397842E-4, -4.45641913851797240494E-3, 1.18139785222060435552E-2,
                       3.58236398605498653373E-2, -2.34591795718243348568E-1, 7.14304917030273074085E-2, 1.00000000000000000320E0 };
#define MAXGAM 171.624376956302725
static double LOGPI = 1.14472988584940017414;
static double SQTPI = 2.50662827463100050242E0;

/* Stirling's formula for the gamma function */
static double STIR[5] = { 7.873113957930936284e-4, -2.2954996161337812638e-4, -2.6813261780578123283e-3, 3.472222216054586673e-3, 8.3333333333348225713e-2 };
#define MAXSTIR 143.01608

static double stirf( double x, nfu_status *status );
static double lgam( double x, int *sgngam, nfu_status *status );
/*
************************************************************
*/
static double stirf( double x, nfu_status * /*status*/ ) {
/* Gamma function computed by Stirling's formula. The polynomial STIR is valid for 33 <= x <= 172.  */

    double y, w, v;

    w = 1.0 / x;
    w = 1.0 + w * nf_polevl( w, STIR, 4 );
    y = G4Exp( x );
    if( x > MAXSTIR ) {                     /* Avoid overflow in pow() */
        v = G4Pow::GetInstance()->powA( x, 0.5 * x - 0.25 );
        y = v * (v / y); }
    else {
        y = G4Pow::GetInstance()->powA( x, x - 0.5 ) / y;
    }
    y = SQTPI * y * w;
    return( y );
}
/*
************************************************************
*/
double nf_gammaFunction( double x, nfu_status *status ) {

    double p, q, z;
    int i, sgngam = 1;

    *status = nfu_badInput;
    if( !isfinite( x ) ) return( x );
    *status = nfu_Okay;

    q = fabs( x );

    if( q > 33.0 ) {
        if( x < 0.0 ) {
            p = floor( q );
            if( p == q ) goto goverf;
            i = (int) p;
            if( ( i & 1 ) == 0 ) sgngam = -1;
            z = q - p;
            if( z > 0.5 ) {
                p += 1.0;
                z = q - p;
            }
            z = q * sin( M_PI * z );
            if( z == 0.0 ) goto goverf;
            z = M_PI / ( fabs( z ) * stirf( q, status ) );
        }
        else {
            z = stirf( x, status );
        }
        return( sgngam * z );
    }

    z = 1.0;
    while( x >= 3.0 ) {
        x -= 1.0;
        z *= x;
    } // Loop checking, 11.06.2015, T. Koi

    while( x < 0.0 ) {
        if( x > -1.E-9 ) goto small;
        z /= x;
        x += 1.0;
    } // Loop checking, 11.06.2015, T. Koi

    while( x < 2.0 ) {
        if( x < 1.e-9 ) goto small;
        z /= x;
        x += 1.0;
    } // Loop checking, 11.06.2015, T. Koi

    if( x == 2.0 ) return( z );

    x -= 2.0;
    p = nf_polevl( x, P, 6 );
    q = nf_polevl( x, Q, 7 );
    return( z * p / q );

small:
    if( x == 0.0 ) goto goverf;
    return( z / ( ( 1.0 + 0.5772156649015329 * x ) * x ) );

goverf:
    return( sgngam * DBL_MAX );
}

/* A[]: Stirling's formula expansion of log gamma
*  B[], C[]: log gamma function between 2 and 3
*/
static double A[] = { 8.11614167470508450300E-4, -5.95061904284301438324E-4, 7.93650340457716943945E-4,
                     -2.77777777730099687205E-3, 8.33333333333331927722E-2 };
static double B[] = { -1.37825152569120859100E3, -3.88016315134637840924E4, -3.31612992738871184744E5,
                      -1.16237097492762307383E6, -1.72173700820839662146E6, -8.53555664245765465627E5 };
static double C[] = { -3.51815701436523470549E2, -1.70642106651881159223E4, -2.20528590553854454839E5,
                      -1.13933444367982507207E6, -2.53252307177582951285E6, -2.01889141433532773231E6 };
static double LS2PI  =  0.91893853320467274178;     /* log( sqrt( 2*pi ) ) */
#define MAXLGM 2.556348e305

/*
************************************************************
*/
double nf_logGammaFunction( double x, nfu_status *status ) {
/* Logarithm of gamma function */

    int sgngam;

    *status = nfu_badInput;
    if( !isfinite( x ) ) return( x );
    *status = nfu_Okay;
    return( lgam( x, &sgngam, status ) );
}
/*
************************************************************
*/
static double lgam( double x, int *sgngam, nfu_status *status ) {

    double p, q, u, w, z;
    int i;

    *sgngam = 1;

    if( x < -34.0 ) {
        q = -x;
        w = lgam( q, sgngam, status );                  /* note this modifies *sgngam! */
        p = floor( q );
        if( p == q ) goto lgsing;
        i = (int) p;
        if( ( i & 1 ) == 0 ) {
            *sgngam = -1; }
        else {
            *sgngam = 1;
        }
        z = q - p;
        if( z > 0.5 ) {
            p += 1.0;
            z = p - q;
        }
        z = q * sin( M_PI * z );
        if( z == 0.0 ) goto lgsing;
        z = LOGPI - G4Log( z ) - w;
        return( z );
    }

    if( x < 13.0 ) {
        z = 1.0;
        p = 0.0;
        u = x;
        while( u >= 3.0 ) {
            p -= 1.0;
            u = x + p;
            z *= u;
        } // Loop checking, 11.06.2015, T. Koi
        while( u < 2.0 ) {
            if( u == 0.0 ) goto lgsing;
            z /= u;
            p += 1.0;
            u = x + p;
        } // Loop checking, 11.06.2015, T. Koi
        if( z < 0.0 ) {
            *sgngam = -1;
            z = -z; }
        else {
            *sgngam = 1;
        }
        if( u == 2.0 ) return( G4Log( z ) );
        p -= 2.0;
        x = x + p;
        p = x * nf_polevl( x, B, 5 ) / nf_p1evl( x, C, 6);
        return( G4Log( z ) + p );
    }

    if( x > MAXLGM ) goto lgsing;
    q = ( x - 0.5 ) * G4Log( x ) - x + LS2PI;
    if( x > 1.0e8 ) return( q );

    p = 1.0 / ( x * x );
    if( x >= 1000.0 ) {
        q += ( ( 7.9365079365079365079365e-4 * p - 2.7777777777777777777778e-3 ) * p + 0.0833333333333333333333 ) / x; }
    else {
        q += nf_polevl( p, A, 4 ) / x;
    }
    return( q );

lgsing:
    return( *sgngam * DBL_MAX );
}

#if defined __cplusplus
}
#endif
