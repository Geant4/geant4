/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#include "nf_Legendre.h"

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

struct nf_Legendre_GaussianQuadrature_degree {
    int n;
    double *weights;
    double *xis;
};

static double sqrt_inv3 = 0.57735026918962576451;        /* sqrt( 1. / 3. ); */

#define n_3 3
static double weights_3[(n_3 + 1)/2] = { 8. / 9., 5. / 9. };
static double xis_3[(n_3 + 1)/2] = { 0., 0.77459666924148337704 };

#define n_4 4
static double weights_4[(n_4 + 1)/2] = { 0.65214515486254614263, 0.34785484513745385737 };
static double xis_4[(n_4 + 1)/2] = { 0.33998104358485626480, 0.86113631159405257522 };

#define n_5 5
static double weights_5[(n_5 + 1)/2] = { 0.568888888888889, 0.478628670499366, 0.236926885056189 };
static double xis_5[(n_5 + 1)/2] = { 0.0, 0.538469310105683, 0.906179845938664 };

#define n_10 10
static double weights_10[(n_10 + 1)/2] = { 0.295524224714752870, 0.269266719309996355, 0.219086362515982044, 0.149451349150580593, 0.066671344308688138 };
static double xis_10[(n_10 + 1)/2] = { 0.148874338981631211, 0.433395394129247191, 0.679409568299024406, 0.865063366688984511, 0.973906528517171720 };

#define n_20 20
static double weights_20[(n_20 + 1)/2] = { 
    0.152753387130725850698, 0.149172986472603746788, 0.142096109318382051329, 0.131688638449176626898, 0.118194531961518417312,
    0.101930119817240435037, 0.083276741576704748725, 0.062672048334109063570, 0.040601429800386941331, 0.017614007139152118312 };
static double xis_20[(n_20 + 1)/2] = { 
    0.076526521133497333755, 0.227785851141645078080, 0.373706088715419560673, 0.510867001950827098004, 0.636053680726515025453,
    0.746331906460150792614, 0.839116971822218823395, 0.912234428251325905868, 0.963971927277913791268, 0.993128599185094924786 };

#define n_40 40
static double weights_40[(n_40 + 1)/2] = {
    0.077505947978424811264, 0.077039818164247965588, 0.076110361900626242372, 0.074723169057968264200, 0.072886582395804059061, 
    0.070611647391286779696, 0.067912045815233903826, 0.064804013456601038075, 0.061306242492928939167, 0.057439769099391551367, 
    0.053227846983936824355, 0.048695807635072232061, 0.043870908185673271992, 0.038782167974472017640, 0.033460195282547847393, 
    0.027937006980023401099, 0.022245849194166957262, 0.016421058381907888713, 0.010498284531152813615, 0.004521277098533191258 };
static double xis_40[(n_40 + 1)/2] = { 
    0.038772417506050821933, 0.116084070675255208483, 0.192697580701371099716, 0.268152185007253681141, 0.341994090825758473007, 
    0.413779204371605001525, 0.483075801686178712909, 0.549467125095128202076, 0.612553889667980237953, 0.671956684614179548379, 
    0.727318255189927103281, 0.778305651426519387695, 0.824612230833311663196, 0.865959503212259503821, 0.902098806968874296728, 
    0.932812808278676533361, 0.957916819213791655805, 0.977259949983774262663, 0.990726238699457006453, 0.998237709710559200350 };

#define nSets 6
static struct nf_Legendre_GaussianQuadrature_degree GaussianQuadrature_degrees[nSets] = { { n_3, weights_3, xis_3 }, { n_4, weights_4, xis_4 }, 
    { n_5, weights_5, xis_5 }, { n_10, weights_10, xis_10 }, { n_20, weights_20, xis_20 }, { n_40, weights_40, xis_40 } };
/*
************************************************************
*/
nfu_status nf_Legendre_GaussianQuadrature( int degree, double x1, double x2, nf_Legendre_GaussianQuadrature_callback func, void *argList, double *integral ) {

    int i, n;
    double x, mu, sum, *weights, *xis;
    nfu_status status = nfu_Okay;

    *integral = 0;
    if( degree < 2 ) {
        status = func( 0.5 * ( x1 + x2 ), integral, argList );
        *integral *= 2.; }
    else if( degree < 4 ) {
        x = 0.5 * ( -sqrt_inv3 * ( x2 - x1 ) + x1 + x2 );
        if( ( status = func( x, integral, argList ) ) == nfu_Okay ) {
            x = 0.5 * (  sqrt_inv3 * ( x2 - x1 ) + x1 + x2 );
            status = func( x, &sum, argList );
            *integral += sum;
        } }
    else {
        for( i = 0; i < nSets - 1; i++ ) {
            if( GaussianQuadrature_degrees[i].n > ( degree + 1 ) / 2 ) break;
        }
        n = ( GaussianQuadrature_degrees[i].n + 1 ) / 2;
        weights = GaussianQuadrature_degrees[i].weights;
        xis = GaussianQuadrature_degrees[i].xis;
        for( i = 0; i < n; i++ ) {
            mu = xis[i];
            x = 0.5 * ( x1 * ( 1 - mu ) + x2 * ( mu + 1 ) );
            if( ( status = func( x, &sum, argList ) ) != nfu_Okay ) break;
            *integral += sum * weights[i];
            if( mu == 0 ) continue;
            x = x1 + x2 - x;
            if( ( status = func( x, &sum, argList ) ) != nfu_Okay ) break;
            *integral += sum * weights[i];
        }
    }
    *integral *= 0.5 * ( x2 - x1 );
    return( status );
}

#if defined __cplusplus
}
#endif
