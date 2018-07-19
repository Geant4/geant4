/*
*   calculate coupling coefficients of angular momenta
*
*   Author:
*                 Kawano, T <kawano@mailaps.org>
*
*   Modified by David Brown <dbrown@bnl.gov>
*       No longer must precompute the logarithm of the factorials.  
*       Also renamed things to make more Python friendly.
*       Finally, fixed a bunch of bugs & confusing conventions
*
*   Functions:
*
*   Note that arguments of those functions must be doubled, namely 1/2 is 1, etc.
*
*   wigner_3j(j1,j2,j3,j4,j5,j6)
*       Wigner's 3J symbol (similar to Clebsh-Gordan)
*               = / j1 j2 j3 \
*                 \ j4 j5 j6 /
*
*   wigner_6j(j1,j2,j3,j4,j5,j6)
*       Wigner's 6J symbol (similar to Racah)
*               = { j1 j2 j3 }
*                 { j4 j5 j6 }
*
*   wigner_9j(j1,j2,j3,j4,j5,j6,j7,j8,j9)
*       Wigner's 9J symbol
*                 / j1 j2 j3 \
*               = | j4 j5 j6 |
*                 \ j7 j8 j9 /
*
*   racah(j1, j2, l2, l1, j3, l3)
*               = W(j1, j2, l2, l1 ; j3, l3)
*               = (-1)^(j1+j2+l1+l2) * { j1 j2 j3 }
*                                      { l1 l2 l3 }
*
*   clebsh_gordan(j1,j2,m1,m2,j3)
*       Clebsh-Gordan coefficient 
*               = <j1,j2,m1,m2|j3,m1+m2>
*               = (-)^(j1-j2+m1+m2) * std::sqrt(2*j3+1) * / j1 j2   j3   \
*                                                    \ m1 m2 -m1-m2 /
*
*   z_coefficient(l1,j1,l2,j2,S,L)
*       Biedenharn's Z-coefficient coefficient
*               =  Z(l1  j1  l2  j2 | S L )
*
*   reduced_matrix_element(L,S,J,l0,j0,l1,j1)
*       Reduced Matrix Element for Tensor Operator
*               = < l1j1 || T(YL,sigma_S)J || l0j0 >
*
* References:
*   A. R. Edmonds, Angular Momentum in Quantum Mechanics, Princeton University Press 1974.
*   E. Condon, and G. Shortley, The Theory of Atomic Spectra, Cambridge, 1935.
*/

#include <stdlib.h>
#include <cmath>

#include "nf_specialFunctions.h"

#if defined __cplusplus
#include <cmath>
#include "G4Exp.hh"
namespace GIDI {
using namespace GIDI;
#endif

static const int    MAX_FACTORIAL    = 200;   // maximal factorial n!  (2 x Lmax)
/*static const double ARRAY_OVER  = 1.0e+300;   // force overflow */
static const double nf_amc_log_fact[] = {0.0, 0.0, 0.69314718056, 1.79175946923, 3.17805383035, 4.78749174278, 6.57925121201, 8.52516136107, 10.6046029027, 12.8018274801, 15.1044125731, 17.5023078459, 19.9872144957, 22.5521638531, 25.1912211827, 27.8992713838, 30.6718601061, 33.5050734501, 36.395445208, 39.3398841872, 42.3356164608, 45.3801388985, 48.4711813518, 51.6066755678, 54.7847293981, 58.003605223, 61.261701761, 64.557538627, 67.8897431372, 71.2570389672, 74.6582363488, 78.0922235533, 81.5579594561, 85.0544670176, 88.5808275422, 92.1361756037, 95.7196945421, 99.3306124548, 102.968198615, 106.631760261, 110.320639715, 114.034211781, 117.7718814, 121.533081515, 125.317271149, 129.123933639, 132.952575036, 136.802722637, 140.673923648, 144.565743946, 148.477766952, 152.409592584, 156.360836303, 160.331128217, 164.320112263, 168.327445448, 172.352797139, 176.395848407, 180.456291418, 184.533828861, 188.628173424, 192.739047288, 196.866181673, 201.009316399, 205.168199483, 209.342586753, 213.532241495, 217.736934114, 221.956441819, 226.190548324, 230.439043566, 234.701723443, 238.978389562, 243.268849003, 247.572914096, 251.89040221, 256.22113555, 260.564940972, 264.921649799, 269.291097651, 273.673124286, 278.06757344, 282.474292688, 286.893133295, 291.323950094, 295.766601351, 300.220948647, 304.686856766, 309.16419358, 313.65282995, 318.15263962, 322.663499127, 327.185287704, 331.717887197, 336.261181979, 340.815058871, 345.379407062, 349.954118041, 354.539085519, 359.13420537, 363.739375556, 368.354496072, 372.979468886, 377.614197874, 382.258588773, 386.912549123, 391.575988217, 396.248817052, 400.930948279, 405.622296161, 410.322776527, 415.032306728, 419.7508056, 424.478193418, 429.214391867, 433.959323995, 438.712914186, 443.475088121, 448.245772745, 453.024896238, 457.812387981, 462.608178527, 467.412199572, 472.224383927, 477.044665493, 481.87297923, 486.709261137, 491.553448223, 496.405478487, 501.265290892, 506.132825342, 511.008022665, 515.890824588, 520.781173716, 525.679013516, 530.584288294, 535.49694318, 540.416924106, 545.344177791, 550.278651724, 555.220294147, 560.169054037, 565.124881095, 570.087725725, 575.057539025, 580.034272767, 585.017879389, 590.008311976, 595.005524249, 600.009470555, 605.020105849, 610.037385686, 615.061266207, 620.091704128, 625.128656731, 630.172081848, 635.221937855, 640.27818366, 645.340778693, 650.409682896, 655.484856711, 660.566261076, 665.653857411, 670.747607612, 675.84747404, 680.953419514, 686.065407302, 691.183401114, 696.307365094, 701.437263809, 706.573062246, 711.714725802, 716.862220279, 722.015511874, 727.174567173, 732.339353147, 737.509837142, 742.685986874, 747.867770425, 753.05515623, 758.248113081, 763.446610113, 768.6506168, 773.860102953, 779.07503871, 784.295394535, 789.521141209, 794.752249826, 799.988691789, 805.230438804, 810.477462876, 815.729736304, 820.987231676, 826.249921865, 831.517780024, 836.790779582, 842.068894242, 847.35209797, 852.640365001, 857.933669826, 863.231987192};

static int parity( int x );
static int max3( int a, int b, int c );
static int max4( int a, int b, int c, int d );
static int min3( int a, int b, int c );
static double w6j0( int, int * );
static double w6j1( int * );
static double cg1( int, int, int );
static double cg2( int, int, int, int, int, int, int, int );
static double cg3( int, int, int, int, int, int );
/*static double triangle( int, int, int );*/
/*
============================================================
*/
double  nf_amc_log_factorial( int n ) {
/*
*       returns ln( n! ).
*/
    if( n > MAX_FACTORIAL ) return( INFINITY );
    if( n < 0 ) return( INFINITY );
    return nf_amc_log_fact[n];
}
/*
============================================================
*/
double  nf_amc_factorial( int n ) {
/*
*       returns n! for pre-computed table. INFINITY is return if n is negative or too large.
*/
    return G4Exp( nf_amc_log_factorial( n ) );
}
/*
============================================================
*/
double nf_amc_wigner_3j( int j1, int j2, int j3, int j4, int j5, int j6 ) {
/*
*       Wigner's 3J symbol (similar to Clebsh-Gordan)
*           = / j1 j2 j3 \
*             \ j4 j5 j6 /
*/
    double cg;

    if( ( j4 + j5 + j6 ) != 0 ) return( 0.0 );
    if( ( cg = nf_amc_clebsh_gordan( j1, j2, j4, j5, j3 ) ) == 0.0 ) return ( 0.0 );
    if( cg == INFINITY ) return( cg );
    return( ( ( ( j1 - j2 - j6 ) % 4 == 0 ) ?  1.0 : -1.0 ) * cg / std::sqrt( j3 + 1.0 ) );          /* BRB j3 + 1 <= 0? */
}
/*
============================================================
*/
double nf_amc_wigner_6j( int j1, int j2, int j3, int j4, int j5, int j6 ) {
/*
*       Wigner's 6J symbol (similar to Racah)
*               = { j1 j2 j3 }
*                 { j4 j5 j6 }
*/
    int i, x[6];

    x[0] = j1; x[1] = j2; x[2] = j3; x[3] = j4; x[4] = j5; x[5] = j6;
    for( i = 0; i < 6; i++ ) if ( x[i] == 0 ) return( w6j0( i, x ) );

    return( w6j1( x ) );
}
/*
============================================================
*/
static double w6j0( int i, int *x ) {

    switch( i ){
       case 0: if ( ( x[1] != x[2] ) || ( x[4] != x[5] ) ) return( 0.0 );
               x[5] = x[3]; x[0] = x[1]; x[3] = x[4];      break;
       case 1: if ( ( x[0] != x[2] ) || ( x[3] != x[5] ) ) return( 0.0 );
               x[5] = x[4];                                break;
       case 2: if ( ( x[0] != x[1] ) || ( x[3] != x[4] ) ) return( 0.0 );
               break;
	       //TK fix bug and add comment on 17-05-23
               //This is the case of 6.3.2 of A. R. Edmonds, Angular Momentum in Quantum Mechanics, Princeton University Press 1974.
       case 3: if ( ( x[1] != x[5] ) || ( x[2] != x[4] ) ) return( 0.0 );
               x[5] = x[0]; x[0] = x[4]; x[3] = x[1];      break;
       case 4: if ( ( x[0] != x[5] ) || ( x[2] != x[3] ) ) return( 0.0 );
               x[5] = x[1];                                break;
       case 5: if ( ( x[0] != x[4] ) || ( x[1] != x[3] ) ) return( 0.0 );
               x[5] = x[2];                                break;
    }

    if( ( x[5] > ( x[0] + x[3] ) ) || ( x[5] < std::abs( x[0] - x[3] ) ) ) return( 0.0 );
    if( x[0] > MAX_FACTORIAL || x[3] > MAX_FACTORIAL ) {        /* BRB Why this test? Why not x[5]? */
        return( INFINITY );
    }

    return( 1.0 / std::sqrt( (double) ( ( x[0] + 1 ) * ( x[3] + 1 ) ) ) * ( ( ( x[0] + x[3] + x[5] ) / 2 ) % 2 != 0 ? -1 : 1 ) );
}
/*
============================================================
*/
static double w6j1( int *x ) {

    double w6j, w;
    int i, k, k1, k2, n, l1, l2, l3, l4, n1, n2, n3, m1, m2, m3, x1, x2, x3, y[4];
    static int a[3][4] = { { 0, 0, 3, 3},
                           { 1, 4, 1, 4},
                           { 2, 5, 5, 2} };

    w6j = 0.0;

    for ( k = 0; k < 4; k++ ){
        x1 = x[ ( a[0][k] ) ];
        x2 = x[ ( a[1][k] ) ];
        x3 = x[ ( a[2][k] ) ];

        n = ( x1 + x2 + x3 ) / 2;
        if( n > MAX_FACTORIAL ) {
            return( INFINITY ); }
        else if( n < 0 ) {
            return( 0.0 );
        }

        if ( ( n1 = n - x3 ) < 0 ) return( 0.0 );
        if ( ( n2 = n - x2 ) < 0 ) return( 0.0 );
        if ( ( n3 = n - x1 ) < 0 ) return( 0.0 );

        y[k] = n + 2;
        w6j += nf_amc_log_fact[n1] + nf_amc_log_fact[n2] + nf_amc_log_fact[n3] - nf_amc_log_fact[n+1];
    }

    n1 = ( x[0] + x[1] + x[3] + x[4] ) / 2;
    n2 = ( x[0] + x[2] + x[3] + x[5] ) / 2;
    n3 = ( x[1] + x[2] + x[4] + x[5] ) / 2;

    k1 = max4( y[0], y[1], y[2], y[3] ) - 1;
    k2 = min3( n1, n2, n3 ) + 1;

    l1 = k1 - y[0] + 1;  m1 = n1 - k1 + 1;
    l2 = k1 - y[1] + 1;  m2 = n2 - k1 + 1;
    l3 = k1 - y[2] + 1;  m3 = n3 - k1 + 1;
    l4 = k1 - y[3] + 1;

    w6j = w = G4Exp( 0.5 * w6j + nf_amc_log_fact[k1] - nf_amc_log_fact[l1] - nf_amc_log_fact[l2] - nf_amc_log_fact[l3] - nf_amc_log_fact[l4]
                             - nf_amc_log_fact[m1] - nf_amc_log_fact[m2] - nf_amc_log_fact[m3] ) * ( ( k1 % 2 ) == 0 ? -1: 1 );
    if( w6j == INFINITY ) return( INFINITY );

    if( k1 != k2 ){
        k = k2 - k1;
        m1 -= k-1; m2 -= k-1; m3 -= k-1;
        l1 += k  ; l2 += k  ; l3 += k  ; l4 += k;

        for ( i = 0; i < k; i++ )
            w6j = w - w6j * ( ( k2 - i ) * ( m1 + i ) * ( m2 + i ) * ( m3 + i ) )
                    / ( ( l1 - i ) * ( l2 - i ) * ( l3 - i ) * ( l4 - i ) );
    }
    return( w6j );
}
/*
============================================================
*/
double nf_amc_wigner_9j( int j1, int j2, int j3, int j4, int j5, int j6, int j7, int j8, int j9 ) {
/*
*       Wigner's 9J symbol
*             / j1 j2 j3 \
*           = | j4 j5 j6 |
*             \ j7 j8 j9 /
*
*/
    int i, i0, i1;
    double rac;

    i0 = max3( std::abs( j1 - j9 ), std::abs( j2 - j6 ), std::abs( j4 - j8 ) );
    i1 = min3(     ( j1 + j9 ),     ( j2 + j6 ),     ( j4 + j8 ) );

    rac = 0.0;
    for ( i = i0; i <= i1; i += 2 ){ 
        rac += nf_amc_racah( j1, j4, j9, j8, j7,  i )
            *  nf_amc_racah( j2, j5,  i, j4, j8, j6 )
            *  nf_amc_racah( j9,  i, j3, j2, j1, j6 ) * ( i + 1 );
        if( rac == INFINITY ) return( INFINITY );
    }

    return( ( ( (int)( ( j1 + j3 + j5 + j8 ) / 2 + j2 + j4 + j9 ) % 4 == 0 ) ?  1.0 : -1.0 ) * rac );
}
/*
============================================================
*/
double nf_amc_racah( int j1, int j2, int l2, int l1, int j3, int l3 ) {
/*
*       Racah coefficient definition in Edmonds (AR Edmonds, "Angular Momentum in Quantum Mechanics", Princeton (1980) is
*       W(j1, j2, l2, l1 ; j3, l3) = (-1)^(j1+j2+l1+l2) * { j1 j2 j3 }
*                                                         { l1 l2 l3 }
*       The call signature of W(...) appears jumbled, but hey, that's the convention.
*
*       This convention is exactly that used by Blatt-Biedenharn (Rev. Mod. Phys. 24, 258 (1952)) too
*/

    double sig;

    sig = ( ( ( j1 + j2 + l1 + l2 ) % 4 == 0 ) ?  1.0 : -1.0 );
    return sig * nf_amc_wigner_6j( j1, j2, j3, l1, l2, l3 );
}

/*
============================================================
*/
/*
static double triangle( int a, int b, int c ) {

   int j1, j2, j3, j4;

   if ( ( j1 = (  a + b - c ) / 2 ) < 0 ) return( 0.0 );
   if ( ( j2 = (  a - b + c ) / 2 ) < 0 ) return( 0.0 );
   if ( ( j3 = ( -a + b + c ) / 2 ) < 0 ) return( 0.0 );
   j4 = ( a + b + c ) / 2 + 1;

   return( std::exp( 0.5 * ( nf_amc_log_fact[j1] + nf_amc_log_fact[j2] + nf_amc_log_fact[j3] - nf_amc_log_fact[j4] ) ) );
}
*/
/*
============================================================
*/
double nf_amc_clebsh_gordan( int j1, int j2, int m1, int m2, int j3 ) {
/*
*       Clebsh-Gordan coefficient 
*           = <j1,j2,m1,m2|j3,m1+m2>
*           = (-)^(j1-j2+m1+m2) * std::sqrt(2*j3+1) * / j1 j2   j3   \
*                                                \ m1 m2 -m1-m2 /
*
*   Note: Last value m3 is preset to m1+m2.  Any other value will evaluate to 0.0.
*/

    int m3, x1, x2, x3, y1, y2, y3;
    double cg = 0.0;

    if ( j1 < 0 || j2 < 0 || j3 < 0) return( 0.0 );
    if ( j1 + j2 + j3 > 2 * MAX_FACTORIAL ) return( INFINITY );

    m3 = m1 + m2;

    if ( ( x1 = ( j1 + m1 ) / 2 + 1 ) <= 0 ) return( 0.0 );
    if ( ( x2 = ( j2 + m2 ) / 2 + 1 ) <= 0 ) return( 0.0 );
    if ( ( x3 = ( j3 - m3 ) / 2 + 1 ) <= 0 ) return( 0.0 );

    if ( ( y1 = x1 - m1 ) <= 0 ) return( 0.0 );
    if ( ( y2 = x2 - m2 ) <= 0 ) return( 0.0 );
    if ( ( y3 = x3 + m3 ) <= 0 ) return( 0.0 );

    if ( j3 == 0 ){
        if ( j1 == j2 ) cg = ( 1.0 / std::sqrt( (double)j1 + 1.0 ) * ( ( y1 % 2 == 0 ) ? -1:1 ) );
    }
    else if ( (j1 == 0 || j2 == 0 ) ){
        if ( ( j1 + j2 ) == j3 ) cg = 1.0;
    }
    else {
        if( m3 == 0 && std::abs( m1 ) <= 1 ){
            if( m1 == 0 ) cg = cg1( x1, x2, x3 );
            else          cg = cg2( x1 + y1 - y2, x3 - 1, x1 + x2 - 2, x1 - y2, j1, j2, j3, m2 );
        }
        else if ( m2 == 0 && std::abs( m1 ) <=1 ){
                          cg = cg2( x1 - y2 + y3, x2 - 1, x1 + x3 - 2, x3 - y1, j1, j3, j3, m1 );
        }
        else if ( m1 == 0 && std::abs( m3 ) <= 1 ){
                          cg = cg2( x1, x1 - 1, x2 + x3 - 2, x2 - y3, j2, j3, j3, -m3 );
        }
        else              cg = cg3( x1, x2, x3, y1, y2, y3 );
    }

    return( cg );
}
/*
============================================================
*/
static double cg1( int x1, int x2, int x3 ) {

    int p1, p2, p3, p4, q1, q2, q3, q4;
    double a;

    p1 = x1 + x2 + x3 - 1; if ( ( p1 % 2 ) != 0 ) return( 0.0 );
    p2 = x1 + x2 - x3;
    p3 =-x1 + x2 + x3;
    p4 = x1 - x2 + x3;
    if ( p2 <= 0 || p3 <= 0 || p4 <= 0 ) return( 0.0 );
    if ( p1 >= MAX_FACTORIAL ) return( INFINITY );

    q1 = ( p1 + 1 ) / 2 - 1; p1--;
    q2 = ( p2 + 1 ) / 2 - 1; p2--;
    q3 = ( p3 + 1 ) / 2 - 1; p3--;
    q4 = ( p4 + 1 ) / 2 - 1; p4--;

    a = nf_amc_log_fact[q1]-( nf_amc_log_fact[q2] + nf_amc_log_fact[q3] + nf_amc_log_fact[q4] )
        + 0.5 * ( nf_amc_log_fact[ 2 * x3 - 1 ] - nf_amc_log_fact[ 2 * x3 - 2 ]
        + nf_amc_log_fact[p2] + nf_amc_log_fact[p3] + nf_amc_log_fact[p4] - nf_amc_log_fact[p1] );

    return( ( ( ( q1 + x1 - x2 ) % 2 == 0 ) ? 1.0 : -1.0 ) * G4Exp( a ) );
}
/*
============================================================
*/
static double cg2( int k, int q0, int z1, int z2, int w1, int w2, int w3, int mm ) {

    int q1, q2, q3, q4, p1, p2, p3, p4;
    double a;

    p1 =  z1 + q0 + 2;
    p2 =  z1 - q0 + 1;
    p3 =  z2 + q0 + 1;
    p4 = -z2 + q0 + 1;
    if ( p2 <= 0 || p3 <= 0 || p4 <= 0) return( 0.0 );
    if ( p1 >= MAX_FACTORIAL ) return( INFINITY );

    q1 = ( p1 + 1 ) / 2 - 1; p1--;
    q2 = ( p2 + 1 ) / 2 - 1; p2--;
    q3 = ( p3 + 1 ) / 2 - 1; p3--;
    q4 = ( p4 + 1 ) / 2 - 1; p4--;

    a = nf_amc_log_fact[q1] - ( nf_amc_log_fact[ q2 ] + nf_amc_log_fact[ q3 ] + nf_amc_log_fact[ q4 ] )
        + 0.5 * ( nf_amc_log_fact[ w3 + 1 ] - nf_amc_log_fact[ w3  ]
                + nf_amc_log_fact[ w1  ]    - nf_amc_log_fact[ w1 + 1 ]
                + nf_amc_log_fact[ w2  ]    - nf_amc_log_fact[ w2 + 1 ]
                + nf_amc_log_fact[ p2 ]     + nf_amc_log_fact[ p3 ] + nf_amc_log_fact[ p4 ] - nf_amc_log_fact[ p1 ] );

    return( ( ( ( q4 + k + ( mm > 0 ) * ( p1 + 2 ) ) % 2 == 0 ) ? -1.0 : 1.0 ) * 2.0 * G4Exp( a ) );
}
/*
============================================================
*/
static double cg3( int x1, int x2, int x3, int y1, int y2, int y3 ) {

    int nx, i, k1, k2, q1, q2, q3, q4, p1, p2, p3, z1, z2, z3;
    double a, cg;

    nx = x1 + x2 + x3 - 1;
    if ( ( z1 = nx - x1 - y1 ) < 0 ) return( 0.0 );
    if ( ( z2 = nx - x2 - y2 ) < 0 ) return( 0.0 );
    if ( ( z3 = nx - x3 - y3 ) < 0 ) return( 0.0 );

    k1 = x2 - y3;
    k2 = y1 - x3;

    q1 = max3( k1, k2, 0 );
    q2 = min3( y1, x2, z3 + 1 ) - 1;
    q3 = q1 - k1;
    q4 = q1 - k2;

    p1 = y1 - q1 - 1;
    p2 = x2 - q1 - 1;
    p3 = z3 - q1;

    a = cg = G4Exp( 0.5 * ( nf_amc_log_fact[ x3 + y3 - 1 ] - nf_amc_log_fact[ x3 + y3 - 2 ] - nf_amc_log_fact[ nx - 1 ]
                 + nf_amc_log_fact[ z1  ]    + nf_amc_log_fact[ z2  ]    + nf_amc_log_fact[ z3  ]
                 + nf_amc_log_fact[ x1 - 1 ] + nf_amc_log_fact[ x2 - 1 ] + nf_amc_log_fact[ x3 - 1 ]
                 + nf_amc_log_fact[ y1 - 1 ] + nf_amc_log_fact[ y2 - 1 ] + nf_amc_log_fact[ y3 - 1 ] )
                 - nf_amc_log_fact[ p1  ]    - nf_amc_log_fact[ p2  ]    - nf_amc_log_fact[ p3  ]
                 - nf_amc_log_fact[ q1  ]    - nf_amc_log_fact[ q3  ]    - nf_amc_log_fact[ q4  ] ) * ( ( ( q1 % 2 ) == 0 ) ? 1 : -1 );
    if( cg == INFINITY ) return( INFINITY );

    if ( q1 != q2 ){
        q3 = q2 - k1;
        q4 = q2 - k2;
        p1 = y1 - q2;
        p2 = x2 - q2;
        p3 = z3 - q2 + 1;
        for( i = 0; i < ( q2 - q1 ); i++ )
            cg = a - cg * ( ( p1 + i ) * ( p2 + i ) * ( p3 + i ) ) / ( ( q2 - i ) * ( q3 - i ) * ( q4 - i ) );
    }
    return( cg );
}
/*
============================================================
*/
double nf_amc_z_coefficient( int l1, int j1, int l2, int j2, int s, int ll ) {
/*
*       Biedenharn's Z-coefficient coefficient
*           =  Z(l1  j1  l2  j2 | S L )
*/
    double z, clebsh_gordan = nf_amc_clebsh_gordan( l1, l2, 0, 0, ll ), racah = nf_amc_racah( l1, j1, l2, j2, s, ll );

    if( ( clebsh_gordan == INFINITY ) || ( racah == INFINITY ) ) return( INFINITY );
    z = ( ( ( -l1 + l2 + ll ) % 8 == 0 ) ? 1.0 : -1.0 )
        * std::sqrt( l1 + 1.0 ) * std::sqrt( l2 + 1.0 ) * std::sqrt( j1 + 1.0 ) * std::sqrt( j2 + 1.0 ) * clebsh_gordan * racah;

   return( z );
}
/*
============================================================
*/
double nf_amc_zbar_coefficient( int l1, int j1, int l2, int j2, int s, int ll ) {
/*
*       Lane & Thomas's Zbar-coefficient coefficient
*           = Zbar(l1  j1  l2  j2 | S L )
*           = (-i)^( -l1 + l2 + ll ) * Z(l1  j1  l2  j2 | S L )
*
*       Lane & Thomas Rev. Mod. Phys. 30, 257-353 (1958).
*       Note, Lane & Thomas define this because they did not like the different phase convention in Blatt & Biedenharn's Z coefficient.  They changed it to get better time-reversal behavior.
*       Froehner uses Lane & Thomas convention as does T. Kawano.
*/
    double zbar, clebsh_gordan = nf_amc_clebsh_gordan( l1, l2, 0, 0, ll ), racah = nf_amc_racah( l1, j1, l2, j2, s, ll );

    if( ( clebsh_gordan == INFINITY ) || ( racah == INFINITY ) ) return( INFINITY );
    zbar = std::sqrt( l1 + 1.0 ) * std::sqrt( l2 + 1.0 ) * std::sqrt( j1 + 1.0 ) * std::sqrt( j2 + 1.0 ) * clebsh_gordan * racah;

   return( zbar );
}
/*
============================================================
*/
double nf_amc_reduced_matrix_element( int lt, int st, int jt, int l0, int j0, int l1, int j1 ) {
/*
*       Reduced Matrix Element for Tensor Operator
*           = < l1j1 || T(YL,sigma_S)J || l0j0 >
*
*   M.B.Johnson, L.W.Owen, G.R.Satchler
*   Phys. Rev. 142, 748 (1966)
*   Note: definition differs from JOS by the factor sqrt(2j1+1)
*/
    int    llt;
    double x1, x2, x3, reduced_mat, clebsh_gordan;

    if ( parity( lt ) != parity( l0 ) * parity( l1 ) ) return( 0.0 );
    if ( std::abs( l0 - l1 ) > lt || ( l0 + l1 ) < lt ) return( 0.0 );
    if ( std::abs( ( j0 - j1 ) / 2 ) > jt || ( ( j0 + j1 ) / 2 ) < jt ) return( 0.0 );

    llt = 2 * lt;
    jt *= 2;
    st *= 2;

    if( ( clebsh_gordan = nf_amc_clebsh_gordan( j1, j0, 1, -1, jt ) ) == INFINITY ) return( INFINITY );

    reduced_mat = 1.0 / std::sqrt( 4 * M_PI ) * clebsh_gordan / std::sqrt( jt + 1.0 )         /* BRB jt + 1 <= 0? */
                * std::sqrt( ( j0 + 1.0 ) * ( j1 + 1.0 ) * ( llt + 1.0 ) )
                * parity( ( j1 - j0 ) / 2 ) * parity( ( -l0 + l1 + lt ) / 2 ) * parity( ( j0 - 1 ) / 2 );

    if( st == 2 ){
        x1 = ( l0 - j0 / 2.0 ) * ( j0 + 1.0 );
        x2 = ( l1 - j1 / 2.0 ) * ( j1 + 1.0 );
        if ( jt == llt ){
            x3 = ( lt == 0 ) ? 0 :  ( x1 - x2 ) / std::sqrt( lt * ( lt + 1.0 ) );
        }
        else if ( jt == ( llt - st ) ){
            x3 = ( lt == 0 ) ? 0 : -( lt + x1 + x2 ) / std::sqrt( lt * ( 2.0 * lt + 1.0 ) );
        }
        else if ( jt == ( llt + st ) ){
            x3 = ( lt + 1 - x1 - x2 ) / std::sqrt( ( 2.0 * lt + 1.0 ) * ( lt + 1.0 ) );
        }
        else{
            x3 = 1.0;
        }
    }
    else x3 = 1.0;
    reduced_mat *= x3;

    return( reduced_mat );
}
/*
============================================================
*/
static int parity( int x ) {

    return( ( ( x / 2 ) % 2 == 0 ) ? 1 : -1 );
}
/*
============================================================
*/
static int max3( int a, int b, int c ) {

    if( a < b ) a = b;
    if( a < c ) a = c;
    return( a );
}
/*
============================================================
*/
static int max4( int a, int b, int c, int d ) {

    if( a < b ) a = b;
    if( a < c ) a = c;
    if( a < d ) a = d;
    return( a );
}
/*
============================================================
*/
static int min3( int a, int b, int c ) {

    if( a > b ) a = b;
    if( a > c ) a = c;
    return( a );
}

#if defined __cplusplus
}
#endif
