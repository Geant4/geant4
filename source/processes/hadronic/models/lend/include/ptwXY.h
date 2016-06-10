/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#ifndef ptwXY_h_included
#define ptwXY_h_included

#include <stdio.h>
#include <stdint.h>

#include <nf_utilities.h>
#include <ptwX.h>

#if defined __cplusplus
    extern "C" {
    namespace GIDI {
#endif

#define ptwXY_minimumSize 10            /* This must be > 0 otherwise some logic will fail. */
#define ptwXY_minimumOverflowSize 4     /* This must be > 0 otherwise some logic will fail. */
#define ptwXY_maxBiSectionMax 20
#define ptwXY_minAccuracy 1e-14
#define ptwXY_sectionSubdivideMax 1 << 16
#define ClosestAllowXFactor 10

typedef enum ptwXY_dataFrom_e { ptwXY_dataFrom_Unknown, ptwXY_dataFrom_Points, ptwXY_dataFrom_Overflow } ptwXY_dataFrom;
typedef enum ptwXY_group_normType_e { ptwXY_group_normType_none, ptwXY_group_normType_dx, ptwXY_group_normType_norm } ptwXY_group_normType;

/* The next macro are used in the routine ptwXY_union. */
#define ptwXY_union_fill 1              /* If filling, union is filled with y value of first ptw. */
#define ptwXY_union_trim 2              /* If trimming, union in only over common domain of ptw1 and ptw2. */
#define ptwXY_union_mergeClosePoints 4  /* If true, union calls ptwXY_mergeClosePoints with eps = 4 * DBL_EPSILON. */
typedef enum ptwXY_sigma_e { ptwXY_sigma_none, ptwXY_sigma_plusMinus, ptwXY_sigma_Minus, ptwXY_sigma_plus } ptwXY_sigma;
typedef enum ptwXY_interpolation_e { ptwXY_interpolationLinLin, ptwXY_interpolationLinLog, ptwXY_interpolationLogLin, ptwXY_interpolationLogLog,
    ptwXY_interpolationFlat, ptwXY_interpolationOther } ptwXY_interpolation;

/*
*  The function ptwXY_getPointsAroundX determines where an x fits into a ptwXY instance. It returns/sets the following.
*
*  if ( some point's x == x )
*      lessThanEqualXPoint is set to point's information (prior, next, index, x, y),
*      greaterThanXPoint is set to a overflowHeader,
*      return( ptwXY_lessEqualGreaterX_equal ).
*   else if ( x < first point's x )
*       lessThanEqualXPoint is set to overflowHeader,
*       greaterThanXPoint is set to first point's information,
*       and greaterThanXPoint.prior points to the overflow which will be before the new point when the new point is inserted into overflowPoints.
*   else if ( x > last point's x )
*       lessThanEqualXPoint is set to last point's information
*       greaterThanXPoint is set to a overflowHeader point
*       and lessThanEqualXPoint.prior points to the overflow which will be before new point when the new point is inserted into overflowPoints.
*   else
*       lessThanEqualXPoint is set to point's information for closes point with point's x <= x
*       greaterThanXPoint is set to point's information for closes point with point's x > x
*/
typedef enum ptwXY_lessEqualGreaterX_e { ptwXY_lessEqualGreaterX_empty, ptwXY_lessEqualGreaterX_lessThan, ptwXY_lessEqualGreaterX_equal,
    ptwXY_lessEqualGreaterX_between, ptwXY_lessEqualGreaterX_greater } ptwXY_lessEqualGreaterX;

typedef
    struct ptwXYPoint_s {
        double x, y;
    } ptwXYPoint;

typedef nfu_status (*ptwXY_createFromFunction_callback)( double x, double *y, void *argList );
typedef nfu_status (*ptwXY_applyFunction_callback)( ptwXYPoint *point, void *argList );
typedef nfu_status (*ptwXY_getValue_callback)( void *argList, double x, double *y, double x1, double y1, double x2, double y2 );

typedef struct {
    char const *interpolationString;
    ptwXY_getValue_callback getValueFunc;
    void *argList;
} ptwXY_interpolationOtherInfo;

typedef
    struct ptwXYOverflowPoint_s {
        struct ptwXYOverflowPoint_s *prior;
        struct ptwXYOverflowPoint_s *next;
        int64_t index;                             /* For overflowHeader set to -1. */
        ptwXYPoint point;
    } ptwXYOverflowPoint;

typedef
    struct ptwXYPoints_s {
        nfu_status status;
        ptwXY_sigma typeX, typeY;
        ptwXY_interpolation interpolation;
        ptwXY_interpolationOtherInfo interpolationOtherInfo;
        int userFlag;
        double biSectionMax;
        double accuracy;
        double minFractional_dx;
        int64_t length;
        int64_t allocatedSize;
        int64_t overflowLength;
        int64_t overflowAllocatedSize;
        int64_t mallocFailedSize;
        ptwXYOverflowPoint overflowHeader;
        ptwXYPoint *points;
        ptwXYOverflowPoint *overflowPoints;
    } ptwXYPoints;

/*
* Routines in ptwXY_core.c
*/
ptwXYPoints *ptwXY_new( ptwXY_interpolation interpolation, ptwXY_interpolationOtherInfo const *interpolationOtherInfo, double biSectionMax,
    double accuracy, int64_t primarySize, int64_t secondarySize, nfu_status *status, int userFlag );
nfu_status ptwXY_setup( ptwXYPoints *ptwXY, ptwXY_interpolation interpolation, ptwXY_interpolationOtherInfo const *interpolationOtherInfo, 
    double biSectionMax, double accuracy, int64_t primarySize, int64_t secondarySize, int userFlag );
ptwXYPoints *ptwXY_create( ptwXY_interpolation interpolation, ptwXY_interpolationOtherInfo const *interpolationOtherInfo, 
    double biSectionMax, double accuracy, int64_t primarySize, int64_t secondarySize, int64_t length, double const *xy, 
    nfu_status *status, int userFlag );
ptwXYPoints *ptwXY_createFrom_Xs_Ys( ptwXY_interpolation interpolation, ptwXY_interpolationOtherInfo const *interpolationOtherInfo, 
    double biSectionMax, double accuracy, int64_t primarySize, int64_t secondarySize, int64_t length, double const *Xs, 
    double const *Ys, nfu_status *status, int userFlag );

nfu_status ptwXY_copy( ptwXYPoints *dest, ptwXYPoints *src );
ptwXYPoints *ptwXY_clone( ptwXYPoints *ptwXY, nfu_status *status );
ptwXYPoints *ptwXY_cloneToInterpolation( ptwXYPoints *ptwXY, ptwXY_interpolation interpolationTo, nfu_status *status );
ptwXYPoints *ptwXY_slice( ptwXYPoints *ptwXY, int64_t index1, int64_t index2, int64_t secondarySize, nfu_status *status );
ptwXYPoints *ptwXY_xSlice( ptwXYPoints *ptwXY, double xMin, double xMax, int64_t secondarySize, int fill, nfu_status *status );
ptwXYPoints *ptwXY_xMinSlice( ptwXYPoints *ptwXY, double xMin, int64_t secondarySize, int fill, nfu_status *status );
ptwXYPoints *ptwXY_xMaxSlice( ptwXYPoints *ptwXY, double xMax, int64_t secondarySize, int fill, nfu_status *status );

ptwXY_interpolation ptwXY_getInterpolation( ptwXYPoints *ptwXY );
char const *ptwXY_getInterpolationString( ptwXYPoints *ptwXY );
nfu_status ptwXY_getStatus( ptwXYPoints *ptwXY );
int ptwXY_getUserFlag( ptwXYPoints *ptwXY );
void ptwXY_setUserFlag( ptwXYPoints *ptwXY, int userFlag );
double ptwXY_getAccuracy( ptwXYPoints *ptwXY );
double ptwXY_setAccuracy( ptwXYPoints *ptwXY, double accuracy );
double ptwXY_getBiSectionMax( ptwXYPoints *ptwXY );
double ptwXY_setBiSectionMax( ptwXYPoints *ptwXY, double biSectionMax );

nfu_status ptwXY_reallocatePoints( ptwXYPoints *ptwXY, int64_t size, int forceSmallerResize );
nfu_status ptwXY_reallocateOverflowPoints( ptwXYPoints *ptwXY, int64_t size );
nfu_status ptwXY_coalescePoints( ptwXYPoints *ptwXY, int64_t size, ptwXYPoint *newPoint, int forceSmallerResize );
nfu_status ptwXY_simpleCoalescePoints( ptwXYPoints *ptwXY );

nfu_status ptwXY_clear( ptwXYPoints *ptwXY );
nfu_status ptwXY_release( ptwXYPoints *ptwXY );
ptwXYPoints *ptwXY_free( ptwXYPoints *ptwXY );

int64_t ptwXY_length( ptwXYPoints *ptwXY );
int64_t ptwXY_getNonOverflowLength( ptwXYPoints const *ptwXY );

nfu_status ptwXY_setXYData( ptwXYPoints *ptwXY, int64_t length, double const *xy );
nfu_status ptwXY_setXYDataFromXsAndYs( ptwXYPoints *ptwXY, int64_t length, double const *x, double const *y );
nfu_status ptwXY_deletePoints( ptwXYPoints *ptwXY, int64_t i1, int64_t i2 );
ptwXYPoint *ptwXY_getPointAtIndex( ptwXYPoints *ptwXY, int64_t index );
ptwXYPoint *ptwXY_getPointAtIndex_Unsafely( ptwXYPoints *ptwXY, int64_t index );
nfu_status ptwXY_getXYPairAtIndex( ptwXYPoints *ptwXY, int64_t index, double *x, double *y );
ptwXY_lessEqualGreaterX ptwXY_getPointsAroundX( ptwXYPoints *ptwXY, double x, ptwXYOverflowPoint *lessThanEqualXPoint, ptwXYOverflowPoint *greaterThanXPoint );
ptwXY_lessEqualGreaterX ptwXY_getPointsAroundX_closeIsEqual( ptwXYPoints *ptwXY, double x, ptwXYOverflowPoint *lessThanEqualXPoint,
        ptwXYOverflowPoint *greaterThanXPoint, double eps, int *closeIsEqual, ptwXYPoint **closePoint );
nfu_status ptwXY_getValueAtX( ptwXYPoints *ptwXY, double x, double *y );
nfu_status ptwXY_setValueAtX( ptwXYPoints *ptwXY, double x, double y );
nfu_status ptwXY_setValueAtX_overrideIfClose( ptwXYPoints *ptwXY, double x, double y, double eps, int override );
nfu_status ptwXY_mergeFromXsAndYs( ptwXYPoints *ptwXY, int length, double *xs, double *ys );
nfu_status ptwXY_mergeFromXYs( ptwXYPoints *ptwXY, int length, double *xys );
nfu_status ptwXY_appendXY( ptwXYPoints *ptwXY, double x, double y );
nfu_status ptwXY_setXYPairAtIndex( ptwXYPoints *ptwXY, int64_t index, double x, double y );

nfu_status ptwXY_getSlopeAtX( ptwXYPoints *ptwXY, double x, const char side, double *slope );

double ptwXY_getXMinAndFrom( ptwXYPoints *ptwXY, ptwXY_dataFrom *dataFrom );
double ptwXY_getXMin( ptwXYPoints *ptwXY );
double ptwXY_getXMaxAndFrom( ptwXYPoints *ptwXY, ptwXY_dataFrom *dataFrom );
double ptwXY_getXMax( ptwXYPoints *ptwXY );
double ptwXY_getYMin( ptwXYPoints *ptwXY );
double ptwXY_getYMax( ptwXYPoints *ptwXY );

/* 
* Methods in ptwXY_methods.c 
*/
nfu_status ptwXY_clip( ptwXYPoints *ptwXY1, double yMin, double yMax );
nfu_status ptwXY_thicken( ptwXYPoints *ptwXY1, int sectionSubdivideMax, double dxMax, double fxMax );
ptwXYPoints *ptwXY_thin( ptwXYPoints *ptwXY1, double accuracy, nfu_status *status );
nfu_status ptwXY_trim( ptwXYPoints *ptwXY );

ptwXYPoints *ptwXY_union( ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, nfu_status *status, int unionOptions );

nfu_status ptwXY_scaleOffsetXAndY( ptwXYPoints *ptwXY, double xScale, double xOffset, double yScale, double yOffset );

/*
* Functions in ptwXY_unitaryOperators.c
*/
nfu_status ptwXY_abs( ptwXYPoints *ptwXY );
nfu_status ptwXY_neg( ptwXYPoints *ptwXY );

/*
* Functions in ptwXY_binaryOperators.c
*/
nfu_status ptwXY_slopeOffset( ptwXYPoints *ptwXY, double slope, double offset );
nfu_status ptwXY_add_double( ptwXYPoints *ptwXY, double value );
nfu_status ptwXY_sub_doubleFrom( ptwXYPoints *ptwXY, double value );
nfu_status ptwXY_sub_fromDouble( ptwXYPoints *ptwXY, double value );
nfu_status ptwXY_mul_double( ptwXYPoints *ptwXY, double value );
nfu_status ptwXY_div_doubleFrom( ptwXYPoints *ptwXY, double value );
nfu_status ptwXY_div_fromDouble( ptwXYPoints *ptwXY, double value );
nfu_status ptwXY_mod( ptwXYPoints *ptwXY, double m, int pythonMod );

ptwXYPoints *ptwXY_binary_ptwXY( ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, double v1, double v2, double v1v2, nfu_status *status );
ptwXYPoints *ptwXY_add_ptwXY( ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, nfu_status *status );
ptwXYPoints *ptwXY_sub_ptwXY( ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, nfu_status *status );
ptwXYPoints *ptwXY_mul_ptwXY( ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, nfu_status *status );
ptwXYPoints *ptwXY_mul2_ptwXY( ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, nfu_status *status );
ptwXYPoints *ptwXY_div_ptwXY( ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, nfu_status *status, int safeDivide );

/* 
* Functions in ptwXY_functions.c 
*/
nfu_status ptwXY_pow( ptwXYPoints *ptwXY, double p );
nfu_status ptwXY_exp( ptwXYPoints *ptwXY, double a );
ptwXYPoints *ptwXY_convolution( ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, nfu_status *status, int mode );

/*
* Functions in ptwXY_interpolation.c
*/
nfu_status ptwXY_interpolatePoint( ptwXY_interpolation interpolation, double x, double *y, double x1, double y1, double x2, double y2 );
ptwXYPoints *ptwXY_flatInterpolationToLinear( ptwXYPoints *ptwXY, double lowerEps, double upperEps, nfu_status *status );
ptwXYPoints *ptwXY_toOtherInterpolation( ptwXYPoints *ptwXY, ptwXY_interpolation interpolation, double accuracy, nfu_status *status );
ptwXYPoints *ptwXY_unitbaseInterpolate( double w, double w1, ptwXYPoints *ptwXY1, double w2, ptwXYPoints *ptwXY2, nfu_status *status );
ptwXYPoints *ptwXY_toUnitbase( ptwXYPoints *ptwXY, nfu_status *status );
ptwXYPoints *ptwXY_fromUnitbase( ptwXYPoints *ptwXY, double xMin, double xMax, nfu_status *status );

/* 
* Functions in ptwXY_convenient.c 
*/
ptwXPoints *ptwXY_getXArray( ptwXYPoints *ptwXY, nfu_status *status );
nfu_status ptwXY_dullEdges( ptwXYPoints *ptwXY, double lowerEps, double upperEps, int positiveXOnly );
nfu_status ptwXY_mergeClosePoints( ptwXYPoints *ptwXY, double epsilon );
ptwXYPoints *ptwXY_intersectionWith_ptwX( ptwXYPoints *ptwXY, ptwXPoints *ptwX, nfu_status *status );
nfu_status ptwXY_areDomainsMutual( ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2 );
nfu_status ptwXY_tweakDomainsToMutualify( ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, int epsilonFactor, double epsilon );
nfu_status ptwXY_mutualifyDomains( ptwXYPoints *ptwXY1, double lowerEps1, double upperEps1, int positiveXOnly1,
                                          ptwXYPoints *ptwXY2, double lowerEps2, double upperEps2, int positiveXOnly2 );
nfu_status ptwXY_copyToC_XY( ptwXYPoints *ptwXY, int64_t index1, int64_t index2, int64_t allocatedSize, int64_t *numberOfPoints, double *xy );
nfu_status ptwXY_valueTo_ptwXAndY( ptwXYPoints *ptwXY, double **xs, double **ys );
ptwXYPoints *ptwXY_valueTo_ptwXY( double x1, double x2, double y, nfu_status *status );
ptwXYPoints *ptwXY_createGaussianCenteredSigma1( double accuracy, nfu_status *status );
ptwXYPoints *ptwXY_createGaussian( double accuracy, double xCenter, double sigma, double amplitude, double xMin, double xMax, 
        double dullEps, nfu_status *status );

/* 
* Functions in ptwXY_misc.c 
*/
void ptwXY_update_biSectionMax( ptwXYPoints *ptwXY1, double oldLength );
ptwXYPoints *ptwXY_createFromFunction( int n, double *xs, ptwXY_createFromFunction_callback func, void *argList, double accuracy, int checkForRoots,
    int biSectionMax, nfu_status *status );
ptwXYPoints *ptwXY_createFromFunction2( ptwXPoints *xs, ptwXY_createFromFunction_callback func, void *argList, double accuracy, int checkForRoots,
    int biSectionMax, nfu_status *status );
nfu_status ptwXY_applyFunction( ptwXYPoints *ptwXY1, ptwXY_applyFunction_callback func, void *argList, int checkForRoots );
ptwXYPoints *ptwXY_fromString( char const *str, ptwXY_interpolation interpolation, ptwXY_interpolationOtherInfo const *interpolationOtherInfo, 
    double biSectionMax, double accuracy, char **endCharacter, nfu_status *status );

void ptwXY_showInteralStructure( ptwXYPoints *ptwXY, FILE *f, int printPointersAsNull );
void ptwXY_simpleWrite( ptwXYPoints *ptwXY, FILE *f, char *format );
void ptwXY_simplePrint( ptwXYPoints *ptwXY, char *format );

/* 
* Functions in ptwXY_integration.c 
*/
nfu_status ptwXY_f_integrate( ptwXY_interpolation interpolation, double x1, double y1, double x2, double y2, double *value );
double ptwXY_integrate( ptwXYPoints *ptwXY, double xMin, double xMax, nfu_status *status );
double ptwXY_integrateDomain( ptwXYPoints *ptwXY, nfu_status *status );
nfu_status ptwXY_normalize( ptwXYPoints *ptwXY1 );
double ptwXY_integrateDomainWithWeight_x( ptwXYPoints *ptwXY, nfu_status *status );
double ptwXY_integrateWithWeight_x( ptwXYPoints *ptwXY, double xMin, double xMax, nfu_status *status );
double ptwXY_integrateDomainWithWeight_sqrt_x( ptwXYPoints *ptwXY, nfu_status *status );
double ptwXY_integrateWithWeight_sqrt_x( ptwXYPoints *ptwXY, double xMin, double xMax, nfu_status *status );
ptwXPoints *ptwXY_groupOneFunction( ptwXYPoints *ptwXY, ptwXPoints *groupBoundaries, ptwXY_group_normType normType, ptwXPoints *ptwX_norm, nfu_status *status );
ptwXPoints *ptwXY_groupTwoFunctions( ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, ptwXPoints *groupBoundaries, ptwXY_group_normType normType, 
        ptwXPoints *ptwX_norm, nfu_status *status );
ptwXPoints *ptwXY_groupThreeFunctions( ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, ptwXYPoints *ptwXY3, ptwXPoints *groupBoundaries,
        ptwXY_group_normType normType, ptwXPoints *ptwX_norm, nfu_status *status );
ptwXPoints *ptwXY_runningIntegral( ptwXYPoints *ptwXY, nfu_status *status );
double ptwXY_integrateWithFunction( ptwXYPoints *ptwXY, ptwXY_createFromFunction_callback func, void *argList,
        double xMin, double xMax, int degree, int recursionLimit, double tolerance, nfu_status *status );

#if defined __cplusplus
    }
    }
#endif

#endif          /* End of ptwXY_h_included. */
