/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
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
typedef enum ptwXY_interpolation_e { ptwXY_interpolationLinLin, ptwXY_interpolationLogLin, ptwXY_interpolationLinLog, 
    ptwXY_interpolationLogLog, ptwXY_interpolationFlat, ptwXY_interpolationOther } ptwXY_interpolation;

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
    ptwXY_lessEqualGreaterX_between, ptwXY_lessEqualGreaterX_greater, ptwXY_lessEqualGreaterX_Error } ptwXY_lessEqualGreaterX;

typedef
    struct ptwXYPoint_s {
        double x, y;
    } ptwXYPoint;

typedef nfu_status (*ptwXY_createFromFunction_callback)( statusMessageReporting *smr, double x, double *y, void *argList );
typedef nfu_status (*ptwXY_applyFunction_callback)( statusMessageReporting *smr, ptwXYPoint *point, void *argList );

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
        ptwXY_interpolation interpolation;
        char const *interpolationString;
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
ptwXYPoints *ptwXY_new( statusMessageReporting *smr, ptwXY_interpolation interpolation, char const *interpolationString, 
        double biSectionMax, double accuracy, int64_t primarySize, int64_t secondarySize, int userFlag );
ptwXYPoints *ptwXY_new2( statusMessageReporting *smr, ptwXY_interpolation interpolation, int64_t primarySize, int64_t secondarySize );
nfu_status ptwXY_initialize( statusMessageReporting *smr, ptwXYPoints *ptwXY, ptwXY_interpolation interpolation, 
        char const *interpolationString, double biSectionMax, double accuracy, int64_t primarySize, int64_t secondarySize, 
        int userFlag );
ptwXYPoints *ptwXY_create( statusMessageReporting *smr, ptwXY_interpolation interpolation, char const *interpolationString, 
        double biSectionMax, double accuracy, int64_t primarySize, int64_t secondarySize, int64_t length, double const *xy, 
        int userFlag );
ptwXYPoints *ptwXY_create2( statusMessageReporting *smr, ptwXY_interpolation interpolation, 
        int64_t primarySize, int64_t secondarySize, int64_t length, double const *xy, int userFlag );
ptwXYPoints *ptwXY_createFrom_Xs_Ys( statusMessageReporting *smr, ptwXY_interpolation interpolation, 
        char const *interpolationString, double biSectionMax, double accuracy, int64_t primarySize, int64_t secondarySize, 
        int64_t length, double const *Xs, double const *Ys, int userFlag );
ptwXYPoints *ptwXY_createFrom_Xs_Ys2( statusMessageReporting *smr, ptwXY_interpolation interpolation, 
        int64_t primarySize, int64_t secondarySize, int64_t length, double const *Xs, double const *Ys, int userFlag );

nfu_status ptwXY_copy( statusMessageReporting *smr, ptwXYPoints *dest, ptwXYPoints *src );
nfu_status ptwXY_copyPointsOnly( statusMessageReporting *smr, ptwXYPoints *dest, ptwXYPoints *src );
ptwXYPoints *ptwXY_clone( statusMessageReporting *smr, ptwXYPoints *ptwXY );
ptwXYPoints *ptwXY_clone2( statusMessageReporting *smr, ptwXYPoints const *ptwXY );
ptwXYPoints *ptwXY_cloneToInterpolation( statusMessageReporting *smr, ptwXYPoints *ptwXY, ptwXY_interpolation interpolationTo );
ptwXYPoints *ptwXY_slice( statusMessageReporting *smr, ptwXYPoints *ptwXY, int64_t index1, int64_t index2, int64_t secondarySize );
ptwXYPoints *ptwXY_domainSlice( statusMessageReporting *smr, ptwXYPoints *ptwXY, double domainMin, double domainMax, 
        int64_t secondarySize, int fill );
ptwXYPoints *ptwXY_domainMinSlice( statusMessageReporting *smr, ptwXYPoints *ptwXY, double domainMin, int64_t secondarySize, int fill );
ptwXYPoints *ptwXY_domainMaxSlice( statusMessageReporting *smr, ptwXYPoints *ptwXY, double domainMax, int64_t secondarySize, int fill );

ptwXY_interpolation ptwXY_getInterpolation( ptwXYPoints *ptwXY );
char const *ptwXY_getInterpolationString( ptwXYPoints *ptwXY );
nfu_status ptwXY_setInterpolationString( ptwXYPoints *ptwXY, char const *interpolationString );
nfu_status ptwXY_getStatus( ptwXYPoints *ptwXY );
int ptwXY_getUserFlag( ptwXYPoints *ptwXY );
void ptwXY_setUserFlag( ptwXYPoints *ptwXY, int userFlag );
double ptwXY_getAccuracy( ptwXYPoints *ptwXY );
double ptwXY_setAccuracy( ptwXYPoints *ptwXY, double accuracy );
double ptwXY_getBiSectionMax( ptwXYPoints *ptwXY );
double ptwXY_setBiSectionMax( ptwXYPoints *ptwXY, double biSectionMax );

nfu_status ptwXY_reallocatePoints( statusMessageReporting *smr, ptwXYPoints *ptwXY, int64_t size, int forceSmallerResize );
nfu_status ptwXY_reallocateOverflowPoints( statusMessageReporting *smr, ptwXYPoints *ptwXY, int64_t size );
nfu_status ptwXY_coalescePoints( statusMessageReporting *smr, ptwXYPoints *ptwXY, int64_t size, ptwXYPoint *newPoint, 
        int forceSmallerResize );
nfu_status ptwXY_simpleCoalescePoints( statusMessageReporting *smr, ptwXYPoints *ptwXY );

nfu_status ptwXY_clear( statusMessageReporting *smr, ptwXYPoints *ptwXY );
nfu_status ptwXY_release( statusMessageReporting *smr, ptwXYPoints *ptwXY );
ptwXYPoints *ptwXY_free( ptwXYPoints *ptwXY );

int64_t ptwXY_length( statusMessageReporting *smr, ptwXYPoints *ptwXY );
int64_t ptwXY_getNonOverflowLength( statusMessageReporting *smr, ptwXYPoints const *ptwXY );

nfu_status ptwXY_startIndex( statusMessageReporting *a_smr, ptwXYPoints *a_ptwXY, double a_x, int64_t *a_startIndex, int64_t *a_length );
nfu_status ptwXY_setXYData( statusMessageReporting *smr, ptwXYPoints *ptwXY, int64_t length, double const *xy );
nfu_status ptwXY_setXYDataFromXsAndYs( statusMessageReporting *smr, ptwXYPoints *ptwXY, int64_t length, double const *x, double const *y );
nfu_status ptwXY_deletePoints( statusMessageReporting *smr, ptwXYPoints *ptwXY, int64_t i1, int64_t i2 );
nfu_status ptwXY_getLowerIndexBoundingX( statusMessageReporting *smr, ptwXYPoints *ptwXY, double x, int64_t *index );
ptwXYPoint *ptwXY_getPointAtIndex( statusMessageReporting *smr, ptwXYPoints *ptwXY, int64_t index );
ptwXYPoint *ptwXY_getPointAtIndex_Unsafely( ptwXYPoints const *ptwXY, int64_t index );
nfu_status ptwXY_getXYPairAtIndex( statusMessageReporting *smr, ptwXYPoints *ptwXY, int64_t index, double *x, double *y );
ptwXY_lessEqualGreaterX ptwXY_getPointsAroundX( statusMessageReporting *smr, ptwXYPoints *ptwXY, double x, 
        ptwXYOverflowPoint *lessThanEqualXPoint, ptwXYOverflowPoint *greaterThanXPoint );
ptwXY_lessEqualGreaterX ptwXY_getPointsAroundX_closeIsEqual( statusMessageReporting *smr, ptwXYPoints *ptwXY, double x, 
        ptwXYOverflowPoint *lessThanEqualXPoint, ptwXYOverflowPoint *greaterThanXPoint, double eps, int *closeIsEqual, 
        ptwXYPoint **closePoint );
nfu_status ptwXY_getValueAtX( statusMessageReporting *smr, ptwXYPoints *ptwXY, double x, double *y );
nfu_status ptwXY_setValueAtX( statusMessageReporting *smr, ptwXYPoints *ptwXY, double x, double y );
nfu_status ptwXY_setValueAtX_overrideIfClose( statusMessageReporting *smr, ptwXYPoints *ptwXY, double x, double y, double eps, int override );
nfu_status ptwXY_mergeFromXsAndYs( statusMessageReporting *smr, ptwXYPoints *ptwXY, int length, double *xs, double *ys );
nfu_status ptwXY_mergeFromXYs( statusMessageReporting *smr, ptwXYPoints *ptwXY, int length, double *xys );
nfu_status ptwXY_appendXY( statusMessageReporting *smr, ptwXYPoints *ptwXY, double x, double y );
nfu_status ptwXY_setXYPairAtIndex( statusMessageReporting *smr, ptwXYPoints *ptwXY, int64_t index, double x, double y );

nfu_status ptwXY_getSlopeAtX( statusMessageReporting *smr, ptwXYPoints *ptwXY, double x, char side, double *slope );

nfu_status ptwXY_domainMinAndFrom( statusMessageReporting *smr, ptwXYPoints *ptwXY, ptwXY_dataFrom *dataFrom, double *value );
nfu_status ptwXY_domainMin( statusMessageReporting *smr, ptwXYPoints *ptwXY, double *value );
nfu_status ptwXY_domainMaxAndFrom( statusMessageReporting *smr, ptwXYPoints *ptwXY, ptwXY_dataFrom *dataFrom, double *value );
nfu_status ptwXY_domainMax( statusMessageReporting *smr, ptwXYPoints *ptwXY, double *value );
nfu_status ptwXY_range( statusMessageReporting *smr, ptwXYPoints *ptwXY, double *rangeMin, double *rangeMax );
nfu_status ptwXY_rangeMin( statusMessageReporting *smr, ptwXYPoints *ptwXY, double *value );
nfu_status ptwXY_rangeMax( statusMessageReporting *smr, ptwXYPoints *ptwXY, double *value );
char const *ptwXY_interpolationToString( ptwXY_interpolation interpolation );
ptwXY_interpolation ptwXY_stringToInterpolation( char const *interpolationString );

/* 
* Methods in ptwXY_methods.c 
*/
nfu_status ptwXY_clip( statusMessageReporting *smr, ptwXYPoints *ptwXY1, double rangeMin, double rangeMax );
nfu_status ptwXY_thicken( statusMessageReporting *smr, ptwXYPoints *ptwXY1, int sectionSubdivideMax, 
        double dDomainMax, double fDomainMax );
ptwXYPoints *ptwXY_thin( statusMessageReporting *smr, ptwXYPoints *ptwXY1, double accuracy );
ptwXYPoints *ptwXY_thinDomain( statusMessageReporting *smr, ptwXYPoints *ptwXY1, double epsilon );
nfu_status ptwXY_trim( statusMessageReporting *smr, ptwXYPoints *ptwXY );

ptwXYPoints *ptwXY_union( statusMessageReporting *smr, ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, int unionOptions );

nfu_status ptwXY_scaleOffsetXAndY( statusMessageReporting *smr, ptwXYPoints *ptwXY, double xScale, double xOffset, 
        double yScale, double yOffset );
nfu_status ptwXY_scaleAndOffsetDomainWith_ptwXYs( statusMessageReporting *smr, ptwXYPoints *ptwXY, ptwXYPoints *offset, ptwXYPoints *slope, int skipLastPoint );

/*
* Functions in ptwXY_unitaryOperators.c
*/
nfu_status ptwXY_abs( statusMessageReporting *smr, ptwXYPoints *ptwXY );
nfu_status ptwXY_neg( statusMessageReporting *smr, ptwXYPoints *ptwXY );

/*
* Functions in ptwXY_binaryOperators.c
*/
nfu_status ptwXY_slopeOffset( statusMessageReporting *smr, ptwXYPoints *ptwXY, double slope, double offset );
nfu_status ptwXY_add_double( statusMessageReporting *smr, ptwXYPoints *ptwXY, double value );
nfu_status ptwXY_sub_doubleFrom( statusMessageReporting *smr, ptwXYPoints *ptwXY, double value );
nfu_status ptwXY_sub_fromDouble( statusMessageReporting *smr, ptwXYPoints *ptwXY, double value );
nfu_status ptwXY_mul_double( statusMessageReporting *smr, ptwXYPoints *ptwXY, double value );
nfu_status ptwXY_div_doubleFrom( statusMessageReporting *smr, ptwXYPoints *ptwXY, double value );
nfu_status ptwXY_div_fromDouble( statusMessageReporting *smr, ptwXYPoints *ptwXY, double value );
nfu_status ptwXY_mod( statusMessageReporting *smr, ptwXYPoints *ptwXY, double m, int pythonMod );

ptwXYPoints *ptwXY_binary_ptwXY( statusMessageReporting *smr, ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, 
        double v1, double v2, double v1v2 );
ptwXYPoints *ptwXY_add_ptwXY( statusMessageReporting *smr, ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2 );
ptwXYPoints *ptwXY_sub_ptwXY( statusMessageReporting *smr, ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2 );
ptwXYPoints *ptwXY_mul_ptwXY( statusMessageReporting *smr, ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2 );
ptwXYPoints *ptwXY_mul2_ptwXY( statusMessageReporting *smr, ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2 );
ptwXYPoints *ptwXY_div_ptwXY( statusMessageReporting *smr, ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, int safeDivide );

/* 
* Functions in ptwXY_functions.c 
*/
nfu_status ptwXY_pow( statusMessageReporting *smr, ptwXYPoints *ptwXY, double p );
nfu_status ptwXY_exp( statusMessageReporting *smr, ptwXYPoints *ptwXY, double a );
ptwXYPoints *ptwXY_convolution( statusMessageReporting *smr, ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, int mode );
ptwXYPoints *ptwXY_inverse( statusMessageReporting *smr, ptwXYPoints *ptwXY );

/*
* Functions in ptwXY_interpolation.c
*/
nfu_status ptwXY_interpolatePoint( statusMessageReporting *smr, ptwXY_interpolation interpolation, double x, double *y, 
        double x1, double y1, double x2, double y2 );
ptwXYPoints *ptwXY_flatInterpolationToLinear( statusMessageReporting *smr, ptwXYPoints *ptwXY, double lowerEps, double upperEps );
ptwXYPoints *ptwXY_toOtherInterpolation( statusMessageReporting *smr, ptwXYPoints *ptwXY, ptwXY_interpolation interpolation, 
        double accuracy );
ptwXYPoints *ptwXY_unitbaseInterpolate( statusMessageReporting *smr, double w, double w1, ptwXYPoints *ptwXY1, 
        double w2, ptwXYPoints *ptwXY2, int scaleRange );
ptwXYPoints *ptwXY_toUnitbase( statusMessageReporting *smr, ptwXYPoints *ptwXY, int scaleRange );
ptwXYPoints *ptwXY_fromUnitbase( statusMessageReporting *smr, ptwXYPoints *ptwXY, double domainMin, double domainMax, 
        int scaleRange );

/* 
* Functions in ptwXY_convenient.c 
*/
ptwXPoints *ptwXY_getXArray( statusMessageReporting *smr, ptwXYPoints *ptwXY );
ptwXPoints *ptwXY_ysMappedToXs( statusMessageReporting *smr, ptwXYPoints *ptwXY, ptwXPoints *Xs, int64_t *offset );
nfu_status ptwXY_mapToXsAndAdd( statusMessageReporting *a_smr, ptwXYPoints *a_ptwXY, int64_t a_offset, int64_t a_length, double const *a_Xs,
        double *a_results, double a_scaleFractor );
nfu_status ptwXY_dullEdges( statusMessageReporting *smr, ptwXYPoints *ptwXY, double lowerEps, double upperEps, int positiveXOnly );
nfu_status ptwXY_mergeClosePoints( statusMessageReporting *smr, ptwXYPoints *ptwXY, double epsilon );
ptwXYPoints *ptwXY_intersectionWith_ptwX( statusMessageReporting *smr, ptwXYPoints *ptwXY, ptwXPoints *ptwX );
nfu_status ptwXY_areDomainsMutual( statusMessageReporting *smr, ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2 );
nfu_status ptwXY_tweakDomainsToMutualify( statusMessageReporting *smr, ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, 
        int epsilonFactor, double epsilon );
nfu_status ptwXY_mutualifyDomains( statusMessageReporting *smr, ptwXYPoints *ptwXY1, double lowerEps1, double upperEps1, 
        int positiveXOnly1, ptwXYPoints *ptwXY2, double lowerEps2, double upperEps2, int positiveXOnly2 );
nfu_status ptwXY_copyToC_XY( statusMessageReporting *smr, ptwXYPoints *ptwXY, int64_t index1, int64_t index2, 
        int64_t allocatedSize, int64_t *numberOfPoints, double *xy );
nfu_status ptwXY_valuesToC_XsAndYs( statusMessageReporting *smr, ptwXYPoints *ptwXY, double **xs, double **ys );
ptwXYPoints *ptwXY_valueTo_ptwXY( statusMessageReporting *smr, double x1, double x2, double y );
ptwXYPoints *ptwXY_createGaussianCenteredSigma1( statusMessageReporting *smr, double accuracy );
ptwXYPoints *ptwXY_createGaussian( statusMessageReporting *smr, double accuracy, double xCenter, double sigma, 
        double amplitude, double domainMin, double domainMax, double dullEps );

/* 
* Functions in ptwXY_misc.c 
*/
double ptwXY_limitAccuracy( double accuracy );
void ptwXY_update_biSectionMax( ptwXYPoints *ptwXY1, double oldLength );
ptwXYPoints *ptwXY_createFromFunction( statusMessageReporting *smr, int n, double *xs, 
        ptwXY_createFromFunction_callback func, void *argList, double accuracy, int checkForRoots, int biSectionMax );
ptwXYPoints *ptwXY_createFromFunction2( statusMessageReporting *smr, ptwXPoints *xs, ptwXY_createFromFunction_callback func, 
        void *argList, double accuracy, int checkForRoots, int biSectionMax );
nfu_status ptwXY_applyFunction( statusMessageReporting *smr, ptwXYPoints *ptwXY1, ptwXY_applyFunction_callback func, 
        void *argList, int checkForRoots );
ptwXYPoints *ptwXY_fromString( statusMessageReporting *smr, char const *str, char sep, ptwXY_interpolation interpolation, char const *interpolationString, 
        double biSectionMax, double accuracy, char **endCharacter, int useSystem_strtod );

void ptwXY_showInteralStructure( ptwXYPoints *ptwXY, FILE *f, int printPointersAsNull );
void ptwXY_simpleWrite( ptwXYPoints *ptwXY, FILE *f, char const *format );
void ptwXY_simplePrint( ptwXYPoints *ptwXY, char const *format );

/* 
* Functions in ptwXY_integration.c 
*/
nfu_status ptwXY_f_integrate( statusMessageReporting *smr, ptwXY_interpolation interpolation, double x1, double y1, 
        double x2, double y2, double *value );
nfu_status ptwXY_integrate( statusMessageReporting *smr, ptwXYPoints *ptwXY, double domainMin, double domainMax, double *value );
nfu_status ptwXY_integrateDomain( statusMessageReporting *smr, ptwXYPoints *ptwXY, double *value );
nfu_status ptwXY_normalize( statusMessageReporting *smr, ptwXYPoints *ptwXY1 );
nfu_status ptwXY_integrateDomainWithWeight_x( statusMessageReporting *smr, ptwXYPoints *ptwXY, double *value );
nfu_status ptwXY_integrateWithWeight_x( statusMessageReporting *smr, ptwXYPoints *ptwXY, double domainMin, double domainMax, 
        double *value );
nfu_status ptwXY_integrateDomainWithWeight_sqrt_x( statusMessageReporting *smr, ptwXYPoints *ptwXY, double *value );
nfu_status ptwXY_integrateWithWeight_sqrt_x( statusMessageReporting *smr, ptwXYPoints *ptwXY, double domainMin, double domainMax, 
        double *value );
ptwXPoints *ptwXY_groupOneFunction( statusMessageReporting *smr, ptwXYPoints *ptwXY, ptwXPoints *groupBoundaries, 
        ptwXY_group_normType normType, ptwXPoints *ptwX_norm );
ptwXPoints *ptwXY_groupTwoFunctions( statusMessageReporting *smr, ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, 
        ptwXPoints *groupBoundaries, ptwXY_group_normType normType, ptwXPoints *ptwX_norm );
ptwXPoints *ptwXY_groupThreeFunctions( statusMessageReporting *smr, ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, 
        ptwXYPoints *ptwXY3, ptwXPoints *groupBoundaries, ptwXY_group_normType normType, ptwXPoints *ptwX_norm );
ptwXPoints *ptwXY_runningIntegral( statusMessageReporting *smr, ptwXYPoints *ptwXY );
nfu_status ptwXY_integrateWithFunction( statusMessageReporting *smr, ptwXYPoints *ptwXY, ptwXY_createFromFunction_callback func, 
        void *argList, double domainMin, double domainMax, int degree, int recursionLimit, double tolerance,
        double *value );
ptwXPoints *ptwXY_equalProbableBins( statusMessageReporting *smr, ptwXYPoints *ptwXY, int numberOfBins );

#if defined __cplusplus
    }
#endif

#endif          /* End of ptwXY_h_included. */
