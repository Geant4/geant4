/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#include <stdio.h>
#include <stdlib.h>
#include <cmath>

#include "nf_utilities.h"

#ifdef WIN32
#include <float.h>
#define is_nan(a) _isnan(a)
/*#define INFINITY (DBL_MAX+DBL_MAX)*/
/*#define NAN (INFINITY-INFINITY)*/
#else
#define is_nan(a) std::isnan(a)
#endif

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

static const char Okay_message[] = "all is okay";
static const char mallocError_message[] = "could not allocate memory";
static const char insufficientMemory_message[] = "user's memory is too small to handle data";
static const char badIndex_message[] = "bad index";
static const char XNotAscending_message[] = "x values are not ascending";
static const char badIndexForX_message[] = "index not correct for x value";
static const char XOutsideDomain_message[] = "x value not in domain";
static const char invalidInterpolation_message[] = "bad x,y values for interpolation";
static const char badSelf_message[] = "source object has bad status value";
static const char divByZero_message[] = "division by zero";
static const char unsupportedInterpolation_message[] = "unsupported interpolation";
static const char unsupportedInterpolationConversion_message[] = "unsupported interpolation conversion";
static const char empty_message[] = "empty instance";
static const char tooFewPoints_message[] = "too few points in instance";
static const char notMutualDomian_message[] = "domains are not mutual";
static const char unknownStatus_message[] = "unknown (i.e., invalid) status value";
static const char badInput_message[] = "bad input to function";
static const char badNorm_message[] = "bad norm";
static const char badIntegrationInput_message[] = "bad integration input";
static const char otherInterpolation_message[] = "other integration not supported";
static const char failedToConverge_message[] = "failed to converge";
static const char oddNumberOfValues_message[] = "odd number of inputted values";

static int nfu_debugging = 0;

/*
************************************************************
*/
double nfu_getNAN( void ) {

    return( NAN );
}
/*
************************************************************
*/
int nfu_isNAN( double d ) {

    return( is_nan( d ) );
}
/*
************************************************************
*/
double nfu_getInfinity( double sign ) {

    if( sign < 0 ) return( -INFINITY );
    return( INFINITY );
}
/*
************************************************************
*/
const char *nfu_statusMessage( nfu_status status ) {

    switch( status ) {
    case nfu_Okay : return( Okay_message );
    case nfu_mallocError : return( mallocError_message );
    case nfu_insufficientMemory : return( insufficientMemory_message );
    case nfu_badIndex : return( badIndex_message );
    case nfu_XNotAscending : return( XNotAscending_message );
    case nfu_badIndexForX : return( badIndexForX_message );
    case nfu_XOutsideDomain : return( XOutsideDomain_message );
    case nfu_invalidInterpolation : return( invalidInterpolation_message );
    case nfu_badSelf : return( badSelf_message );
    case nfu_divByZero : return( divByZero_message );
    case nfu_unsupportedInterpolation : return( unsupportedInterpolation_message );
    case nfu_unsupportedInterpolationConversion : return( unsupportedInterpolationConversion_message );
    case nfu_empty : return( empty_message );
    case nfu_tooFewPoints : return( tooFewPoints_message );
    case nfu_domainsNotMutual : return( notMutualDomian_message );
    case nfu_badInput : return( badInput_message );
    case nfu_badNorm : return( badNorm_message );
    case nfu_badIntegrationInput : return( badIntegrationInput_message );
    case nfu_otherInterpolation : return( otherInterpolation_message );
    case nfu_failedToConverge : return( failedToConverge_message );
    case nfu_oddNumberOfValues : return( oddNumberOfValues_message );
    }
    return( unknownStatus_message );
}
/*
************************************************************
*/
void nfu_setMemoryDebugMode( int mode ) {

    nfu_debugging = mode;
}
/*
************************************************************
*/
void *nfu_malloc( size_t size ) {

    void *p = malloc( size );

    if( nfu_debugging ) printf( "nfu_malloc  %12p size = %8llu\n", p, (long long unsigned) size );
    return( p );
}
/*
************************************************************
*/
void *nfu_calloc( size_t size, size_t n ) {

    void *p = calloc( size, n );

    if( nfu_debugging ) printf( "nfu_calloc  %12p size = %8llu, n = %8llu\n", p, (long long unsigned) size, (long long unsigned) n );
    return( p );
}
/*
************************************************************
*/
void *nfu_realloc( size_t size, void *old ) {

    void *p = realloc( old, size );

    if( nfu_debugging ) printf( "nfu_realloc %12p size = %8llu", p, (long long unsigned) size );
    return( p );
}
/*
************************************************************
*/
void *nfu_free( void *p ) {

    if( p != NULL ) {
        if( nfu_debugging ) printf( "nfu_free    %12p\n", p );
        free( p );
    }
    return( NULL );
}
/*
********************************************************
*/
void nfu_printMsg( char *fmt, ... ) {

    va_list args;

    va_start( args, fmt );
    vfprintf( stderr, fmt, args );
    fprintf( stderr, "\n" );
    va_end( args );
}
/*
********************************************************
*/
void nfu_printErrorMsg( char *fmt, ... ) {

    va_list args;

    va_start( args, fmt );
    vfprintf( stderr, fmt, args );
    fprintf( stderr, "\n" );
    va_end( args );

    exit( EXIT_FAILURE );
}

#if defined __cplusplus
}
#endif
