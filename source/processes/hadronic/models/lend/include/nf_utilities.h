/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#ifndef nf_utilities_h_included
#define nf_utilities_h_included

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdarg.h>

#define NUMERICALFUNCTIONS_SVN_VERSION 110+

#define nf_floatToShortestString_trimZeros   ( 1 << 0 )
#define nf_floatToShortestString_keepPeriod  ( 1 << 1 )
#define nf_floatToShortestString_includeSign ( 1 << 2 )

#if defined __cplusplus
    extern "C" {
    namespace GIDI {
#endif

typedef enum nfu_status_e {         nfu_Okay,           nfu_mallocError,            nfu_insufficientMemory,     
    nfu_badIndex,                   nfu_XNotAscending,  nfu_badIndexForX,           nfu_XOutsideDomain,             
    nfu_invalidInterpolation,       nfu_badSelf,        nfu_divByZero,              nfu_unsupportedInterpolationConversion, 
    nfu_unsupportedInterpolation,   nfu_empty,          nfu_tooFewPoints,           nfu_domainsNotMutual,                   
    nfu_badInput,                   nfu_badNorm,        nfu_badIntegrationInput,    nfu_otherInterpolation,
    nfu_failedToConverge,           nfu_oddNumberOfValues } nfu_status;

/*
* Functions in nf_utilities.c
*/
double nfu_getNAN( void );
int nfu_isNAN( double d );
double nfu_getInfinity( double sign );
const char *nfu_statusMessage( nfu_status status );
void nfu_setMemoryDebugMode( int mode );
void *nfu_malloc( size_t size );
void *nfu_calloc( size_t size, size_t n );
void *nfu_realloc( size_t size, void *old );
void *nfu_free( void *p );
void nfu_printMsg( char *fmt, ... );
void nfu_printErrorMsg( char *fmt, ... );

/*
* Functions in nf_stringToDoubles.c
*/
nfu_status nfu_stringToListOfDoubles( char const *str, int64_t *numberConverted, double **doublePtr, char **endCharacter );
char *nf_floatToShortestString( double value, int significantDigits, int favorEFormBy, int flags );

#if defined __cplusplus
    }
    }
#endif

#endif          /* End of nf_utilities_h_included. */
