/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#ifndef nf_utilities_h_included
#define nf_utilities_h_included

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>

#include <statusMessageReporting.h>

extern int nfu_SMR_libraryID;

#ifdef __APPLE__
#ifndef NFU_USEHEAP
#define NFU_USEHEAP 1
#endif
#endif

#define nf_floatToShortestString_trimZeros   ( 1 << 0 )
#define nf_floatToShortestString_keepPeriod  ( 1 << 1 )
#define nf_floatToShortestString_includeSign ( 1 << 2 )

#if defined __cplusplus
    extern "C" {
#endif

typedef enum nfu_status_e {         
    nfu_Okay,               
    nfu_Error,               
    nfu_badSelf,            
    nfu_mallocError,            
    nfu_insufficientMemory,     
    nfu_badIndex,                   
    nfu_XNotAscending,      
    nfu_badIndexForX,           
    nfu_XOutsideDomain,             
    nfu_invalidInterpolation,       
    nfu_divByZero,              
    nfu_unsupportedInterpolationConversion, 
    nfu_unsupportedInterpolation,   
    nfu_empty,              
    nfu_tooFewPoints,           
    nfu_domainsNotMutual,                   
    nfu_badInput,                   
    nfu_badNorm,            
    nfu_badIntegrationInput,    
    nfu_otherInterpolation,
    nfu_flatInterpolation,
    nfu_failedToConverge,           
    nfu_oddNumberOfValues,  
    nfu_badLogValue 
} nfu_status;

/*
* Functions in nf_utilities.c
*/
int nfu_setup( void );
double nfu_getNAN( void );
int nfu_isNAN( double d );
double nfu_getInfinity( double sign );
const char *nfu_statusMessage( nfu_status status );
void nfu_setMemoryDebugMode( int mode );
void nfu_printMsg( char const *fmt, ... );
void nfu_printErrorMsg( char const *fmt, ... );

/*
* These function are to be deleted when conversion to statusMessageReporting is completed.
*/
void *nfu_malloc( size_t size );
void *nfu_calloc( size_t size, size_t n );
void *nfu_realloc( size_t size, void *old );
void *nfu_free( void *p );
/*
* Functions in nf_stringToInt32s.c
*/
int32_t *nfu_stringToListOfInt32s( statusMessageReporting *smr, char const *str, char sep, int64_t *numberConverted, char **endCharacter );
int nfu_stringToInt32( statusMessageReporting *smr, char const *str, char **endCharacter, int32_t *value );
/*
* Functions in nf_stringToDoubles.c
*/
double *nfu_stringToListOfDoubles( statusMessageReporting *smr, char const *str, char sep, int64_t *numberConverted, 
        char **endCharacter, int useSystem_strtod );
double nf_strtod( char const *ptr, char **endCharacter );
char *nf_floatToShortestString( double value, int significantDigits, int favorEFormBy, int flags );

#if defined __cplusplus
    }
#endif

#endif          /* End of nf_utilities_h_included. */
