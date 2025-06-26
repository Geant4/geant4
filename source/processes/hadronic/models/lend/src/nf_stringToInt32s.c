/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include <stdlib.h>
#include <stdint.h>
#include <ctype.h>

#include "nf_utilities.h"

#define numberOfStaticInt32s ( 100 * 1000 )

#ifndef INT32_MIN
#define INT32_MIN -2147483648
#define INT32_MAX 2147483647
#endif

static int32_t *nfu_stringToListOfInt32s_2( statusMessageReporting *smr, char const *str, char sep, int64_t *numberConverted, char **endCharacter );
/*
========================================================================
*/
int32_t *nfu_stringToListOfInt32s( statusMessageReporting *smr, char const *str, char sep, int64_t *numberConverted, 
        char **endCharacter ) {

    if( strchr( "0123456789.+-eE", sep ) != NULL ) {
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_badInput, "Invalid sep ='%c'.", sep );
        return( NULL );
    }

    *numberConverted = 0;
    *endCharacter = (char *) str;
    if( isspace( sep ) ) sep = ' ';             /* Make it the space character if any white space as it simplifies logic below. */
    return( nfu_stringToListOfInt32s_2( smr, str, sep, numberConverted, endCharacter ) );
}
/*
========================================================================
*/
static int32_t *nfu_stringToListOfInt32s_2( statusMessageReporting *smr, char const *str, char sep, int64_t *numberConverted, 
        char **endCharacter ) {

    int64_t i1, i2, numberConverted_initial = *numberConverted;
    int32_t *Int32Ptr = NULL;
#if NFU_USEHEAP
    int32_t *staticInt32s = (int32_t *) smr_malloc2( smr, (size_t) numberOfStaticInt32s * sizeof( int32_t ), 0, "staticInt32s" );
    if( staticInt32s == NULL ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( NULL );
    }
#else
    int32_t staticInt32s[numberOfStaticInt32s];
#endif

    for( i1 = 0; i1 < numberOfStaticInt32s; i1++, (*numberConverted)++ ) {
        if(  *numberConverted == 0 ) {
            if( nfu_stringToInt32( smr, str, endCharacter, &staticInt32s[i1] ) != 0 ) {
                *endCharacter = (char *) str;
                smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
#if NFU_USEHEAP
                free( staticInt32s );
#endif
                return( NULL );
            } }
        else {                                  /* Check that there is one sep character and allow for arbitrary number of white spaces. */
            char const *str2 = str;

            while( isspace( *str2 ) ) ++str2;   /* Only need to check for white spaces before sep character as strtol will ignore after. */
            if( sep != ' ' ) {
                if( *str2 == sep ) {
                    ++str2; }
                else {
                    str2 = str;
                }
            }
            if( str < str2 ) {
                if( nfu_stringToInt32( smr, str2, endCharacter, &staticInt32s[i1] ) != 0 ) {
#if NFU_USEHEAP
                    free( staticInt32s );
#endif
                    *endCharacter = (char *) str;
                    smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
                    return( NULL );
                }
            }
            if( str2 == (char const *) *endCharacter ) *endCharacter = (char *) str;
        }
        if( str == (char const *) *endCharacter ) {
            int64_t number = *numberConverted;

            if( *numberConverted == 0 ) number = 1;
            if( ( Int32Ptr = (int32_t *) smr_malloc2( smr, (size_t) number * sizeof( int32_t ), 0, "Int32Ptr" ) ) == NULL ) {
#if NFU_USEHEAP
                free( staticInt32s );
#endif
                smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
                return( NULL );
            }
            break;
        }
        str = (char const *) *endCharacter;
    }

    if( Int32Ptr == NULL ) Int32Ptr = nfu_stringToListOfInt32s_2( smr, str, sep, numberConverted, endCharacter );
    if( Int32Ptr != NULL ) {
        int32_t *Int32Ptr2 = &(Int32Ptr[numberConverted_initial]);
        char *end = *endCharacter;

        for( i2 = 0; i2 < i1; i2++, Int32Ptr2++ ) *Int32Ptr2 = staticInt32s[i2];
        while( isspace( *end ) ) ++end;
        if( *end == 0 ) *endCharacter = end;
    }

#if NFU_USEHEAP
    free( staticInt32s );
#endif

    return( Int32Ptr );
}
/*
========================================================================
*/
int nfu_stringToInt32( statusMessageReporting *smr, char const *str, char **endCharacter, int32_t *value ) {

    long lValue = strtol( str, endCharacter, 10 );

    if( lValue < INT32_MIN ) {
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_badInput, "int32_t underflow: %l", lValue );
        return( -1 ); }
    else if( lValue > INT32_MAX ) {
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_badInput, "int32_t overflow: %l", lValue );
        return( 1 );
    }
    *value = (int) lValue;
    return( 0 );
}
