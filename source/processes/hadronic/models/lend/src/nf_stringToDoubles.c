/*
# <<BEGIN-copyright>>
# Copyright 2019, Lawrence Livermore National Security, LLC.
# This file is part of the gidiplus package (https://github.com/LLNL/gidiplus).
# gidiplus is licensed under the MIT license (see https://opensource.org/licenses/MIT).
# SPDX-License-Identifier: MIT
# <<END-copyright>>
*/

#include <stdlib.h>
#include <ctype.h>

#include "nf_utilities.h"

#ifdef _WIN32
#include <float.h>
#endif

#define numberOfStaticDoubles ( 100 * 1000 )

static double *nfu_stringToListOfDoubles2( statusMessageReporting *smr, char const *str, char sep, int64_t *numberConverted, 
        char **endCharacter, int useSystem_strtod );
/*
========================================================================
*/
double *nfu_stringToListOfDoubles( statusMessageReporting *smr, char const *str, char sep, int64_t *numberConverted, 
        char **endCharacter, int useSystem_strtod ) {

    if( strchr( "0123456789.+-eE", sep ) != NULL ) {
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_badInput, "Invalid sep ='%c'.", sep );
        return( NULL );
    }

    *numberConverted = 0;
    *endCharacter = (char *) str;
    if( isspace( sep ) ) sep = ' ';             /* Make it the space character if any white space as it simplifies logic below. */
    return( nfu_stringToListOfDoubles2( smr, str, sep, numberConverted, endCharacter, useSystem_strtod ) );
}
/*
========================================================================
*/
static double *nfu_stringToListOfDoubles2( statusMessageReporting *smr, char const *str, char sep, int64_t *numberConverted, 
        char **endCharacter, int useSystem_strtod ) {

    int64_t i1, i2, numberConverted_initial = *numberConverted;
    double *doublePtr = NULL;
    nfu_status status = nfu_Okay;
    double (*_strtod)( char const *str, char **endCharacter );

#if NFU_USEHEAP
    double *staticDoubles = (double *) smr_malloc2( smr, (size_t) numberOfStaticDoubles * sizeof( double ), 0, "staticDoubles" );
    if( staticDoubles == NULL ) {
        smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
        return( NULL );
    }
#else
    double staticDoubles[numberOfStaticDoubles];
#endif

    _strtod = useSystem_strtod ? strtod : nf_strtod;

    for( i1 = 0; i1 < numberOfStaticDoubles; i1++, (*numberConverted)++ ) {
        if(  *numberConverted == 0 ) {
            staticDoubles[i1] = _strtod( str, endCharacter ); }
        else {                                  /* Check that there is one sep character and allow for arbitrary number of white spaces. */
            char const *str2 = str;

            while( isspace( *str2 ) ) ++str2;   /* Only need to check for white spaces before sep character as strtod will ignore after. */
            if( sep != ' ' ) {
                if( *str2 == sep ) {
                    ++str2; }
                else {
                    str2 = str;
                }
            }
            if( str < str2 ) staticDoubles[i1] = _strtod( str2, endCharacter );
            if( str2 == (char const *) *endCharacter ) *endCharacter = (char *) str;
        }
        if( str == (char const *) *endCharacter ) {
            int64_t number = *numberConverted;
            if( *numberConverted == 0 ) number = 1;
            if( ( doublePtr = (double *) smr_malloc2( smr, (size_t) number * sizeof( double ), 0, "doublePtr" ) ) == NULL ) {
                smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
#if NFU_USEHEAP
                free( staticDoubles );
#endif
                return( NULL );
            }
            break;
        }
        str = (char const *) *endCharacter;
    }

    if( ( status == nfu_Okay ) && ( doublePtr == NULL ) )
        doublePtr = nfu_stringToListOfDoubles2( smr, str, sep, numberConverted, endCharacter, useSystem_strtod );
    if( doublePtr != NULL ) {
        double *doublePtr2 = &(doublePtr[numberConverted_initial]);
        char *end = *endCharacter;

        for( i2 = 0; i2 < i1; i2++, doublePtr2++ ) *doublePtr2 = staticDoubles[i2];
        while( isspace( *end ) ) ++end;
        if( *end == 0 ) *endCharacter = end;
    }

#if NFU_USEHEAP
    free( staticDoubles );
#endif
    return( doublePtr );
}

#define valid_digit( c ) ( ( c ) >= '0' && ( c ) <= '9')
/*
============================================================
*/
double nf_strtod( char const *str, char **endCharacter ) {

    *endCharacter = (char *) str;
    char *ptr = *endCharacter;

    while( isspace( *ptr ) ) ++ptr;                         /* Skip leading white space, if any. */

    double sign = 1.0;                                      /* Get sign, if any. */
    if( *ptr == '-' ) {
        sign = -1.0;
        ++ptr; } 
    else if( *ptr == '+' ) {
        ++ptr;
    }

    double value = 0.0;                                     /* Get digits before decimal point or exponent, if any. */
    for( ; valid_digit( *ptr ); ++ptr ) value = value * 10.0 + ( *ptr - '0' );

    if( *ptr == '.' ) {                                     /* Get digits after decimal point, if any. */
        double invPow10 = 0.1;
        ++ptr;
        while( valid_digit( *ptr ) ) {
            value += ( *ptr - '0' ) * invPow10;
            invPow10 *= 0.1;
            ++ptr;
        }
    }

    if( ( *ptr == 'e' ) || ( *ptr == 'E' ) ) {              /* Handle exponent, if any. */
        int negativeExponent = 0;

        ++ptr;                                            /* Get sign of exponent, if any. */
        if( *ptr == '-' ) {
            negativeExponent = 1;
            ++ptr; }
        else if( *ptr == '+' ) {
            ++ptr;
        }

        unsigned int exponent = 0;                      /* Get digits of exponent. There must be at least 1. */
        for( ; valid_digit( *ptr ); ++ptr ) exponent = exponent * 10 + ( *ptr - '0' );
        if( exponent > 308 ) {
            return( strtod( str, endCharacter ) );
        }

        double scale = 1.0;                             /* Calculate scaling factor. */
        if( exponent == 0 ) negativeExponent = 0;
        while( exponent >= 50 ) { scale *= 1E50; exponent -= 50; }
        while( exponent >=  8 ) { scale *= 1E8;  exponent -=  8; }
        while( exponent >   0 ) { scale *= 10.0; exponent -=  1; }

        if( negativeExponent ) scale = 1.0 / scale;
        value *= scale;
    }

    *endCharacter = ptr;
    return( sign * value );
}
/*
============================================================
*/
char *nf_floatToShortestString( double value, int significantDigits, int favorEFormBy, int flags ) {

    int n1, ne, nf, digitsRightOfPeriod_f, exponent;
    char Str_e[512], Str_f[512], *Str_r = Str_e, Fmt[32], *e1, *e2;
    const char *sign = "";

    if( flags & nf_floatToShortestString_includeSign ) sign = "+";

    if( !isfinite( value ) ) {
        sprintf( Fmt, "%%%sf", sign );
        sprintf( Str_e, Fmt, value );
        return( strdup( Str_e ) );
    }

    significantDigits--;
    if( significantDigits < 0 ) significantDigits = 0;
    if( significantDigits > 24 ) significantDigits = 24;

    sprintf( Fmt, "%%%s.%de", sign, significantDigits );
    sprintf( Str_e, Fmt, value );

    e1 = strchr( Str_e, 'e' );
    if( significantDigits == 0 ) {
        if( *(e1 - 1) != '.' ) {
            char *e3;

            e2 = strchr( e1, 0 );
            e3 = e2 + 1;
            for( ; e2 != e1; e2--, e3-- ) *e3 = *e2;
            *(e1++) = '.';
        }
    }
    *e1 = 0;
    n1 = (int) strlen( Str_e ) - 1;
    if( flags & nf_floatToShortestString_trimZeros ) while( Str_e[n1] == '0' ) n1--;
    ne = flags & nf_floatToShortestString_keepPeriod;
    if( !( flags & nf_floatToShortestString_keepPeriod ) ) if( Str_e[n1] == '.' ) n1--;
    n1++;
    Str_e[n1] = 0;

    e1++;
    exponent = (int) strtol( e1, &e2, 10 );
    if( exponent != 0 ) {               /* If 0, the exponent was "e+00". */
        for( e1 = Str_e; *e1 != 0; e1++ ) ;
        sprintf( e1, "e%d", exponent );

        digitsRightOfPeriod_f = significantDigits - exponent;
        if( ( digitsRightOfPeriod_f > 25 ) || ( exponent > 50 ) ) return( strdup( Str_r ) );
        if( digitsRightOfPeriod_f < 0 ) digitsRightOfPeriod_f = 0;

        sprintf( Fmt, "%%%s.%df", sign, digitsRightOfPeriod_f );
        sprintf( Str_f, Fmt, value );

        ne = (int) strlen( Str_e );
        nf = (int) strlen( Str_f );
        if( strchr( Str_f, '.' ) != NULL ) {        /* '.' in string. */
            if( flags & nf_floatToShortestString_trimZeros ) while( Str_f[nf-1] == '0' ) nf--;
            if( Str_f[nf-1] == '.' ) {
                if( !( flags & nf_floatToShortestString_keepPeriod ) ) nf--;
            } }
        else {      /* Maybe we want a '.' else it looks like an integer, "12345." vs "12345". */
            if( flags & nf_floatToShortestString_keepPeriod ) {
                Str_f[nf] = '.';
                nf++;
            }
        }
        Str_f[nf] = 0;

        if( ( nf + favorEFormBy ) < ne ) Str_r = Str_f;
    }
    return( strdup( Str_r ) );
}
