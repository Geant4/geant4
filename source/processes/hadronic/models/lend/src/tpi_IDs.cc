/*
# <<BEGIN-copyright>>
# Copyright (c) 2010, Lawrence Livermore National Security, LLC. 
# Produced at the Lawrence Livermore National Laboratory 
# Written by Bret R. Beck, beck6@llnl.gov. 
# CODE-461393
# All rights reserved. 
#  
# This file is part of GIDI. For details, see nuclear.llnl.gov. 
# Please also read the "Additional BSD Notice" at nuclear.llnl.gov. 
# 
# Redistribution and use in source and binary forms, with or without modification, 
# are permitted provided that the following conditions are met: 
#
#      1) Redistributions of source code must retain the above copyright notice, 
#         this list of conditions and the disclaimer below.
#      2) Redistributions in binary form must reproduce the above copyright notice, 
#         this list of conditions and the disclaimer (as noted below) in the 
#          documentation and/or other materials provided with the distribution.
#      3) Neither the name of the LLNS/LLNL nor the names of its contributors may be 
#         used to endorse or promote products derived from this software without 
#         specific prior written permission. 
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT 
# SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC, THE U.S. DEPARTMENT OF ENERGY OR 
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS 
# OR SERVICES;  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED 
# AND ON  ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, 
# EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
# <<END-copyright>>
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include <tpi_IDs.h>

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

static tpi_channelID *_tpi_channelID_parse2( statusMessageReporting *smr, char const *str, char const *origStr, int isDecayChannel, char **EOP );
static tpi_spectralID *_tpi_spectralID_parse2( statusMessageReporting *smr, char const *str, char const *origStr, char **EOP );
/*
***************************************************
*/
tpi_channelID *tpi_channelID_allocate( statusMessageReporting *smr ) {

    tpi_channelID *channelID;

    //if( ( channelID = xData_malloc2( smr, sizeof( tpi_channelID ), 0, "channelID" ) ) == NULL ) return( NULL );
    if( ( channelID = (tpi_channelID*) xData_malloc2( smr, sizeof( tpi_channelID ), 0, "channelID" ) ) == NULL ) return( NULL );
    tpi_channelID_initialize( smr, channelID );
    return( channelID );
}
/*
***************************************************
*/
//int tpi_channelID_initialize( statusMessageReporting *smr, tpi_channelID *channelID ) {
int tpi_channelID_initialize( statusMessageReporting *, tpi_channelID *channelID ) {

    memset( channelID, 0, sizeof( tpi_channelID ) );
    return( 0 );
}
/*
***************************************************
*/
void *tpi_channelID_free( statusMessageReporting *smr, tpi_channelID *channelID ) {

    if( channelID != NULL ) {
        tpi_channelID_release( smr, channelID );
        xData_free( smr, channelID );
    }
    return( NULL );
}
/*
***************************************************
*/
int tpi_channelID_release( statusMessageReporting *smr, tpi_channelID *channelID ) {

    tpi_spectralID *spectralID, *next;

    for( spectralID = channelID->spectralIDs; spectralID != NULL; spectralID = next ) {
        next = spectralID->next;
        tpi_spectralID_free( smr, spectralID );
    }
    tpi_channelID_initialize( smr, channelID );
    return( 0 );
}
/*
***************************************************
*/
tpi_channelID *tpi_channelID_parse( statusMessageReporting *smr, char const *str, char **EOP ) {

    return( _tpi_channelID_parse2( smr, str, str, 0, EOP ) );
}
/*
***************************************************
*/
static tpi_channelID *_tpi_channelID_parse2( statusMessageReporting *smr, char const *str, char const *origStr, int isDecayChannel, char **EOP ) {

    tpi_channelID *channelID;
    tpi_spectralID *spectralID, *priorSpectral;

    *EOP = (char *) tpi_misc_firstNonWhiteSpace( smr, (char *) str );
    if( **EOP == 0 ) {
        smr_setMessageError( smr, NULL, __FILE__, __LINE__, 1, "Empty channel string to parse for '%s'", origStr );
        return( NULL );
    }
    if( ( channelID = tpi_channelID_allocate( smr ) ) == NULL ) return( NULL );
    priorSpectral= (tpi_spectralID *) &(channelID->spectralIDs);
    while( 1 ) {
        //if( ( spectralID = _tpi_spectralID_parse2( smr, *EOP, origStr, EOP ) ) == NULL ) return( tpi_channelID_free( smr, channelID ) );
        if( ( spectralID = _tpi_spectralID_parse2( smr, *EOP, origStr, EOP ) ) == NULL ) return( (tpi_channelID*)tpi_channelID_free( smr, channelID ) );
        priorSpectral->next = spectralID;
        priorSpectral = spectralID;
        *EOP = (char *) tpi_misc_firstNonWhiteSpace( smr, *EOP );
        if( **EOP == 0 ) break;
        if( isDecayChannel && ( **EOP == ')' ) ) break;
        if( **EOP != '+' ) {
            smr_setMessageError( smr, NULL, __FILE__, __LINE__, 1, "Missing '+' (or maybe ')') in channelID at '%s' of '%s'", *EOP, origStr );
            //return( tpi_channelID_free( smr, channelID ) );
            return( (tpi_channelID*)tpi_channelID_free( smr, channelID ) );
        }
        (*EOP)++;
    }
    return( channelID );
} 
/*
***************************************************
*/
int tpi_channelID_toString( statusMessageReporting *smr, tpi_channelID *channelID, gString *gStr ) {

    return( tpi_channelID_toStringSans( smr, channelID, gStr, NULL ) );
}
/*
***************************************************
*/
int tpi_channelID_toStringSanRevision( statusMessageReporting *smr, tpi_channelID *channelID, gString *gStr ) {

    //char *sans[] = { "revision", NULL };
    char *sans[] = { (char*)"revision", NULL };

    return( tpi_channelID_toStringSans( smr, channelID, gStr, sans ) );
}
/*
***************************************************
*/
int tpi_channelID_toStringSans( statusMessageReporting *smr, tpi_channelID *channelID, gString *gStr, char *sans[] ) {

    tpi_spectralID *spectralID;

    for( spectralID = channelID->spectralIDs; spectralID != NULL; spectralID = spectralID->next ) {
        if( ( tpi_spectralID_toStringSans( smr, spectralID, gStr, sans ) ) != 0 ) return( 1 );
        if( spectralID->next != NULL ) {
            if( ( gString_addTo( smr, gStr, " + " ) ) != 0 ) return( 1 );
        }
    }
    return( 0 );
}
/*
***************************************************
*/
tpi_spectralID *tpi_spectralID_allocate( statusMessageReporting *smr ) {

    tpi_spectralID *spectralID;

    //if( ( spectralID = xData_malloc2( smr, sizeof( tpi_spectralID ), 1, "spectralID" ) ) == NULL ) return( NULL );
    if( ( spectralID = (tpi_spectralID*) xData_malloc2( smr, sizeof( tpi_spectralID ), 1, "spectralID" ) ) == NULL ) return( NULL );
    tpi_spectralID_initialize( smr, spectralID );
    return( spectralID );
}
/*
***************************************************
*/
//int tpi_spectralID_initialize( statusMessageReporting *smr, tpi_spectralID *spectralID ) {
int tpi_spectralID_initialize( statusMessageReporting *, tpi_spectralID *spectralID ) {

    memset( spectralID, 0, sizeof( tpi_spectralID ) );
    return( 0 );
}
/*
***************************************************
*/
void *tpi_spectralID_free( statusMessageReporting *smr, tpi_spectralID *spectralID ) {

    if( spectralID != NULL ) {
        tpi_spectralID_release( smr, spectralID );
        xData_free( smr, spectralID );
    }
    return( NULL );
}
/*
***************************************************
*/
int tpi_spectralID_release( statusMessageReporting *smr, tpi_spectralID *spectralID ) {

    tpi_spectralIDQualifier *qualifier, *next;

    if( spectralID->name != NULL ) free( spectralID->name );
    for( qualifier = spectralID->qualifiers; qualifier != NULL; qualifier = next ) {
        next = qualifier->next;
        xData_free( smr, qualifier );
    }
    if( spectralID->decayChannel != NULL ) tpi_channelID_free( smr, spectralID->decayChannel );
    tpi_spectralID_initialize( smr, spectralID );
    return( 0 );
}
/*
***************************************************
*/
tpi_spectralID *tpi_spectralID_parse( statusMessageReporting *smr, char const *str, char **EOP ) {

    return( _tpi_spectralID_parse2( smr, str, str, EOP ) );
}
/*
***************************************************
*/
static tpi_spectralID *_tpi_spectralID_parse2( statusMessageReporting *smr, char const *str, char const *origStr, char **EOP ) {

    int breakup = 0, i1, i2, m;
    double d;
    char const *s, *q;
    char *e;
    char c, bOrC;
    tpi_spectralID *spectralID;
    tpi_spectralIDQualifier *qualifier, *priorQualifier = NULL;

    *EOP = (char *) tpi_misc_firstNonWhiteSpace( smr, str );
    if( **EOP == 0 ) {
        smr_setMessageError( smr, NULL, __FILE__, __LINE__, 1, "Empty spectralID string to parse for '%s'", origStr );
        return( NULL );
    }
    if( **EOP == '(' ) {               /* Breakup spectralID, like '(Be_8 -> He_4[multipliticy:"2"])' */
        breakup = 1;
        (*EOP)++;
        *EOP = (char *) tpi_misc_firstNonWhiteSpace( smr, *EOP );
    }
    if( !isalpha( **EOP ) ) {
        smr_setMessageError( smr, NULL, __FILE__, __LINE__, 1, "Invalid spectralID name '%s' in '%s'", *EOP, origStr );
        return( NULL );
    }
    for( s = *EOP, i1 = 0; ( isalnum( **EOP ) ) || ( **EOP == '_' ); (*EOP)++ ) i1++;
    if( ( spectralID = tpi_spectralID_allocate( smr ) ) == NULL ) return( NULL );
    //if( ( spectralID->name = tpi_misc_allocateAndCopyNCharacters( smr, s, i1 ) ) == NULL ) return( tpi_spectralID_free( smr, spectralID ) );
    if( ( spectralID->name = tpi_misc_allocateAndCopyNCharacters( smr, s, i1 ) ) == NULL ) return( (tpi_spectralID*) tpi_spectralID_free( smr, spectralID ) );
    *EOP = (char *) tpi_misc_firstNonWhiteSpace( smr, *EOP );
    if( **EOP == '[' ) {               /* Spectral has qualifiers. */
        priorQualifier = (tpi_spectralIDQualifier *) &(spectralID->qualifiers);
        bOrC = '[';
        while( **EOP == bOrC ) {
            bOrC = ',';
            (*EOP)++;                                                   /* Get qualifier's name. */
            *EOP = (char *) tpi_misc_firstNonWhiteSpace( smr, *EOP );
            if( !isalpha( **EOP ) ) {
                smr_setMessageError( smr, NULL, __FILE__, __LINE__, 1, "Invalid qualifier name '%s' in '%s'", *EOP, origStr );
                //return( tpi_spectralID_free( smr, spectralID ) );
                return( (tpi_spectralID*) tpi_spectralID_free( smr, spectralID ) );
            }
            for( s = *EOP, i1 = 0; ( isalnum( **EOP ) ) || ( **EOP == '_' ); (*EOP)++ ) i1++;

            *EOP = (char *) tpi_misc_firstNonWhiteSpace( smr, *EOP );                  /* Skip qualifier's separator. */
            if( **EOP != ':' ) {
                smr_setMessageError( smr, NULL, __FILE__, __LINE__, 1, "Missing ':' in qualifier defintion at '%s' in '%s'", *EOP, origStr );
                //return( tpi_spectralID_free( smr, spectralID ) );
                return( (tpi_spectralID*) tpi_spectralID_free( smr, spectralID ) );
            }

            (*EOP)++;                                                    /* Get qualifier's value. */
            *EOP = (char *) tpi_misc_firstNonWhiteSpace( smr, *EOP );
            if( ( **EOP != '"' ) && ( **EOP != '\'' ) ) {
                smr_setMessageError( smr, NULL, __FILE__, __LINE__, 1, "Missing start quote in qualifier defintion at '%s' in '%s'", *EOP, origStr );
                //return( tpi_spectralID_free( smr, spectralID ) );
                return( (tpi_spectralID*) tpi_spectralID_free( smr, spectralID ) );
            }
            c = **EOP;
            (*EOP)++;
            for( q = *EOP, i2 = 0; ( **EOP != c ) && ( **EOP != 0 ); (*EOP)++ ) i2++;
            if( **EOP != c ) {
                smr_setMessageError( smr, NULL, __FILE__, __LINE__, 1, "Missing end quote in qualifier defintion at '%s' in '%s'", *EOP, origStr );
                //return( tpi_spectralID_free( smr, spectralID ) );
                return( (tpi_spectralID*) tpi_spectralID_free( smr, spectralID ) );
            }
            (*EOP)++;
            *EOP = (char *) tpi_misc_firstNonWhiteSpace( smr, *EOP );

            //if( ( qualifier = xData_malloc2( smr, sizeof( tpi_spectralIDQualifier ) + i1 + i2 + 2, 1, "qualifier" ) ) == NULL ) 
            if( ( qualifier = (tpi_spectralIDQualifier*) xData_malloc2( smr, sizeof( tpi_spectralIDQualifier ) + i1 + i2 + 2, 1, "qualifier" ) ) == NULL ) 
                //return( tpi_spectralID_free( smr, spectralID ) );
                return( (tpi_spectralID*) tpi_spectralID_free( smr, spectralID ) );
            qualifier->next = NULL;
            qualifier->name = (char *) &(qualifier[1]);
            qualifier->value = &(qualifier->name[i1+1]);
            strncpy( qualifier->name, s, i1 );
            qualifier->name[i1] = 0;
            strncpy( qualifier->value, q, i2 );
            qualifier->value[i2] = 0;

            if( strcmp( qualifier->name, "revision" ) == 0 ) {
                spectralID->revision = qualifier->value; }
            else if( strcmp( qualifier->name, "multiplicity" ) == 0 ) {
                spectralID->multiplicityStr = qualifier->value; }
            else if( strcmp( qualifier->name, "level" ) == 0 ) {
                spectralID->levelStr = qualifier->value;
            }
            priorQualifier->next = qualifier;
            priorQualifier = qualifier;
        }
        if( **EOP != ']' ) {
            smr_setMessageError( smr, NULL, __FILE__, __LINE__, 1, "Missing ']' for qualifier at '%s' in '%s'", *EOP, origStr );
            //return( tpi_spectralID_free( smr, spectralID ) );
            return( (tpi_spectralID*) tpi_spectralID_free( smr, spectralID ) );
        }
        (*EOP)++;
        if( spectralID->multiplicityStr != NULL ) {
            m = strtol( spectralID->multiplicityStr, &e, 10 );
            if( ( *e == 0 ) && ( e != spectralID->multiplicityStr ) ) spectralID->multiplicity = m;
        }
        if( spectralID->levelStr != NULL ) {
            d = strtod( spectralID->levelStr, &e );
            if( ( *e == 0 ) && ( e != spectralID->levelStr ) ) spectralID->level = d;
        }
    }
    if( breakup ) {
        *EOP = (char *) tpi_misc_firstNonWhiteSpace( smr, *EOP );
        if( ( **EOP != '-' ) || ( (*EOP)[1] != '>' ) ) {
            smr_setMessageError( smr, NULL, __FILE__, __LINE__, 1, "Missing '->' for breakup at '%s' in '%s'", *EOP, origStr );
            //return( tpi_spectralID_free( smr, spectralID ) );
            return( (tpi_spectralID*) tpi_spectralID_free( smr, spectralID ) );
        }
        *EOP += 2;
        if( (spectralID->decayChannel = _tpi_channelID_parse2( smr, *EOP, origStr, breakup, EOP )) == NULL ) 
            //return( tpi_spectralID_free(smr, spectralID) );
            return( (tpi_spectralID*) tpi_spectralID_free(smr, spectralID) );
        if( **EOP != ')' ) {
            smr_setMessageError( smr, NULL, __FILE__, __LINE__, 1, "Missing ')' for breakup at '%s' in '%s'", *EOP, origStr );
            //return( tpi_spectralID_free( smr, spectralID ) );
            return( (tpi_spectralID*) tpi_spectralID_free( smr, spectralID ) );
        }
        (*EOP)++;
    }
    return( spectralID );
}
/*
***************************************************
*/
int tpi_spectralID_toString( statusMessageReporting *smr, tpi_spectralID *spectralID, gString *gStr ) {

    return( tpi_spectralID_toStringSans( smr, spectralID, gStr, NULL ) );
}
/*
***************************************************
*/
int tpi_spectralID_toStringSanRevision( statusMessageReporting *smr, tpi_spectralID *spectralID, gString *gStr ) {

    //char *sans[] = { "revision", NULL };
    char *sans[] = { (char*)"revision", NULL };

    return( tpi_spectralID_toStringSans( smr, spectralID, gStr, sans ) );
}
/*
***************************************************
*/

int tpi_spectralID_toStringSans( statusMessageReporting *smr, tpi_spectralID *spectralID, gString *gStr, char *sans[] ) {

    tpi_spectralIDQualifier *qualifier;
    int i;
    char **san, *sSan[] = { NULL };

    if( sans == NULL ) sans = sSan;
    if( spectralID->decayChannel != NULL ) if( gString_addTo( smr, gStr, "(" ) != 0 ) return( 1 );
    if( ( gString_addTo( smr, gStr, spectralID->name ) ) != 0 ) return( 1 );
    if( spectralID->qualifiers != NULL ) {
        for( qualifier = spectralID->qualifiers, i = 0; qualifier != NULL; qualifier = qualifier->next ) i++;
        for( qualifier = spectralID->qualifiers; qualifier != NULL; qualifier = qualifier->next ) {
            for( san = (char **) sans; *san != NULL; san++ ) {
                if( strcmp( *san, qualifier->name ) == 0 ) {
                    i--;
                    break;
                }
            }
        }
        if( i > 0 ) {
            if( gString_addTo( smr, gStr, "[" ) != 0 ) return( 1 );
            for( qualifier = spectralID->qualifiers; qualifier != NULL; qualifier = qualifier->next ) {
                for( san = (char **) sans; *san != NULL; san++ ) if( strcmp( *san, qualifier->name ) == 0 ) break;
                if( *san != NULL ) continue;
                if( gString_addTo( smr, gStr, qualifier->name ) != 0 ) return( 1 );
                if( gString_addTo( smr, gStr, ":'" ) != 0 ) return( 1 );
                if( gString_addTo( smr, gStr, qualifier->value ) != 0 ) return( 1 );
                if( gString_addTo( smr, gStr, "'") != 0 ) return( 1 );
                if( i == 1 ) {
                    if( gString_addTo( smr, gStr, "]" ) != 0 ) return( 1 ); }
                else {
                    if( gString_addTo( smr, gStr, ", " ) != 0 ) return( 1 );
                }
                i--;
            }
        }
    }
    if( spectralID->decayChannel != NULL ) {
        if( ( gString_addTo( smr, gStr, " -> " ) ) != 0 ) return( 1 );
        if( ( tpi_channelID_toStringSans( smr, spectralID->decayChannel, gStr, sans ) ) != 0 ) return( 1 );
        if( ( gString_addTo( smr, gStr, ")" ) ) != 0 ) return( 1 );
    }
    return( 0 );
}
/*
***************************************************
*/
//char const *tpi_misc_firstNonWhiteSpace( statusMessageReporting *smr, char const *str ) {
char const *tpi_misc_firstNonWhiteSpace( statusMessageReporting *, char const *str ) {

    char const *s;

    for( s = str; ( *s != 0 ) && isspace( *s ); s++ ) ;
    return( s );
}
/*
***************************************************
*/
char *tpi_misc_allocateAndCopyNCharacters( statusMessageReporting *smr, char const *str, int n ) {

    char *s;

    //if( ( s = xData_malloc2( smr, n + 1, 0, "s" ) ) == NULL ) return( NULL );
    if( ( s = (char*) xData_malloc2( smr, n + 1, 0, "s" ) ) == NULL ) return( NULL );
    strncpy( s, str, n );
    s[n] = 0;
    return( s );
}

#if defined __cplusplus
}
#endif
