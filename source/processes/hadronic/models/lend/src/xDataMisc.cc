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
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#if !defined(WIN32) || defined(__GNUC__)
   #include <unistd.h>
#else
   #include <direct.h>
   #define getcwd _getcwd
#endif
#include "xData.h"

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif


/*
************************************************************
*/
void *xData_malloc( statusMessageReporting *smr, size_t size, int zero, const char *forItem, const char *file, int line ) {

    void *p = xData_realloc( smr, NULL, size, forItem, file, line );
    int i;
    char *c;
    long long *l;

    if( ( p != NULL ) && zero ) {
        //for( i = 0, l = (long long *) p; i < size / sizeof( long long ); i++, l++ ) *l = 0;
        for( i = 0, l = (long long *) p; i < (int)( size / sizeof( long long ) ); i++, l++ ) *l = 0;
        //for( i = sizeof( long long ) * i, c = (char *) l; i < size; i++, c++ ) *c = 0;
        for( i = sizeof( long long ) * i, c = (char *) l; i < (int) size; i++, c++ ) *c = 0;
    }

    return( p );
}
/*
************************************************************
*/
void *xData_realloc( statusMessageReporting *smr, void *pOld, size_t size, const char *forItem, const char *file, int line ) {

    void *p = realloc( pOld, size );

    if( ( p == NULL ) && ( smr != NULL ) ) {
        smr_setMessageError( smr, NULL, file, line, -1, " xData_realloc: failed to realloc size = %llu for variable %s\n", 
            (unsigned long long) size, forItem );
    }
    return( p );
}
/*
************************************************************
*/
//void *xData_free( statusMessageReporting *smr, void *p ) {
void *xData_free( statusMessageReporting *, void *p ) {

    if( p != NULL ) free( p );
    return( NULL );
}
/*
************************************************************
*/
char *xDataMisc_allocateCopyString( statusMessageReporting *smr, const char *s, const char *forItem, const char *file, int line ) {
/*
*   User must free returned string.
*/
    char *c;

    //if( ( c = xData_malloc( smr, strlen( s ) + 1, 0, forItem, file, line ) ) != NULL ) {
    if( ( c = (char*) xData_malloc( smr, strlen( s ) + 1, 0, forItem, file, line ) ) != NULL ) {
        strcpy( c, s );
    }
    return( c );
}
/*
************************************************************
*/
char *xDataMisc_getAbsPath( statusMessageReporting *smr, const char *fileName ) {
/*
*   User must free returned string.
*/
    int n = strlen( fileName ) + 1, nCwd = 0;
    //char *absPath, cwd[4 * 1024] = "", *p, *needle;
    char *absPath, cwd[4 * 1024 + 1] = "", *p, *needle;

    if( fileName[0] != '/' ) {
        //if( getcwd( cwd, sizeof( cwd ) + 1 ) == NULL ) {
        if( getcwd( cwd, sizeof( cwd ) ) == NULL ) {
            smr_setMessageError( smr, NULL, __FILE__, __LINE__, -1, "hardwired cwd too small" );
            return( NULL );
        }
        nCwd = strlen( cwd );
        n += nCwd + 1;                                  /* cwd + '/'. */
    }
    //if( ( absPath = xData_malloc2( smr, n, 0, "absPath" ) ) == NULL ) return( NULL );
    if( ( absPath = (char*) xData_malloc2( smr, n, 0, "absPath" ) ) == NULL ) return( NULL );
    if( fileName[0] != '/' ) {
        strcpy( absPath, cwd );
        strcat( absPath, "/" );
        strcat( absPath, fileName ); }
    else {
        strcpy( absPath, fileName );
    }

    while( 1 ) {                                        /* Remove all ./ from path. */
        if( ( needle = strstr( absPath, "/./" ) ) == NULL ) break;
        p = needle;
        for( needle += 2; *needle; p++, needle++ ) *p = *needle;
        *p = 0;
    }

    while( 1 ) {                                        /* Remove all ../ from path. */
        if( ( needle = strstr( absPath, "/../" ) ) == NULL ) break;
        p = needle - 1;
        while( ( p > absPath ) && ( *p != '/' ) ) p--;
        if( *p != '/' ) break;                           /* This should not happen if path is legit, I think, and I do not know what to do so will leave it. */
        if( p == absPath ) break;                       /* Ditto. */
        for( needle += 3; *needle; p++, needle++ ) *p = *needle;
        *p = 0;
    }
    return( absPath );
}
/*
************************************************************
*/
int xData_setMessageError_ReturnInt( int value, statusMessageReporting *smr, void *userInterface, const char *packageName, int lineNumber, int code, 
    const char *fmt, ... ) {

    va_list args;

    va_start( args, fmt );
    smr_setMessageError( smr, userInterface, packageName, lineNumber, code, fmt, args );
    va_end( args );
    return( value );
}

#if defined __cplusplus
}
#endif
