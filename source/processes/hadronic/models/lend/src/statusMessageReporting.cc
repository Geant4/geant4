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
#include <stdarg.h>

#include "statusMessageReporting.h"

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

static const char smr_mallocFailed[] = "statusMessageReporting could not allocate memory for message";

static int smr_setAllocationFailure( statusMessageReporting *smr, const char *fmt, va_list *args );
static int smr_setMessage( statusMessageReporting *smr, void *userInterface, const char *file, int line, int code, 
    enum smr_status status, const char *fmt, va_list *args );
static char *smr_getFullMessage2( char const *fmt, ... );
/*
************************************************************
*/
int smr_initialize( statusMessageReporting *smr ) {

    smr->status = smr_status_Ok;
    smr->packageName[0] = 0;
    smr->line= -1;
    smr->code = 0;
    smr->message = NULL;
    return( 0 );
}
/*
************************************************************
*/
int smr_release( statusMessageReporting *smr ) {

    if( smr->message != NULL ) {
        if( smr->message != smr_mallocFailed ) free( smr->message );
    }
    return( smr_initialize( smr ) );
}
/*
************************************************************
*/
int smr_setMessageInfo( statusMessageReporting *smr, void *userInterface, const char *file, int line, int code, const char *fmt, ... ) {

    int status;
    va_list args;

    va_start( args, fmt );
    status = smr_setMessage( smr, userInterface, file, line, code, smr_status_Info, fmt, &args );
    va_end( args );
    return( status );
}
/*
************************************************************
*/
int smr_vsetMessageInfo( statusMessageReporting *smr, void *userInterface, const char *file, int line, int code, const char *fmt, va_list *args ) {

    int status = smr_setMessage( smr, userInterface, file, line, code, smr_status_Info, fmt, args );
    return( status );
}
/*
************************************************************
*/
int smr_setMessageError( statusMessageReporting *smr, void *userInterface, const char *file, int line, int code, const char *fmt, ... ) {

    int status;
    va_list args;

    va_start( args, fmt );
    status = smr_setMessage( smr, userInterface, file, line, code, smr_status_Error, fmt, &args );
    va_end( args );
    return( status );
}
/*
************************************************************
*/
int smr_vsetMessageError( statusMessageReporting *smr, void *userInterface, const char *file, int line, int code, const char *fmt, va_list *args ) {

    int status = smr_setMessage( smr, userInterface, file, line, code, smr_status_Error, fmt, args );
    return( status );
}
/*
************************************************************
*/
char *smr_allocateFormatMessage( const char *fmt, ... ) {

    char *s;
    va_list args;

    va_start( args, fmt );
    s = smr_vallocateFormatMessage( fmt, &args );
    va_end( args );
    return( s );
}
/*
************************************************************
*/
char *smr_vallocateFormatMessage( const char *fmt, va_list *args ) {

    int n, size = 128;
    char *message = NULL;
    va_list args_;

    while( 1 ) {
        if( ( message = (char *) realloc( message, size ) ) == NULL ) return( NULL );
        //TK110426
#if defined WIN32
        args_ = *args;
#elif defined __IBMCPP__
        va_copy( args_, *args );
#else
        __va_copy( args_, *args );
#endif
        n = vsnprintf( message, size, fmt, args_ );
        va_end( args_ );
        if( ( n > -1 ) && ( n < size ) ) break;
        if( n > -1 ) {      /* glibc 2.1 */
            size = n + 3; }
        else {              /* glibc 2.0 */
            size += 128;
        }
    }
    return( message );
}
/*
************************************************************
*/
static int smr_setMessage( statusMessageReporting *smr, void *userInterface, const char *file, int line, int code, 
    enum smr_status status, const char *fmt, va_list *args ) {

    char *userMsg = NULL;
    int userSize;

    if( smr == NULL ) return( 0 );
    smr_release( smr );
    smr->status = status;
    if( file != NULL ) strncpy( smr->packageName, file, smr_maximumPackageNameSize );
    smr->packageName[smr_maximumPackageNameSize-1] = 0;
    smr->line= line;
    smr->code = code;

    if( ( smr->message = smr_vallocateFormatMessage( fmt, args ) ) == NULL ) return( smr_setAllocationFailure( smr, fmt, args ) );
    if( userInterface != NULL ) {
        if( ( userSize = (*(smr_userInterface *) userInterface)( (void *) userInterface, NULL ) ) > 0 ) {
            //if( (smr->message = realloc(smr->message, strlen( smr->message ) + userSize + 2)) == NULL ) return( smr_setAllocationFailure( smr, fmt, args ) );
            if( (smr->message = (char*) realloc(smr->message, strlen( smr->message ) + userSize + 2)) == NULL ) return( smr_setAllocationFailure( smr, fmt, args ) );
            strcat( smr->message, "\n" );
            userSize = (*(smr_userInterface *) userInterface)( (void *) userInterface, &userMsg );
            if( userSize < 0 ) return( smr_setAllocationFailure( smr, fmt, args ) );
            if( userSize > 0 ) {
                strcat( smr->message, userMsg );
                free( userMsg );
            }
        }
    }
    return( 0 );
}
/*
************************************************************
*/
static int smr_setAllocationFailure( statusMessageReporting *smr, const char *fmt, va_list *args ) {

    vfprintf( stderr, fmt, *args );      /* Assume calling routine calls va_end( args ). */
    fprintf( stderr, "\nAt line %d of %s\n", smr->line, smr->packageName );
    smr->status = smr_status_Fatal;
    smr->message = (char *) smr_mallocFailed;
    return( 1 );
}
/*
************************************************************
*/
int smr_isOk( statusMessageReporting *smr ) { 

    if( smr == NULL ) return( 1 );
    return( smr->status == smr_status_Ok );
}
/*
************************************************************
*/
int smr_isInfo( statusMessageReporting *smr ) { 

    if( smr == NULL ) return( 1 );
    return( smr->status == smr_status_Info );
}
/*
************************************************************
*/
int smr_isError( statusMessageReporting *smr ) { 

    if( smr == NULL ) return( 1 );
    return( ( smr->status == smr_status_Error ) || ( smr->status == smr_status_Fatal ) );
}
/*
************************************************************
*/
int smr_isFatal( statusMessageReporting *smr ) { 

    if( smr == NULL ) return( 1 );
    return( smr->status == smr_status_Fatal );
}
/*
************************************************************
*/
const char *smr_getMessage( statusMessageReporting *smr ) {

    return( smr->message );
}
/*
************************************************************
*/
char *smr_getFullMessage( statusMessageReporting *smr ) {

    return( smr_getFullMessage2( "%s\nAt line %d of %s", smr->message, smr->line, smr->packageName ) );
}
/*
************************************************************
*/
static char *smr_getFullMessage2( char const *fmt, ... ) {

    va_list args;
    char *message;

    va_start( args, fmt );
    message = smr_vallocateFormatMessage( fmt, &args );
    va_end( args );
    return( message );
}
/*
************************************************************
*/
void smr_print( statusMessageReporting *smr, FILE *f, int clear ) {

    if( smr->message != NULL ) fprintf( f, "%s\nAt line %d of %s\n", smr->message, smr->line, smr->packageName );
    if( clear ) smr_release( smr );
}

#if defined __cplusplus
}
#endif
