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
#include "gString.h"
#include <xData.h>
#include <string.h>

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

/*
***************************************************
*/
int gString_initialize( statusMessageReporting *smr, gString *gStr, int size, int increment ) {

    if( size > 0 ) {
        //if( ( gStr->gStr = xData_malloc2( smr, size + 1, 0, "gStr->gStr" ) ) == NULL ) return( 1 );
        if( ( gStr->gStr = (char*) xData_malloc2( smr, size + 1, 0, "gStr->gStr" ) ) == NULL ) return( 1 );
        gStr->length = 1;
        gStr->gStr[0] = 0; }
    else {
        size = 0;
        gStr->length = 0;
        gStr->gStr = NULL;
    }
    gStr->allocated = size;
    //if( increment < gString_minIncrement ) increment = gString_minIncrement;
    if( increment < (int) gString_minIncrement ) increment = gString_minIncrement;
    gStr->increment = increment;
    return( 0 );
}
/*
***************************************************
*/
int gString_release( statusMessageReporting *smr, gString *gStr ) {

    if( gStr->gStr != NULL ) free( gStr->gStr );
    gString_initialize( smr, gStr, 0, gStr->increment );
    return( 0 );
}
/*
***************************************************
*/
//int gString_clear( statusMessageReporting *smr, gString *gStr ) {
int gString_clear( statusMessageReporting *, gString *gStr ) {

    if( gStr->gStr != NULL ) {
        gStr->length = 1;
        gStr->gStr[0] = 0;
    }
    return( 0 );
}
/*
***************************************************
*/
int gString_addTo( statusMessageReporting *smr, gString *gStr, char const *str ) {

    int n, size = strlen( str );

    if( gStr->gStr == NULL ) {
        if( gString_initialize( smr, gStr, size + 1, gStr->increment ) != 0 ) return( 1 ); }
    else if( ( gStr->length + size ) > gStr->allocated ) {
        n = gStr->increment;
        if( n < size ) n = size;
        //if( ( gStr->gStr = xData_realloc2( smr, gStr->gStr, gStr->allocated + n, "gStr->gStr" ) ) == NULL ) return( 1 );
        if( ( gStr->gStr = (char*) xData_realloc2( smr, gStr->gStr, gStr->allocated + n, "gStr->gStr" ) ) == NULL ) return( 1 );
        gStr->allocated += n;
    }
    strcpy( &(gStr->gStr[gStr->length - 1]), str );
    gStr->length = gStr->length + size;
    return( 0 );
}
/*
***************************************************
*/
//char const *gString_string( statusMessageReporting *smr, gString *gStr ) {
char const *gString_string( statusMessageReporting *, gString *gStr ) {

    return( gStr->gStr );
}
/*
***************************************************
*/
//int gString_length( statusMessageReporting *smr, gString *gStr ) {
int gString_length( statusMessageReporting *, gString *gStr ) {

    return( gStr->length );
}
/*
***************************************************
*/
//int gString_allocated( statusMessageReporting *smr, gString *gStr ) {
int gString_allocated( statusMessageReporting *, gString *gStr ) {

    return( gStr->allocated );
}
/*
***************************************************
*/
//int gString_increment( statusMessageReporting *smr, gString *gStr ) {
int gString_increment( statusMessageReporting *, gString *gStr ) {

    return( gStr->increment );
}

#if defined __cplusplus
}
#endif
