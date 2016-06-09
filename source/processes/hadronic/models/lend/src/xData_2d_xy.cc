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
#include <limits.h>
#include <ctype.h>
#include "xData.h"

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif 

//char const * const xData_twod_xy_ID = "2d.xy";

static int toData( statusMessageReporting *smr, xDataType *xDT, xData_attributionList *attributes, const char *text );
static char *toString( statusMessageReporting *smr, xDataType *xDT );
static int release( statusMessageReporting *smr, xDataType *xDT );
/*
************************************************************
*/
int xData_init_2d_xy( statusMessageReporting *smr, xData_element *element ) {

    xDataType *xDT = &(element->xDataTypeInfo);

    xDT->status = xData_xDataType_Ok;
    xDT->typeString = xData_twod_xy_ID;
    xDT->element = element;
    xDT->toData = toData;
    xDT->toString = toString;
    xDT->release = release;
    xDT->data = NULL;
    return( xData_xDataTypeConvertAttributes( smr, element ) );
}
/*
************************************************************
*/
int xData_is_2d_xy( statusMessageReporting *smr, xDataType *xDT, int setMsg ) {

    return( xData_is_xDataType( smr, xDT, xData_twod_xy_ID, setMsg ) );
}
/*
************************************************************
*/
int xData_isElement_2d_xy( statusMessageReporting *smr, xData_element *element, int setMsg ) {

    return( xData_is_2d_xy( smr, &(element->xDataTypeInfo), setMsg ) );
}
/*
************************************************************
*/
double *xData_2d_xy_allocateCopyData( statusMessageReporting *smr, xData_element *element, xData_Int *length ) {

    xData_Int i;
    xDataType *xDT = &(element->xDataTypeInfo);
    double *data;

    if( !xData_isElement_2d_xy( smr, element, 1 ) ) return( NULL );
    *length = xDT->length;
    //if( ( data = xData_malloc2( smr, 2 * xDT->length * sizeof( double ), 0, "data" ) ) ) {
    if( ( data = (double*)xData_malloc2( smr, 2 * xDT->length * sizeof( double ), 0, "data" ) ) ) {
        for( i = 0; i < 2 * xDT->length; i++ ) data[i] = ((double *) xDT->data)[i];
    }
    return( data );
}
/*
************************************************************
*/
int xData_2d_xy_free_copyData( statusMessageReporting *smr, void *data ) {

    xData_free( smr, data );
    return( 0 );
}
/*
************************************************************
*/
//static int toData( statusMessageReporting *smr, xDataType *xDT, xData_attributionList *attributes, const char *text ) {
static int toData( statusMessageReporting *smr, xDataType *xDT, xData_attributionList *, const char *text ) {

    xData_Int i, status = 0;
    char *e;
    const char *s;
    double *p;
    void *smrUser = xData_get_smrUserInterfaceFromElement( xDT->element );

    if( xDT->status != xData_xDataType_Ok ) return( xData_setMessageError_ReturnInt( 1, smr, smrUser, __FILE__, __LINE__,
        1, "bad xDataType instance" ) );
    release( smr, xDT );
    if( ( xDT->data = xData_malloc2( smr, 2 * xDT->length * sizeof( double ), 0, "data" ) ) == NULL ) return( 1 );
    for( i = 0, s = text, p = (double *) xDT->data; i < 2 * xDT->length; i++, p++, s = e ) {
        if( xData_stringTo_double( smr, smrUser, s, p, " \n", &e ) ) { status = 1; break; }
    }
    if( status == 0 ) {
        while( isspace( *e ) ) e++;
        if( *e != 0 ) {
            smr_setMessageError( smr, smrUser, __FILE__, __LINE__, 1, "text contains extra data = %s", e );
            status = 1;
        }
    }
    if( status != 0 ) release( smr, xDT );
    return( status );
}
/*
************************************************************
*/
//static char *toString( statusMessageReporting *smr, xDataType *xDT ) {
static char *toString( statusMessageReporting *, xDataType *xDT ) {

    xData_Int i, n = xDT->length;
    char *str, *p;
    double *data = (double *) xDT->data;

    if( n < 0 ) n = 0;
    if( ( str = (char *) malloc( 2 * ( n + 1 ) * 17 ) ) == NULL ) return( NULL );
    for( i = 0, p = str; i < n; i++, p += 2 * 16 + 1, data += 2 ) {
        sprintf( p, " %15.7e %15.7e\n", *data, data[1] );
    }
    *p = 0;
    return( str );
}
/*
************************************************************
*/
static int release( statusMessageReporting *smr, xDataType *xDT ) {

    if( xDT->data != NULL ) xDT->data = xData_free( smr, xDT->data );
    return( xDT->status = xData_xDataType_Ok );
}

#if defined __cplusplus
}
#endif
