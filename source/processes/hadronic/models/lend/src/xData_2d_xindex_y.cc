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

//char const * const xData_twod_xindex_y_ID = "2d.xindex_y";

static int toData( statusMessageReporting *smr, xDataType *xDT, xData_attributionList *attributes, const char *text );
static char *toString( statusMessageReporting *smr, xDataType *xDT );
static int release( statusMessageReporting *smr, xDataType *xDT );
static double *xData_2d_xindex_y_toFilled( statusMessageReporting *smr, xData_element *element, double *Xs, int size );
/*
************************************************************
*/
int xData_init_2d_xindex_y( statusMessageReporting *smr, xData_element *element ) {

    xDataType *xDT = &(element->xDataTypeInfo);

    xDT->status = xData_xDataType_Ok;
    xDT->typeString = xData_twod_xindex_y_ID;
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
int xData_is_2d_xindex_y( statusMessageReporting *smr, xDataType *xDT, int setMsg ) {

    return( xData_is_xDataType( smr, xDT, xData_twod_xindex_y_ID, setMsg ) );
}
/*
************************************************************
*/
int xData_isElement_2d_xindex_y( statusMessageReporting *smr, xData_element *element, int setMsg ) {

    return( xData_is_2d_xindex_y( smr, &(element->xDataTypeInfo), setMsg ) );
}
/*
************************************************************
*/
//static int toData( statusMessageReporting *smr, xDataType *xDT, xData_attributionList *attributes, const char *text ) {
static int toData( statusMessageReporting *smr, xDataType *xDT, xData_attributionList *, const char *text ) {

    xData_Int i, status = 0;
    char *e;
    const char *s;
    xData_2d_xindex_y *p;
    void *smrUser = xData_get_smrUserInterfaceFromElement( xDT->element );

    if( xDT->status != xData_xDataType_Ok ) return( xData_setMessageError_ReturnInt( 1, smr, smrUser, __FILE__, __LINE__, 1, "bad xDataType instance" ) );
    release( smr, xDT );
    if( ( xDT->data = xData_malloc2( smr, 2 * xDT->length * sizeof( xData_2d_xindex_y ), 0, "data" ) ) == NULL ) return( 1 );
    for( i = 0, s = text, p = (xData_2d_xindex_y *) xDT->data; i < xDT->length; i++, p++, s = e ) {
        if( xData_stringTo_xData_Int( smr, smrUser, s, &(p->index), " \n", &e ) ) { status = 1; break; }
        s = e;
        if( xData_stringTo_double( smr, smrUser, s, &(p->value), " \n", &e ) ) { status = 1; break; }
    }
    if( status == 0 ) {
        while( isspace( *e ) ) e++;
        if( *e != 0 ) {
            smr_setMessageError( smr, smrUser, __FILE__, __LINE__, 1, "2d_xindex_y contains extra data = %s", e );
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

    xData_Int i, n = xDT->length, recordSize = 5 + 16 + 1, indexMax = 9;
    int iFmt = 4;
    char *str, *p, fmt[32] = " %99d %15.7e\n";
    xData_2d_xindex_y *data = (xData_2d_xindex_y *) xDT->data;

    if( n < 0 ) n = 0;
    for( i = 0; i < n; i++, data++ ) {
        while( ( data->index > indexMax ) && ( indexMax > 0 ) ) {
            indexMax = 10 * indexMax + 9;
            recordSize++;
            iFmt++;
        }
    }
    sprintf( fmt, " %%%dld %%15.7e\n", iFmt );
    if( ( str = (char *) malloc( recordSize * ( n + 1 ) ) ) == NULL ) return( NULL );
    for( i = 0, p = str, data = (xData_2d_xindex_y *) xDT->data; i < n; i++, p += recordSize, data++ ) {
        sprintf( p, fmt, data->index, data->value );
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
/*
************************************************************
*/
xData_Int *xData_2d_xindex_y_rawIndices( statusMessageReporting *smr, xData_element *element ) {
/*
*   Returns NULL if length is 0 or memory could not be allocated.
*/
    xDataType *xDT = &(element->xDataTypeInfo);
    xData_Int i, index = xDT->start, *values;
    xData_2d_xindex_y *data = (xData_2d_xindex_y *) xDT->data;

    if( xDT->length == 0 ) return( NULL );
    //if( ( values = xData_malloc2( smr, xDT->length * sizeof( xData_Int ), 0, "values" ) ) == NULL ) return( NULL );
    if( ( values = (xData_Int*) xData_malloc2( smr, xDT->length * sizeof( xData_Int ), 0, "values" ) ) == NULL ) return( NULL );
    for( i = 0; i < xDT->length; i++ ) {
        values[i] = index;
        index += data[i].index;
    }
    return( values );
}
/*
************************************************************
*/
int xData_2d_xindex_y_free_rawIndices( statusMessageReporting *smr, void *data ) {

    xData_free( smr, data );
    return( 0 );
}
/*
************************************************************
*/
double *xData_2d_xindex_y_toXYs( statusMessageReporting *smr, xData_element *element, double *Xs ) {
/*
*   Returns NULL if length is 0 or memory could not be allocated.
*/
    xDataType *xDT = &(element->xDataTypeInfo);
    xData_Int i, index = xDT->start;
    double *values = NULL, *p;
    xData_2d_xindex_y *data = (xData_2d_xindex_y *) xDT->data;

    if( xDT->length == 0 ) return( NULL );
    //if( ( values = xData_malloc2( smr, 2 * xDT->length * sizeof( double ), 0, "values" ) ) == NULL ) return( NULL );
    if( ( values = (double*) xData_malloc2( smr, 2 * xDT->length * sizeof( double ), 0, "values" ) ) == NULL ) return( NULL );
    p = values;
    for( i = 0; i < xDT->length; i++, p++ ) {
        index += data[i].index;
        *p = Xs[index];
        p++;
        *p = data[i].value;
    }
    return( values );
}
/*
************************************************************
*/
double *xData_2d_xindex_y_toFilledYs( statusMessageReporting *smr, xData_element *element, double *Xs ) {

    return( xData_2d_xindex_y_toFilled( smr, element, Xs, 1 ) );
}
/*
************************************************************
*/
int xData_2d_xindex_y_free_toFilledYs( statusMessageReporting *smr, void *data ) {

    xData_free( smr, data );
    return( 0 );
}
/*
************************************************************
*/
double *xData_2d_xindex_y_toFilledXYs( statusMessageReporting *smr, xData_element *element, double *Xs ) {

    return( xData_2d_xindex_y_toFilled( smr, element, Xs, 2 ) );
}
/*
************************************************************
*/
static double *xData_2d_xindex_y_toFilled( statusMessageReporting *smr, xData_element *element, double *Xs, int size ) {
/*
*   Returns NULL if length is 0 or memory could not be allocated.
*/
    xDataType *xDT = &(element->xDataTypeInfo);
    xData_Int i, j, index = xDT->start, length = xDT->end - xDT->start;
    double x1, x2, *x, y1, y2, *values = NULL, *p;
    xData_2d_xindex_y *data = (xData_2d_xindex_y *) xDT->data;

    if( xDT->length == 0 ) return( NULL );
    //if( ( values = xData_malloc2( smr, size * length * sizeof( double ), 0, "values" ) ) == NULL ) return( NULL );
    if( ( values = (double*) xData_malloc2( smr, size * length * sizeof( double ), 0, "values" ) ) == NULL ) return( NULL );
    p = values;
    x = &(Xs[xDT->start]);
    x2 = 0.;                                        /* Dummy initializations, as x1 and y1 (set by x2 and y2) are not used first time thru loop. */
    y2 = 0.;
    for( i = 0; i < xDT->length; i++, x++, p++ ) {
        index += data[i].index;                     /* Note, data[0].index is 0; otherwise, following logic would not work. */
        x1 = x2;
        x2 = Xs[index];
        y1 = y2;
        y2 = data[i].value;
        for( j = data[i].index; j > 1; j--, x++, p++ ) {
            if( size == 2 ) {
                *p = *x;
                p++;
            }
            *p = ( y1 * ( x2 - *x ) + y2 * ( *x - x1 ) ) / ( x2 - x1 );
        }
        if( size == 2 ) *(p++) = *x;
        *p = y2;
    }
    return( values );
}

#if defined __cplusplus
}
#endif
