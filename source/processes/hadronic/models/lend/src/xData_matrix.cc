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
#include <limits.h>
#include <ctype.h>
#include "xData.h"

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

//char const * const xData_matrix_ID = "matrix";

static int toData( statusMessageReporting *smr, xDataType *xDT, xData_attributionList *attributes, const char *text );
static char *toString( statusMessageReporting *smr, xDataType *xDT );
static int release( statusMessageReporting *smr, xDataType *xDT );
/*
************************************************************
*/
int xData_init_matrix( statusMessageReporting *smr, xData_element *element ) {

    xDataType *xDT = &(element->xDataTypeInfo);

    xDT->status = xData_xDataType_Ok;
    xDT->typeString = xData_matrix_ID;
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
int xData_is_matrix( statusMessageReporting *smr, xDataType *xDT, int setMsg ) {

    return( xData_is_xDataType( smr, xDT, xData_matrix_ID, setMsg ) );
}
/*
************************************************************
*/
int xData_isElement_matrix( statusMessageReporting *smr, xData_element *element, int setMsg ) {

    return( xData_is_matrix( smr, &(element->xDataTypeInfo), setMsg ) );
}
/*
************************************************************
*/
xData_matrix *xData_matrix_copyData( statusMessageReporting *smr, xData_element *element ) {

    xData_Int i, n;
    xDataType *xDT = &(element->xDataTypeInfo);
    xData_matrix *oldMatrix = (xData_matrix *) xDT->data, *newMatrix;
    double *oldP, *newP;

    if( !xData_isElement_matrix( smr, element, 1 ) ) return( NULL );
    n = oldMatrix->rows * oldMatrix->columns;
    if( ( newMatrix = (xData_matrix *) xData_malloc2( smr, sizeof( xData_matrix ) + xDT->length * sizeof( xData_matrix_rowStartEnd ) +
        n * sizeof( double ), 0, "data" ) ) == NULL ) return( NULL );
    newMatrix->rows = oldMatrix->rows;
    newMatrix->columns = oldMatrix->columns;
    newMatrix->rowStartEnds = (xData_matrix_rowStartEnd *) &(newMatrix[1]);
    newMatrix->values = (double *) &(newMatrix->rowStartEnds[xDT->length]);
    for( i = 0; i < xDT->length; i++ ) newMatrix->rowStartEnds[i] = oldMatrix->rowStartEnds[i];
    for( i = 0, oldP = oldMatrix->values, newP = newMatrix->values; i < n; i++, oldP++, newP++ ) *newP = *oldP;
    return( newMatrix );
}
/*
************************************************************
*/
int xData_matrix_free_copyData( statusMessageReporting *smr, void *data ) {

    xData_free( smr, data );
    return( 0 );
}
/*
************************************************************
*/
static int toData( statusMessageReporting *smr, xDataType *xDT, xData_attributionList *attributes, const char *text ) {

    xData_Int i, j, row, start, end, rows, columns, status = 0;
    char *e;
    const char *s, *size;
    double *p;
    xData_matrix *matrix;
    void *smrUser = xData_get_smrUserInterfaceFromElement( xDT->element );

    if( xDT->status != xData_xDataType_Ok ) return( xData_setMessageError_ReturnInt( 1, smr, smrUser, __FILE__, __LINE__, 1, "bad xDataType instance" ) );
    release( smr, xDT );

    if( ( size = xData_getAttributesValue( attributes, "size" ) ) == NULL ) return( xData_setMessageError_ReturnInt( 1, smr, 
        smrUser, __FILE__, __LINE__, 1, "xData missing \"size\" attribute" ) );
    if( xData_stringTo_xData_Int( smr, smrUser, size, &rows, " ,", &e ) ) return( 1 );
    while( isspace( *e ) ) e++;
    if( *e != ',' ) {
        smr_setMessageError( smr, smrUser, __FILE__, __LINE__, 1, "matrix size attribute missing \",\" separator" );
        return( 1 );
    }
    s = (const char *) ++e;
    if( xData_stringTo_xData_Int( smr, smrUser, s, &columns, "", &e ) ) return( 1 );
    if( ( xDT->data = (xData_matrix *) xData_malloc2( NULL, sizeof( xData_matrix ) + xDT->length * sizeof( xData_matrix_rowStartEnd ) + 
        rows * columns * sizeof( double ), 0, "xDT->data" ) ) == NULL ) return( 1 );
    matrix = (xData_matrix *) xDT->data;
    matrix->rows = rows;
    matrix->columns = columns;
    matrix->rowStartEnds = (xData_matrix_rowStartEnd *) &(matrix[1]);
    matrix->values = (double *) &(matrix->rowStartEnds[xDT->length]);

    for( i = 0, s = text; i < xDT->length; i++ ) {
        if( xData_stringTo_xData_Int( smr, smrUser, s, &row, " \n", &e ) ) { status = 1; break; }
        if( ( row < 0 ) || ( row >= rows ) ) {
            status = 1;
            smr_setMessageError( smr, smrUser, __FILE__, __LINE__, 1, "row = %lld out-of-range (valid range is 0 <= row <= %lld)", row, rows );
            break;
        }
        if( xData_stringTo_xData_Int( smr, smrUser, e, &start, " \n", &e ) ) { status = 1; break; }
        if( ( start < 0 ) || ( start > columns ) ) {
            status = 1;
            smr_setMessageError( smr, smrUser, __FILE__, __LINE__, 1, "start = %lld out-of-range (valid range is 0 <= start <= %lld)", start, columns );
            break;
        }
        if( xData_stringTo_xData_Int( smr, smrUser, e, &end, " \n", &e ) ) { status = 1; break; }
        if( ( end < start ) || ( end > columns ) ) {
            status = 1;
            smr_setMessageError( smr, smrUser, __FILE__, __LINE__, 1, "end = %lld out-of-range (valid range is %lld <= end <= %lld)", end, start, columns );
            break;
        }
        matrix->rowStartEnds[i].row = row;
        matrix->rowStartEnds[i].start = start;
        matrix->rowStartEnds[i].end = end;
        p = &(matrix->values[row * columns]);
        for( j = 0; j < start; j++ ) *(p++) = 0.;
        for( s = e; j < end; j++, p++, s = e ) {
            if( xData_stringTo_double( smr, smrUser, s, p, " \n", &e ) ) { status = 1; break; }
        }
        if( status != 0 ) break;
        for( ; j < columns; j++ ) *(p++) = 0.;
    }
    if( status == 0 ) {
        while( isspace( *e ) ) e++;
        if( *e != 0 ) {
            smr_setMessageError( smr, smrUser, __FILE__, __LINE__, 1, "matrix contains extra data = %s", e );
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

    xData_Int i, n = 1, start, end, iRow, iColumn;
    char *str, *p;
    xData_matrix *matrix = (xData_matrix *) xDT->data;
    double *row;

    for( iRow = 0; iRow < matrix->rows; iRow++ ) {
        row = &(matrix->values[iRow * matrix->columns]);
        for( iColumn = 0; iColumn < matrix->columns; iColumn++ ) if( row[iColumn] != 0. ) break;
        start = iColumn;
        for( end = iColumn; iColumn < matrix->columns; iColumn++ ) if( row[iColumn] != 0. ) end = iColumn + 1;
        if( start < end ) n += 10 * 3 + 17 * ( end - start + 1 );
    }
    if( ( str = (char *) xData_malloc2( NULL, n, 0, "str" ) ) == NULL ) return( NULL );
    for( iRow = 0, p = str; iRow < matrix->rows; iRow++ ) {
        row = &(matrix->values[iRow * matrix->columns]);
        for( iColumn = 0; iColumn < matrix->columns; iColumn++ ) if( row[iColumn] != 0. ) break;
        start = iColumn;
        for( end = iColumn; iColumn < matrix->columns; iColumn++ ) if( row[iColumn] != 0. ) end = iColumn + 1;
        if( start < end ) {
            sprintf( p, "%3d %3d %3d", iRow, start, end );
            p += strlen( p );
            for( i = start; i < end; i++, p += 17 ) sprintf( p, " %16.9e", row[i] );
            sprintf( p, "\n" );
            p++;
        }
        *p = 0;
    }
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
int getRowStartEndAtIndex( statusMessageReporting *smr, xDataType *xDT, xData_Int index, xData_Int *row, xData_Int *start, xData_Int *end ) {

    int status = 0;
    xData_matrix *matrix = (xData_matrix *) xDT->data;

    if( !xData_is_matrix( smr, xDT, 1 ) ) return( 1 );
    if( ( index < 0 ) || ( index >= xDT->length ) ) {
        smr_setMessageError( smr, xData_get_smrUserInterfaceFromElement( xDT->element ), __FILE__, __LINE__, 1, 
            "index = %lld out of range (valid range 0 <= index < %lld)", index, xDT->length );
        status = 1; }
    else {
        *row = matrix->rowStartEnds[index].row;
        *start = matrix->rowStartEnds[index].start;
        *end = matrix->rowStartEnds[index].end;
    }
    return( status );
}

#if defined __cplusplus
}
#endif
