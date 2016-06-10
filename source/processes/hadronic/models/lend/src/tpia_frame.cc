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
#include <tpia_target.h>
#include <tpia_misc.h>
#include "G4Types.hh"

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

static G4ThreadLocal int nlab_Str = 0, nCOM_Str = 0;
static const char lab_Str[] = "lab", COM_Str[] = "centerOfMass";
/*
************************************************************
*/
//int tpia_frame_clear( statusMessageReporting *smr, tpia_data_frame *frame ) {
int tpia_frame_clear( statusMessageReporting *, tpia_data_frame *frame ) {

    frame->frames = 0;
    return( 0 );
}
/*
************************************************************
*/
int tpia_frame_setFromElement( statusMessageReporting *smr, xData_element *element, int dimension, tpia_data_frame *frame ) {

    char const *value = xData_getAttributesValueInElement( element, "frame" );

    if( value == NULL ) {
        tpia_misc_setMessageError_Element( smr, NULL, element, __FILE__, __LINE__, 1, "element is missing frame attribute" );
        return( 1 );
    }
    return( tpia_frame_setFromString( smr, element->fullName, value, dimension, frame ) );
}
/*
************************************************************
*/
int tpia_frame_setFromString( statusMessageReporting *smr, const char *forItem, const char *value, int dimension, tpia_data_frame *frame ) {

    char const *e;
    int status = 1, i;

    for( i = 0; i < tpia_maxNumberOfFrames; i++ ) tpia_frame_setColumn( smr, frame, i, tpia_referenceFrame_None );
    if( dimension > tpia_maxNumberOfFrames ) {
        smr_setMessageError( smr, NULL, __FILE__, __LINE__, 1, "dimension argument = %d is greater than tpia_maxNumberOfFrames = %d", dimension,
            tpia_maxNumberOfFrames );
        return( status );
    }
    for( i = 0, e = value; ( i < dimension ) && ( *e != 0 ); i++ ) {
        if( strstr( e, lab_Str ) == e ) {
            tpia_frame_setColumn( smr, frame, i, tpia_referenceFrame_lab );
            e += 3; }
        else if( strstr( e, COM_Str ) == e ) {
            tpia_frame_setColumn( smr, frame, i, tpia_referenceFrame_COM );
            e += 12; }
        else {
            smr_setMessageError( smr, NULL, __FILE__, __LINE__, 1, "bad frame '%s' for %s", value, forItem );
            break;
        }
        if( *e != 0 ) {
            if( *e != ',' ) {
                smr_setMessageError( smr, NULL, __FILE__, __LINE__, 1, "bad separater for frame '%s' for %s", value, forItem );
                break;
            }
            e++;
        }
    }
    if( smr_isOk( smr ) ) {
        if( i == dimension ) {
            if( *e == 0 ) {
                status = 0; }
            else {
                smr_setMessageError( smr, NULL, __FILE__, __LINE__, 1, "extra values for frame '%s' for %s", value, forItem );
            } }
        else {
            smr_setMessageError( smr, NULL, __FILE__, __LINE__, 1, "missing values for frame '%s' for %s", value, forItem );
        }
    }
    return( status );
}
/*
************************************************************
*/
//int tpia_frame_getDimensions( statusMessageReporting *smr, tpia_data_frame *frame ) {
int tpia_frame_getDimensions( statusMessageReporting *, tpia_data_frame *frame ) {

    int i, dimension = 0;
    unsigned value = frame->frames;

    for( i = 0; i < tpia_maxNumberOfFrames; i++ ) {
        if( ( value & 3 ) == 0 ) break;
        dimension++;
        value = value >> 2;
    }
    return( dimension );
}
/*
************************************************************
*/
//char *tpia_frame_toString( statusMessageReporting *smr, const char *fromItem, tpia_data_frame *frame ) {
char *tpia_frame_toString( statusMessageReporting *smr, const char *, tpia_data_frame *frame ) {

    int i, n = tpia_frame_getDimensions( smr, frame ), value, nStr = 0;
    char *str = NULL, *p;

    if( nlab_Str == 0 ) {
        nlab_Str = strlen( lab_Str );
        nCOM_Str = strlen( COM_Str );
    }
    for( i = 0; i < n; i++ ) {
        value = tpia_frame_getColumn( smr, frame, i );
        if( value == tpia_referenceFrame_COM ) {
            nStr += nCOM_Str + 1; }
        else if( value == tpia_referenceFrame_lab ) {
            nStr += nlab_Str + 1; }
        else {
            smr_setMessageError( smr, NULL, __FILE__, __LINE__, 1, "bad frame value = %d for column = %d", value, i );
            return( NULL );
        }
    }
    if( nStr == 0 ) nStr = 1;
    //if( ( str = xData_malloc2( smr, nStr, 1, "str" ) ) == NULL ) return( NULL );
    if( ( str = (char*) xData_malloc2( smr, nStr, 1, "str" ) ) == NULL ) return( NULL );
    for( i = 0, p = str - 1; i < n; i++ ) {
        p++;
        value = tpia_frame_getColumn( smr, frame, i );
        if( value == tpia_referenceFrame_COM ) {
            strcpy( p, COM_Str );
            p += nCOM_Str; }
        else {
            strcpy( p, lab_Str );
            p += nlab_Str;
        }
        *p = ',';
    }
    *p = 0;
    return( str );
}
/*
************************************************************
*/
int tpia_frame_setColumns( statusMessageReporting *smr, tpia_data_frame *frame, int nColumns, int *values ) {

    int i;

    tpia_frame_clear( smr, frame );
    for( i = 0; i < nColumns; i++ ) if( tpia_frame_setColumn( smr, frame, i, values[i] ) != 0 ) return( 1 );

    return( 0 );
}
/*
************************************************************
*/
int tpia_frame_setColumn( statusMessageReporting *smr, tpia_data_frame *frame, int column, int value ) {

    int i;

    if( ( column < 0 ) || ( column >= tpia_maxNumberOfFrames ) ) {
        smr_setMessageError( smr, NULL, __FILE__, __LINE__, 1, "bad column = %d value for setting frame (0 <= column < %d)", column, tpia_maxNumberOfFrames );
        return( 1 ); }
    else {
        if( ( value < tpia_referenceFrame_None ) || ( value > tpia_referenceFrame_Max ) ) {
            smr_setMessageError( smr, NULL, __FILE__, __LINE__, 1, "bad frame value = %d for column = %d", value, column );
            return( 1 ); }
        else {
            i = 3 << ( 2 * column );
            frame->frames = frame->frames & (~i);
            value = value << ( 2 * column );
            frame->frames = frame->frames | value;
        }
    }
    return( 0 );
}
/*
************************************************************
*/
int tpia_frame_getColumn( statusMessageReporting *smr, tpia_data_frame *frame, int column ) {

    unsigned int value;

    if( ( column < 0 ) || ( column >= tpia_maxNumberOfFrames ) ) {
        smr_setMessageError( smr, NULL, __FILE__, __LINE__, 1, "bad column = %d value for getting frame (0 <= column < %d)", column, tpia_maxNumberOfFrames );
        return( -1 ); }
    else {
        value = frame->frames >> ( 2 * column );
    }
    return( (int) ( value & 3 ) );
}

#if defined __cplusplus
}
#endif
