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
#include <string.h>
#include <ctype.h>

#include <tpia_target.h>
#include <tpia_misc.h>

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

/*
************************************************************
*/
tpia_multiplicity *tpia_multiplicity_create( statusMessageReporting *smr ) {

    tpia_multiplicity *multiplicity;

    //if( ( multiplicity = xData_malloc2( smr, sizeof( tpia_multiplicity ), 0, "multiplicity" ) ) == NULL ) return( NULL );
    if( ( multiplicity = (tpia_multiplicity*) xData_malloc2( smr, sizeof( tpia_multiplicity ), 0, "multiplicity" ) ) == NULL ) return( NULL );
    if( tpia_multiplicity_initialize( smr, multiplicity ) ) multiplicity = tpia_multiplicity_free( smr, multiplicity );
    return( multiplicity );
}
/*
************************************************************
*/
int tpia_multiplicity_initialize( statusMessageReporting *smr, tpia_multiplicity *multiplicity ) {

    memset( multiplicity, 0, sizeof( tpia_multiplicity ) );
    if( tpia_frame_setFromString( smr, "", "", 0, &(multiplicity->frame) ) ) return( 1 );
    return( 0 );
}
/*
************************************************************
*/
tpia_multiplicity *tpia_multiplicity_free( statusMessageReporting *smr, tpia_multiplicity *multiplicity ) {

    tpia_multiplicity_release( smr, multiplicity );
    xData_free( smr, multiplicity );
    return( NULL );
}
/*
************************************************************
*/
int tpia_multiplicity_release( statusMessageReporting *smr, tpia_multiplicity *multiplicity ) {

    //multiplicity->grouped.data = xData_free( smr, multiplicity->grouped.data );
    multiplicity->grouped.data = (double*) xData_free( smr, multiplicity->grouped.data );
    //multiplicity->pointwise = xData_free( smr, multiplicity->pointwise );
    multiplicity->pointwise = (double*) xData_free( smr, multiplicity->pointwise );
    tpia_multiplicity_initialize( smr, multiplicity );
    return( 0 );
}
/*
************************************************************
*/
tpia_multiplicity *tpia_multiplicity_createGetFromElement( statusMessageReporting *smr, xData_element *multiplicityElement, int nGroups ) {

    tpia_multiplicity *multiplicity;

    if( ( multiplicity = tpia_multiplicity_create( smr ) ) == NULL ) return( NULL );
    if( tpia_multiplicity_getFromElement( smr, multiplicityElement, multiplicity, nGroups ) != 0 ) multiplicity = tpia_multiplicity_free( smr, multiplicity  );
    return( multiplicity );
}
/*
************************************************************
*/
//int tpia_multiplicity_getFromElement( statusMessageReporting *smr, xData_element *multiplicityElement, tpia_multiplicity *multiplicity, int nGroups ) {
int tpia_multiplicity_getFromElement( statusMessageReporting *smr, xData_element *multiplicityElement, tpia_multiplicity *multiplicity, int ) {

    const char *timeScale;
    int isDelayedNeutrons;
    xData_element *data;

    xData_addToAccessed( smr, multiplicityElement, 1 );
    if( ( tpia_frame_setFromElement( smr, multiplicityElement, 2, &(multiplicity->frame) ) ) != 0 ) return( 1 );
    if( tpia_multiplicity_getTimeScaleFromElement( smr, multiplicityElement, &timeScale, &isDelayedNeutrons, &(multiplicity->timeScale) ) ) return( 1 );
    for( data = xData_getFirstElement( multiplicityElement ); data != NULL; data = xData_getNextElement( data ) ) {
        if( strcmp( data->name, "grouped" ) == 0 ) {
            if( tpia_misc_get2d_xShared_yHistogram_data_Grouped( smr, data, &(multiplicity->grouped) ) ) return( 1 ); }
        else if( strcmp( data->name, "pointwise" ) == 0 ) {
            if( ( multiplicity->pointwise = tpia_misc_get2dx_y_data( smr, data, &multiplicity->numberOfPointwise ) ) == NULL ) return( 1 ); }
        else {
            tpia_misc_setMessageError_Element( smr, NULL, multiplicityElement, __FILE__, __LINE__, 1, "unsupported multiplicity type = %s", data->name );
            return( 1 );
        }
    }
    return( 0 );
}
/*
************************************************************
*/
int tpia_multiplicity_getTimeScaleFromElement( statusMessageReporting *smr, xData_element *element, const char **timeScale, int *isDelayedNeutrons, 
        double *dTimeScale ) {

    const char *p;
    char *e;

    *isDelayedNeutrons = 0;
    *dTimeScale = 0.;
    *timeScale = xData_getAttributesValue( &(element->attributes), "timeScale" );
    if( *timeScale != NULL ) {
        if( strcmp( *timeScale, "prompt" ) ) {
            if( strncmp( *timeScale, "delayed", 7 ) ) {
                tpia_misc_setMessageError_Element( smr, NULL, element, __FILE__, __LINE__, 1, "invalid timeScale attribute = %s", *timeScale );
                return( 1 ); }
            else {
                for( p = &((*timeScale)[7]); *p; p++ ) if( !isspace( *p ) ) break;
                if( *p != '=' ) {
                    tpia_misc_setMessageError_Element( smr, NULL, element, __FILE__, __LINE__, 1, "timeScale attribute '%s' missing '='", *timeScale );
                    return( 1 );
                }
                p++;
                *dTimeScale = strtod( p, &e );
                if( *e != 0 ) {
                    tpia_misc_setMessageError_Element( smr, NULL, element, __FILE__, __LINE__, 1, "could not convert timeScale attribute '%s' to double",
                        *timeScale );
                    return( 1 );
                }
                *isDelayedNeutrons = 1;
            }
        }
    }
    return( 0 );
}

#if defined __cplusplus
}
#endif
