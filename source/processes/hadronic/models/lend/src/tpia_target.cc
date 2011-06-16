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
#include <tpia_target.h>
#include <tpia_misc.h>

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

static int _tpia_target_releaseAndReturnOne( statusMessageReporting *smr, tpia_target *target );
/*
************************************************************
*/
tpia_target *tpia_target_create( statusMessageReporting *smr ) {

    tpia_target *target;

    //if( ( target = xData_malloc2( smr, sizeof( tpia_target ), 0, "target" ) ) == NULL ) return( NULL );
    if( ( target = (tpia_target*) xData_malloc2( smr, sizeof( tpia_target ), 0, "target" ) ) == NULL ) return( NULL );
    if( tpia_target_initialize( smr, target ) ) target = tpia_target_free( smr, target );
    return( target );
}
/*
************************************************************
*/
int tpia_target_initialize( statusMessageReporting *smr, tpia_target *target ) {

    memset( target, 0, sizeof( tpia_target ) );
    tpia_samplingMethods_initialize( smr, &(target->samplingMethods) );
    return( 0 );
}
/*
************************************************************
*/
tpia_target *tpia_target_createRead( statusMessageReporting *smr, const char *fileName ) {

    tpia_target *target;

    if( ( target = tpia_target_create( smr ) ) == NULL ) return( NULL );
    //if( tpia_target_read( smr, target, fileName ) != 0 ) target = xData_free( smr, target );
    if( tpia_target_read( smr, target, fileName ) != 0 ) target = (tpia_target*) xData_free( smr, target );
    return( target );
}
/*
************************************************************
*/
int tpia_target_readFromMap( statusMessageReporting *smr, tpia_target *target, tpia_map *map, const char *evaluation, const char *projectileName, 
        const char *targetName ) {

    char *targetPath;

    if( ( targetPath = tpia_map_findTarget( smr, map, evaluation, projectileName, targetName ) ) == NULL ) return( 1 );
    return( tpia_target_read( smr, target, targetPath ) );
}
/*
************************************************************
*/
tpia_target *tpia_target_createReadFromMap( statusMessageReporting *smr, tpia_map *map, const char *evaluation, const char *projectileName, 
        const char *targetName ) {

    char *targetPath;
    tpia_target *target;

    targetPath = tpia_map_findTarget( smr, map, evaluation, projectileName, targetName );
    if( targetPath == NULL ) return( NULL );
    target = tpia_target_createRead( smr, targetPath );
    xData_free( smr, targetPath );
    return( target );
}
/*
************************************************************
*/
tpia_target *tpia_target_free( statusMessageReporting *smr, tpia_target *target ) {

    tpia_target_release( smr, target );
    xData_free( smr, target );
    return( NULL );
}
/*
************************************************************
*/
int tpia_target_release( statusMessageReporting *smr, tpia_target *target ) {

    int i;

    //target->path = xData_free( smr, target->path );
    target->path = (char*) xData_free( smr, target->path );
    //target->absPath = xData_free( smr, target->absPath );
    target->absPath = (char*) xData_free( smr, target->absPath );
    xData_releaseAttributionList( smr, &(target->attributes) );
    for( i = 0; i < target->nHeatedTargets; i++ ) {
        //target->heatedTargets[i].path = xData_free( smr, target->heatedTargets[i].path );
        target->heatedTargets[i].path = (char*) xData_free( smr, target->heatedTargets[i].path );
        //target->heatedTargets[i].contents = xData_free( smr, target->heatedTargets[i].contents );
        target->heatedTargets[i].contents = (char*) xData_free( smr, target->heatedTargets[i].contents );
            if( target->heatedTargets[i].heatedTarget != NULL ) tpia_target_heated_free( smr, target->heatedTargets[i].heatedTarget );
    }
    //target->heatedTargets = xData_free( smr, target->heatedTargets );
    target->heatedTargets = (tpia_target_heated_info*) xData_free( smr, target->heatedTargets );
    //target->readHeatedTargets = xData_free( smr, target->readHeatedTargets );
    target->readHeatedTargets = (tpia_target_heated_info**) xData_free( smr, target->readHeatedTargets );
    tpia_target_initialize( smr, target );
    return( 0 );
}
/*
************************************************************
*/
int tpia_target_read( statusMessageReporting *smr, tpia_target *target, const char *fileName ) {
/*
*   If a target has already been read into this target, user must have called tpia_target_release before calling this routine.
*   Otherwise, there will be memory leaks.
*/
    xData_document *doc;
    xData_element *element, *child;
    int i, iHeated, nHeated = 0, status = 1;
    double temperature;
    //fix for gcc4.6 warings 110602
    //char *pReturnValue, *name;
    char *name;
    char const *contents;

    tpia_target_initialize( smr, target );
    if( ( target->path = xDataMisc_allocateCopyString2( smr, fileName, "path" ) ) == NULL ) return( status );
    if( ( target->absPath = xDataMisc_getAbsPath( smr, fileName ) ) == NULL ) return( _tpia_target_releaseAndReturnOne( smr, target ) );
    if( ( doc = xData_parseReadFile( smr, fileName, NULL, NULL ) ) == NULL ) return( _tpia_target_releaseAndReturnOne( smr, target ) );
    element = xData_getDocumentsElement( doc );
    if( strcmp( element->name, "xTarget" ) != 0 ) {
        tpia_misc_setMessageError_Element( smr, NULL, element, __FILE__, __LINE__, 1, "input file's top element must be xTarget and not %s", element->name ); }
    else {
        //pReturnValue = ( xData_copyAttributionList( smr, &(target->attributes), &(element->attributes) ) ) ? NULL : target->path;
        //fix for gcc4.6 warings 110602
        xData_copyAttributionList( smr, &(target->attributes),&(element->attributes) );
        name = tpia_misc_pointerToAttributeIfAllOk2( smr, element, 1, &(target->attributes), "projectile" );
        if( smr_isOk( smr ) ) target->projectileID = tpia_particle_getInternalID( smr, name );
        if( smr_isOk( smr ) && ( name = tpia_misc_pointerToAttributeIfAllOk2( smr, element, 1, &(target->attributes), "target" ) ) != NULL ) {
            if( smr_isOk( smr ) && ( target->targetID = tpia_particle_getInternalID( smr, name ) ) != NULL ) {
                status = 0;
                for( nHeated = 0, child = xData_getFirstElement( element ); child != NULL; nHeated++, child = xData_getNextElement( child ) ) {
                    if( strcmp( child->name, "target" ) != 0 ) {
                        tpia_misc_setMessageError_Element( smr, NULL, element, __FILE__, __LINE__, 1, "element can only have target sub-elements%s", 
                                element->name );
                        status = 1;
                        break;
                    }
                }
                if( status == 0 ) {
                    //if( ( target->heatedTargets = xData_malloc2( smr, nHeated * sizeof( tpia_target_heated_info ), 1, "heatedTargets" ) ) == NULL ) {
                    if( ( target->heatedTargets = (tpia_target_heated_info*) xData_malloc2( smr, nHeated * sizeof( tpia_target_heated_info ), 1, "heatedTargets" ) ) == NULL ) {
                        status = 1; }
                    else {
                        //if( ( target->readHeatedTargets = xData_malloc2( smr, nHeated * sizeof( tpia_target_heated_info * ), 1, "heatedTargets" ) ) == NULL ) 
                        if( ( target->readHeatedTargets = (tpia_target_heated_info**) xData_malloc2( smr, nHeated * sizeof( tpia_target_heated_info * ), 1, "heatedTargets" ) ) == NULL ) 
                                status = 1;
                    }
                    for( nHeated = 0, child = xData_getFirstElement( element ); ( status == 0 ) && ( child != NULL ); 
                            nHeated++, child = xData_getNextElement( child ) ) {
                        if( ( i = xData_convertAttributeToDouble( smr, child, "temperature", &temperature ) ) != 0 ) {
                            if( i > 0 ) smr_setMessageError( smr, NULL, __FILE__, __LINE__, 1, "target does not have a temperature attribute" );
                            status = 1;
                            break;
                        }
                        for( iHeated = 0; iHeated < nHeated; iHeated++ ) if( target->heatedTargets[iHeated].temperature > temperature ) break;
                        if( iHeated < nHeated ) for( i = nHeated; i >= iHeated; i-- ) target->heatedTargets[i+1] = target->heatedTargets[i];
                        target->heatedTargets[iHeated].temperature = temperature;
                        target->heatedTargets[iHeated].path = NULL;
                        target->heatedTargets[iHeated].contents = NULL;
                        target->heatedTargets[iHeated].heatedTarget = NULL;
                        if( ( contents = xData_getAttributesValueInElement( child, "contents" ) ) != NULL ) {
                            if( ( target->heatedTargets[iHeated].contents = xDataMisc_allocateCopyString2( smr, contents, "contents" ) ) == NULL ) {
                                status = 1;
                                break;
                            }
                        }
                        if( ( contents = xData_getAttributesValueInElement( child, "file" ) ) == NULL ) {
                            status = 1;
                            break;
                        }
                        //if((target->heatedTargets[iHeated].path = xData_malloc2(smr, strlen( target->absPath ) + strlen( contents ) + 2, 0, "path")) == NULL) {
                        if((target->heatedTargets[iHeated].path = (char*) xData_malloc2(smr, strlen( target->absPath ) + strlen( contents ) + 2, 0, "path")) == NULL) {
                            status = 1;
                            break;
                        }
                        strcpy( target->heatedTargets[iHeated].path, target->absPath );
                        *strrchr( target->heatedTargets[iHeated].path, '/' ) = 0;
                        strcat( target->heatedTargets[iHeated].path, "/" );
                        strcat( target->heatedTargets[iHeated].path, contents );
                        target->nHeatedTargets++;
                    }
                }
            }
        }
    }
    xData_parseFree( smr, doc );
    if( status == 0 ) {
        for( i = 0; i < nHeated; i++ ) target->heatedTargets[i].ordinal = i;
        for( i = 0; i < nHeated; i++ ) if( target->heatedTargets[i].contents == NULL ) break;
        if( i == nHeated ) i = 0;                                           /* All heated targets are crossSection only. */
        if( tpia_target_readHeatedTarget( smr, target, i, 0 ) == 0 ) {
            target->baseHeatedTarget = target->heatedTargets[i].heatedTarget; }
        else {
            tpia_target_release( NULL, target );
            status = 1;
        } }
    else {
        tpia_target_release( smr, target );
    }
    return( status );
}
/*
************************************************************
*/
//char *tpia_target_getAttributesValue( statusMessageReporting *smr, tpia_target *target, char const *name ) {
char *tpia_target_getAttributesValue( statusMessageReporting *, tpia_target *target, char const *name ) {

    return( xData_getAttributesValue( &(target->attributes), name ) );
}
/*
************************************************************
*/
//int tpia_target_getTemperatures( statusMessageReporting *smr, tpia_target *target, double *temperatures ) {
int tpia_target_getTemperatures( statusMessageReporting *, tpia_target *target, double *temperatures ) {

    int i;

    if( temperatures != NULL ) for( i = 0; i < target->nHeatedTargets; i++ ) temperatures[i] = target->heatedTargets[i].temperature;
    return( target->nHeatedTargets );
}
/*
************************************************************
*/
int tpia_target_readHeatedTarget( statusMessageReporting *smr, tpia_target *target, int index, int checkElememtsForAccess ) {

    int i;

    if( ( index < 0 ) || ( index >= target->nHeatedTargets ) ) {
        smr_setMessageError( smr, NULL, __FILE__, __LINE__, 1, "temperature index = %d out of range (0 <= index < %d", index, target->nHeatedTargets );
        return( -1 );
    }
    if( target->heatedTargets[index].heatedTarget != NULL ) return( 1 );
    target->heatedTargets[index].heatedTarget = tpia_target_heated_createRead( smr, target->heatedTargets[index].path, checkElememtsForAccess );
    if( target->heatedTargets[index].heatedTarget != NULL ) {
        target->heatedTargets[index].heatedTarget->ordinal = target->heatedTargets[index].ordinal;
        for( i = target->nReadHeatedTargets; i > 0; i-- ) {
            if( target->readHeatedTargets[i-1]->temperature < target->heatedTargets[index].temperature ) break;
            target->readHeatedTargets[i] = target->readHeatedTargets[i-1];
        }
        target->readHeatedTargets[i] = &(target->heatedTargets[i]);
        target->nReadHeatedTargets++;
    }
    return( ( target->heatedTargets[index].heatedTarget == NULL ? -1 : 0 ) );
}
/*
************************************************************
*/
tpia_target_heated *tpia_target_getHeatedTargetAtIndex_ReadIfNeeded( statusMessageReporting *smr, tpia_target *target, int index ) {

    if( ( index < 0 ) || ( index >= target->nHeatedTargets ) ) {
        smr_setMessageError( smr, NULL, __FILE__, __LINE__, 1, "temperature index = %d out of range (0 <= index < %d", index, target->nHeatedTargets );
        return( NULL );
    }
    if( target->heatedTargets[index].heatedTarget == NULL ) tpia_target_readHeatedTarget( smr, target, index, 0 );
    return( target->heatedTargets[index].heatedTarget );
}
/*
************************************************************
*/
tpia_target_heated *tpia_target_getHeatedTargetAtTIndex( statusMessageReporting *smr, tpia_target *target, int index ) {

    if( ( index < 0 ) || ( index >= target->nHeatedTargets ) ) {
        smr_setMessageError( smr, NULL, __FILE__, __LINE__, 1, "temperature index = %d out of range (0 <= index < %d", index, target->nHeatedTargets );
        return( NULL );
    }
    if( target->heatedTargets[index].heatedTarget == NULL ) {
        smr_setMessageError( smr, NULL, __FILE__, __LINE__, 1, "temperature index = %d not read in", index );
        return( NULL );
    }
    return( target->heatedTargets[index].heatedTarget );
}
/*
************************************************************
*/
int tpia_target_numberOfChannels( statusMessageReporting *smr, tpia_target *target ) {

    return( tpia_target_heated_numberOfChannels( smr, target->baseHeatedTarget ) );
}
/*
************************************************************
*/
int tpia_target_numberOfProductionChannels( statusMessageReporting *smr, tpia_target *target ) {

    return( tpia_target_heated_numberOfProductionChannels( smr, target->baseHeatedTarget ) );
}
/*
************************************************************
*/
xData_Int tpia_target_getEnergyGridAtTIndex( statusMessageReporting *smr, tpia_target *target, int index, double **energyGrid ) {

    tpia_target_heated *heated = tpia_target_getHeatedTargetAtTIndex( smr, target, index );

    if( !smr_isOk( smr ) ) return( -1 );
    return( tpia_target_heated_getEnergyGrid( smr, heated, energyGrid ) );
}
/*
************************************************************
*/
tpia_1dData *tpia_target_getTotalCrossSectionAtTIndex( statusMessageReporting *smr, tpia_target *target, int index, int crossSectionType ) {

    tpia_target_heated *heated = tpia_target_getHeatedTargetAtTIndex( smr, target, index );

    if( !smr_isOk( smr ) ) return( NULL );
    if( crossSectionType == tpia_crossSectionType_grouped ) {
        return( &(heated->totalCrossSectionGrouped) ); }
    else if( crossSectionType == tpia_crossSectionType_pointwise ) {
        return( &(heated->totalCrossSectionPointwise) );
    }
    smr_setMessageError( smr, NULL, __FILE__, __LINE__, 1, "Invalue crossSectionType = %d", crossSectionType );
    return( NULL );
}
/*
************************************************************
*/
double tpia_target_getTotalCrossSectionAtTAndE( statusMessageReporting *smr, tpia_target *target, double T, xData_Int iEg, double e_in,
    int crossSectionType ) {

    int i;
    double xsec = 0., xsec1, xsec2;

    for( i = 0; i < target->nReadHeatedTargets; i++ ) if( target->readHeatedTargets[i]->temperature > T ) break;
    if( i == 0 ) {
        xsec = tpia_target_heated_getTotalCrossSectionAtE( smr, target->readHeatedTargets[0]->heatedTarget, iEg, e_in, crossSectionType ); }
    else if( i == target->nReadHeatedTargets ) {
        xsec = tpia_target_heated_getTotalCrossSectionAtE( smr, target->readHeatedTargets[i-1]->heatedTarget, iEg, e_in, crossSectionType ); }
    else {
        xsec1 = tpia_target_heated_getTotalCrossSectionAtE( smr, target->readHeatedTargets[i-1]->heatedTarget, iEg, e_in, crossSectionType );
        xsec2 = tpia_target_heated_getTotalCrossSectionAtE( smr, target->readHeatedTargets[i  ]->heatedTarget, iEg, e_in, crossSectionType );
        xsec = ( ( target->readHeatedTargets[i]->temperature - T ) * xsec1 + ( T - target->readHeatedTargets[i-1]->temperature ) * xsec2 ) / 
            ( target->readHeatedTargets[i]->temperature - target->readHeatedTargets[i-1]->temperature );
    }

    return( xsec );
}
/*
************************************************************
*/
double tpia_target_getIndexChannelCrossSectionAtE( statusMessageReporting *smr, tpia_target *target, int index, double T, xData_Int iEg, double e_in,
    int crossSectionType ) {

    int i;
    double xsec = 0., xsec1, xsec2;

    for( i = 0; i < target->nReadHeatedTargets; i++ ) if( target->readHeatedTargets[i]->temperature > T ) break;
    if( i == 0 ) {
        xsec = tpia_target_heated_getIndexChannelCrossSectionAtE( smr, target->readHeatedTargets[0]->heatedTarget, index, iEg, e_in, crossSectionType ); }
    else if( i == target->nReadHeatedTargets ) {
        xsec = tpia_target_heated_getIndexChannelCrossSectionAtE( smr, target->readHeatedTargets[i-1]->heatedTarget, index, iEg, e_in, crossSectionType ); }
    else {
        xsec1 = tpia_target_heated_getIndexChannelCrossSectionAtE(smr, target->readHeatedTargets[i-1]->heatedTarget, index, iEg, e_in, crossSectionType);
        xsec2 = tpia_target_heated_getIndexChannelCrossSectionAtE(smr, target->readHeatedTargets[i  ]->heatedTarget, index, iEg, e_in, crossSectionType);
        xsec = ( ( target->readHeatedTargets[i]->temperature - T ) * xsec1 + ( T - target->readHeatedTargets[i-1]->temperature ) * xsec2 ) / 
            ( target->readHeatedTargets[i]->temperature - target->readHeatedTargets[i-1]->temperature );
    }

    return( xsec );
}
/*
************************************************************
*/
//int tpia_target_sampleIndexChannelProductsAtE( statusMessageReporting *smr, tpia_target *target, int index, double T, 
int tpia_target_sampleIndexChannelProductsAtE( statusMessageReporting *smr, tpia_target *target, int index, double , 
        tpia_decaySamplingInfo *decaySamplingInfo, int nProductData, tpia_productOutgoingData *productData ) {

    return( tpia_target_heated_sampleIndexChannelProductsAtE( smr, target->baseHeatedTarget, index, decaySamplingInfo,
        nProductData, productData ) );
}
/*
************************************************************
*/
static int _tpia_target_releaseAndReturnOne( statusMessageReporting *smr, tpia_target *target ) {

    tpia_target_release( smr, target );
    return( 1 );
}

#if defined __cplusplus
}
#endif
