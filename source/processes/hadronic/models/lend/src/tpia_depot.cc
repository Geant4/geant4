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
#include <tpia_depot.h>

#if defined __cplusplus
namespace GIDI{
using namespace GIDI;
#endif

/*
************************************************************
*/
tpia_depot *tpia_depot_create( statusMessageReporting *smr, const char *projectileName ) {

    tpia_depot *depot;

    //if( ( depot = xData_malloc2( smr, sizeof( tpia_depot ), 0, "map" ) ) == NULL ) return( NULL );
    if( ( depot = (tpia_depot*) xData_malloc2( smr, sizeof( tpia_depot ), 0, "map" ) ) == NULL ) return( NULL );
    if( tpia_depot_initialize( smr, depot, projectileName ) ) depot = tpia_depot_free( depot, 0 );
    return( depot );
}
/*
************************************************************
*/
int tpia_depot_initialize( statusMessageReporting *smr, tpia_depot *depot, const char *projectileName ) {

    memset( depot, 0, sizeof( tpia_depot ) );
    depot->status = 0;
    depot->projectileName = NULL;
    depot->numberOfTargets = 0;
    depot->targets = NULL;
    depot->map = NULL;
    //if( ( depot->projectileName = xData_malloc2( smr, strlen( projectileName ) + 1, 0, "projectile" ) ) == NULL ) return( 1 );
    if( ( depot->projectileName = (char*) xData_malloc2( smr, strlen( projectileName ) + 1, 0, "projectile" ) ) == NULL ) return( 1 );
    return( 0 );
}
/*
************************************************************
*/
tpia_depot *tpia_depot_free( tpia_depot *depot, int freeMap ) {

    tpia_depot_release( depot, freeMap );
    xData_free( NULL, depot );
    return( NULL );
}
/*
************************************************************
*/
int tpia_depot_release( tpia_depot *depot, int freeMap ) {

    tpia_targetEntry *next, *targetEntry;

    if( depot->projectileName != NULL ) xData_free( NULL, depot->projectileName );
    for( targetEntry = depot->targets; targetEntry != NULL; targetEntry = next ) {
        next = targetEntry->next;
        tpia_target_free( NULL, targetEntry->target );
        xData_free( NULL, targetEntry );
    }
    depot->numberOfTargets = 0;
    depot->targets = NULL;
    //if( freeMap && ( depot->map != NULL ) ) depot->map = tpia_map_free( NULL, depot->map );
    if( freeMap && ( depot->map != NULL ) ) depot->map = (tpia_map*) tpia_map_free( NULL, depot->map );
    return( depot->status = 0 );
}
/*
************************************************************
*/
//int tpia_depot_setMap( statusMessageReporting *smr, tpia_depot *depot, tpia_map *map ) {
int tpia_depot_setMap( statusMessageReporting *, tpia_depot *depot, tpia_map *map ) {

    depot->map = map;
    return( 0 );
}
/*
************************************************************
*/
int tpia_depot_setMapFromFilename( statusMessageReporting *smr, tpia_depot *depot, const char *basePath, const char *mapFileName ) {

    if( ( depot->map = tpia_map_readFile( smr, basePath, mapFileName ) ) == NULL ) return( 1 );
    return( 0 );
}
/*
************************************************************
*/
tpia_map *tpia_depot_releaseMap( tpia_depot *depot ) {

    tpia_map *map = depot->map;

    depot->map = NULL;
    return( map );
}
/*
************************************************************
*/
int tpia_depot_freeMap( tpia_depot *depot ) {

    tpia_map *map = tpia_depot_releaseMap( depot );

    if( map != NULL ) tpia_map_free( NULL, map );
    return( 0 );
}
/*
************************************************************
*/
tpia_targetEntry *tpia_depot_getFirstTargetEntry( tpia_depot *depot ) {

    return( depot->targets );
}
/*
************************************************************
*/
tpia_targetEntry *tpia_depot_getNextTargetEntry( tpia_targetEntry *targetEntry ) {

    return( targetEntry->next );
}
/*
************************************************************
*/
tpia_target *tpia_depot_addTarget( statusMessageReporting *smr, tpia_depot *depot, const char *evaluation, const char *targetName ) {

    return( tpia_depot_addTargetFromMap( smr, depot, depot->map, evaluation, targetName ) );
}
/*
************************************************************
*/
tpia_target *tpia_depot_addTargetFromMap( statusMessageReporting *smr, tpia_depot *depot, tpia_map *map, const char *evaluation, const char *targetName ) {

    tpia_targetEntry *targetEntry;
    tpia_target *target;

    for( targetEntry = tpia_depot_getFirstTargetEntry( depot ); targetEntry != NULL; targetEntry = tpia_depot_getNextTargetEntry( targetEntry ) ) {
        if( !strcmp( targetEntry->target->targetID->name, targetName ) ) {
            smr_setMessageError( smr, NULL, __FILE__, __LINE__, 1, "depot already contains target = %s ", targetName );
            return( NULL );
        }
    }
    target = tpia_target_createReadFromMap( smr, map, evaluation, depot->projectileName, targetName );
    return( target );
}

#if defined __cplusplus
}
#endif
