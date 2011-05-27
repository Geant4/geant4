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
#include <limits.h>
#ifdef WIN32
   #include <Shlwapi.h>
#endif
#include <tpia_map.h>

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

static tpia_mapEntry *_tpia_map_addEntry( statusMessageReporting *smr, tpia_map *map, enum tpia_mapEntry_type type, const char *schema, const char *path,
    const char *evaluation, const char *projectile, const char *target );
static char *_tpia_map_findTarget2( statusMessageReporting *smr, tpia_map *map, const char *evaluation, const char *projectile, const char *target );
static int _tpia_map_findAllOfTarget2( statusMessageReporting *smr, tpia_map *mapAllOfTarget, tpia_map *map, const char *projectile, const char *targetName );
static int _tpia_map_walkTree2( statusMessageReporting *smr, tpia_map *map, int level, int (*handler)( tpia_mapEntry *entry, int level, void *userData), 
    void *userData );
static void _tpia_map_simpleWrite2( FILE *f, tpia_map *map, int level );
static int _tpia_map_smrUserInterface( void *userData, char **str );
/*
************************************************************
*/
tpia_map *tpia_map_create( statusMessageReporting *smr ) {

    tpia_map *map;
    
    //if( ( map = xData_malloc2( smr, sizeof( tpia_map ), 0, "map" ) ) == NULL ) return( NULL );
    if( ( map = (tpia_map*) xData_malloc2( smr, sizeof( tpia_map ), 0, "map" ) ) == NULL ) return( NULL );
    //if( tpia_map_initialize( smr, map ) ) map = tpia_map_free( NULL, map );
    if( tpia_map_initialize( smr, map ) ) map = (tpia_map*) tpia_map_free( NULL, map );
    return( map );
}
/*
************************************************************
*/
//int tpia_map_initialize( statusMessageReporting *smr, tpia_map *map ) {
int tpia_map_initialize( statusMessageReporting *, tpia_map *map ) {

    memset( map, 0, sizeof( tpia_map ) );
    map->status = tpia_map_status_Ok;
    map->smrUserInterface.smrUserInterface = _tpia_map_smrUserInterface;
    map->smrUserInterface.map = map;
    map->path = NULL;
    map->mapFileName = NULL;
    map->numberOfEntries = 0;
    map->mapEntries = NULL;
    return( 0 );
}
/*
************************************************************
*/
tpia_map *tpia_map_readFile( statusMessageReporting *smr, const char *basePath, const char *mapFileName ) {
/*
*   If an error occurrs, map is freed and NULL is returned.
*/
    int n = 0;
    xData_document *doc;
    xData_element *element;
    xData_element *child;
    tpia_map *map;
    const char *evaluation, *projectile, *targetName, *path, *schema;
#ifndef WIN32
    char realPath[2 * ( PATH_MAX + 1 )], *p = &(realPath[PATH_MAX+1]);
#endif
#ifdef WIN32
    char realPath[2 * ( _MAX_PATH + 1 )], *p = &(realPath[_MAX_PATH+1]);
#endif

    if( ( map = tpia_map_create( smr ) ) == NULL ) return( NULL );

    if( ( basePath == NULL ) || ( mapFileName[0] == '/' ) ) {
        strcpy( realPath, mapFileName ); }
    else {
        strcpy( realPath, basePath );
        strcat( realPath, "/" );
        strcat( realPath, mapFileName );
    }
#ifndef WIN32
    if( realpath( realPath, p ) == NULL ) {
        smr_setMessageError( smr, NULL, __FILE__, __LINE__, tpia_map_status_mapParsing, "No map file %s\n", mapFileName );
        //return( tpia_map_free( NULL, map ) );
        return( (tpia_map*) tpia_map_free( NULL, map ) );
    }
#endif
    n = strlen( p ) + 2;
    //if( ( map->path = xData_malloc2( smr, 2 * n, 0, "map->path" ) ) == NULL ) return( tpia_map_free( NULL, map ) );
    //if( ( map->path = (char*) xData_malloc2( smr, 2 * n, 0, "map->path" ) ) == NULL ) return( tpia_map_free( NULL, map ) );
    if( ( map->path = (char*) xData_malloc2( smr, 2 * n, 0, "map->path" ) ) == NULL ) return( (tpia_map*) tpia_map_free( NULL, map ) );
    map->mapFileName = &(map->path[n + 1]);
    strcpy( map->mapFileName, p );
    strcpy( map->path, p );
    if( ( p = strrchr( map->path, '/' ) ) != NULL ) {
        *p = 0; }
    else {
        strcpy( map->path, "." );
    }

    //if( ( doc = xData_parseReadFile( smr, map->mapFileName, NULL, NULL ) ) == NULL ) return( tpia_map_free( NULL, map ) );
    if( ( doc = xData_parseReadFile( smr, map->mapFileName, NULL, NULL ) ) == NULL ) return( (tpia_map*) tpia_map_free( NULL, map ) );

    element = xData_getDocumentsElement( doc );
    for( child = xData_getFirstElement( element ); child != NULL; child = xData_getNextElement( child ) ) {
        if( !strcmp( child->name, "path" ) ) {
            if( ( path = xData_getAttributesValueInElement( child , "path" ) ) == NULL ) {
                smr_setMessageError( smr, &(map->smrUserInterface), __FILE__, __LINE__, tpia_map_status_mapParsing, "path missing path attribute" );
                break;
            }
            if( ( projectile = xData_getAttributesValueInElement( child , "projectile" ) ) == NULL ) {
                smr_setMessageError( smr, &(map->smrUserInterface), __FILE__, __LINE__, tpia_map_status_mapParsing, "path missing projectile attribute" );
                break;
            }
            tpia_map_addPath( smr, map, path, projectile ); }
        else if( !strcmp( child->name, "target" ) ) {
            if( ( schema = xData_getAttributesValueInElement( child , "schema" ) ) == NULL ) {
                smr_setMessageError( smr, &(map->smrUserInterface), __FILE__, __LINE__, tpia_map_status_mapParsing, "target missing 'schema' attribute" );
                break;
            }
            if( ( path = xData_getAttributesValueInElement( child , "path" ) ) == NULL ) {
                smr_setMessageError( smr, &(map->smrUserInterface), __FILE__, __LINE__, tpia_map_status_mapParsing, "target missing 'path' attribute" );
                break;
            }
            if( ( evaluation = xData_getAttributesValueInElement( child , "evaluation" ) ) == NULL ) {
                smr_setMessageError( smr, &(map->smrUserInterface), __FILE__, __LINE__, tpia_map_status_mapParsing, "target missing 'evaluation' attribute" );
                break;
            }
            if( ( projectile = xData_getAttributesValueInElement( child , "projectile" ) ) == NULL ) {
                smr_setMessageError( smr, &(map->smrUserInterface), __FILE__, __LINE__, tpia_map_status_mapParsing, "target missing 'projectile' attribute" );
                break;
            }
            if( ( targetName = xData_getAttributesValueInElement( child , "target" ) ) == NULL ) {
                smr_setMessageError( smr, &(map->smrUserInterface), __FILE__, __LINE__, tpia_map_status_mapParsing, "target missing 'target' attribute" );
                break;
            }
            tpia_map_addTarget( smr, map, schema, path, evaluation, projectile, targetName ); }
        else {
            smr_setMessageError( smr, &(map->smrUserInterface), __FILE__, __LINE__, tpia_map_status_mapParsing, "invalid element = %s", child->name );
        }
        if( !smr_isOk( smr ) ) break;
    }
    xData_parseFree( smr, doc );
    //if( !smr_isOk( smr ) ) map = tpia_map_free( NULL, map );
    if( !smr_isOk( smr ) ) map = (tpia_map*) tpia_map_free( NULL, map );
    return( map );
}
/*
************************************************************
*/
void *tpia_map_free( statusMessageReporting *smr, tpia_map *map ) {

    tpia_map_release( smr, map );
    xData_free( smr, map );
    return( NULL );
}
/*
************************************************************
*/
void tpia_map_release( statusMessageReporting *smr, tpia_map *map ) {

    tpia_mapEntry *entry, *next;

    if( map->path != NULL ) xData_free( NULL, map->path );
    for( entry = map->mapEntries; entry != NULL; entry = next ) {
        next = entry->next;
        if( entry->schema != NULL ) xData_free( NULL, entry->schema );
        if( entry->path != NULL ) xData_free( NULL, entry->path );
        if( entry->evaluation != NULL ) xData_free( NULL, entry->evaluation );
        if( entry->projectile != NULL ) xData_free( NULL, entry->projectile );
        if( entry->targetName != NULL ) xData_free( NULL, entry->targetName );
        if( entry->map != NULL ) tpia_map_free( smr, entry->map );
        xData_free( NULL, entry );
    }
    map->numberOfEntries = 0;
    map->mapEntries = NULL;
    map->status = tpia_map_status_Ok;
}
/*
************************************************************
*/
tpia_mapEntry *tpia_map_getFirstEntry( tpia_map *map ) {

    return( map->mapEntries );
}
/*
************************************************************
*/
tpia_mapEntry *tpia_map_getNextEntry( tpia_mapEntry *entry ) {

    return( entry->next );
}
/*
************************************************************
*/
int tpia_map_addTarget( statusMessageReporting *smr, tpia_map *map, const char *schema, const char *path, const char *evaluation, const char *projectile, const char *target ) {

    return( _tpia_map_addEntry( smr, map, tpia_mapEntry_type_target, schema, path, evaluation, projectile, target ) != NULL );
}
/*
************************************************************
*/
int tpia_map_addPath( statusMessageReporting *smr, tpia_map *map, const char *path, const char *projectile ) {

    tpia_mapEntry *entry = _tpia_map_addEntry( smr, map, tpia_mapEntry_type_path, NULL, path, NULL, projectile, NULL );

    if( entry != NULL ) {
        if( ( entry->map = tpia_map_readFile( smr, map->path, entry->path ) ) == NULL ) entry = NULL;
    }
    return( entry != NULL );
}
/*
************************************************************
*/
static tpia_mapEntry *_tpia_map_addEntry( statusMessageReporting *smr, tpia_map *map, enum tpia_mapEntry_type type, const char *schema, const char *path, 
    const char *evaluation, const char *projectile, const char *targetName ) {

    tpia_mapEntry *p;
    tpia_mapEntry *entry;

    //if( ( entry = xData_malloc2( smr, sizeof( tpia_mapEntry ), 1, "entry" ) ) == NULL ) return( NULL );
    if( ( entry = (tpia_mapEntry*) xData_malloc2( smr, sizeof( tpia_mapEntry ), 1, "entry" ) ) == NULL ) return( NULL );
    entry->next = NULL;
    entry->type = type;
    entry->path = NULL;
    entry->map = NULL;
    if( path != NULL ) {
        //if( ( entry->path = xData_malloc2( smr, strlen( path ) + 1, 0, "path" ) ) == NULL ) {
        if( ( entry->path = (char*) xData_malloc2( smr, strlen( path ) + 1, 0, "path" ) ) == NULL ) {
            xData_free( smr, entry );
            return( NULL );
        }
        strcpy( entry->path, path );
    }
    entry->evaluation = NULL;
    if( evaluation != NULL ) {
        //if( ( entry->evaluation = xData_malloc2( smr, strlen( evaluation ) + 1, 0, "evaluation" ) ) == NULL ) {
        if( ( entry->evaluation = (char*) xData_malloc2( smr, strlen( evaluation ) + 1, 0, "evaluation" ) ) == NULL ) {
            xData_free( smr, entry->path );
            xData_free( smr, entry );
            return( NULL );
        }
        strcpy( entry->evaluation, evaluation );
    }
    entry->projectile = NULL;
    if( projectile != NULL ) {
        //if( ( entry->projectile = xData_malloc2( smr, strlen( projectile ) + 1, 0, "projectile" ) ) == NULL ) {
        if( ( entry->projectile = (char*) xData_malloc2( smr, strlen( projectile ) + 1, 0, "projectile" ) ) == NULL ) {
            xData_free( smr, entry->evaluation );
            xData_free( smr, entry->path );
            xData_free( smr, entry );
            return( NULL );
        }
        strcpy( entry->projectile, projectile );
    }
    entry->targetName = NULL;
    if( targetName != NULL ) {
        //if( ( entry->targetName = xData_malloc2( smr, strlen( targetName ) + 1, 0, "target" ) ) == NULL ) {
        if( ( entry->targetName = (char*) xData_malloc2( smr, strlen( targetName ) + 1, 0, "target" ) ) == NULL ) {
            xData_free( smr, entry->path );
            xData_free( smr, entry->evaluation );
            xData_free( smr, entry->projectile );
            xData_free( smr, entry );
            return( NULL );
        }
        strcpy( entry->targetName, targetName );
    }
    entry->schema = NULL;
    if( schema != NULL ) {
        //if( ( entry->schema = xData_malloc2( smr, strlen( schema ) + 1, 0, "schema" ) ) == NULL ) {
        if( ( entry->schema = (char*) xData_malloc2( smr, strlen( schema ) + 1, 0, "schema" ) ) == NULL ) {
            xData_free( smr, entry->path );
            xData_free( smr, entry->evaluation );
            xData_free( smr, entry->projectile );
            xData_free( smr, entry->targetName );
            xData_free( smr, entry );
            return( NULL );
        }
        strcpy( entry->schema, schema );
    }

    if( map->mapEntries == NULL ) {
        map->mapEntries = entry; }
    else {
        for( p = map->mapEntries; p->next != NULL; p = p->next ){;}
        p->next = entry;
    }
    map->numberOfEntries++;
    return( entry );
}
/*
************************************************************
*/
char *tpia_map_findTarget( statusMessageReporting *smr, tpia_map *map, const char *evaluation, const char *projectile, const char *targetName ) {
/*
* Calling routine must free returned pointer.
*/
    char *path;

    if( map->status != tpia_map_status_Ok ) return( NULL );

    path = _tpia_map_findTarget2( smr, map, evaluation, projectile, targetName );
    if( ( path == NULL ) && smr_isOk( smr ) ) {
        if( evaluation == NULL ) {
            smr_setMessageInfo( smr, &(map->smrUserInterface), __FILE__, __LINE__, 1, "target %s for projectile %s not found", targetName, projectile ); }
        else {
            smr_setMessageInfo( smr, &(map->smrUserInterface), __FILE__, __LINE__, 1, "target %s for projectile %s and evaluation %s not found", targetName, projectile, evaluation );
        }
    }
    return( path );
}
/*
************************************************************
*/
static char *_tpia_map_findTarget2( statusMessageReporting *smr, tpia_map *map, const char *evaluation, const char *projectile, const char *targetName ) {

    tpia_mapEntry *entry;
    char *path = NULL;
    int n, status;

    for( entry = map->mapEntries; entry != NULL; entry = entry->next ) {
        switch( entry->type ) {
        case tpia_mapEntry_type_target :
            if( !strcmp( projectile, entry->projectile ) && ( !strcmp( targetName, entry->targetName ) ) ) {
                if( evaluation == NULL ) {
                    status = 1; }
                else {
                    status = !strcmp( evaluation,  entry->evaluation );
                }
                if( status ) {
                    n = strlen( map->path ) + 1 + strlen( entry->path ) + 1;
                    //if( ( path = xData_malloc2( smr, n, 0, "path" ) ) == NULL ) return( NULL );
                    if( ( path = (char*) xData_malloc2( smr, n, 0, "path" ) ) == NULL ) return( NULL );
                    strcpy( path, map->path );
                    strcat( path, "/" );
                    if( entry->path[0] == '/' ) {
                        strcpy( path, entry->path ); }
                    else {
                        strcat( path, entry->path );
                    }
                    return( path );
                }
            }
            break;
        case tpia_mapEntry_type_path :
            if( !strcmp( projectile, entry->projectile ) ) {
                if( ( path = _tpia_map_findTarget2( smr, entry->map, evaluation, projectile, targetName ) ) != NULL ) return( path );
            }
            break;
        default :
            smr_setMessageInfo( smr, &(map->smrUserInterface), __FILE__, __LINE__, tpia_map_status_UnknownType, "unknown type = %d", entry->type );
            return( NULL );
        }
    }
    return( NULL );
}
/*
************************************************************
*/
tpia_map *tpia_map_findAllOfTarget( statusMessageReporting *smr, tpia_map *map, const char *projectile, const char *targetName ) {
/*
* Calling routine must free returned pointer.
*/
    int status;
    tpia_map *mapAllOfTarget;
    
    if( map->status != tpia_map_status_Ok ) return( NULL );
    if( ( mapAllOfTarget = tpia_map_create( smr ) ) == NULL ) return( NULL );
    status = _tpia_map_findAllOfTarget2( smr, mapAllOfTarget, map, projectile, targetName );
    //if( ( status != 0 ) ) mapAllOfTarget = tpia_map_free( smr, mapAllOfTarget );
    if( ( status != 0 ) ) mapAllOfTarget = (tpia_map*) tpia_map_free( smr, mapAllOfTarget );
    return( mapAllOfTarget );
}
/*
************************************************************
*/
static int _tpia_map_findAllOfTarget2( statusMessageReporting *smr, tpia_map *mapAllOfTarget, tpia_map *map, const char *projectile, const char *targetName ) {

    tpia_mapEntry *entry;

    for( entry = map->mapEntries; entry != NULL; entry = entry->next ) {
        switch( entry->type ) {
        case tpia_mapEntry_type_target :
            if( !strcmp( projectile, entry->projectile ) && ( !strcmp( targetName, entry->targetName ) ) ) {
                if( _tpia_map_addEntry( smr, mapAllOfTarget, entry->type, entry->schema, entry->path, entry->evaluation, entry->projectile, 
                    entry->targetName ) == NULL ) return( 1 );
            }
            break;
        case tpia_mapEntry_type_path :
            if( !strcmp( projectile, entry->projectile ) ) {
                if( _tpia_map_findAllOfTarget2( smr, mapAllOfTarget, entry->map, projectile, targetName ) != 0 ) return( 1 );
            }
            break;
        default :
            smr_setMessageInfo( smr, &(map->smrUserInterface), __FILE__, __LINE__, tpia_map_status_UnknownType, "unknown type = %d", entry->type );
            return( 1 );
        }
    }
    return( 0 );
}
/*
************************************************************
*/
char *tpia_map_getFullPath( statusMessageReporting *smr, tpia_map *map, const char *endPath ) {

    char *path;

    if( endPath[0] == '/' ) {
        //if( ( path = xData_malloc2( smr, strlen( endPath ) + 1, 0, "path" ) ) == NULL ) return( NULL );
        if( ( path = (char*) xData_malloc2( smr, strlen( endPath ) + 1, 0, "path" ) ) == NULL ) return( NULL );
        path[0] = 0; }
    else {
        //if( ( path = xData_malloc2( smr, strlen( map->path ) + strlen( endPath ) + 2, 0, "path" ) ) == NULL ) return( NULL );
        if( ( path = (char*) xData_malloc2( smr, strlen( map->path ) + strlen( endPath ) + 2, 0, "path" ) ) == NULL ) return( NULL );
        strcpy( path, map->path );
        strcat( path, "/" );
    }
    strcat( path, endPath );
    return( path );
}
/*
************************************************************
*/
int tpia_map_walkTree( statusMessageReporting *smr, tpia_map *map, int (*handler)( tpia_mapEntry *entry, int level, void *userData), void *userData ) {

    return( _tpia_map_walkTree2( smr, map, 0, handler, userData ) );
}
/*
************************************************************
*/
static int _tpia_map_walkTree2( statusMessageReporting *smr, tpia_map *map, int level, int (*handler)( tpia_mapEntry *entry, int level, void *userData), 
    void *userData ) {
    
    tpia_mapEntry *entry;

    for( entry = map->mapEntries; entry != NULL; entry = entry->next ) {
        if( handler( entry, level, userData ) != 0 ) return( 1 );
        if( entry->type == tpia_mapEntry_type_path ) if( _tpia_map_walkTree2( smr, entry->map, level + 1, handler, userData ) != 0 ) return( 1 );
    }
    return( 0 );
}
/*
************************************************************
*/
char *tpia_map_toXMLString( statusMessageReporting *smr, tpia_map *map ) {

    tpia_mapEntry *entry;
    char *s, *p;
    char targetFormat[] = "<target schema=\"%s\" evaluation=\"%s\" projectile=\"%s\" target=\"%s\" path=\"%s\"/>\n";
    char pathFormat[] = "<path projectile=\"%s\" path=\"%s\"/>\n";
    char start[] = "<map>\n";
    char end[] = "</map>";
    int n = 0, nStart = strlen( start ), nEnd = strlen( end );
    int nTarget = strlen( targetFormat ) - 10, nPath = strlen( pathFormat ) - 4;

    if( map->status != tpia_map_status_Ok ) return( NULL );

    n = nStart + nEnd + 1;
    for( entry = map->mapEntries; entry != NULL; entry = entry->next ) {
        switch( entry->type ) {
        case tpia_mapEntry_type_target :
            n += strlen( entry->schema ) + strlen( entry->path ) + strlen( entry->evaluation ) + strlen( entry->projectile ) + strlen( entry->targetName ) + nTarget;
            break;
        case tpia_mapEntry_type_path :
            n += strlen( entry->path ) + strlen( entry->projectile ) + nPath;
            break;
        default :
            smr_setMessageInfo( smr, &(map->smrUserInterface), __FILE__, __LINE__, tpia_map_status_UnknownType, "unknown type = %d", entry->type );
            return( NULL );
        }
    }

    if( ( s = (char *) xData_malloc2( smr, n, 0, "xml string" ) ) == NULL ) return( NULL );
    p = s;
    strcpy( p, start );
    while( *p ) p++;
    for( entry = map->mapEntries; entry != NULL; entry = entry->next ) {
        switch( entry->type ) {
        case tpia_mapEntry_type_target :
            sprintf( p, targetFormat, entry->schema, entry->evaluation, entry->projectile, entry->targetName, entry->path );
            break;
        case tpia_mapEntry_type_path :
            sprintf( p, pathFormat, entry->projectile, entry->path );
            break;
        }
        while( *p ) p++;
    }
    strcpy( p, end );
    return( s );
}
/*
************************************************************
*/
void tpia_map_simpleWrite( FILE *f, tpia_map *map ) { _tpia_map_simpleWrite2( f, map, 0 ); }
/*
************************************************************
*/
static void _tpia_map_simpleWrite2( FILE *f, tpia_map *map, int level ) {

    tpia_mapEntry *entry;
    char sLevel[] = "                ";
    int n = strlen( sLevel ) / 4;

    if( map->status != tpia_map_status_Ok ) {
        fprintf( f, "Bad map status = %d\n", map->status );
        return;
    }
    if( level < n ) sLevel[4 * level] = 0;
    fprintf( f, "%smap->path = %s\n", sLevel, map->path );
    fprintf( f, "%smap->mapFileName = %s\n", sLevel, map->mapFileName );
    for( entry = map->mapEntries; entry != NULL; entry = entry->next ) {
        switch( entry->type ) {
        case tpia_mapEntry_type_target :
            fprintf( f, "%sType = target: schema = %s: evaluation = %s: projectile = %s: target = %s: path = %s\n", sLevel, entry->schema, 
                entry->evaluation, entry->projectile, entry->targetName, entry->path );
            break;
        case tpia_mapEntry_type_path :
            fprintf( f, "%sType =   path: projectile = %s: path = %s\n", sLevel, entry->projectile, entry->path );
            _tpia_map_simpleWrite2( f, entry->map, level + 1 );
            break;
        default :
            fprintf( f, "%sUnknown type = %d\n", sLevel, entry->type );
        }
    }
}
/*
************************************************************
*/
static int _tpia_map_smrUserInterface( void *userData, char **str ) {

    tpia_map_smr *smrUserInterface = (tpia_map_smr *) userData;
    char fnl[] = "map file = ";
    int size = strlen( fnl ) + strlen( smrUserInterface->map->mapFileName ) + 1;

    if( str != NULL ) {
        //if( ( *str = xData_malloc2( NULL, size, 0, "mapFileName" ) ) == NULL ) return( -1 );
        if( ( *str = (char*) xData_malloc2( NULL, size, 0, "mapFileName" ) ) == NULL ) return( -1 );
        strcpy( *str, fnl );
        strcat( *str, smrUserInterface->map->mapFileName );
    }
    return( size );
}

#if defined __cplusplus
}
#endif
