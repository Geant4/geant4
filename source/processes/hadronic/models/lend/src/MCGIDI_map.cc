/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/
#include <string.h>
#include <limits.h>
#include <PoPs.h>
#include "MCGIDI_map.h"
#include <xDataTOM_importXML_private.h>

#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

#if defined(WIN32) || defined(__MINGW32__)
#include <windows.h>
#define realpath( a, b ) GetFullPathName( a, PATH_MAX, b, NULL )
#endif

static int aliasesNeeded = 1;

#if defined __cplusplus
    extern "C" {
    namespace GIDI {
    using namespace GIDI;
#endif
static MCGIDI_mapEntry *_MCGIDI_map_addEntry( statusMessageReporting *smr, MCGIDI_map *map, enum MCGIDI_mapEntry_type type, const char *schema, const char *path,
    const char *evaluation, const char *projectile, const char *target );
static char *_MCGIDI_map_findTargetViaPoPIDs2( statusMessageReporting *smr, MCGIDI_map *map, const char *evaluation, 
            int projectile_PoPID, int target_PoPID );
static int _MCGIDI_map_findAllOfTargetViaPoPIDs2( statusMessageReporting *smr, MCGIDI_map *mapAllOfTarget, MCGIDI_map *map, 
        int projectile_PoPID, int target_PoPID );
static int _MCGIDI_map_walkTree2( statusMessageReporting *smr, MCGIDI_map *map, int level, int (*handler)( MCGIDI_mapEntry *entry, int level, void *userData), 
    void *userData );
static void _MCGIDI_map_simpleWrite2( FILE *f, MCGIDI_map *map, int level );
static char *_MCGIDI_map_smrUserInterface( void *userData );
#if defined __cplusplus
    }
    }
#endif

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif
/*
************************************************************
*/
MCGIDI_map *MCGIDI_map_new( statusMessageReporting *smr ) {

    MCGIDI_map *map;
    
    if( ( map = (MCGIDI_map *) smr_malloc2( smr, sizeof( MCGIDI_map ), 0, "map" ) ) == NULL ) return( NULL );
    if( MCGIDI_map_initialize( smr, map ) ) map = (MCGIDI_map *) MCGIDI_map_free( NULL, map );
    return( map );
}
/*
************************************************************
*/
int MCGIDI_map_initialize( statusMessageReporting *smr, MCGIDI_map *map ) {

    memset( map, 0, sizeof( MCGIDI_map ) );
    map->status = MCGIDI_map_status_Ok;
    map->smrUserInterface.smrUserInterface = _MCGIDI_map_smrUserInterface;
    map->smrUserInterface.map = map;
    map->path = NULL;
    map->mapFileName = NULL;
    map->numberOfEntries = 0;
    map->mapEntries = NULL;

/*
*   Add some default aliases. This is a kludge until aliases are fully supported.
*/
if( aliasesNeeded ) { /* Support all meta-stables in ENDF/B-VII.1 */
    char const *aliases[] = { "Co58m1",  "Ag110m1",  "Cd115m1",  "Te127m1",  "Te129m1",  "Pm148m1",  "Ho166m1",  "Am242m1",  "Am244m1",  "Es254m1" };
    char const *targets[] = { "Co58_e1", "Ag110_e2", "Cd115_e1", "Te127_e2", "Te129_e1", "Pm148_e2", "Ho166_e1", "Am242_e2", "Am244_e1", "Es254_e2" };
    int i1, n1 = sizeof( aliases ) / sizeof( aliases[1] );
    

    for( i1 = 0; i1 < n1; i1++ ) {
        lPoPs_addParticleIfNeeded( smr, targets[i1], NULL );
        if( !smr_isOk( smr ) ) return( 1 );
        PoPs_addAliasIfNeeded( smr, targets[i1], aliases[i1] );
        if( !smr_isOk( smr ) ) return( 1 );
    }
    aliasesNeeded = 0;
}
    return( 0 );
}
/*
************************************************************
*/
MCGIDI_map *MCGIDI_map_readFile( statusMessageReporting *smr, const char *basePath, const char *mapFileName ) {
/*
*   If an error occurrs, map is freed and NULL is returned.
*/
    int n = 0;
    xDataXML_document *doc;
    xDataXML_element *element;
    xDataXML_element *child;
    MCGIDI_map *map;
    const char *evaluation, *projectile, *targetName, *path, *schema;
    char realPath[2 * ( PATH_MAX + 1 )], *p = &(realPath[PATH_MAX+1]);

    if( ( map = MCGIDI_map_new( smr ) ) == NULL ) return( NULL );

    if( ( basePath == NULL ) || ( mapFileName[0] == '/' ) ) {
        strcpy( realPath, mapFileName ); }
    else {
        strcpy( realPath, basePath );
        strcat( realPath, "/" );
        strcat( realPath, mapFileName );
    }
    if( realpath( realPath, p ) == NULL ) {
        smr_setReportError2( smr, smr_unknownID, MCGIDI_map_status_mapParsing, "No map file %s\n", mapFileName );
        return( (MCGIDI_map *) MCGIDI_map_free( NULL, map ) );
    }
    n = (int) strlen( p ) + 2;
    if( ( map->path = (char *) smr_malloc2( smr, 2 * n, 0, "map->path" ) ) == NULL ) return( (MCGIDI_map *) MCGIDI_map_free( NULL, map ) );
    map->mapFileName = &(map->path[n + 1]);
    strcpy( map->mapFileName, p );
    strcpy( map->path, p );
    if( ( p = strrchr( map->path, '/' ) ) != NULL ) {
        *p = 0; }
    else {
        strcpy( map->path, "." );
    }

    if( ( doc = xDataXML_importFile2( smr, map->mapFileName ) ) == NULL ) return( (MCGIDI_map *) MCGIDI_map_free( NULL, map ) );

    element = xDataXML_getDocumentsElement( doc );
    for( child = xDataXML_getFirstElement( element ); child != NULL; child = xDataXML_getNextElement( child ) ) {
        if( strcmp( child->name, "path" ) == 0 ) {
            if( ( path = xDataXML_getAttributesValueInElement( child , "path" ) ) == NULL ) {
                smr_setReportError3p( smr, &(map->smrUserInterface), smr_unknownID, MCGIDI_map_status_mapParsing, "path missing path attribute" );
                break;
            }
            MCGIDI_map_addPath( smr, map, path ); }
        else if( strcmp( child->name, "target" ) == 0 ) {
            if( ( schema = xDataXML_getAttributesValueInElement( child , "schema" ) ) == NULL ) {
                smr_setReportError3p( smr, &(map->smrUserInterface), smr_unknownID, MCGIDI_map_status_mapParsing, "target missing 'schema' attribute" );
                break;
            }
            if( ( path = xDataXML_getAttributesValueInElement( child , "path" ) ) == NULL ) {
                smr_setReportError3p( smr, &(map->smrUserInterface), smr_unknownID, MCGIDI_map_status_mapParsing, "target missing 'path' attribute" );
                break;
            }
            if( ( evaluation = xDataXML_getAttributesValueInElement( child , "evaluation" ) ) == NULL ) {
                smr_setReportError3p( smr, &(map->smrUserInterface), smr_unknownID, MCGIDI_map_status_mapParsing, "target missing 'evaluation' attribute" );
                break;
            }
            if( ( projectile = xDataXML_getAttributesValueInElement( child , "projectile" ) ) == NULL ) {
                smr_setReportError3p( smr, &(map->smrUserInterface), smr_unknownID, MCGIDI_map_status_mapParsing, "target missing 'projectile' attribute" );
                break;
            }
            if( ( targetName = xDataXML_getAttributesValueInElement( child , "target" ) ) == NULL ) {
                smr_setReportError3p( smr, &(map->smrUserInterface), smr_unknownID, MCGIDI_map_status_mapParsing, "target missing 'target' attribute" );
                break;
            }
            MCGIDI_map_addTarget( smr, map, schema, path, evaluation, projectile, targetName ); }
        else {
            smr_setReportError3( smr, &(map->smrUserInterface), smr_unknownID, MCGIDI_map_status_mapParsing, "invalid element = %s", child->name );
        }
        if( !smr_isOk( smr ) ) break;
    }
    xDataXML_freeDoc( smr, doc );
    if( !smr_isOk( smr ) ) map = (MCGIDI_map *) MCGIDI_map_free( NULL, map );
    return( map );
}
/*
************************************************************
*/
void *MCGIDI_map_free( statusMessageReporting *smr, MCGIDI_map *map ) {

    MCGIDI_map_release( smr, map );
    smr_freeMemory( (void **) &map );
    return( NULL );
}
/*
************************************************************
*/
void MCGIDI_map_release( statusMessageReporting *smr, MCGIDI_map *map ) {

    MCGIDI_mapEntry *entry, *next;

    if( map->path != NULL ) smr_freeMemory( (void **) &(map->path) );
    for( entry = map->mapEntries; entry != NULL; entry = next ) {
        next = entry->next;
        if( entry->schema != NULL ) smr_freeMemory( (void **) &(entry->schema) );
        if( entry->path != NULL ) smr_freeMemory( (void **) &(entry->path) );
        if( entry->evaluation != NULL ) smr_freeMemory( (void **) &(entry->evaluation) );
        if( entry->projectile != NULL ) smr_freeMemory( (void **) &(entry->projectile) );
        if( entry->targetName != NULL ) smr_freeMemory( (void **) &(entry->targetName) );
        if( entry->map != NULL ) MCGIDI_map_free( smr, entry->map );
        smr_freeMemory( (void **) &entry );
    }
    map->numberOfEntries = 0;
    map->mapEntries = NULL;
    map->status = MCGIDI_map_status_Ok;
}
/*
************************************************************
*/
MCGIDI_mapEntry *MCGIDI_map_getFirstEntry( MCGIDI_map *map ) {

    return( map->mapEntries );
}
/*
************************************************************
*/
MCGIDI_mapEntry *MCGIDI_map_getNextEntry( MCGIDI_mapEntry *entry ) {

    return( entry->next );
}
/*
************************************************************
*/
int MCGIDI_map_addTarget( statusMessageReporting *smr, MCGIDI_map *map, const char *schema, const char *path, const char *evaluation, const char *projectile, const char *target ) {

    return( _MCGIDI_map_addEntry( smr, map, MCGIDI_mapEntry_type_target, schema, path, evaluation, projectile, target ) != NULL );
}
/*
************************************************************
*/
int MCGIDI_map_addPath( statusMessageReporting *smr, MCGIDI_map *map, const char *path ) {

    MCGIDI_mapEntry *entry = _MCGIDI_map_addEntry( smr, map, MCGIDI_mapEntry_type_path, NULL, path, NULL, NULL, NULL );

    if( entry != NULL ) {
        if( ( entry->map = MCGIDI_map_readFile( smr, map->path, entry->path ) ) == NULL ) entry = NULL;
    }
    return( entry != NULL );
}
/*
************************************************************
*/
static MCGIDI_mapEntry *_MCGIDI_map_addEntry( statusMessageReporting *smr, MCGIDI_map *map, enum MCGIDI_mapEntry_type type, const char *schema, 
    const char *path, const char *evaluation, const char *projectile, const char *targetName ) {

    MCGIDI_mapEntry *p;
    MCGIDI_mapEntry *entry;

    if( ( entry = (MCGIDI_mapEntry * ) smr_malloc2( smr, sizeof( MCGIDI_mapEntry ), 1, "entry" ) ) == NULL ) return( NULL );
    entry->next = NULL;
    entry->type = type;
    entry->parent = map;
    entry->schema = NULL;
    entry->path = NULL;
    entry->evaluation = NULL;
    entry->projectile = NULL;
    entry->targetName = NULL;
    entry->globalPoPsIndexProjectile = entry->globalPoPsIndexTarget = -1;
    entry->map = NULL;

    if( path != NULL ) {
        if( ( entry->path = (char *) smr_malloc2( smr, strlen( path ) + 1, 0, "path" ) ) == NULL ) goto err;
        strcpy( entry->path, path );
    }

    if( evaluation != NULL ) {
        if( ( entry->evaluation = (char *) smr_malloc2( smr, strlen( evaluation ) + 1, 0, "evaluation" ) ) == NULL ) goto err;
        strcpy( entry->evaluation, evaluation );
    }

    if( projectile != NULL ) {
        if( ( entry->globalPoPsIndexProjectile = lPoPs_addParticleIfNeeded( smr, projectile, "LLNL" ) ) < 0 ) goto err;
        if( ( entry->projectile = (char *) smr_malloc2( smr, strlen( projectile ) + 1, 0, "projectile" ) ) == NULL ) goto err;
        strcpy( entry->projectile, projectile );
    }

    if( targetName != NULL ) {
        if( ( entry->globalPoPsIndexTarget = lPoPs_addParticleIfNeeded( smr, targetName, "LLNL" ) ) < 0 ) goto err;
        if( ( entry->targetName = (char *) smr_malloc2( smr, strlen( targetName ) + 1, 0, "target" ) ) == NULL ) goto err;
        strcpy( entry->targetName, targetName );
    }

    if( schema != NULL ) {
        if( ( entry->schema = (char *) smr_malloc2( smr, strlen( schema ) + 1, 0, "schema" ) ) == NULL ) goto err;
        strcpy( entry->schema, schema );
    }

    if( map->mapEntries == NULL ) {
        map->mapEntries = entry; }
    else {
        for( p = map->mapEntries; p->next != NULL; p = p->next );
        p->next = entry;
    }
    map->numberOfEntries++;
    return( entry );

err:
    smr_freeMemory( (void **) &(entry->path) );
    smr_freeMemory( (void **) &(entry->evaluation) );
    smr_freeMemory( (void **) &(entry->projectile) );
    smr_freeMemory( (void **) &(entry->targetName) );
    smr_freeMemory( (void **) &entry );
    return( NULL );
}
/*
************************************************************
*/
char *MCGIDI_map_findTargetViaPoPIDs( statusMessageReporting *smr, MCGIDI_map *map, const char *evaluation, 
    int projectile_PoPID, int target_PoPID ) {
/*
* Calling routine must free returned pointer.
*/
    char *path;
    char const *projectileName = PoPs_getName_atIndex( smr, projectile_PoPID );
    char const *targetName = PoPs_getName_atIndex( smr, target_PoPID );

    if( !smr_isOk( smr ) ) return( NULL );
    if( map->status != MCGIDI_map_status_Ok ) return( NULL );

    path = _MCGIDI_map_findTargetViaPoPIDs2( smr, map, evaluation, projectile_PoPID, target_PoPID );
    if( ( path == NULL ) && smr_isOk( smr ) ) {
        if( evaluation == NULL ) {
            smr_setReportInfo3( smr, &(map->smrUserInterface), smr_unknownID, 1, "target %s for projectile %s not found", 
                targetName, projectileName ); }
        else {
            smr_setReportInfo3( smr, &(map->smrUserInterface), smr_unknownID, 1, "target %s for projectile %s and evaluation %s not found", 
                targetName, projectileName, evaluation );
        }
    }
    return( path );
}
/*
************************************************************
*/
static char *_MCGIDI_map_findTargetViaPoPIDs2( statusMessageReporting *smr, MCGIDI_map *map, const char *evaluation, 
            int projectile_PoPID, int target_PoPID ) {

    MCGIDI_mapEntry *entry;
    char *path = NULL;
    int n, status;

    if( evaluation != NULL ) {
        if( strlen( evaluation ) == 0 ) evaluation = NULL;
    }

    for( entry = map->mapEntries; entry != NULL; entry = entry->next ) {
        switch( entry->type ) {
        case MCGIDI_mapEntry_type_target :
            if( ( projectile_PoPID == entry->globalPoPsIndexProjectile ) && ( target_PoPID == entry->globalPoPsIndexTarget ) ) {
                if( evaluation == NULL ) {
                    status = 1; }
                else {
                    status = strcmp( evaluation,  entry->evaluation ) == 0;
                }
                if( status ) {
                    n = (int) strlen( map->path ) + 1 + (int) strlen( entry->path ) + 1;
                    if( ( path = (char * ) smr_malloc2( smr, n, 0, "path" ) ) == NULL ) return( NULL );
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
        case MCGIDI_mapEntry_type_path :
            if( ( path = _MCGIDI_map_findTargetViaPoPIDs2( smr, entry->map, evaluation, projectile_PoPID, target_PoPID ) ) != NULL ) return( path );
            break;
        default :
            smr_setReportInfo3( smr, &(map->smrUserInterface), smr_unknownID, MCGIDI_map_status_UnknownType, "unknown type = %d", entry->type );
            return( NULL );
        }
    }
    return( NULL );
}
/*
************************************************************
*/
char *MCGIDI_map_findTarget( statusMessageReporting *smr, MCGIDI_map *map, const char *evaluation, const char *projectile, const char *targetName ) {

    int projectile_PoPID, target_PoPID;

    if( ( projectile_PoPID = lPoPs_addParticleIfNeeded( smr, projectile, "LLNL" ) ) < 0 ) return( NULL );
    if( ( target_PoPID     = lPoPs_addParticleIfNeeded( smr, targetName, "LLNL" ) ) < 0 ) return( NULL );
    return( MCGIDI_map_findTargetViaPoPIDs( smr, map, evaluation, projectile_PoPID, target_PoPID ) );
}
/*
************************************************************
*/
MCGIDI_map *MCGIDI_map_findAllOfTargetViaPoPIDs( statusMessageReporting *smr, MCGIDI_map *map, int projectile_PoPID, 
    int target_PoPID ) {
/*
* Calling routine must free returned pointer.
*/
    int status;
    MCGIDI_map *mapAllOfTarget;
    
    if( map->status != MCGIDI_map_status_Ok ) return( NULL );
    if( ( mapAllOfTarget = MCGIDI_map_new( smr ) ) == NULL ) return( NULL );
    status = _MCGIDI_map_findAllOfTargetViaPoPIDs2( smr, mapAllOfTarget, map, projectile_PoPID, target_PoPID );
    if( ( status != 0 ) ) mapAllOfTarget = (MCGIDI_map *) MCGIDI_map_free( smr, mapAllOfTarget );
    return( mapAllOfTarget );
}
/*
************************************************************
*/
static int _MCGIDI_map_findAllOfTargetViaPoPIDs2( statusMessageReporting *smr, MCGIDI_map *mapAllOfTarget, MCGIDI_map *map, 
        int projectile_PoPID, int target_PoPID ) {

    MCGIDI_mapEntry *entry;

    for( entry = map->mapEntries; entry != NULL; entry = entry->next ) {
        switch( entry->type ) {
        case MCGIDI_mapEntry_type_target :
            if( ( projectile_PoPID == entry->globalPoPsIndexProjectile ) && ( target_PoPID == entry->globalPoPsIndexTarget ) ) {
                if( _MCGIDI_map_addEntry( smr, mapAllOfTarget, entry->type, entry->schema, entry->path, entry->evaluation, entry->projectile, 
                    entry->targetName ) == NULL ) return( 1 );
            }
            break;
        case MCGIDI_mapEntry_type_path :
            if( _MCGIDI_map_findAllOfTargetViaPoPIDs2( smr, mapAllOfTarget, entry->map, projectile_PoPID, target_PoPID ) != 0 ) return( 1 );
            break;
        default :
            smr_setReportInfo3( smr, &(map->smrUserInterface), smr_unknownID, MCGIDI_map_status_UnknownType, "unknown type = %d", entry->type );
            return( 1 );
        }
    }
    return( 0 );
}
/*
************************************************************
*/
MCGIDI_map *MCGIDI_map_findAllOfTarget( statusMessageReporting *smr, MCGIDI_map *map, const char *projectile, const char *targetName ) {

    int projectile_PoPID, target_PoPID;
    
    if( ( projectile_PoPID = lPoPs_addParticleIfNeeded( smr, projectile, "LLNL" ) ) < 0 ) return( NULL );
    if( ( target_PoPID     = lPoPs_addParticleIfNeeded( smr, targetName, "LLNL" ) ) < 0 ) return( NULL );
    return( MCGIDI_map_findAllOfTargetViaPoPIDs( smr, map, projectile_PoPID, target_PoPID ) );
}
/*
************************************************************
*/
char *MCGIDI_map_getFullPath( statusMessageReporting *smr, MCGIDI_map *map, const char *endPath ) {

    char *path;

    if( endPath[0] == '/' ) {
        if( ( path = (char *) smr_malloc2( smr, strlen( endPath ) + 1, 0, "path" ) ) == NULL ) return( NULL );
        path[0] = 0; }
    else {
        if( ( path = (char *) smr_malloc2( smr, strlen( map->path ) + strlen( endPath ) + 2, 0, "path" ) ) == NULL ) return( NULL );
        strcpy( path, map->path );
        strcat( path, "/" );
    }
    strcat( path, endPath );
    return( path );
}
/*
************************************************************
*/
char *MCGIDI_map_getTargetsFullPath( statusMessageReporting *smr, MCGIDI_mapEntry *target ) {

    char *path = NULL;
    MCGIDI_map *map = target->parent;

    switch( target->type ) {
    case MCGIDI_mapEntry_type_target :
        path = MCGIDI_map_getFullPath( smr, map, target->path );
        break;
    case MCGIDI_mapEntry_type_path :
        smr_setReportInfo3p( smr, &(map->smrUserInterface), smr_unknownID, MCGIDI_map_status_UnknownType, "path type not allowed" );
        break;
    default :
        smr_setReportInfo3( smr, &(map->smrUserInterface), smr_unknownID, MCGIDI_map_status_UnknownType, "unknown type = %d", target->type );
        break;
    }
    return( path );
}
/*
************************************************************
*/
static int _MCGIDI_map_walkTree2( statusMessageReporting *smr, MCGIDI_map *map, int level, int (*handler)( MCGIDI_mapEntry *entry, int level, void *userData), 
    void *userData ) {
    
    MCGIDI_mapEntry *entry;

    for( entry = map->mapEntries; entry != NULL; entry = entry->next ) {
        if( handler( entry, level, userData ) != 0 ) return( 1 );
        if( entry->type == MCGIDI_mapEntry_type_path ) if( _MCGIDI_map_walkTree2( smr, entry->map, level + 1, handler, userData ) != 0 ) return( 1 );
    }
    return( 0 );
}
/*
************************************************************
*/
int MCGIDI_map_walkTree( statusMessageReporting *smr, MCGIDI_map *map, int (*handler)( MCGIDI_mapEntry *entry, int level, void *userData), void *userData ) {

    return( _MCGIDI_map_walkTree2( smr, map, 0, handler, userData ) );
}
/*
************************************************************
*/
char *MCGIDI_map_toXMLString( statusMessageReporting *smr, MCGIDI_map *map ) {

    MCGIDI_mapEntry *entry;
    char *s, *p;
    char targetFormat[] = "<target schema=\"%s\" evaluation=\"%s\" projectile=\"%s\" target=\"%s\" path=\"%s\"/>\n";
    char pathFormat[] = "<path projectile=\"%s\" path=\"%s\"/>\n";
    char start[] = "<map>\n";
    char end[] = "</map>";
    int n = 0, nStart = (int) strlen( start ), nEnd = (int) strlen( end );
    int nTarget = (int) strlen( targetFormat ) - 10, nPath = (int) strlen( pathFormat ) - 4;

    if( map->status != MCGIDI_map_status_Ok ) return( NULL );

    n = nStart + nEnd + 1;
    for( entry = map->mapEntries; entry != NULL; entry = entry->next ) {
        switch( entry->type ) {
        case MCGIDI_mapEntry_type_target :
            n += (int) ( strlen( entry->schema ) + strlen( entry->path ) + strlen( entry->evaluation ) + strlen( entry->projectile ) + strlen( entry->targetName ) + nTarget );
            break;
        case MCGIDI_mapEntry_type_path :
            n += (int ) strlen( entry->path ) + (int ) strlen( entry->projectile ) + nPath;
            break;
        default :
            smr_setReportInfo3( smr, &(map->smrUserInterface), smr_unknownID, MCGIDI_map_status_UnknownType, "unknown type = %d", entry->type );
            return( NULL );
        }
    }

    if( ( s = (char *) smr_malloc2( smr, n, 0, "xml string" ) ) == NULL ) return( NULL );
    p = s;
    strcpy( p, start );
    while( *p ) p++; // Loop checking, 11.06.2015, T. Koi
    for( entry = map->mapEntries; entry != NULL; entry = entry->next ) {
        switch( entry->type ) {
        case MCGIDI_mapEntry_type_target :
            snprintf( p, sizeof start, targetFormat, entry->schema, entry->evaluation, entry->projectile, entry->targetName, entry->path );
            break;
        case MCGIDI_mapEntry_type_path :
            snprintf( p, sizeof start, pathFormat, entry->projectile, entry->path );
            break;
        }
        while( *p ) p++; // Loop checking, 11.06.2015, T. Koi
    }
    strcpy( p, end );
    return( s );
}
/*
************************************************************
*/
void MCGIDI_map_simpleWrite( FILE *f, MCGIDI_map *map ) { _MCGIDI_map_simpleWrite2( f, map, 0 ); }
/*
************************************************************
*/
static void _MCGIDI_map_simpleWrite2( FILE *f, MCGIDI_map *map, int level ) {

    MCGIDI_mapEntry *entry;
    char sLevel[] = "                        ";
    int n = (int ) strlen( sLevel ) / 4;

    if( map->status != MCGIDI_map_status_Ok ) {
        fprintf( f, "Bad map status = %d\n", map->status );
        return;
    }
    if( level < n ) sLevel[4 * level] = 0;
    fprintf( f, "%smap->path = %s\n", sLevel, map->path );
    fprintf( f, "%smap->mapFileName = %s\n", sLevel, map->mapFileName );
    for( entry = map->mapEntries; entry != NULL; entry = entry->next ) {
        switch( entry->type ) {
        case MCGIDI_mapEntry_type_target :
            fprintf( f, "%sType = target: schema = %s: evaluation = %s: projectile = %s: target = %s: path = %s\n", sLevel, entry->schema, 
                entry->evaluation, entry->projectile, entry->targetName, entry->path );
            break;
        case MCGIDI_mapEntry_type_path :
            fprintf( f, "%sType =   path: path = %s\n", sLevel, entry->path );
            _MCGIDI_map_simpleWrite2( f, entry->map, level + 1 );
            break;
        default :
            fprintf( f, "%sUnknown type = %d\n", sLevel, entry->type );
        }
    }
}
/*
************************************************************
*/
static char *_MCGIDI_map_smrUserInterface( void *userData ) {

    MCGIDI_map_smr *smrUserInterface = (MCGIDI_map_smr *) userData;

    return( smr_allocateFormatMessage( "map file = %s", smrUserInterface->map->mapFileName ) );
}

#if defined __cplusplus
}
#endif

