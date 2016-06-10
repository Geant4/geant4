/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/
#ifndef MCGIDI_map_h_included
#define MCGIDI_map_h_included

#include <statusMessageReporting.h>

#if defined __cplusplus
    extern "C" {
    namespace GIDI {
#endif

enum MCGIDI_map_status { MCGIDI_map_status_Ok, MCGIDI_map_status_memory, MCGIDI_map_status_mapParsing,
    MCGIDI_map_status_UnknownType };
enum MCGIDI_mapEntry_type { MCGIDI_mapEntry_type_target, MCGIDI_mapEntry_type_path };

typedef struct MCGIDI_map_s MCGIDI_map;
typedef struct MCGIDI_mapEntry_s MCGIDI_mapEntry;
typedef struct MCGIDI_map_smr_s MCGIDI_map_smr;

struct MCGIDI_map_smr_s {
    smr_userInterface smrUserInterface;
    MCGIDI_map *map;
};

struct MCGIDI_mapEntry_s {
    MCGIDI_mapEntry *next;
    enum MCGIDI_mapEntry_type type;
    MCGIDI_map *parent;
    char *schema;
    char *path;
    char *evaluation;
    char *projectile;
    char *targetName;
    int globalPoPsIndexProjectile, globalPoPsIndexTarget;
    MCGIDI_map *map;
};

struct MCGIDI_map_s {
    enum MCGIDI_map_status status;
    MCGIDI_map_smr smrUserInterface;
    char *path;
    char *mapFileName;
    int numberOfEntries;
    MCGIDI_mapEntry *mapEntries;
};

MCGIDI_map *MCGIDI_map_new( statusMessageReporting *smr );
int MCGIDI_map_initialize( statusMessageReporting *smr, MCGIDI_map *map );
MCGIDI_map *MCGIDI_map_readFile( statusMessageReporting *smr, const char *basePath, const char *mapFileName );
void *MCGIDI_map_free( statusMessageReporting *smr, MCGIDI_map *map );
void MCGIDI_map_release( statusMessageReporting *smr, MCGIDI_map *map );
MCGIDI_mapEntry *MCGIDI_map_getFirstEntry( MCGIDI_map *map );
MCGIDI_mapEntry *MCGIDI_map_getNextEntry( MCGIDI_mapEntry *entry );
int MCGIDI_map_addTarget( statusMessageReporting *smr, MCGIDI_map *map, const char *method, const char *path, const char *evaluation, const char *projectile, const char *targetName );
int MCGIDI_map_addPath( statusMessageReporting *smr, MCGIDI_map *map, const char *path );
char *MCGIDI_map_findTargetViaPoPIDs( statusMessageReporting *smr, MCGIDI_map *map, const char *evaluation, int projectile_PoPID, int target_PoPID );
char *MCGIDI_map_findTarget( statusMessageReporting *smr, MCGIDI_map *map, const char *evaluation, const char *projectile, const char *targetName );
MCGIDI_map *MCGIDI_map_findAllOfTargetViaPoPIDs( statusMessageReporting *smr, MCGIDI_map *map, int projectile_PoPID, int target_PoPID );
MCGIDI_map *MCGIDI_map_findAllOfTarget( statusMessageReporting *smr, MCGIDI_map *map, const char *projectile, const char *targetName );
char *MCGIDI_map_getFullPath( statusMessageReporting *smr, MCGIDI_map *map, const char *endPath );
char *MCGIDI_map_getTargetsFullPath( statusMessageReporting *smr, MCGIDI_mapEntry *target );
int MCGIDI_map_walkTree( statusMessageReporting *smr, MCGIDI_map *map, int (*handler)( MCGIDI_mapEntry *entry, int level, void *userData), void *userData );
char *MCGIDI_map_toXMLString( statusMessageReporting *smr, MCGIDI_map *map );
void MCGIDI_map_simpleWrite( FILE *f, MCGIDI_map *map );

#if defined __cplusplus
    }
    }
#endif

#endif          /* End of MCGIDI_map_h_included. */
