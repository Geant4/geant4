/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#ifndef PoPs_private_h_included
#define PoPs_private_h_included

#if defined __cplusplus
    extern "C" {
    namespace GIDI {
#endif

typedef struct unitsDB_s unitsDB;
typedef struct PoPs_s PoPs;

struct unitsDB_s {
    int numberOfUnits;
    int allocated;
    char const **unsorted;
};

struct PoPs_s {
    int numberOfParticles;
    int allocated;
    PoP **pops;
    PoP **sorted;
};

int PoPs_releasePrivate( statusMessageReporting *smr );

char const *unitsDB_addUnitIfNeeded( statusMessageReporting *smr, char const *unit );
int unitsDB_index( statusMessageReporting *smr, char const *unit );
char const *unitsDB_stringFromIndex( statusMessageReporting *smr, int index );

#if defined __cplusplus
    }
    }
#endif

#endif          /* End of PoPs_private_h_included. */
