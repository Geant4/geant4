#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "PoPs.h"
#include "PoPs_private.h"

/*
    In PoPs_addParticleIfNeeded and unitsDB_addUnitIfNeeded, smr_malloc2 and not smr_realloc2 is used so that the current database is not
    lost if more memory cannot be allocated (not sure that this is needed, maybe should crash).
*/

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

#define incrementalSize 1000

#define MeV2eV 1e6
#define MeV2keV 1e3
#define AMU2MeV 931.494028
#define AMU2eV ( MeV2eV * 931.494028 )
#define K2MeV 8.6173856922566752e-11
#define K2eV ( MeV2eV * K2MeV )

typedef struct unitConversions_s unitConversions;

struct unitConversions_s {
    char const *_from;
    char const *_to;
    double ratio;
};

int PoPs_smr_ID = smr_unknownID;
static int referenceCount = 0;
static char versionStr[64] = "";

/*
*   For MPI the following need to be broadcasted.
*/
static unitsDB unitsRoot = { 0, 0, NULL };
static PoPs popsRoot = { 0, 0, NULL, NULL };
/*
*   End need to MPI broadcasted.
*/

static unitConversions conversions[] = { { "amu", "eV/c**2", AMU2eV }, { "amu", "MeV/c**2", AMU2MeV }, { "MeV/c**2", "eV/c**2", MeV2eV }, 
    { "MeV", "eV", MeV2eV }, { "MeV", "keV", MeV2keV }, { "K", "MeV", K2MeV }, { "K", "eV", K2eV } };

static char const *PoPs_genreStrings[] = { "invalid", "unknown", "alias", "photon", "lepton", "quark", "meson", "baryon", "nucleus", "atom" };

static int PoPs_particleProperIndex( int index );
static int PoPs_sortedParticleIndex( char const *name );
static int unitsDB_release( void );
/*
========================================================================
*/
const char *PoPs_version( void ) {

    if( versionStr[0] == 0 ) snprintf( versionStr, sizeof versionStr, "PoPs version %d.%d.%d", POPS_VERSION_MAJOR, POPS_VERSION_MINOR, POPS_VERSION_PATCHLEVEL );
    return( versionStr );
}
/*
========================================================================
*/
int PoPs_versionMajor( void ) { return( POPS_VERSION_MAJOR ); }
int PoPs_versionMinor( void ) { return( POPS_VERSION_MINOR ); }
int PoPs_versionPatchLevel( void ) { return( POPS_VERSION_PATCHLEVEL ); }
/*
========================================================================
*/
int PoPs_register( void ) {

    if( referenceCount < 0 ) return( -1 );
    return( ++referenceCount );
}
/*
========================================================================
*/
int PoPs_readDatabase( statusMessageReporting *smr, char const *fileName ) {

    return( PoPs_particleReadDatabase( smr, fileName ) );
}
/*
========================================================================
*/
int PoPs_release( statusMessageReporting *smr ) {

    referenceCount--;
    if( referenceCount != 0 ) return( referenceCount );
    PoPs_releasePrivate( smr );
    return( 0 );
}
/*
========================================================================
*/
int PoPs_releasePrivate( statusMessageReporting * /*smr*/ ) {

    int i;

    for( i = 0; i < popsRoot.numberOfParticles; i++ ) PoP_free( popsRoot.pops[i] );
    smr_freeMemory( (void **) &(popsRoot.pops) );
    popsRoot.sorted = NULL;
    popsRoot.numberOfParticles = 0;
    popsRoot.allocated = 0;
    unitsDB_release( );
    return( 0 );
}
/*
========================================================================
*/
PoP *PoPs_addParticleIfNeeded( statusMessageReporting *smr, PoP *pop ) {
/*
    If particle with name pop->name is already in popsRoot, returns the pointer to the existing particle.
    A NULL is returned if adding particle to popsRoot fails.
*/
    int i, index = PoPs_sortedParticleIndex( pop->name );

    if( index >= 0 ) return( popsRoot.pops[PoPs_particleProperIndex( popsRoot.sorted[index]->index )] );
    if( popsRoot.numberOfParticles == popsRoot.allocated ) {
        int size = popsRoot.allocated + incrementalSize;
        PoP **sorted, **pops = (PoP **) smr_malloc2( smr, 2 * size * sizeof( PoPs * ), 0, "pops" );

        if( pops == NULL ) return( NULL );
        sorted = &(pops[size]);
        for( i = 0; i < popsRoot.numberOfParticles; i++ ) {
            pops[i] = popsRoot.pops[i];
            sorted[i] = popsRoot.sorted[i];
        }
        smr_freeMemory( (void **) &(popsRoot.pops) );
        popsRoot.pops = pops;
        popsRoot.sorted = sorted;
        popsRoot.allocated = size;
    }
    popsRoot.pops[popsRoot.numberOfParticles] = pop;
    index = -index - 1;
    for( i = popsRoot.numberOfParticles; i > index; i-- ) popsRoot.sorted[i] = popsRoot.sorted[i-1];
    popsRoot.sorted[index] = pop;
    pop->index = popsRoot.numberOfParticles;
    popsRoot.numberOfParticles++;
    if( pop->genre == PoPs_genre_alias ) {  /* Add pop->index to end of list of particles aliased by pop->properIndex. */
        PoP *pop2;

        for( pop2 = popsRoot.pops[pop->properIndex]; pop2->aliasIndex >= 0; pop2 = popsRoot.pops[pop2->aliasIndex] ) ;
        pop2->aliasIndex = pop->index;
    }
    return( pop );
}
/*
========================================================================
*/
PoP *PoPs_copyAddParticleIfNeeded( statusMessageReporting *smr, PoP *pop ) {
/*
    If particle with name pop->name is already in popsRoot, return the address of the existing particle.
    If particle is not in popsRoot then copy particle to a new 'PoP *', add the copied PoP to popsRoot and return its address.
    A NULL is return if particle coping fails or adding particle to popsRoot fails.
*/

    int index = PoPs_particleIndex( pop->name );
    PoP *newPoP;

    if( index >= 0 ) return( popsRoot.pops[index] );

    if( ( newPoP = (PoP *) smr_malloc2( smr, sizeof( PoP ), 0, "newPoP" ) ) == NULL ) return( NULL );
    if( PoP_copyParticle( smr, newPoP, pop ) ) {
        smr_freeMemory( (void **) &newPoP );
        return( NULL );
    }
    if( PoPs_addParticleIfNeeded( smr, newPoP ) == NULL ) {
        PoP_free( newPoP );
        return( NULL );
    }
    return( newPoP );
}
/*
========================================================================
*/
PoP *PoPs_addAliasIfNeeded( statusMessageReporting *smr, char const *name, char const *alias ) {

    PoP *pop = PoP_makeAlias( smr, name, alias );

    if( pop != NULL ) {
        if( pop->index < 0 ) {
            if( PoPs_addParticleIfNeeded( smr, pop ) == NULL ) {
                PoP_free( pop );
                return( NULL );
            }
        }
    }
    
    return( pop );
}
/*
========================================================================
*/
int PoPs_numberOfParticle( void ) {

    return( popsRoot.numberOfParticles );
}
/*
========================================================================
*/
int PoPs_particleIndex( char const *name ) {
/*
    A negative number is return if particle is not in popsRoot. Else, the Id of the real (not aliased) particle is returned.
*/
    int index = PoPs_sortedParticleIndex( name );

    if( index >= 0 ) index = PoPs_particleProperIndex( popsRoot.sorted[index]->index );
    return( index );
}
/*
========================================================================
*/
int PoPs_particleIndex_smr( statusMessageReporting *smr, char const *name, char const *file, int line, char const *func ) {

    int index = PoPs_particleIndex( name );

    if( index < 0 )
        smr_setReportError( smr, NULL, file, line, func, PoPs_smr_ID, PoPs_errorToken_badName, "particle '%s' not in PoPs", name );
    return( index );
}
/*
========================================================================
*/
static int PoPs_particleProperIndex( int index ) {

    while( popsRoot.pops[index]->properIndex >= 0 ) index = popsRoot.pops[index]->properIndex;  /* For alias particles. */  // Loop checking, 11.05.2015, T. Koi
    return( index );
}
/*
========================================================================
*/
static int PoPs_sortedParticleIndex( char const *name ) {
/*
    If name is a particle in popsRoot, its index in the sorted list is returned; otherwise,
    a negative number is returned. For a particle not found, its index would be -returnValue + 1 if added;
*/
    int low = 0, mid, high = popsRoot.numberOfParticles, iCmp;

    if( high == 0 ) return( -1 );
    while( ( high - low ) > 1 ) {
        mid = ( low + high ) >> 1;
        iCmp = strcmp( name, popsRoot.sorted[mid]->name );
        if( iCmp == 0 ) return( mid );
        if( iCmp > 0 ) {
            low = mid; }
        else {
            high = mid;
        }
    }  // Loop checking, 11.05.2015, T. Koi
    if( high == 1 ) {           /* First point is not checked as loop exits when ( high = 1 ) - ( low = 0 ) <= 1 ). */
        if( !strcmp( name, popsRoot.sorted[0]->name ) ) return( 0 );            /* First name is a match. */
        if( strcmp( name, popsRoot.sorted[0]->name ) < 0 ) return( -1 );        /* name is less than first name. */
    }
    if( high < popsRoot.numberOfParticles ) {
        if( strcmp( name, popsRoot.sorted[high]->name ) == 0 ) return( high );
    }
    return( -high - 1 );
}
/*
========================================================================
*/
double PoPs_getMassInUnitOf( statusMessageReporting *smr, char const *name, char const *unit ) {

    int index = PoPs_particleIndex_smr( smr, name, __FILE__, __LINE__, __func__ );

    if( index < 0 ) return( -1. );
    return( PoPs_getMassInUnitOf_atIndex( smr, index, unit ) );
}
/*
========================================================================
*/
char const *PoPs_getName_atIndex( statusMessageReporting *smr, int index ) {

    if( ( index < 0 ) || ( index >= popsRoot.numberOfParticles ) ) {
        smr_setReportError2( smr, PoPs_smr_ID, PoPs_errorToken_badIndex, "index %d not in PoPs", index );
        return( NULL );
    }
    return( popsRoot.pops[index]->name );
}
/*
========================================================================
*/
double PoPs_getMassInUnitOf_atIndex( statusMessageReporting *smr, int index, char const *unit ) {

    double mass = -1.;

    if( ( index < 0 ) || ( index >= popsRoot.numberOfParticles ) ) {
        smr_setReportError2( smr, PoPs_smr_ID, PoPs_errorToken_badIndex, "index %d not in PoPs", index ); }
    else {
        mass = PoP_getMassInUnitOf( smr, popsRoot.pops[index], unit );
    }

    return( mass );
}
/*
========================================================================
*/
enum PoPs_genre PoPs_getGenre( statusMessageReporting *smr, char const *name ) {

    int index = PoPs_particleIndex_smr( smr, name, __FILE__, __LINE__, __func__ );

    if( index < 0 ) return( PoPs_genre_invalid );
    return( popsRoot.pops[index]->genre );
}
/*
========================================================================
*/
enum PoPs_genre PoPs_getGenre_atIndex( statusMessageReporting *smr, int index ) {

    enum PoPs_genre genre = PoPs_genre_invalid;

    if( ( index < 0 ) || ( index >= popsRoot.numberOfParticles ) ) {
        smr_setReportError2( smr, PoPs_smr_ID, PoPs_errorToken_badIndex, "index %d not in PoPs", index ); }
    else {
        genre = popsRoot.pops[index]->genre;
    }
    return( genre );
}
/*
========================================================================
*/
int PoPs_getZ_A_l( statusMessageReporting *smr, char const *name, int *Z, int *A, int *l ) {

    int index = PoPs_particleIndex_smr( smr, name, __FILE__, __LINE__, __func__ );

    if( index < 0 ) return( -1 );
    return( PoPs_getZ_A_l_atIndex( smr, index, Z, A, l ) );
}
/*
========================================================================
*/
int PoPs_getZ_A_l_atIndex( statusMessageReporting *smr, int index, int *Z, int *A, int *l ) {

    if( ( index < 0 ) || ( index >= popsRoot.numberOfParticles ) ) {
        smr_setReportError2( smr, PoPs_smr_ID, PoPs_errorToken_badIndex, "index %d not in PoPs", index );
        return( -1 );
    }
    *Z = popsRoot.pops[index]->Z;
    *A = popsRoot.pops[index]->A;
    *l = 0;
    return( 0 );
}
/*
========================================================================
*/
int PoPs_hasNucleus( statusMessageReporting *smr, char const *name, int protonIsNucleus ) {

    int index = PoPs_particleIndex_smr( smr, name, __FILE__, __LINE__, __func__ );

    if( index < 0 ) return( -1 );
    return( PoPs_hasNucleus_atIndex( smr, index, protonIsNucleus ) );
}
/*
========================================================================
*/
int PoPs_hasNucleus_atIndex( statusMessageReporting *smr, int index, int protonIsNucleus ) {
/*
*   If an error is encountered, a negative value is returned. A value greater than 0 means the particle
*   contains a nucleus (is an atom, ion or nucleus). Otherwise, a 0 is returned.
*/
    if( ( index < 0 ) || ( index >= popsRoot.numberOfParticles ) ) {
        smr_setReportError2( smr, PoPs_smr_ID, PoPs_errorToken_badIndex, "index %d not in PoPs", index );
        return( -1 );
    }
    if( ( popsRoot.pops[index]->genre == PoPs_genre_nucleus ) || ( popsRoot.pops[index]->genre == PoPs_genre_atom ) ) return( 1 );
    if( protonIsNucleus ) {
        if( strcmp( "p", popsRoot.pops[index]->name ) == 0 ) return( 1 );
    }
    return( 0 );
}
/*
========================================================================
*/
char const *PoPs_getAtomsName( statusMessageReporting *smr, char const *name ) {

    int index = PoPs_particleIndex_smr( smr, name, __FILE__, __LINE__, __func__ );

    if( index < 0 ) return( NULL );
    return( PoPs_getAtomsName_atIndex( smr, index ) );
}
/*
========================================================================
*/
char const *PoPs_getAtomsName_atIndex( statusMessageReporting *smr, int index ) {

    int atomIndex = PoPs_getAtomsIndex_atIndex( smr, index );

    if( atomIndex < 0 ) return( NULL );
    return( popsRoot.pops[atomIndex]->name );
}
/*
========================================================================
*/
int PoPs_getAtomsIndex( statusMessageReporting *smr, char const *name ) {

    int index = PoPs_particleIndex_smr( smr, name, __FILE__, __LINE__, __func__ );

    if( index < 0 ) return( index );
    return( PoPs_getAtomsIndex_atIndex( smr, index ) );
}
/*
========================================================================
*/
int PoPs_getAtomsIndex_atIndex( statusMessageReporting *smr, int index ) {

    char const *p = NULL;

    if( ( index < 0 ) || ( index >= popsRoot.numberOfParticles ) ) {
        smr_setReportError2( smr, PoPs_smr_ID, PoPs_errorToken_badIndex, "index %d not in PoPs", index );
        return( -1 );
    }

    if( popsRoot.pops[index]->genre == PoPs_genre_atom ) return( index );

    if(      strcmp( "p", popsRoot.pops[index]->name ) == 0 ) {
        p = "H1"; }
    else {
        if( popsRoot.pops[index]->genre != PoPs_genre_nucleus ) return( -1 );
        else if( strcmp( "h2", popsRoot.pops[index]->name ) == 0 ) {
            p = "H2"; }
        else if( strcmp( "h3", popsRoot.pops[index]->name ) == 0 ) {
            p = "H3"; }
        else if( strcmp( "he3", popsRoot.pops[index]->name ) == 0 ) {
            p = "He3"; }
        else if( strcmp( "he4", popsRoot.pops[index]->name ) == 0 ) {
            p = "He4";
        }
    }
    if( p != NULL ) return( PoPs_particleIndex_smr( smr, p, __FILE__, __LINE__, __func__ ) );
    return( -1 );
}
/*
========================================================================
*/
PoP *PoPs_getParticle_atIndex( int index ) {

    if( ( index < 0 ) || ( index >= popsRoot.numberOfParticles ) ) return( NULL );
    return( popsRoot.pops[index] );
}
/*
========================================================================
*/
char const *PoPs_genreTokenToString( enum PoPs_genre genre ) {

    if( genre < PoPs_genre_invalid ) return( NULL );
    if( genre > PoPs_genre_atom ) return( NULL );
    return( PoPs_genreStrings[genre] );
}
/*
========================================================================
*/
void PoPs_print( int sorted ) {

    PoPs_write( stdout, sorted );
}
/*
========================================================================
*/
void PoPs_write( FILE *f, int sorted ) {

    int i1, properIndex;
    PoP *pop;

    fprintf( f, "Mass units: number of units = %d\n", unitsRoot.numberOfUnits );
    for( i1 = 0; i1 < unitsRoot.numberOfUnits; i1++ ) {
        fprintf( f, " %s", unitsRoot.unsorted[i1] );
    }
    fprintf( f, "\n\n" );
    fprintf( f, "Particles: number of particles = %d\n", popsRoot.numberOfParticles );
    fprintf( f, " name                      index   genre            mass             hasNucleus    alias info\n" );
    fprintf( f, "                                                                           Z   A l\n" );
    fprintf( f, " --------------------------------------------------------------------------------------------\n" );
    for( i1 = 0; i1 < popsRoot.numberOfParticles; i1++ ) {
        if( sorted ) {
            pop = popsRoot.sorted[i1]; }
        else {
            pop = popsRoot.pops[i1];
        }
        properIndex = PoPs_particleProperIndex( pop->index );
        fprintf( f, " %-24s %6d   %-10s %15.8e %-6s", pop->name, pop->index, PoPs_genreTokenToString( pop->genre ), 
            popsRoot.pops[properIndex]->mass, popsRoot.pops[properIndex]->massUnit );
        if( PoPs_hasNucleus( NULL, pop->name, 0 ) ) {
            fprintf( f, " T" ); }
        else {
            fprintf( f, "  " );
        }
        if( PoPs_hasNucleus( NULL, pop->name, 1 ) ) {
            fprintf( f, " T" ); }
        else {
            fprintf( f, "  " );
        }
        if( pop->Z + pop->A > 0 ) {
            fprintf( f, " %3d %3d", pop->Z, pop->A );
            if( pop->l > 0 ) {
                fprintf( f, " %d", pop->l ); }
            else {
                fprintf( f, "  " );
            } }
        else {
            fprintf( f, "          " );
        }
        if( pop->genre == PoPs_genre_alias ) {
            fprintf( f, " %s (%d)", popsRoot.pops[properIndex]->name, popsRoot.pops[properIndex]->index ); }
        else {
            int aliasIndex;

            for( aliasIndex = pop->aliasIndex; aliasIndex >= 0; aliasIndex = popsRoot.pops[aliasIndex]->aliasIndex ) fprintf( f, " %d", aliasIndex );
        }
        fprintf( f, "\n" );
    }
}

/*
==========================   PoP functions    ==========================
*/
/*
========================================================================
*/
PoP *PoP_new( statusMessageReporting *smr ) {

    PoP *pop;

    if( ( pop = (PoP *) smr_malloc2( smr, sizeof( PoP ), 0, "pop" ) ) == NULL ) return( NULL );
    if( PoP_initialize( smr, pop ) != 0 ) pop = PoP_free( pop );
    return( pop );
}
/*
========================================================================
*/
int PoP_initialize( statusMessageReporting * /*smr*/, PoP *pop ) {

    pop->index = -1;
    pop->properIndex = -1;
    pop->aliasIndex = -1;
    pop->genre = PoPs_genre_unknown;
    pop->name = NULL;
    pop->Z = 0;
    pop->A = 0;
    pop->mass = 0.0;
    pop->massUnit = NULL;
    return( 0 );
}
/*
========================================================================
*/
int PoP_release( PoP *pop ) {

    if( pop->name != NULL ) smr_freeMemory( (void **) &(pop->name ) );
    PoP_initialize( NULL, pop );                                /* Make it clean in case someone trys to use if. */
    return( 0 );
}
/*
========================================================================
*/
PoP *PoP_free( PoP *pop ) {

    PoP_release( pop );
    smr_freeMemory( (void **) &pop );
    return( NULL );
}
/*
========================================================================
*/
int PoP_copyParticle( statusMessageReporting *smr, PoP *desc, PoP *src ) {

    desc->index = -1;
    desc->properIndex = src->properIndex;
    desc->aliasIndex = src->aliasIndex;
    desc->genre = src->genre;
    if( ( desc->name = smr_allocateCopyString2( smr, src->name, "desc->name" ) ) == NULL ) return( 1 );
    desc->Z = src->Z;
    desc->A = src->A;
    desc->l = src->l;
    desc->mass = src->mass;
    desc->massUnit = src->massUnit;

    return( 0 );
}
/*
========================================================================
*/
PoP *PoP_makeParticle( statusMessageReporting *smr, enum PoPs_genre genre, char const *name, double mass, char const *massUnit ) {

    PoP *pop;
    
    if( ( pop = PoP_new( smr ) ) == NULL ) return( NULL );
    if( ( pop->name = smr_allocateCopyString2( smr, name, "name" ) ) == NULL ) {
        PoP_free( pop );
        return( NULL );
    }
    pop->genre = genre;
    pop->mass = mass;
    if( ( pop->massUnit = unitsDB_addUnitIfNeeded( smr, massUnit ) ) == NULL ) pop = PoP_free( pop );
    return( pop );
}
/*
========================================================================
*/
int PoP_setZ_A_l( statusMessageReporting * /*smr*/, PoP *pop, int Z, int A, int l ) {

    pop->Z = Z;
    pop->A = A;
    pop->l = l;
    return( 0 );
}
/*
========================================================================
*/
int PoP_getIndex( PoP *pop ) {

    return( pop->index );
}
/*
========================================================================
*/
char const *PoP_getName( PoP *pop ) {

    return( pop->name );
}
/*
========================================================================
*/
double PoP_getMassInUnitOf( statusMessageReporting *smr, PoP *pop, char const *unit ) {

    double mass = -1., ratio;
    /* PoP *pop_ = pop;*/

    /*if( pop->genre == PoPs_genre_alias ) pop_ = popsRoot.pops[PoPs_particleProperIndex( pop->index )];*/
    if( PoPs_unitConversionRatio( pop->massUnit, unit, &ratio ) != 0 ) {
        smr_setReportError2( smr, PoPs_smr_ID, PoPs_errorToken_badUnitConversion, "could not convert unit '%s' to '%s'", pop->massUnit, unit ); }
    else {
        mass = pop->mass * ratio;
    }

    return( mass );
}

/*
==========================  alias functions   ==========================
*/
/*
========================================================================
*/
PoP *PoP_makeAlias( statusMessageReporting *smr, char const *name, char const *alias ) {

    int properIndex = PoPs_particleIndex( name ), aliasIndex = PoPs_particleIndex( alias );
    PoP *pop;

    if( properIndex < 0 ) {
        smr_setReportError2( smr, PoPs_smr_ID, PoPs_errorToken_badName, "proper particle '%s' not in PoPs for alias '%s'", name, alias );
        return( NULL );
    }
    if( aliasIndex >= 0 ) {     /* alias has already been defined. */
        PoP *truePop = popsRoot.pops[aliasIndex];

        for( pop = truePop; strcmp( alias, pop->name ); pop = popsRoot.pops[aliasIndex] ) aliasIndex = pop->aliasIndex;
        if( pop->genre != PoPs_genre_alias ) {
            smr_setReportError2( smr, PoPs_smr_ID, PoPs_errorToken_badName, "particle '%s' already in PoPs and not an alias", alias );
            return( NULL );
        }
        if( pop->properIndex != properIndex ) {
            smr_setReportError2( smr, PoPs_smr_ID, PoPs_errorToken_badName, "particle '%s' already an alias for '%s', cannot re-alias to '%s'", 
                alias, truePop->name, name );
            return( NULL );
        } }
    else {
        if( ( pop = PoP_new( smr ) ) == NULL ) return( NULL );
        if( ( pop->name = smr_allocateCopyString2( smr, alias, "name" ) ) == NULL ) {
            PoP_free( pop );
            return( NULL );
        }
        pop->properIndex = properIndex;
        pop->genre = PoPs_genre_alias;
    }
    return( pop );
}

/*
==========================  unitsDB functions  =========================
*/
/*
========================================================================
*/
static int unitsDB_release( void ) {

    int i;

    for( i = 0; i < unitsRoot.numberOfUnits; i++ ) smr_freeMemory( (void **) &(unitsRoot.unsorted[i]) );
    smr_freeMemory( (void **) &(unitsRoot.unsorted) );
    unitsRoot.numberOfUnits = 0;
    unitsRoot.allocated = 0;
    return( 0 );
}
/*
========================================================================
*/
char const *unitsDB_addUnitIfNeeded( statusMessageReporting *smr, char const *unit ) {

    int i;

    for( i = 0; i < unitsRoot.numberOfUnits; i++ ) {
        if( strcmp(  unit, unitsRoot.unsorted[i] ) == 0 ) return( unitsRoot.unsorted[i] );
    }
    if( unitsRoot.numberOfUnits == unitsRoot.allocated ) {
        int size = unitsRoot.allocated + 20;
        char const **unsorted = (char const **) smr_malloc2( smr, size * sizeof( char * ), 0, "unsorted" );

        if( unsorted == NULL ) return( NULL );
        for( i = 0; i < unitsRoot.numberOfUnits; i++ ) unsorted[i] = unitsRoot.unsorted[i];
        smr_freeMemory( (void **) &(unitsRoot.unsorted) );
        unitsRoot.unsorted = unsorted;
        unitsRoot.allocated = size;
    }
    if( ( unitsRoot.unsorted[unitsRoot.numberOfUnits] = smr_allocateCopyString2( smr, unit, "unitsRoot.unsorted[unitsRoot.numberOfUnits]" ) ) == NULL ) 
        return( NULL );
    unitsRoot.numberOfUnits++;
    return( unitsRoot.unsorted[unitsRoot.numberOfUnits - 1] );
}
/*
========================================================================
*/
int unitsDB_index( statusMessageReporting * /*smr*/, char const *unit ) {

    int i;

    for( i = 0; i < unitsRoot.numberOfUnits; i++ ) {
        if( !strcmp( unit, unitsRoot.unsorted[i] ) ) return( i );
    }
    return( -1 );
}
/*
========================================================================
*/
char const *unitsDB_stringFromIndex( statusMessageReporting *smr, int index ) {

    if( ( index < 0 ) || ( index >= unitsRoot.numberOfUnits ) ) {
        smr_setReportError2( smr, PoPs_smr_ID, 1, "index = %d out of baounds [0 to %d)", index, unitsRoot.numberOfUnits );
        return( NULL );
    }
    return( unitsRoot.unsorted[index] );
}
/*
========================================================================
*/
int PoPs_unitConversionRatio( char const *_from, char const *_to, double *ratio ) {

    int i, n = sizeof( conversions ) / sizeof( conversions[0] );

    *ratio = 1.;
    if( strcmp( _from, _to ) == 0 ) return( 0 );
    for( i = 0; i < n; i++ ) {
        if( strcmp( conversions[i]._from, _from ) == 0 ) {
            if( strcmp( conversions[i]._to, _to ) == 0 ) {
                *ratio = conversions[i].ratio;
                return( 0 );
            } }
        else if( strcmp( conversions[i]._to, _from ) == 0 ) {
            if( strcmp( conversions[i]._from, _to ) == 0 ) {
                *ratio = 1. / conversions[i].ratio;
                return( 0 );
            }
        }
    }
    return( 1 );
}
#ifdef PoPs_MPI
#include "PoPs_Bcast_private.h"
/*
========================================================================
*/
int PoPs_Bcast( statusMessageReporting *smr, MPI_Comm comm, int rank ) {

    return( PoPs_Bcast2( smr, comm, rank, &unitsRoot, &popsRoot ) );
}
#endif

#if defined __cplusplus
}
#endif
