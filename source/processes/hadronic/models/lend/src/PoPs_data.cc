/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/
#include <string.h>

#include "PoPs.h"
#include "PoPs_private.h"
#include "PoPs_data.h"

#ifdef POPS_BDFLS
#include <cbdfls.h>
#endif

#if defined __cplusplus
    extern "C" {
namespace GIDI {
using namespace GIDI;
#endif

static double PoPs_getBDFLS_mass( char const *name, PoP *pop, double mass );

#if defined __cplusplus
    }
    }
#endif
/*
========================================================================
*/

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

int PoPs_particleReadDatabase( statusMessageReporting *smr, char const * /*name*/ ) {

    int i1, n1 = sizeof( PoPDatas ) / sizeof( PoPDatas[0] );
    PoP *pop;
    char ZAName[32];

    for( i1 = 0; i1 < n1; ++i1 ) {
        if( ( pop = PoPs_particleCreateLoadInfo( smr, PoPDatas[i1].name ) ) == NULL ) return( 1 );
        if( PoPs_addParticleIfNeeded( smr, pop ) == pop ) {
            if( ( pop->genre == PoPs_genre_atom ) && ( pop->Z < 110 ) ) {
                snprintf( ZAName, sizeof ZAName, "%d%.3d", pop->Z, pop->A );
                if( lPoPs_addParticleIfNeeded( smr, ZAName, "LLNL" ) < 0 ) return( 1 );
            } }
        else {
            PoP_free( pop );
        }
        if( smr_isOk( smr ) == 0 ) return( 1 );
    }
    if( lPoPs_addParticleIfNeeded( smr, "gamma", "LLNL" ) < 0 ) return( 1 );
    if( lPoPs_addParticleIfNeeded( smr, "g", "LLNL" ) < 0 ) return( 1 );
    return( 0 );
}
/*
========================================================================
*/
PoP *PoPs_particleCreateLoadInfo( statusMessageReporting *smr, const char *name ) {

    PoP *pop;

    if( ( pop = PoP_new( smr ) ) != NULL ) {
        if( PoPs_particleLoadInfo( smr, name, pop ) != 0 ) pop = PoP_free( pop );
    }
    return( pop );
}
/*
========================================================================
*/
int PoPs_particleLoadInfo( statusMessageReporting *smr, const char *name, PoP *pop ) {

    int i, n = sizeof( PoPDatas ) / sizeof( PoPDatas[0] );

    if( ( pop->name = smr_allocateCopyString2( smr, name, "name" ) ) == NULL ) return( -1 );
    for( i = 0; i < n; i++ ) {
        if( strcmp( PoPDatas[i].name, name ) == 0 ) {
            pop->genre = PoPDatas[i].genre;
            pop->Z = PoPDatas[i].Z;
            pop->A = 0;
            if( PoPDatas[i].N >= 0 ) pop->A = pop->Z + PoPDatas[i].N;
            pop->l = PoPDatas[i].nuclearLevel;
            pop->mass = PoPs_getBDFLS_mass( name, pop, PoPDatas[i].mass );
            pop->massUnit = unitsDB_addUnitIfNeeded( smr, "amu" );
            return( 0 );
        }
    }
    smr_freeMemory( (void **) &(pop->name) );
    smr_setReportError2( smr, smr_unknownID, 1, "particle %s not in database", name );
    return( -1 );
}

static void *BDFLS_Data = NULL;

/*
========================================================================
*/
static double PoPs_getBDFLS_mass( char const * /*name*/, PoP * /*pop*/, double mass ) {

#ifdef POPS_BDFLS

    int ZA = 1000 * pop->Z + pop->A;
    double mass_ = -1;

    if( BDFLS_Data == NULL ) return( mass );
    if( ZA > 0 ) {
        mass_ = cbdflsGetMass( (cbdfls_file *) BDFLS_Data, ZA ); }
    else if( pop->genre == PoPs_genre_lepton ) {
        if( pop->name[0] == 'e' ) mass_ = cbdflsGetMass( (cbdfls_file *) BDFLS_Data, 8 );
    }
    if( mass_ < 0 ) mass_ = mass;
    mass = mass_;
#endif
    return( mass );
}
/*
========================================================================
*/
int PoPs_setBDFLS_File( char const *name ) {

#ifdef POPS_BDFLS

    cbdfls_file *p;
    cbdflsErrors Error;

    if( BDFLS_Data != NULL ) cbdflsRelease( (cbdfls_file *) BDFLS_Data );
    BDFLS_Data = NULL;
    if( name != NULL ) {
        if( ( p = cbdflsOpen( name, &Error ) ) == NULL ) return( 1 );
        BDFLS_Data = (void *) p;
    }
#else
    if( name == NULL ) BDFLS_Data = NULL;   /* Do something with name so compilers do not complain. */
#endif
    return( 0 );
}

#if defined __cplusplus
}
#endif
