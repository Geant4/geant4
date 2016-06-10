/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#include <string.h>
#include <PoPs.h>

#include "MCGIDI.h"

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

/*
************************************************************
*/
MCGIDI_POPs *MCGIDI_POPs_new( statusMessageReporting *smr, int size ) {

    MCGIDI_POPs *pops;

    if( ( pops = (MCGIDI_POPs *) smr_malloc2( smr, sizeof( MCGIDI_POPs ), 0, "pops->sorted" ) ) == NULL ) return( NULL );
    if( MCGIDI_POPs_initial( smr, pops, size ) ) smr_freeMemory( (void **) &pops );
    return( pops );
}
/*
************************************************************
*/
int MCGIDI_POPs_initial( statusMessageReporting * /*smr*/, MCGIDI_POPs *pops, int size ) {

    memset( pops, 0, sizeof( MCGIDI_POPs ) );
    if( size < 10 ) size = 10;
    pops->increment = size;
    
    return( 0 );
}
/*
************************************************************
*/
void *MCGIDI_POPs_free( MCGIDI_POPs *pops ) {

    if( pops == NULL ) return( NULL );
    MCGIDI_POPs_release( pops );
    smr_freeMemory( (void **) &pops );
    return( NULL );
}
/*
************************************************************
*/
int MCGIDI_POPs_release( MCGIDI_POPs *pops ) {

    MCGIDI_POP *pop, *next;

    if( pops == NULL ) return( 0 );
    for( pop = pops->first; pop != NULL; pop = next ) {
        next = pop->next;
        MCGIDI_POP_free( pop );
    }
    smr_freeMemory( (void **) &(pops->sorted) );
    MCGIDI_POPs_initial( NULL, pops, 0 );
    return( 0 );
}
/*
************************************************************
*/
MCGIDI_POP *MCGIDI_POPs_addParticleIfNeeded( statusMessageReporting *smr, MCGIDI_POPs *pops, char const *name, double mass_MeV, 
        double level_MeV, MCGIDI_POP *parent, int globalParticle ) {

    int i, index;
    MCGIDI_POP *pop;

    if( ( index = MCGIDI_POPs_findParticleIndex( pops, name ) ) >= 0 ) return( pops->sorted[index] );
    if( pops->size == pops->numberOfPOPs ) {
        int size = pops->size + pops->increment;
        MCGIDI_POP **sorted = (MCGIDI_POP **) smr_malloc2( smr, size * sizeof( MCGIDI_POP * ), 0, "sorted" );

        if( sorted == NULL ) return( NULL );
        for( i = 0; i < pops->numberOfPOPs; i++ ) sorted[i] = pops->sorted[i];
        smr_freeMemory( (void **) &(pops->sorted) );
        pops->sorted = sorted;
        pops->size = size;
    }
    index = -index - 1;
    if( ( pop = MCGIDI_POP_new( smr, name, mass_MeV, level_MeV, parent ) ) == NULL ) return( NULL );
    for( i = pops->numberOfPOPs; i > index; i-- ) pops->sorted[i] = pops->sorted[i-1];
    pops->sorted[index] = pop;
    if( pops->first == NULL ) {
        pops->first = pop; }
    else {
        pops->last->next = pop;
    }
    pops->last = pop;
    pops->numberOfPOPs++;
    pop->globalPoPsIndex = -1;
    if( globalParticle ) {
        if( ( pop->globalPoPsIndex = lPoPs_addParticleIfNeeded( smr, name, "LLNL" ) ) < 0 ) return( NULL );
    }
    return( pop );
}
/*
************************************************************
*/
int MCGIDI_POPs_findParticleIndex( MCGIDI_POPs *pops, char const *name ) {

    int iCmp = 0, min = 0, mid = 0, max = pops->numberOfPOPs;

    if( max == 0 ) return( -1 );
    while( ( max - min ) > 1 ) {
        mid = ( min + max ) / 2;
        iCmp = strcmp( name, pops->sorted[mid]->name );
        if( iCmp == 0 ) return( mid );
        if( iCmp < 0 ) {
            max = mid; }
        else {
            min = mid;
        }
    }  // Loop checking, 11.05.2015, T. Koi
    if( max == 1 ) {            /* First point is not checked as loop exits when ( max = 1 ) - ( min = 0 ) !> 1 ). */
        if( strcmp( name, pops->sorted[0]->name ) == 0 ) return( 0 );
    }
    if( max < pops->numberOfPOPs ) {
        if( strcmp( name, pops->sorted[max]->name ) == 0 ) return( max );
    }
    if( max == 1 ) {
        if( strcmp( name, pops->sorted[0]->name ) < 0 ) return( -1 );
    }
    return( -max - 1 );
}
/*
************************************************************
*/
MCGIDI_POP *MCGIDI_POPs_findParticle( MCGIDI_POPs *pops, char const *name ) {

    int index = MCGIDI_POPs_findParticleIndex( pops, name );

    if( index < 0 ) return( NULL );
    return( pops->sorted[index] );
}
/*
************************************************************
*/
void MCGIDI_POPs_writeSortedList( MCGIDI_POPs *pops, FILE *f ) {

    int i;

    fprintf( f, "POPs Information: n = %d\n", pops->numberOfPOPs );
    for( i = 0; i < pops->numberOfPOPs; i++ ) fprintf( f, "    %-20s  %e\n", pops->sorted[i]->name, pops->sorted[i]->mass_MeV );
}
/*
************************************************************
*/
void MCGIDI_POPs_printSortedList( MCGIDI_POPs *pops ) {

    MCGIDI_POPs_writeSortedList( pops, stdout );
}


/*
********* MCGIDI_POP routines *********
*/
/*
************************************************************
*/
MCGIDI_POP *MCGIDI_POP_new( statusMessageReporting *smr, char const *name, double mass_MeV, double level_MeV, MCGIDI_POP *parent ) {

    int Z, A, m, level;
    MCGIDI_POP *pop = (MCGIDI_POP *) smr_malloc2( smr, sizeof( MCGIDI_POP ), 0, "pop" );

    if( pop == NULL ) return( NULL );
    pop->next = NULL;
    pop->parent = parent;
    if( ( pop->name = smr_allocateCopyString2( smr, name, "pop->name" ) ) == NULL ) {
        smr_freeMemory( (void **) &pop );
        return( NULL );
    }
    MCGIDI_miscNameToZAm( smr, name, &Z, &A, &m, &level );
    pop->Z = Z;
    pop->A = A;
    pop->level = level;
    pop->m = m;
    pop->mass_MeV = mass_MeV;
    pop->level_MeV = level_MeV;
    pop->numberOfGammaBranchs = 0;
    pop->gammas = NULL;
    return( pop );
}
/*
************************************************************
*/
MCGIDI_POP *MCGIDI_POP_free( MCGIDI_POP *pop ) {

    if( pop == NULL ) return( NULL );
    MCGIDI_POP_release( pop );
    smr_freeMemory( (void **) &pop );
    return( NULL );
}
/*
************************************************************
*/
MCGIDI_POP *MCGIDI_POP_release( MCGIDI_POP *pop ) {

    if( pop == NULL ) return( NULL );
    smr_freeMemory( (void **) &(pop->name) );
    pop->numberOfGammaBranchs = 0;
    if( pop->gammas != NULL ) smr_freeMemory( (void **) &(pop->gammas) );
    return( NULL );
}
/*
************************************************************
*/
double MCGIDI_POP_getMass_MeV( MCGIDI_POP *pop ) {

    return( pop->mass_MeV );
}

#if defined __cplusplus
}
#endif

