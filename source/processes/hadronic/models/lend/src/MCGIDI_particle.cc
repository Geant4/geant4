/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/
#include <string.h>
#include "MCGIDI.h"

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

static int initialSizeOfList = 1000, incrementalSizeOfList = 1000;
static int numberOfParticles = 0, sizeOfParticleSortedList = 0;
static MCGIDI_particle **particleSortedList = NULL;
static MCGIDI_particle *particleList = NULL, *particleListEnd = NULL;
/*
************************************************************
*/
MCGIDI_particle *MCGIDI_particle_new( statusMessageReporting *smr ) {

    MCGIDI_particle *particle = (MCGIDI_particle *) smr_malloc2( smr, sizeof( MCGIDI_particle ), 0, "particle" );

    if( particle == NULL ) return( NULL );
    MCGIDI_particle_initialize( smr, particle );
    return( particle );
}
/*
************************************************************
*/
int MCGIDI_particle_initialize( statusMessageReporting * /*smr*/, MCGIDI_particle *particle ) {

    memset( particle, 0, sizeof( MCGIDI_particle ) );
    return( 0 );
}
/*
************************************************************
*/
MCGIDI_particle *MCGIDI_particle_free( statusMessageReporting *smr, MCGIDI_particle *particle ) {

    int i, j;
    MCGIDI_particle **p;

    for( i = 0, p = particleSortedList; i < numberOfParticles; i++, p++ ) {
        if( *p == particle ) {
            numberOfParticles--;
            for( j = i; j < numberOfParticles; j++, p++ ) *p = p[1];
            break;
        }
    }
    if( particle == particleListEnd ) particleListEnd = particle->prior;
    if( particle == particleList ) particleList = particle->next;
    if( particle->prior != NULL ) particle->prior->next = particle->next;
    if( particle->next != NULL ) particle->next->prior = particle->prior;
    MCGIDI_particle_release( smr, particle );
    smr_freeMemory( (void **) &particle );
    return( NULL );
}
/*
************************************************************
*/
int MCGIDI_particle_release( statusMessageReporting * /*smr*/, MCGIDI_particle *particle ) {

    smr_freeMemory( (void **) &(particle->name) );
    return( 0 );
}
/*
************************************************************
*/
int MCGIDI_particle_freeInternalList( statusMessageReporting *smr ) {

    while( particleList != NULL ) MCGIDI_particle_free( smr, particleList );  // Loop checking, 11.06.2015, T. Koi
    particleSortedList = (MCGIDI_particle **) smr_freeMemory( (void **) &particleSortedList );
    return( 0 );
}
/*
************************************************************
*/
MCGIDI_particle *MCGIDI_particle_getInternalID( statusMessageReporting *smr, const char * const name, MCGIDI_POPs *pops ) {

    int i, iCmp, min, mid, max, Z, A, m, level;
    MCGIDI_particle *particle;
    MCGIDI_POP *pop;

    iCmp = 0;
    min = mid = 0;
    max = numberOfParticles;
    while( min != max ) {  // Loop checking, 11.06.2015, T. Koi
        mid = ( min + max ) / 2;
        iCmp = strcmp( name, particleSortedList[mid]->name );
        if( iCmp == 0 ) return( particleSortedList[mid] );
        if( iCmp < 0 ) {
            max = mid - 1;
            if( mid == 0 ) max = 0; }
        else {
            min = mid + 1;
            if( min > max ) min = max;
        }
    }
    mid = min;
    if( numberOfParticles > 0 ) {
        iCmp = strcmp( name, particleSortedList[mid]->name );
        if( iCmp == 0 ) return( particleSortedList[mid] );
        if( ( iCmp < 0 ) && ( mid != 0 ) ) {
            mid--;
            iCmp = strcmp( name, particleSortedList[mid]->name );
        }
    }

    if( ( particle = MCGIDI_particle_new( smr ) ) == NULL ) return( NULL );
    if( ( particle->name = smr_allocateCopyString( smr, name, "particle->name", __FILE__, __LINE__, __func__ ) ) == NULL ) goto err;
    if( MCGIDI_miscNameToZAm( smr, name, &Z, &A, &m, &level ) != 0 ) goto err;
    particle->prior = NULL;
    particle->next = NULL;
    particle->Z = Z;
    particle->A = A;
    particle->m = m;
    if( ( pop = MCGIDI_POPs_findParticle( pops, name ) ) == NULL ) {    /* This should not happend. */
        particle->mass_MeV = MCGIDI_AMU2MeV * MCGIDI_particleMass_AMU( smr, name ); }
    else {
        particle->mass_MeV = pop->mass_MeV;
    }
    if( !smr_isOk( smr ) ) goto err;

    if( sizeOfParticleSortedList < ( numberOfParticles + 1 ) ) {
        if( sizeOfParticleSortedList == 0 ) {
            sizeOfParticleSortedList = initialSizeOfList; }
        else {
            sizeOfParticleSortedList += incrementalSizeOfList;
        }
        if( ( particleSortedList = (MCGIDI_particle **) smr_realloc2( smr, particleSortedList, sizeOfParticleSortedList * sizeof( MCGIDI_particle * ), 
            "particleSortedList" ) ) == NULL ) goto err;
    }

    if( particleList == NULL ) {
        particle->ordinal = 0;
        particleListEnd = particleList = particle; }
    else {
        particle->ordinal = particleListEnd->ordinal + 1;
        particle->prior = particleListEnd;
        particleListEnd->next = particle;
        particleListEnd = particle;
    }

    if( ( mid != 0 ) || ( iCmp > 0 ) ) mid++;
    for( i = numberOfParticles; i > mid; i-- ) particleSortedList[i] = particleSortedList[i-1];
    particleSortedList[mid] = particle;
    numberOfParticles++;

    return( particle );

err:
    MCGIDI_particle_free( smr, particle );
    return( NULL );
}
/*
************************************************************
*/
int MCGIDI_particle_printInternalSortedList( statusMessageReporting * /*smr*/ ) {

    int i;

    for( i = 0; i < numberOfParticles; i++ ) printf( "%s\n", particleSortedList[i]->name );
    return( 0 );
}

#if defined __cplusplus
}
#endif

