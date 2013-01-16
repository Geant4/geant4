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
#include "G4Types.hh"

#if defined __cplusplus
namespace GIDI {
using namespace GIDI;
#endif

static G4ThreadLocal int initialSizeOfList = 100, incrementalSizeOfList = 100;
static G4ThreadLocal int numberOfParticles = 0, sizeOfParticleSortedList = 0;
static G4ThreadLocal tpia_particle **particleSortedList = NULL;
static G4ThreadLocal tpia_particle *particleList = NULL, *particleListEnd = NULL;
/*
************************************************************
*/
tpia_particle *tpia_particle_create( statusMessageReporting *smr ) {

    //tpia_particle *particle = xData_malloc2( smr, sizeof( tpia_particle ), 0, "particle" );
    tpia_particle *particle = (tpia_particle*) xData_malloc2( smr, sizeof( tpia_particle ), 0, "particle" );

    if( particle == NULL ) return( NULL );
    tpia_particle_initialize( smr, particle );
    return( particle );
}
/*
************************************************************
*/
//int tpia_particle_initialize( statusMessageReporting *smr, tpia_particle *particle ) {
int tpia_particle_initialize( statusMessageReporting *, tpia_particle *particle ) {

    memset( particle, 0, sizeof( tpia_particle ) );
    return( 0 );
}
/*
************************************************************
*/
tpia_particle *tpia_particle_free( statusMessageReporting *smr, tpia_particle *particle ) {

    int i, j;
    tpia_particle **p;

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
    tpia_particle_release( smr, particle );
    xData_free( smr, particle );
    return( NULL );
}
/*
************************************************************
*/
int tpia_particle_release( statusMessageReporting *smr, tpia_particle *particle ) {

    if( particle->spectralID != NULL ) tpi_spectralID_free( smr, particle->spectralID );
    return( 0 );
}
/*
************************************************************
*/
int tpia_particle_freeInternalList( statusMessageReporting *smr ) {

    while( particleList != NULL ) tpia_particle_free( smr, particleList );
    //particleSortedList = xData_free( smr, particleSortedList );
    particleSortedList = (tpia_particle**) xData_free( smr, particleSortedList );
    return( 0 );
}
/*
************************************************************
*/
tpia_particle *tpia_particle_getInternalID( statusMessageReporting *smr, const char * const name ) {

    int i, iCmp, min, mid, max, Z, A, m;
    tpia_particle *particle;
    char *EOP;

    iCmp = 0;
    min = mid = 0;
    max = numberOfParticles;
    while( min != max ) {
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

    if( ( particle = tpia_particle_create( smr ) ) == NULL ) return( NULL );
    if( ( particle->spectralID = tpi_spectralID_parse( smr, name, &(EOP) ) ) == NULL ) return( tpia_particle_free( smr, particle ) );
    particle->name = particle->spectralID->name;
    if( tpia_miscNameToZAm( smr, particle->name, &Z, &A, &m ) != 0 ) return( tpia_particle_free( smr, particle ) );
    particle->prior = NULL;
    particle->next = NULL;
    particle->Z = Z;
    particle->A = A;
    particle->m = m;
    particle->mass = tpia_particleMass_AMU( smr, particle->name );
    if( !smr_isOk( smr ) ) return( tpia_particle_free( smr, particle ) );
    particle->fullMass_MeV = tpia_AMU2MeV * particle->mass + particle->spectralID->level;

    if( sizeOfParticleSortedList < ( numberOfParticles + 1 ) ) {
        if( sizeOfParticleSortedList == 0 ) {
            sizeOfParticleSortedList = initialSizeOfList; }
        else {
            sizeOfParticleSortedList += incrementalSizeOfList;
        }
        //if( ( particleSortedList = xData_realloc2( smr, particleSortedList, sizeOfParticleSortedList * sizeof( tpia_particle * ), 
        if( ( particleSortedList = (tpia_particle** ) xData_realloc2( smr, particleSortedList, sizeOfParticleSortedList * sizeof( tpia_particle * ), 
            "particleSortedList" ) ) == NULL ) return( tpia_particle_free( smr, particle ) );
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
}
/*
************************************************************
*/
//int tpia_particle_printInternalSortedList( statusMessageReporting *smr ) {
int tpia_particle_printInternalSortedList( statusMessageReporting * ) {

    int i;

    for( i = 0; i < numberOfParticles; i++ ) printf( "%s\n", particleSortedList[i]->name );
    return( 0 );
}

#if defined __cplusplus
}
#endif
