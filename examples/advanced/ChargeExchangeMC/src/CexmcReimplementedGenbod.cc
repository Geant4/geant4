//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/*
 * =============================================================================
 *
 *       Filename:  CexmcReimplementedGenbod.cc
 *
 *    Description:  reimplemented GENBOD
 *                  (mostly adopted from ROOT TGenPhaseSpace)
 *
 *        Version:  1.0
 *        Created:  08.09.2010 18:52:39
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#include <cmath>
#include <Randomize.hh>
#include <G4PhysicalConstants.hh>
#include <G4SystemOfUnits.hh>
#include "CexmcReimplementedGenbod.hh"
#include "CexmcException.hh"


namespace
{
    G4int  DoubleMax( const void *  a, const void *  b ) 
    {
       G4double  aa( *( ( G4double * )a ) );
       G4double  bb( *( ( G4double * )b ) ); 

       if ( aa > bb )
           return  1;

       if ( aa < bb )
           return -1;

       return 0;
    }


    G4double  PDK( G4double  a, G4double  b, G4double  c )
    {
        G4double  x( ( a - b - c ) * ( a + b + c ) * ( a - b + c ) *
                     ( a + b - c ) );
        x = std::sqrt( x ) / ( 2 * a );

        return x;
    }
}


CexmcReimplementedGenbod::CexmcReimplementedGenbod() : maxWeight( 0. ),
    nmbOfOutputParticles( 0 )
{
}


G4double  CexmcReimplementedGenbod::Generate( void )
{
    // Generate a random final state.
    // The function returns the weigth of the current event.
    // Note that Momentum, Energy units are Gev/C, GeV

    G4double  te_minus_tm( totalEnergy - totalMass );
    G4double  rno[ maxParticles ];
    rno[ 0 ] = 0;

    if ( nmbOfOutputParticles > 2 )
    {
        for ( G4int  i( 1 ); i < nmbOfOutputParticles - 1; ++i )
        {
            rno[ i ] = G4UniformRand();
        }
        qsort( rno + 1, nmbOfOutputParticles - 2, sizeof( G4double ),
               DoubleMax );
    }
    rno[ nmbOfOutputParticles - 1 ] = 1;

    G4double  invMas[ maxParticles ];
    G4double  sum( 0 );

    for ( int  i( 0 ); i < nmbOfOutputParticles; ++i )
    {
        sum += outVec[ i ].mass / GeV;
        invMas[ i ] = rno[ i ] * te_minus_tm / GeV + sum;
    }

    //
    //-----> compute the weight of the current event
    //
    G4double  wt( maxWeight );
    G4double  pd[ maxParticles ];

    for ( int  i( 0 ); i < nmbOfOutputParticles - 1; ++i )
    {
        pd[ i ] = PDK( invMas[ i + 1 ], invMas[ i ],
                       outVec[ i + 1 ].mass / GeV );
        wt *= pd[ i ];
    }

    //
    //-----> complete specification of event (Raubold-Lynch method)
    //
    outVec[ 0 ].lVec->setPx( 0. );
    outVec[ 0 ].lVec->setPy( pd[ 0 ] );
    outVec[ 0 ].lVec->setPz( 0. );
    outVec[ 0 ].lVec->setE( std::sqrt( pd[ 0 ] * pd[ 0 ] +
                                       outVec[ 0 ].mass / GeV *
                                       outVec[ 0 ].mass / GeV ) );

    G4int  i( 1 );

    while ( true )
    {
        outVec[ i ].lVec->setPx( 0. );
        outVec[ i ].lVec->setPy( -pd[ i - 1 ] );
        outVec[ i ].lVec->setPz( 0. );
        outVec[ i ].lVec->setE( std::sqrt( pd[ i - 1 ] * pd[ i - 1 ] +
                                           outVec[ i ].mass / GeV *
                                           outVec[ i ].mass / GeV ) );

        G4double  cZ( 2 * G4UniformRand() - 1 );
        G4double  sZ( std::sqrt( 1 - cZ * cZ ) );
        G4double  angY( 2 * pi * G4UniformRand() );
        G4double  cY( std::cos( angY ) );
        G4double  sY( std::sin( angY ) );

        for ( int  j( 0 ); j <= i; ++j )
        {
            G4LorentzVector *  v( outVec[ j ].lVec );
            G4double           x( v->px() );
            G4double           y( v->py() );
            v->setPx( cZ * x - sZ * y );
            v->setPy( sZ * x + cZ * y );   // rotation around Z
            x = v->px();
            G4double           z( v->pz() );
            v->setPx( cY * x - sY * z );
            v->setPz( sY * x + cY * z );   // rotation around Y
        }

        if ( i == nmbOfOutputParticles - 1 )
            break;

        G4double  beta( pd[ i ] / std::sqrt( pd[ i ] * pd[ i ] +
                                             invMas[ i ] * invMas[ i ] ) );
        for ( int  j( 0 ); j <= i; ++j )
            outVec[ j ].lVec->boost( 0, beta, 0 );

        ++i;
    }

    for ( int  j( 0 ); j < nmbOfOutputParticles; ++j )
        *outVec[ j ].lVec *= GeV;

    //
    //---> return the weigth of event
    //
    return wt;
}


void  CexmcReimplementedGenbod::ParticleChangeHook( void )
{
    nmbOfOutputParticles = outVec.size();

    if ( nmbOfOutputParticles < 2 || nmbOfOutputParticles > maxParticles )
        throw CexmcException( CexmcKinematicsException );

    SetMaxWeight();
}


void  CexmcReimplementedGenbod::FermiEnergyDepStatusChangeHook( void )
{
    SetMaxWeight();
}


void  CexmcReimplementedGenbod::SetMaxWeight( void )
{
    G4double  te_minus_tm( totalEnergy - totalMass );

    if ( fermiEnergyDepIsOn )
    {
        // ffq[] = pi * (2*pi)**(N-2) / (N-2)!
        G4double  ffq[] = { 0 
                     ,3.141592, 19.73921, 62.01255, 129.8788, 204.0131
                     ,256.3704, 268.4705, 240.9780, 189.2637
                     ,132.1308,  83.0202,  47.4210,  24.8295
                     ,12.0006,   5.3858,   2.2560,   0.8859 };
        maxWeight =
                std::pow( te_minus_tm / GeV, nmbOfOutputParticles - 2 ) *
                ffq[ nmbOfOutputParticles - 1 ] / ( totalEnergy / GeV );
    }
    else
    {
        G4double  emmax( ( te_minus_tm + outVec[ 0 ].mass ) / GeV );
        G4double  emmin( 0. );
        G4double  wtmax( 1. );

        for ( G4int  i( 1 ); i < nmbOfOutputParticles; ++i )
        {
            emmin += outVec[ i - 1 ].mass / GeV;
            emmax += outVec[ i ].mass / GeV;
            wtmax *= PDK( emmax, emmin, outVec[ i ].mass / GeV );
        }
        maxWeight = 1 / wtmax;
    }
}

