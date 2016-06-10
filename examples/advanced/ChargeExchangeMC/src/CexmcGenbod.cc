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
 *       Filename:  CexmcGenbod.cc
 *
 *    Description:  original fortran routine GENBOD wrapper
 *
 *        Version:  1.0
 *        Created:  08.09.2010 14:36:23
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifdef CEXMC_USE_GENBOD

#include <G4SystemOfUnits.hh>
#include "CexmcGenbod.hh"
#include "CexmcException.hh"


extern "C"
{
    extern int  genbod_( void );
}


extern struct  genbod_in_data
{
    int    np;
    float  tecm;
    float  amass[ 18 ];
    int    kgenev;
}  genin_;


extern struct  genbod_out_data
{
    float  pcm[ 18 ][ 5 ];
    float  wt;
}  genout_;


CexmcGenbod::CexmcGenbod()
{
    genin_.np = 0;
    genin_.kgenev = 1;
}


G4bool  CexmcGenbod::CheckKinematics( void )
{
    totalEnergy = 0.;
    for ( CexmcPhaseSpaceInVector::const_iterator  k( inVec.begin() );
                                                        k != inVec.end(); ++k )
    {
        totalEnergy += ( *k )->e();
    }

    /* epsilon is needed to compensate float to fortran real cast accuracy,
     * the value 3E-6 was found experimentally, maybe has to be made bigger,
     * but not smaller */
    const float  epsilon( 3E-6 );

    return totalEnergy - totalMass > 0.0f + epsilon;
}


G4double  CexmcGenbod::Generate( void )
{
    genin_.tecm = totalEnergy / GeV;

    genbod_();

    float ( *pcm )[ 5 ]( &genout_.pcm[ 0 ] );

    for ( CexmcPhaseSpaceOutVector::iterator  k( outVec.begin() );
                                                        k != outVec.end(); ++k )
    {
        k->lVec->setPx( ( *pcm )[ 0 ] * GeV );
        k->lVec->setPy( ( *pcm )[ 1 ] * GeV );
        k->lVec->setPz( ( *pcm )[ 2 ] * GeV );
        k->lVec->setE( ( *pcm++ )[ 3 ] * GeV );
    }

    return genout_.wt;
}


void  CexmcGenbod::ParticleChangeHook( void )
{
    size_t  nmbOfOutputParticles( outVec.size() );

    if ( nmbOfOutputParticles < 2 || nmbOfOutputParticles > 18 )
        throw CexmcException( CexmcKinematicsException );

    genin_.np = nmbOfOutputParticles;

    float *  amass( genin_.amass );

    for ( CexmcPhaseSpaceOutVector::const_iterator  k( outVec.begin() );
                                                        k != outVec.end(); ++k )
    {
        *amass++ = k->mass / GeV;
    }
}


void  CexmcGenbod::FermiEnergyDepStatusChangeHook( void )
{
    genin_.kgenev = fermiEnergyDepIsOn ? 2 : 1;
}

#endif

