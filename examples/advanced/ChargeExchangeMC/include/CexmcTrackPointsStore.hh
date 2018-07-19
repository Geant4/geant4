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
 *       Filename:  CexmcTrackPointsStore.hh
 *
 *    Description:  store const references of track points of interest
 *
 *        Version:  1.0
 *        Created:  25.11.2009 13:41:43
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_TRACK_POINTS_STORE_HH
#define CEXMC_TRACK_POINTS_STORE_HH

#include <G4Allocator.hh>
#include "CexmcTrackPointInfo.hh"


struct  CexmcTrackPointsStore
{
    CexmcTrackPointsStore( const  CexmcTrackPointInfo &  monitorTP_,
    const  CexmcTrackPointInfo &  targetTPBeamParticle_,
    const  CexmcTrackPointInfo &  targetTPOutputParticle_,
    const  CexmcTrackPointInfo &  targetTPNucleusParticle_,
    const  CexmcTrackPointInfo &  targetTPOutputParticleDecayProductParticle1_,
    const  CexmcTrackPointInfo &  targetTPOutputParticleDecayProductParticle2_,
    const  CexmcTrackPointInfo &  vetoCounterTPLeft_,
    const  CexmcTrackPointInfo &  vetoCounterTPRight_,
    const  CexmcTrackPointInfo &  calorimeterTPLeft_,
    const  CexmcTrackPointInfo &  calorimeterTPRight_ ) :
        monitorTP( monitorTP_ ), targetTPBeamParticle( targetTPBeamParticle_ ),
        targetTPOutputParticle( targetTPOutputParticle_ ),
        targetTPNucleusParticle( targetTPNucleusParticle_ ),
        targetTPOutputParticleDecayProductParticle1(
                                targetTPOutputParticleDecayProductParticle1_ ),
        targetTPOutputParticleDecayProductParticle2(
                                targetTPOutputParticleDecayProductParticle2_ ),
        vetoCounterTPLeft( vetoCounterTPLeft_ ),
        vetoCounterTPRight( vetoCounterTPRight_ ),
        calorimeterTPLeft( calorimeterTPLeft_ ),
        calorimeterTPRight( calorimeterTPRight_ )
    {}

    void *  operator new( size_t  size );

    void    operator delete( void *  obj );

    const CexmcTrackPointInfo &  monitorTP;
    
    const CexmcTrackPointInfo &  targetTPBeamParticle;
    
    const CexmcTrackPointInfo &  targetTPOutputParticle;
    
    const CexmcTrackPointInfo &  targetTPNucleusParticle;
    
    const CexmcTrackPointInfo &  targetTPOutputParticleDecayProductParticle1;
    
    const CexmcTrackPointInfo &  targetTPOutputParticleDecayProductParticle2;

    const CexmcTrackPointInfo &  vetoCounterTPLeft;

    const CexmcTrackPointInfo &  vetoCounterTPRight;

    const CexmcTrackPointInfo &  calorimeterTPLeft;

    const CexmcTrackPointInfo &  calorimeterTPRight;
};


extern G4Allocator< CexmcTrackPointsStore >  trackPointsStoreAllocator;


inline void *  CexmcTrackPointsStore::operator new( size_t )
{
    return trackPointsStoreAllocator.MallocSingle();
}


inline void  CexmcTrackPointsStore::operator delete( void *  obj )
{
    trackPointsStoreAllocator.FreeSingle(
                        reinterpret_cast< CexmcTrackPointsStore * >( obj ) );
}


#endif

