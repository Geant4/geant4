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
    CexmcTrackPointsStore( const  CexmcTrackPointInfo &  monitorTP,
    const  CexmcTrackPointInfo &  targetTPBeamParticle,
    const  CexmcTrackPointInfo &  targetTPOutputParticle,
    const  CexmcTrackPointInfo &  targetTPNucleusParticle,
    const  CexmcTrackPointInfo &  targetTPOutputParticleDecayProductParticle1,
    const  CexmcTrackPointInfo &  targetTPOutputParticleDecayProductParticle2,
    const  CexmcTrackPointInfo &  vetoCounterTPLeft,
    const  CexmcTrackPointInfo &  vetoCounterTPRight,
    const  CexmcTrackPointInfo &  calorimeterTPLeft,
    const  CexmcTrackPointInfo &  calorimeterTPRight ) :
        monitorTP( monitorTP ), targetTPBeamParticle( targetTPBeamParticle ),
        targetTPOutputParticle( targetTPOutputParticle ),
        targetTPNucleusParticle( targetTPNucleusParticle ),
        targetTPOutputParticleDecayProductParticle1 (
                                targetTPOutputParticleDecayProductParticle1 ),
        targetTPOutputParticleDecayProductParticle2(
                                targetTPOutputParticleDecayProductParticle2 ),
        vetoCounterTPLeft( vetoCounterTPLeft ),
        vetoCounterTPRight( vetoCounterTPRight ),
        calorimeterTPLeft( calorimeterTPLeft ),
        calorimeterTPRight( calorimeterTPRight )
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

