/*
 * ============================================================================
 *
 *       Filename:  CexmcSimpleTrackPointInfoStore.cc
 *
 *    Description:  serialization helper for track point info objects
 *
 *        Version:  1.0
 *        Created:  31.12.2009 14:07:04
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * ============================================================================
 */

#ifdef CEXMC_USE_PERSISTENCY

#include <G4ParticleTable.hh>
#include "CexmcSimpleTrackPointInfoStore.hh"
#include "CexmcTrackPointInfo.hh"
#include "CexmcException.hh"


CexmcSimpleTrackPointInfoStore::CexmcSimpleTrackPointInfoStore()
{
}


CexmcSimpleTrackPointInfoStore::CexmcSimpleTrackPointInfoStore(
                                        const CexmcTrackPointInfo &  tpInfo )
{
    positionLocal = tpInfo.positionLocal;
    positionWorld = tpInfo.positionWorld;
    directionLocal = tpInfo.directionLocal;
    directionWorld = tpInfo.directionWorld;
    momentumAmp = tpInfo.momentumAmp;
    if ( tpInfo.trackId == CexmcInvalidTrackId )
        particlePDGEncoding = 0;
    else
        particlePDGEncoding = tpInfo.particle->GetPDGEncoding();
    trackId = tpInfo.trackId;
    trackType = tpInfo.trackType;
}


CexmcSimpleTrackPointInfoStore::operator CexmcTrackPointInfo() const
{
    G4ParticleDefinition *  particleDefinition(
                    G4ParticleTable::GetParticleTable()->FindParticle(
                                                        particlePDGEncoding ) );
    if ( ! particleDefinition && trackId != CexmcInvalidTrackId )
        throw CexmcException( CexmcWeirdException );

    return CexmcTrackPointInfo( positionLocal, positionWorld, directionLocal,
                        directionWorld, momentumAmp, particleDefinition,
                        trackId, trackType );
}

#endif

