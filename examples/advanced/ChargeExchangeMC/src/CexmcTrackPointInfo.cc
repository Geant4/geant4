/*
 * ============================================================================
 *
 *       Filename:  CexmcTrackPointInfo.cc
 *
 *    Description:  single track point information
 *
 *        Version:  1.0
 *        Created:  28.12.2009 22:55:18
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * ============================================================================
 */

#include <iostream>
#include "CexmcTrackPointInfo.hh"


std::ostream &  operator<<( std::ostream &  out,
                            const CexmcTrackPointInfo &  trackPointInfo )
{
    if ( ! trackPointInfo.IsValid() )
        return out << "tp is not valid";

    const char *  trackTypeInfo = "???";

    switch ( trackPointInfo.trackType )
    {
    case CexmcBeamParticleTrack :
        trackTypeInfo = "bp";
        break;
    case CexmcOutputParticleTrack :
        trackTypeInfo = "op";
        break;
    case CexmcNucleusParticleTrack :
        trackTypeInfo = "np";
        break;
    case CexmcOutputParticleDecayProductTrack :
        trackTypeInfo = "opdp";
        break;
    default :
        break;
    }

    std::ostream::fmtflags  savedFlags( out.flags() );
    std::streamsize         prec( out.precision() );

    out.precision( 4 );

    out << trackPointInfo.particle->GetParticleName() << " [" <<
           trackPointInfo.trackId << "," << trackTypeInfo << "] " <<
           G4BestUnit( trackPointInfo.momentumAmp, "Energy" ) << " :  " <<
           G4BestUnit( trackPointInfo.positionLocal, "Length" ) << " [" <<
           trackPointInfo.directionLocal.x() << ", " <<
           trackPointInfo.directionLocal.y() << ", " <<
           trackPointInfo.directionLocal.z() << "]";
#ifdef CEXMC_DEBUG_TP
    out << std::endl << "                            < in world: " <<
           G4BestUnit( trackPointInfo.positionWorld, "Length" ) << " [" <<
           trackPointInfo.directionWorld.x() << ", " <<
           trackPointInfo.directionWorld.y() << ", " <<
           trackPointInfo.directionWorld.z() << "] >";
#endif

    out.precision( prec );
    out.flags( savedFlags );

    return out;
}

