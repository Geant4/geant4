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

