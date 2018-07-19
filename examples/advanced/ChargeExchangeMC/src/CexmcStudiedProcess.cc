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
 *       Filename:  CexmcStudiedProcess.cc
 *
 *    Description:  studied process in the target
 *
 *        Version:  1.0
 *        Created:  14.09.2010 18:37:01
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#include <G4VParticleChange.hh>
#include "CexmcStudiedProcess.hh"
#include "CexmcPhysicsManager.hh"
#include "CexmcIncidentParticleTrackInfo.hh"
#include "CexmcCommon.hh"


CexmcStudiedProcess::CexmcStudiedProcess(
        CexmcPhysicsManager *  physicsManager_, G4ProcessType  processType ) :
    G4WrapperProcess( CexmcStudiedProcessFirstName, processType ),
    physicsManager( physicsManager_ )
{
}


G4double  CexmcStudiedProcess::PostStepGetPhysicalInteractionLength(
            const G4Track &  track, G4double , G4ForceCondition *  condition )
{
    *condition = NotForced;

    if ( ! physicsManager->IsStudiedProcessAllowed() )
        return CexmcDblMax;

    CexmcTrackInfo *  trackInfo( static_cast< CexmcTrackInfo * >(
                                                track.GetUserInformation() ) );

    if ( ! trackInfo ||
         trackInfo->GetTypeInfo() != CexmcIncidentParticleTrackType )
        return CexmcDblMax;

    CexmcIncidentParticleTrackInfo *  theTrackInfo(
                static_cast< CexmcIncidentParticleTrackInfo * >( trackInfo ) );

    if ( ! theTrackInfo->IsStudiedProcessActivated() )
        return CexmcDblMax;

    return theTrackInfo->GetFinalTrackLengthInTarget() -
            theTrackInfo->GetCurrentTrackLengthInTarget();
}


G4VParticleChange *  CexmcStudiedProcess::PostStepDoIt( const G4Track &  track,
                                                        const G4Step &  step )
{
    G4VParticleChange *  particleChange( pRegProcess->PostStepDoIt( track,
                                                                    step ) );

    if ( particleChange && particleChange->GetTrackStatus() == fStopAndKill )
        physicsManager->IncrementNumberOfTriggeredStudiedInteractions();

    return particleChange;
}

