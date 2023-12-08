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
//
/// \file PhysChemIO.cc
/// \brief Implementation of the PhysChemIO class

#include "PhysChemIO.hh"
#include "PhysSteppingAction.hh"
#include "PhysAnalysis.hh"
#include "G4Track.hh"
#include "G4NavigationHistory.hh"
#include "G4RunManager.hh"

#ifdef USE_MPI
#include "G4MPImanager.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysChemIO::PhysChemIO(PhysSteppingAction* stepAction) : G4VPhysChemIO(),
fSteppingAction(stepAction)
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysChemIO::CreateWaterMolecule(G4int electronicModif, G4int electronicLevel,
                                    G4double /*energy*/,
                                    const G4Track* theIncomingTrack)
{
    //L.T. Anh:  to correct electronicLevel in G4DNAChemistryManager, 
    //see in G4DNAChemistryManager::CreateWaterMolecule
    electronicLevel = 4 - electronicLevel; 
    //Rel pos
    G4ThreeVector relPos;
    auto touchable = theIncomingTrack->GetStep()->GetPreStepPoint()->GetTouchable();
    relPos = touchable->GetHistory()->GetTopTransform().TransformPoint(theIncomingTrack->GetPosition());
    
    // Get the flag of the current volume
    G4int volumeFlag =(G4int)fSteppingAction->SetupVolumeFlag(
        theIncomingTrack->GetStep()->GetPreStepPoint()->GetPhysicalVolume()->GetName());
    
    if(    volumeFlag == 161  // voxelStraight
        || volumeFlag == 162 // voxelRight
        || volumeFlag == 163 // voxelLeft
        || volumeFlag == 164 // voxelUp
        || volumeFlag == 165 // voxelDown
        || volumeFlag == 261  // voxelStraight2
        || volumeFlag == 262 // voxelRight2
        || volumeFlag == 263 // voxelLeft2
        || volumeFlag == 264 // voxelUp2
        || volumeFlag == 265) // voxelDown2
    {
        // Get the volume copy number
        G4int volumeCpNum = touchable->GetCopyNumber();
        //theIncomingTrack->GetStep()->GetPreStepPoint()->GetPhysicalVolume()->GetUserID();

        // Default flag values
        G4String motherVolumeName = "";
        G4int motherVolumeFlag = -1;
        G4int motherVolumeCpNum = -1;

        // Mother volume informations

        // Be sure there is a mother volume to ask for
        if(theIncomingTrack->GetStep()->GetPreStepPoint()->GetTouchable()->GetHistoryDepth() >0)
        {
            G4VPhysicalVolume* motherVol = theIncomingTrack->GetStep()->GetPreStepPoint()->GetTouchable()->GetVolume(1);

            // General infos
            motherVolumeName = motherVol->GetName();
            motherVolumeFlag = (G4int)fSteppingAction->SetupVolumeFlag(motherVolumeName);
            motherVolumeCpNum = motherVol->GetCopyNo();
        }
        G4int eventId = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
#ifdef USE_MPI
        auto g4MPI = G4MPImanager::GetManager();
        if (g4MPI->IsSlave()) { // update eventID only for slave, cause rank_master=0
            G4int rank = g4MPI->GetRank();
            eventId += g4MPI->GetEventsInMaster() + (rank-1)*g4MPI->GetEventsInSlave();
        }
#endif
        
        InfoForChemGeo aInfo;
        aInfo.fType = 1; // water=1
        aInfo.fState = G4double( electronicModif );
        aInfo.fElectronicLevel = G4double( electronicLevel );
        aInfo.fX = theIncomingTrack->GetPosition().x()/nm;
        aInfo.fY = theIncomingTrack->GetPosition().y()/nm;
        aInfo.fZ = theIncomingTrack->GetPosition().z()/nm;
        aInfo.fParentTrackID = G4double( theIncomingTrack->GetTrackID() );
        aInfo.fEventNumber = G4double(eventId);
        aInfo.fVolume = G4double( volumeFlag );
        aInfo.fVolumeCopyNumber = G4double( volumeCpNum);
        aInfo.fMotherVolume = G4double( motherVolumeFlag );
        aInfo.fMotherVolumeCopyNumber = G4double( motherVolumeCpNum );
        aInfo.fRelX = relPos.x()/nm;
        aInfo.fRelY = relPos.y()/nm;
        aInfo.fRelZ = relPos.z()/nm;
        PhysAnalysis::GetAnalysis()->AddInfoForChemGeo(aInfo);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysChemIO::CreateSolvatedElectron(const G4Track* theIncomingTrack, G4ThreeVector* finalPosition)
{
    G4ThreeVector pos;
    if(finalPosition) pos = *finalPosition;
    else pos = theIncomingTrack->GetPosition();

    // Rel pos
    G4ThreeVector relPos;
    const G4VTouchable* touchable = theIncomingTrack->GetStep()->GetPreStepPoint()->GetTouchable();
    relPos = touchable->GetHistory()->GetTopTransform().TransformPoint(pos);

    // Current volume infos
    G4int volumeFlag = (G4int)fSteppingAction->SetupVolumeFlag(
        theIncomingTrack->GetStep()->GetPreStepPoint()->GetPhysicalVolume()->GetName());
    if(    volumeFlag == 161  // voxelStraight
        || volumeFlag == 162 // voxelRight
        || volumeFlag == 163 // voxelLeft
        || volumeFlag == 164 // voxelUp
        || volumeFlag == 165 // voxelDown
        || volumeFlag == 261  // voxelStraight2
        || volumeFlag == 262 // voxelRight2
        || volumeFlag == 263 // voxelLeft2
        || volumeFlag == 264 // voxelUp2
        || volumeFlag == 265) // voxelDown2
    {
        G4int volumeCpNum = touchable->GetCopyNumber();

        G4String motherVolumeName = "";
        G4int motherVolumeFlag = -1;
        G4int motherVolumeCpNum = -1;

        // Mother volume informations
        // Be sure there is a mother volume to ask for
        if(theIncomingTrack->GetStep()->GetPreStepPoint()->GetTouchable()->GetHistoryDepth() >0)
        {
            G4VPhysicalVolume* motherVol = theIncomingTrack->GetStep()->GetPreStepPoint()->GetTouchable()->GetVolume(1);
            // General infos
            motherVolumeName = motherVol->GetName();
            motherVolumeFlag = (G4int)fSteppingAction->SetupVolumeFlag(motherVolumeName);
            motherVolumeCpNum = motherVol->GetCopyNo();
        }
        G4int eventId = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
#ifdef USE_MPI
        auto g4MPI = G4MPImanager::GetManager();
        if (g4MPI->IsSlave()) { // update eventID only for slave, cause rank_master=0
            G4int rank = g4MPI->GetRank();
            eventId += g4MPI->GetEventsInMaster() + (rank-1)*g4MPI->GetEventsInSlave();
        }
#endif
        
        InfoForChemGeo aInfo;
        aInfo.fType = 2; // / solvated electron=2
        aInfo.fState = -1; // no state for solvated electron
        aInfo.fElectronicLevel = -1; // no electronic level for solvated electron
        aInfo.fX = pos.x()/nm;
        aInfo.fY = pos.y()/nm;
        aInfo.fZ = pos.z()/nm;
        aInfo.fParentTrackID = G4double( theIncomingTrack->GetTrackID() );
        aInfo.fEventNumber = G4double(eventId);
        aInfo.fVolume = G4double( volumeFlag );
        aInfo.fVolumeCopyNumber = G4double( volumeCpNum );
        aInfo.fMotherVolume = G4double( motherVolumeFlag );
        aInfo.fMotherVolumeCopyNumber = G4double( motherVolumeCpNum );
        aInfo.fRelX = relPos.x()/nm;
        aInfo.fRelY = relPos.y()/nm;
        aInfo.fRelZ = relPos.z()/nm;
        PhysAnalysis::GetAnalysis()->AddInfoForChemGeo(aInfo);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......