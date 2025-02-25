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
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class

#include "SteppingAction.hh"
#include "Analysis.hh"
#include "G4SteppingManager.hh"
#include "G4VTouchable.hh"
#include "G4VPhysicalVolume.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SteppingManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4RunManager.hh"
#include "G4Proton.hh"
#include "G4Electron.hh"
#include "G4Alpha.hh"
#include "G4DNAGenericIonsManager.hh"

#ifdef USE_MPI 
#include "G4MPImanager.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SteppingAction::SteppingAction(EventAction* pEvent)
{
    fEventAction = pEvent;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SteppingAction::UserSteppingAction(const G4Step* step)
{
    SetupFlags(step);
    
    SetupVoxelCopyNumber(step);

    if(fFlagVolume>0)
    {
        fEventAction->AddEdep(step->GetTotalEnergyDeposit()/eV);
    }
    //G4cout<<step->GetPreStepPoint()->GetTouchable()->GetVolume()->GetName()<<" ; "<<step->GetPreStepPoint()->GetTouchable()->GetVolume()->GetName()<<G4endl;
    // Check we are in a DNA volume
    if(fFlagVolume == 1 // d1
            || fFlagVolume == 11 // p1
            || fFlagVolume == 2 // d2
            || fFlagVolume == 22 // p2
            || fFlagVolume == 3 // cyto
            || fFlagVolume == 4 // gua
            || fFlagVolume == 5 // thy
            || fFlagVolume == 6 // ade
            || fFlagVolume == 7 // d1_w
            || fFlagVolume == 71 // p1_w
            || fFlagVolume == 8 // d2_w
            || fFlagVolume == 81 // p2_w
            || fFlagVolume == 9 // ade_w
            || fFlagVolume == 10 // gua_w
            || fFlagVolume == 13 // cyto_w
            || fFlagVolume == 12) // thy_w
    {
        // *****************************************************
        // Saving physical stage informations
        // *****************************************************
        if(step->GetPostStepPoint()->GetTouchable()->GetVolume() ) // avoid asking non existing information in postStep
        {
            G4double x=step->GetPreStepPoint()->GetPosition().x()/nanometer;
            G4double y=step->GetPreStepPoint()->GetPosition().y()/nanometer;
            G4double z=step->GetPreStepPoint()->GetPosition().z()/nanometer;

            G4double dE = step->GetTotalEnergyDeposit()/eV;
            G4double copyNo = G4double (step->GetPreStepPoint()->GetTouchable()->GetVolume()->GetCopyNo() );
                                        
            G4int eventId = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
#ifdef USE_MPI
            auto g4MPI = G4MPImanager::GetManager();
            if (g4MPI->IsSlave()) { // update eventID only for slave, cause rank_master=0
                G4int rank = g4MPI->GetRank();
                eventId += g4MPI->GetEventsInMaster() + (rank-1)*g4MPI->GetEventsInSlave();
            }
#endif
            // ***************************
            // put to vector
            // ***************************
            InfoInPhysStage aInfo;
            aInfo.fFlagParticle = fFlagParticle;
            aInfo.fFlagParentID = fFlagParentID;
            aInfo.fFlagProcess = fFlagProcess;
            aInfo.fX = x;
            aInfo.fY = y;
            aInfo.fZ = z;
            aInfo.fEdep = dE;
            aInfo.fEventNumber = eventId;
            aInfo.fVolumeName = fFlagVolume;
            aInfo.fCopyNumber = copyNo;
            aInfo.fLastMetVoxelCopyNum = fLastMetVoxelCopyNumber;
            Analysis::GetAnalysis()->AddInfoInPhysStage(aInfo);
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SteppingAction::SetupFlags(const G4Step* step)
{
    fFlagParticle=0;
    fFlagProcess=0;
    fFlagParentID=0;
    fFlagVolume=0;

    fFlagParentID = step->GetTrack()->GetParentID();

    SetupParticleAndProcessFlags(step);

    fFlagVolume = SetupVolumeFlag(step->GetPreStepPoint()->GetTouchable()->GetVolume()->GetName() );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SteppingAction::SetupParticleAndProcessFlags(const G4Step *step)
{
    fFlagParticle = 0;
    fFlagProcess = 0;
    auto partDef = step->GetTrack()->GetDynamicParticle()->GetDefinition();
    auto* instance = G4DNAGenericIonsManager::Instance(); 

    if (partDef == G4Electron::ElectronDefinition())       fFlagParticle = 1;
    if (partDef == G4Proton::ProtonDefinition())   fFlagParticle = 2;
    if (partDef == instance->GetIon("hydrogen")) fFlagParticle = 3;
    if (partDef == G4Alpha::AlphaDefinition())    fFlagParticle = 4;
    if (partDef == instance->GetIon("alpha+"))   fFlagParticle = 5;
    if (partDef == instance->GetIon("helium"))   fFlagParticle = 6;

    G4int procSubtype = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessSubType();


    if (fFlagParticle == 1) { // e-
        if (procSubtype == 58) fFlagProcess = 10; //"e-_G4DNAElectronSolvation"
        if (procSubtype == 51) fFlagProcess = 11; //"e-_G4DNAElastic"
        if (procSubtype == 52) fFlagProcess = 12; //"e-_G4DNAExcitation"
        if (procSubtype == 53) fFlagProcess = 13; //"e-_G4DNAIonisation"
        if (procSubtype == 55) fFlagProcess = 14; //"e-_G4DNAAttachment"
        if (procSubtype == 54) fFlagProcess = 15; //"e-_G4DNAVibExcitation"
    }

    if (fFlagParticle == 2) { // Proton
        if (procSubtype == 52) fFlagProcess = 17; //"proton_G4DNAExcitation"
        if (procSubtype == 53) fFlagProcess = 18; //"proton_G4DNAIonisation"
        if (procSubtype == 56) fFlagProcess = 19; //"proton_G4DNAChargeDecrease"
    }

    if (fFlagParticle == 3) { // hydrogen
        if (procSubtype == 52) fFlagProcess = 20; //"hydrogen_G4DNAExcitation"
        if (procSubtype == 53) fFlagProcess = 21; //"hydrogen_G4DNAIonisation"
        if (procSubtype == 57) fFlagProcess = 22; //"hydrogen_G4DNAChargeIncrease"
    }

    if (fFlagParticle == 4) { // alpha
        if (procSubtype == 52) fFlagProcess = 23; //"alpha_G4DNAExcitation"
        if (procSubtype == 53) fFlagProcess = 24; //"alpha_G4DNAIonisation"
        if (procSubtype == 56) fFlagProcess = 25; //"alpha_G4DNAChargeDecrease"
    }

    if (fFlagParticle == 5) { // alpha+
        if (procSubtype == 52) fFlagProcess = 26; //"alpha+_G4DNAExcitation"
        if (procSubtype == 53) fFlagProcess = 27; //"alpha+_G4DNAIonisation"
        if (procSubtype == 56) fFlagProcess = 28; //"alpha+_G4DNAChargeDecrease"
        if (procSubtype == 57) fFlagProcess = 29; //"alpha+_G4DNAChargeIncrease"
    }

    if (fFlagParticle == 6) { // helium
        if (procSubtype == 52) fFlagProcess = 30; //"helium__G4DNAExcitation"
        if (procSubtype == 53) fFlagProcess = 31; //"helium__G4DNAIonisation"
        if (procSubtype == 57) fFlagProcess = 32; //"helium__G4DNAChargeIncrease"
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int SteppingAction::SetupVolumeFlag(const G4String& volumeName)
{
    G4int flagVolume(-1);

    if(volumeName=="deoxyribose1_phys") flagVolume = 1;
    else if(volumeName=="phosphate1_phys") flagVolume = 11;
    else if(volumeName=="deoxyribose2_phys") flagVolume = 2;
    else if(volumeName=="phosphate2_phys") flagVolume = 22;
    else if(volumeName=="base_cytosine_phys") flagVolume = 3;
    else if(volumeName=="base_guanine_phys") flagVolume = 4;
    else if(volumeName=="base_thymine_phys") flagVolume = 5;
    else if(volumeName=="base_adenine_phys") flagVolume = 6;
    else if(volumeName=="deoxyribose1_water_phys") flagVolume = 7;
    else if(volumeName=="phosphate1_water_phys") flagVolume = 71;
    else if(volumeName=="deoxyribose2_water_phys") flagVolume = 8;
    else if(volumeName=="phosphate2_water_phys") flagVolume = 81;
    else if(volumeName=="base_adenine_water_phys") flagVolume = 9;
    else if(volumeName=="base_guanine_water_phys") flagVolume = 10;
    else if(volumeName=="base_cytosine_water_phys") flagVolume = 13;
    else if(volumeName=="base_thymine_water_phys") flagVolume = 12;
    else if(volumeName=="fiber") flagVolume = 110;
    else if(volumeName=="voxelStraight" || volumeName=="VoxelStraight") flagVolume = 161;
    else if(volumeName=="voxelRight" || volumeName=="VoxelRight") flagVolume = 162;
    else if(volumeName=="voxelLeft" || volumeName=="VoxelLeft") flagVolume = 163;
    else if(volumeName=="voxelUp" || volumeName=="VoxelUp") flagVolume = 164;
    else if(volumeName=="voxelDown" || volumeName=="VoxelDown") flagVolume = 165;
    else if(volumeName=="voxelStraight2" || volumeName=="VoxelStraight2") flagVolume = 261;
    else if(volumeName=="voxelRight2" || volumeName=="VoxelRight2") flagVolume = 262;
    else if(volumeName=="voxelLeft2" || volumeName=="VoxelLeft2") flagVolume = 263;
    else if(volumeName=="voxelUp2" || volumeName=="VoxelUp2") flagVolume = 264;
    else if(volumeName=="voxelDown2" || volumeName=="VoxelDown2") flagVolume = 265;
    else if(volumeName=="physWorld") flagVolume = -1.;
    else if(volumeName=="wrapper") flagVolume = 130;
    else if(volumeName=="histone_phys") flagVolume = 140;
    else if(volumeName=="nucleus_pl") flagVolume = 200;
    return flagVolume;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SteppingAction::SetupVoxelCopyNumber(const G4Step* step)
{
    // Each time we will do a step in a voxel, we will setup the lastMetVoxelCpNum flag.
    // This way, this number will always be the one of the last voxel met by the particle.
    // Therefore, in a DNA volume, the number will be the one of the container voxel.

    if(step->GetPostStepPoint()->GetTouchable()->GetVolume() ) // avoid asking non existing information in postStep
    {
        const G4String& volPost = step->GetPostStepPoint()->GetTouchable()->GetVolume()->GetName();

        // If we enter a voxel or stay in it
        if(    volPost == "VoxelStraight" || volPost=="voxelStraight"
            || volPost == "VoxelRight" || volPost=="voxelRight"
            || volPost == "VoxelLeft" || volPost=="voxelLeft"
            || volPost == "VoxelUp" || volPost=="voxelUp"
            || volPost == "VoxelDown" || volPost=="voxelDown"
            || volPost == "VoxelStraight2" || volPost=="voxelStraight2"
            || volPost == "VoxelRight2" || volPost=="voxelRight2"
            || volPost == "VoxelLeft2" || volPost=="voxelLeft2"
            || volPost == "VoxelUp2" || volPost=="voxelUp2"
            || volPost == "VoxelDown2" || volPost=="voxelDown2")
        {
            fLastMetVoxelCopyNumber = step->GetPostStepPoint()->GetTouchable()->GetCopyNumber(); 
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

