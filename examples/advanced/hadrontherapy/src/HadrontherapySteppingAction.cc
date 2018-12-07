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
// Hadrontherapy advanced example for Geant4
// See more at: https://twiki.cern.ch/twiki/bin/view/Geant4/AdvancedExamplesHadrontherapy

#include "G4SteppingManager.hh"
#include "G4TrackVector.hh"
#include "HadrontherapySteppingAction.hh"
#include "G4ios.hh"
#include "G4SteppingManager.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4TrackVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4UserEventAction.hh"
#include "G4TransportationManager.hh"
#include "G4VSensitiveDetector.hh"
#include "HadrontherapyRunAction.hh"
#include "HadrontherapyMatrix.hh"
#include "G4SystemOfUnits.hh"

/////////////////////////////////////////////////////////////////////////////
HadrontherapySteppingAction::HadrontherapySteppingAction( HadrontherapyRunAction *run)
{
    runAction = run;
}

/////////////////////////////////////////////////////////////////////////////
HadrontherapySteppingAction::~HadrontherapySteppingAction()
{
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapySteppingAction::UserSteppingAction(const G4Step* aStep)
{
    
    // The followings are calls to usefuls information retrieved at the step level
    // Please, comment out them if want to use
    
   // G4Track* theTrack = aStep->GetTrack();
    
    G4StepPoint* PreStep = aStep->GetPreStepPoint();
    G4StepPoint* PostStep = aStep->GetPostStepPoint();
    
    G4TouchableHandle touchPreStep = PreStep->GetTouchableHandle();
    G4TouchableHandle touchPostStep = PostStep->GetTouchableHandle();
    
    //G4double PreStepX =PreStep->GetPosition().x();
    //G4double PreStepY =PreStep->GetPosition().y();
    //G4double PreStepZ =PreStep->GetPosition().z();
    
    //G4double PostStepX =PostStep->GetPosition().x();
    //G4double PostStepY =PostStep->GetPosition().y();
    //G4double PostStepZ  =PostStep->GetPosition().z();
    
    //To get the current volume:
    G4VPhysicalVolume* volumePre = touchPreStep->GetVolume();
    //G4VPhysicalVolume* volumePost =touchPostStep->GetVolume();
    
    //To get its name:
    G4String namePre = volumePre->GetName();

    
    // positions in the global coordinate system:
    //G4ThreeVector posPreStep  = PreStep->GetPosition();
    //G4ThreeVector posPostStep = PostStep->GetPosition();
    
    //G4int eventNum = G4RunManager::GetRunManager() -> GetCurrentEvent() -> GetEventID();
    
    //G4double parentID =aStep->GetTrack()->GetParentID();
    //G4double trackID =aStep->GetTrack()->GetTrackID();
    
    G4double eKin = aStep -> GetPreStepPoint() -> GetKineticEnergy();
    
    G4double PosX = aStep->GetTrack()->GetPosition().x();
    G4double PosY = aStep->GetTrack()->GetPosition().y();
    G4double PosZ = aStep->GetTrack()->GetPosition().z();
    
    G4String volume=  aStep->GetTrack()->GetVolume()->GetName();
    G4Track* theTrack = aStep->GetTrack();
    
    //G4String material= aStep -> GetTrack() -> GetMaterial() -> GetName();
    //G4cout << "material   " << material << G4endl;
    //G4String volume=  aStep->GetTrack()->GetVolume()->GetName();
    //G4String pvname= pv-> GetName();
    
    G4String particleName = aStep->GetTrack()->GetDefinition()->GetParticleName();

    G4double momentumX =  aStep->GetTrack()->GetMomentumDirection().x();
    G4double momentumY =  aStep->GetTrack()->GetMomentumDirection().y();
    G4double momentumZ =  aStep->GetTrack()->GetMomentumDirection().z();
    
    
    G4ParticleDefinition *particleDef = theTrack -> GetDefinition();
    G4int pdg = particleDef ->GetPDGEncoding();
    
    if(namePre == "VirtualLayer")
    {
        std::ofstream WriteDataIn("Virtual_Layer.txt", std::ios::app);
        WriteDataIn
        
        <<   eKin             <<" " //  1
        <<   PosX             <<" " //  2
        <<   PosY             <<" " //  3
        <<   PosZ             <<" " //  4
        <<   momentumX        <<" " //  5
        <<   momentumY        <<" " //  6
        <<   momentumZ        <<" " //  7
        <<   pdg
        //<<   theTrack         << '\t' << "   "
        
        <<   G4endl;
        
        theTrack -> SetTrackStatus(fKillTrackAndSecondaries);
        
        
    }
    
    
    
  
}



