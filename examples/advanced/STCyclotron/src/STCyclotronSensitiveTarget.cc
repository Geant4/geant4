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
//  Author: F. Poignant, floriane.poignant@gmail.com
//
// file STCyclotronSensitiveTarget.cc
//
#include "STCyclotronAnalysis.hh"
#include "STCyclotronRun.hh"
#include "STCyclotronSensitiveTarget.hh"

#include "G4RunManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"
#include "G4Step.hh"
#include "G4SteppingManager.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4DecayTable.hh"
#include "G4VDecayChannel.hh"
#include "G4RadioactiveDecay.hh"
#include "G4TrackVector.hh"
#include "G4VProcess.hh"
#include "G4Tubs.hh"
#include <map>

STCyclotronSensitiveTarget::STCyclotronSensitiveTarget(G4String name,
						       STCyclotronDetectorConstruction* det)
  : G4VSensitiveDetector(name), 
    fDet(det)
{ 
  fTempTrack = 0;
  fTempTrack1 = 0;
  fTempEnergy = 0.;
  fTempVector = G4ThreeVector(0.,0.,0.);
  fTrack=0;
}

STCyclotronSensitiveTarget::~STCyclotronSensitiveTarget()
{
 delete fTrack;
 }

G4bool STCyclotronSensitiveTarget::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  
  STCyclotronRun* fRun = static_cast<STCyclotronRun*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  fTrack = aStep->GetTrack();

  auto analysisManager = G4AnalysisManager::Instance();
  
  //----------------------------------------------
  //        Volume info
  //----------------------------------------------
  G4double targetHalfDiameter= (fDet->GetTargetDiameter())/2.;
  
  //----------------------------------------------
  //       Step information
  //----------------------------------------------
  
  G4double edep = aStep->GetTotalEnergyDeposit();
  G4double energy = aStep->GetPreStepPoint()->GetKineticEnergy();
  G4ThreeVector momentumDirection = aStep->GetPreStepPoint()->GetMomentumDirection();
  G4ThreeVector vectorPosition = aStep->GetPreStepPoint()->GetPosition();
  
  //----------------------------------------------
  //        Track
  //----------------------------------------------

  G4ParticleDefinition* thePartDef = fTrack->GetDefinition();
  G4String partType= fTrack->GetDefinition()->GetParticleType();
  G4String name = fTrack->GetDefinition()->GetParticleName();
  G4double timeLife = fTrack->GetDefinition()->GetPDGLifeTime(); //<---
  const G4VProcess* process = fTrack->GetCreatorProcess();
 
  //----------------------------------------------
  //    Collect general information concerning all of the particles
  //    Collect energy deposition ; separe decay case to beam case
  //----------------------------------------------

  fRun->EnergyDepositionTarget(edep);

  
  //----------------------------------------------
  //Collect information about protons and deuterons
  //----------------------------------------------

  if(name == "proton" || name == "deuteron")
    {
   
      if(fTrack->GetTrackID()!=fTempTrack && (momentumDirection.getZ()>0.) &&
	 vectorPosition.getX()<targetHalfDiameter &&  
	 vectorPosition.getX()>-targetHalfDiameter &&
	 vectorPosition.getY()<targetHalfDiameter &&  
	 vectorPosition.getY()>-targetHalfDiameter)
	{
	  analysisManager->FillH2(0,vectorPosition.getX(),vectorPosition.getY());
	  analysisManager->FillH1(0,energy);
	  fRun->CountParticlesTarget(); 
	  fTempTrack = fTrack->GetTrackID();
	}
    
      if(fTempTrack1 == 0)
	{
	  fTempTrack1 = fTrack->GetTrackID();
	}

      if(fTrack->GetTrackID()!=fTempTrack1 && (momentumDirection.getZ()>0.) &&
	 fTempVector.getX()<targetHalfDiameter &&  
	 fTempVector.getX()>-targetHalfDiameter &&
	 fTempVector.getY()<targetHalfDiameter &&  
	 fTempVector.getY()>-targetHalfDiameter )
	{

	  analysisManager->FillH2(4,fTempVector.getX(),fTempVector.getY());
	  analysisManager->FillH1(2,fTempEnergy);
	  fTempTrack1 = fTrack->GetTrackID(); 
	}
      
      fTempVector =  aStep->GetPostStepPoint()->GetPosition();      //vectorPosition;
      fTempEnergy = aStep->GetPostStepPoint()->GetKineticEnergy();  //energy;

      analysisManager->FillH2(3,vectorPosition.getZ(),energy);  
    
    }
 
  //----------------------------------------------
  //    Store ID for particles that are 
  //    not protons/electrons or deuterons
  //----------------------------------------------

  if((name != "proton") && (name != "e-") && (name != "deuteron"))
    {
      fRun->StoreIsotopeID(fTrack->GetTrackID(),name);
    }

  //----------------------------------------------
  //  Collect of information for unstable isotopes
  //  generated from an interaction with the target
  //----------------------------------------------

  if (name!="deuteron")
    {
      if (( partType == "nucleus") &&  !(thePartDef->GetPDGStable()) && (fTrack->GetCurrentStepNumber()==1) && timeLife!=0.)
	{
	  //G4cout << "Saving unstable particles ..." << G4endl;
	  G4int Z=thePartDef->GetAtomicNumber();
	  G4int A=thePartDef->GetAtomicMass();
	  analysisManager->FillH2(2,Z,A);
	  
	  //----------------------------------------------
	  //   isotopes count
	  //----------------------------------------------
	  
	  fRun->PrimaryIsotopeCountTarget(name,timeLife); 
	  analysisManager->FillH1(4,fTrack->GetPosition().getZ());
	  //particle that created the nucleus
	  std::map<G4int,G4String> parentID = fRun->GetIsotopeID();
	  G4String nameParent = parentID[fTrack->GetParentID()];
	  fRun->ParticleParent(name, process->GetProcessName()); 
	  
	  //G4cout << name << " : " << process->GetProcessName() << " with track ID " << fTrack->GetTrackID() << " and step ID " << fTrack->GetCurrentStepNumber() << G4endl;

	}
    }
  
  //----------------------------------------------
  //   Collect of information for stable isotopes
  //   generated from an interaction with the target
  //----------------------------------------------

  if (name!="deuteron")
    {
      if (( partType == "nucleus") &&  (thePartDef->GetPDGStable()) && (process->GetProcessName() != "RadioactiveDecay")  && (fTrack->GetCurrentStepNumber()==1) )
	{
	  //----------------------------------------------
	  //    isotopes count
	  //----------------------------------------------
	  fRun->CountStableIsotopes(name); 
	}
    }

  
  //----------------------------------------------
  //   Collect unstable isotopes from decay
  //----------------------------------------------

  if (( partType == "nucleus") &&  !(thePartDef->GetPDGStable()) && (process->GetProcessName() == "RadioactiveDecay")  && (fTrack->GetCurrentStepNumber()==1) && timeLife!=0)
    {
      std::map<G4int,G4String>::iterator itbis;
      std::map<G4int,G4String> parentID = fRun->GetIsotopeID();
      G4String nameParent = parentID[fTrack->GetParentID()];
      fRun->DecayIsotopeCountTarget(name,nameParent,timeLife);
    }
  
  //----------------------------------------------
  //   Collect any other particles emitted
  //----------------------------------------------

  if((partType!="nucleus")&&(name!="proton")&&(name!="deuteron"))
    {
      fRun->ParticleCountTarget(name); 
      //Condition so the particle will be counted for only one step
      if((fTrack->GetCurrentStepNumber()==1))
	{
	  if(process->GetProcessName() != "RadioactiveDecay")
	    {
	      if(name=="e+"){
		analysisManager->FillH1(5,energy);
	      }
	      if(name=="e-"){
		analysisManager->FillH1(6,energy);
	      }
	      if(name=="gamma"){
		analysisManager->FillH1(7,energy);
	      }
	      if(name=="neutron"){
		analysisManager->FillH1(8,energy);
	      }
	    }
      
	  if(process->GetProcessName() == "RadioactiveDecay")
	    {
	      if(name=="e+"){
		analysisManager->FillH1(9,energy);
	      }
	      if(name=="e-"){
		analysisManager->FillH1(10,energy);
	      }
	      if(name=="gamma"){
		analysisManager->FillH1(11,energy);
	      }
	      if(name=="neutron"){
		analysisManager->FillH1(12,energy);
	      }
	      if(name=="nu_e"){
		analysisManager->FillH1(13,energy);
	      }
	      if(name=="anti_nu_e"){
		analysisManager->FillH1(14,energy);
	      }
	    }
	}
    }


  fRun->SetTargetVolume(fDet->GetTargetVolume());
  fRun->SetTargetThickness(fDet->GetTargetThickness());
  fRun->SetTargetDiameter(fDet->GetTargetDiameter());

  return true;

}


