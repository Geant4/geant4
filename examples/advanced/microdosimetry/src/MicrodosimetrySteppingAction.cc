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
// -------------------------------------------------------------------
// $Id: MicrodosimetrySteppingAction.cc,v 1.1 2007-10-09 08:00:29 sincerti Exp $
// -------------------------------------------------------------------

#include "G4SteppingManager.hh"

#include "MicrodosimetrySteppingAction.hh"
#include "MicrodosimetryRunAction.hh"
#include "MicrodosimetryDetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

MicrodosimetrySteppingAction::MicrodosimetrySteppingAction(MicrodosimetryRunAction* run,MicrodosimetryDetectorConstruction* det)
:Run(run),Detector(det)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

MicrodosimetrySteppingAction::~MicrodosimetrySteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MicrodosimetrySteppingAction::UserSteppingAction(const G4Step* aStep)
  
{ 


// TOTAL DOSE DEPOSIT AND DOSE DEPOSIT WITHIN A PHANTOM VOXEL

if ( 	   (aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName()  == "physicalNucleus")
  	&& (aStep->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "physicalNucleus")
  	&& (
	
	aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName() == "e-"
	||
	aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName() == "alpha"
	||
	aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName() == "alpha+" 
	||
	aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName() == "helium") 

	)
	
	{ 
   	 G4double dose = (e_SI*(aStep->GetTotalEnergyDeposit()/eV))/(Run->GetMassNucleus());
   	 Run->AddDoseN(dose);

	 G4ThreeVector v;
    	 Run->AddDoseBox(aStep->GetPreStepPoint()->GetTouchableHandle()->GetReplicaNumber(),
	  aStep->GetTotalEnergyDeposit()/eV);
	}

if ( 	   (aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName()  == "physicalCytoplasm")
  	&& (aStep->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "physicalCytoplasm")
  	&& (
	
	aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName() == "e-"
	||
	aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName() == "alpha" 
	||
	aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName() == "alpha+" 
	||
	aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName() == "helium") 

	)
	{ 
   	 G4double dose = (e_SI*(aStep->GetTotalEnergyDeposit()/eV))/(Run->GetMassCytoplasm());
   	 Run->AddDoseC(dose);

	 G4ThreeVector v;
    	 Run->AddDoseBox(aStep->GetPreStepPoint()->GetTouchableHandle()->GetReplicaNumber(),
	  aStep->GetTotalEnergyDeposit()/eV);
 	}
}
