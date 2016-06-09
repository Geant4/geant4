//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: RemSimSteppingAction.cc,v 1.6 2004/11/23 15:43:41 guatelli Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
//
// Author: Susanna Guatelli (guatelli@ge.infn.it)
//

#include "G4ios.hh"
#include "G4SteppingManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4ParticleDefinition.hh"
#include "G4VPhysicalVolume.hh"
#include "RemSimSteppingAction.hh"
#include "RemSimPrimaryGeneratorAction.hh"
#include "G4TrackStatus.hh"
#include "RemSimSteppingActionMessenger.hh"
#ifdef G4ANALYSIS_USE
#include "RemSimAnalysisManager.hh"
#endif

RemSimSteppingAction::RemSimSteppingAction(RemSimPrimaryGeneratorAction* primary):
  primaryAction(primary)
{ 
  hadronic = "Off";
  messenger = new RemSimSteppingActionMessenger(this);
}

RemSimSteppingAction::~RemSimSteppingAction()
{ 
  delete messenger;
}

void RemSimSteppingAction::UserSteppingAction(const G4Step* aStep)
{ 
  G4String oldVolumeName = aStep -> GetPreStepPoint()->  
                                      GetPhysicalVolume()-> GetName();  
 
  // Store in histograms useful information concerning primary particles

#ifdef G4ANALYSIS_USE    
  RemSimAnalysisManager* analysis = RemSimAnalysisManager::getInstance();

  G4VPhysicalVolume* volume = aStep -> GetPostStepPoint() -> GetPhysicalVolume();
 
  G4String particleName = (aStep -> GetTrack() -> GetDynamicParticle()
			   -> GetDefinition() -> GetParticleName());

  if(oldVolumeName != "phantom") 
    { 
      if (volume) 
	{
	  G4String volumeName = volume -> GetName();
	  if (volumeName == "phantom") 
	    {  
	      G4double initialEnergy = primaryAction -> GetInitialEnergy(); 
	      G4double particleEnergy = aStep -> GetTrack() -> GetKineticEnergy();
	      if ( particleName == "proton" || 
                   particleName == "alpha"  ||
		   particleName == "IonO16" || 
		   particleName == "IonC12" || 
                   particleName == "IonFe52"|| 
		   particleName == "IonSi28")
		{
		 analysis -> PrimaryInitialEnergyIn(initialEnergy);
	         analysis -> PrimaryEnergyIn(particleEnergy);
		}
	    }    
	}
    }

  if(oldVolumeName == "phantom") 
    { 
      if (volume) 
	{
	  G4String volumeName = volume -> GetName();
             
	  if (volumeName != "phantom") 
	    {
	     G4double initialEnergy = primaryAction -> GetInitialEnergy(); 
	     G4double particleEnergy = aStep->GetTrack() -> GetKineticEnergy();
             if ( particleName == "proton" || 
                  particleName == "alpha"  || 
                  particleName == "IonO16" || 
		  particleName == "IonC12" || 
                  particleName == "IonFe52"|| 
		  particleName == "IonSi28" )
		{
	         analysis -> PrimaryInitialEnergyOut(initialEnergy); 
	         analysis -> PrimaryEnergyOut(particleEnergy);
		}
	    }
	}
    }
#endif

  //analysis of hadronic physics
  if (hadronic == "On")
    {
      if(aStep -> GetPostStepPoint() -> GetProcessDefinedStep() != NULL)
	{
	  G4String process = aStep -> GetPostStepPoint() ->
                           GetProcessDefinedStep() ->GetProcessName();

	  if ((process != "Transportation") &&
	      (process != "msc") && 
              (process != "LowEnergyIoni") &&
	      (process != "LowEnergyBrem") && 
              (process != "eIoni") &&
	      (process != "hIoni") &&
              (process != "eBrem") && 
              (process != "compt") &&
	      (process != "phot")  && 
              (process != "conv")  &&
              (process != "annihil") &&
              (process != "hLowEIoni") &&
	      (process != "LowEnBrem") && 
              (process != "LowEnCompton") && 
              (process != "LowEnPhotoElec") && 
              (process != "LowEnRayleigh") && 
              (process != "LowEnConversion"))
	      G4cout << "Hadronic Process:" << process << G4endl;
	}
    }
}
void RemSimSteppingAction::SetHadronicAnalysis(G4String value)    
{
  hadronic = value;
}






