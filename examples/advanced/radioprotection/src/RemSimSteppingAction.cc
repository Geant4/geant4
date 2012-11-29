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
// $Id$
//
//
// Author: Susanna Guatelli (susanna@uow.edu.au)
//

#include "RemSimSteppingAction.hh"
#include "RemSimPrimaryGeneratorAction.hh"
#include "RemSimSteppingActionMessenger.hh"
#include "RemSimAnalysis.hh"

#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4SteppingManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4ParticleDefinition.hh"
#include "G4VPhysicalVolume.hh"
#include "G4TrackStatus.hh"

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
 
  // Store in histograms useful information concerning secondary particles

  G4VPhysicalVolume* volume = aStep -> GetPostStepPoint() -> GetPhysicalVolume();
 
 // Get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
 
  // Secondary particles reaching the astronaut
  if(oldVolumeName != "phantom") 
    { 
      if (volume) 
	{
	  G4String volumeName = volume -> GetName();
	 if (volumeName == "phantom") 
	    { 
               G4String particleName = (aStep -> GetTrack() -> GetDynamicParticle()
			   -> GetDefinition() -> GetParticleName());

	       G4double particleEnergy = aStep -> GetTrack() -> GetKineticEnergy();
	     
	       if(aStep -> GetTrack() -> GetTrackID()!= 1) 
               // analysis of secondary particles reaching the astronaut
		{
		  if (particleName == "proton") analysisManager -> FillH1(1, particleEnergy/MeV); 
		
	 	 else{
	      	       if (particleName == "neutron") analysisManager -> FillH1(2, particleEnergy/MeV);
	
 			  else{
 	                    	if (particleName == "pi+" || particleName == "pi-" ||particleName == "pi0" ) 	       
			 	    analysisManager -> FillH1(3,particleEnergy/MeV);		   
		   	  else{
		       		if (particleName == "alpha") analysisManager -> FillH1(4, particleEnergy/MeV); 	
		               }
			       }
			}
		}
    	}
   }
}
  

  // Analysis of the secondary particles generated in the phantom

  G4SteppingManager*  steppingManager = fpSteppingManager;
  G4Track* theTrack = aStep -> GetTrack();

  // check if it is alive
  if(theTrack-> GetTrackStatus() == fAlive) { return; }
	
  // Retrieve the secondary particles
   G4TrackVector* fSecondary = steppingManager -> GetfSecondary();
     
   for(size_t lp1=0;lp1<(*fSecondary).size(); lp1++)
     { 
       // Retrieve the info about the generation of secondary particles in the phantom and
       // in the vehicle
       G4String volumeName = (*fSecondary)[lp1] -> GetVolume() -> GetName(); 
       G4String secondaryParticleName =  (*fSecondary)[lp1]->GetDefinition() -> GetParticleName();  
       G4double secondaryParticleKineticEnergy =  (*fSecondary)[lp1] -> GetKineticEnergy(); 
       G4String process = (*fSecondary)[lp1]-> GetCreatorProcess()-> GetProcessName();   
      
       // If the secondaries are originated in the phantom....
       if (volumeName == "phantom")
	 {	   
	   
	   if (secondaryParticleName == "proton") analysisManager -> FillH1(5, secondaryParticleKineticEnergy/MeV); 
	      
	   else
	     {
	       if (secondaryParticleName == "neutron") analysisManager -> FillH1(6, secondaryParticleKineticEnergy/MeV);
		   
	       else{
			 if (secondaryParticleName == "pi+" ||
			     secondaryParticleName == "pi-"||secondaryParticleName == "pi0")   		
		     	    analysisManager -> FillH1(7, secondaryParticleKineticEnergy/MeV);
		 		
			else{
		  	      if (secondaryParticleName == "alpha") 
				 analysisManager -> FillH1(8, secondaryParticleKineticEnergy/MeV);
		       	    }
		     }
		   }

		 }

	       }

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
