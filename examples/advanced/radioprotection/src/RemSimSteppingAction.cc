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
// $Id: RemSimSteppingAction.cc,v 1.2 2004-03-12 10:55:55 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Author: Susanna Guatelli (guatelli@ge.infn.it)
//
// History:
// -----------
//27 May  2003   S.Guatelli    first code review 
//17 May  2003   S. Guatelli   1st implementation
//
// -------------------------------------------------------------------
#include "G4ios.hh"
#include <math.h> // standard c math library
#include "G4SteppingManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4ParticleDefinition.hh"
#include "G4VPhysicalVolume.hh"
#include "RemSimSteppingAction.hh"
#include "RemSimEventAction.hh"
#include "RemSimPrimaryGeneratorAction.hh"
#include "RemSimDetectorConstruction.hh"
#include "G4TrackStatus.hh"
#ifdef G4ANALYSIS_USE
#include "RemSimAnalysisManager.hh"
#endif

RemSimSteppingAction::RemSimSteppingAction(RemSimPrimaryGeneratorAction* primary,RemSimEventAction* event, RemSimDetectorConstruction* det):
  primaryAction(primary),
  eventAction(event),IDnow(-2),IDold(-1),
  detector(det)
{ }

RemSimSteppingAction::~RemSimSteppingAction()
{ }

void RemSimSteppingAction::UserSteppingAction(const G4Step* aStep)
{ 
  G4int evno = eventAction->GetEventNo(); 
  G4int trno = aStep->GetTrack()->GetTrackID(); 
  IDnow = trno + evno*10000; 
  
#ifdef G4ANALYSIS_USE    
 RemSimAnalysisManager* analysis = RemSimAnalysisManager::getInstance();

 G4String oldVolumeName = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();  
 G4VPhysicalVolume* volume = aStep->GetPostStepPoint()->GetPhysicalVolume();
 
 G4String particleName = (aStep->GetTrack()->GetDynamicParticle()
      ->GetDefinition()-> GetParticleName());

 if(aStep->GetTrack()->GetTrackID()!=1)
   { 
     G4int stepNumber = aStep->GetTrack()->GetCurrentStepNumber();
     //Initial Kinetic Energy
      if (stepNumber == 1)
       {
	 G4double particleEnergy = aStep->GetTrack()->GetKineticEnergy();
	 if (particleName =="neutron")
	   analysis->neutronEnergyDistribution(particleEnergy/MeV); 

	 if (particleName =="gamma" )
	   analysis->photonEnergyDistribution(particleEnergy/MeV); 
 
	 if (particleName =="e-" )
	   analysis->electronEnergyDistribution(particleEnergy/MeV); 
    
	 G4double charge = aStep->GetTrack()->GetDynamicParticle()->
	   GetDefinition()->GetPDGCharge();

	 if (charge!=0)
	   {if ((particleName !="alpha") && (particleName !="proton")
		&& (particleName !="e+") && (particleName !="e-") )
	   analysis->hadronEnergyDistribution(particleEnergy/MeV); 
	   }
       }
   }

 //Analysis of other particles of interest
   if(oldVolumeName == "world") 
   { 
    if (volume) 
      {
        G4String volumeName = volume->GetName();
       if (volumeName == "phantom") 
     {   
      G4double particleEnergy = aStep->GetTrack()->GetKineticEnergy();
     
     if (particleName == "proton" || particleName == "alpha" )
	{
	  analysis -> hadronEnergySpectrum1(particleEnergy/MeV);         
	}
      if (particleName == "e-" || particleName == "e+" )
	{
	  analysis -> leptonsEnergySpectrum1(particleEnergy/MeV);         
	}
     if (particleName == "gamma" )
	{
	  analysis -> gammaEnergySpectrum1(particleEnergy/MeV);         
	}
     }
      }
   }

   if(oldVolumeName == "phantom") 
     { 
       if (volume) 
	 {
	   G4String volumeName = volume->GetName();
             
	   if (volumeName == "detector") 
	     {  
	       G4double particleEnergy = aStep->GetTrack()->GetKineticEnergy();
	       if (particleName == "proton" || particleName == "alpha" )
		 {
		   analysis -> hadronEnergySpectrum3(particleEnergy/MeV);
		 }
	       if (particleName == "e-" ||particleName == "e+" )
		 {
		   analysis -> leptonsEnergySpectrum3(particleEnergy/MeV);         
		 }
	       if (particleName == "gamma" )
		 {
		   analysis -> gammaEnergySpectrum3(particleEnergy/MeV);         		 }
	     }
	 }
     }     


   if(oldVolumeName == "phantom") 
     { 
       if (volume) 
	 {
	   G4String volumeName = volume->GetName();
             
	   if (volumeName == "phantom") 
	     {  
	       G4double zIn=aStep->GetPreStepPoint()->GetPosition().z();
	       G4double zOut=aStep->GetPostStepPoint()->GetPosition().z();
    
	       G4double particleEnergy = aStep->GetTrack()->GetKineticEnergy();
  
	       if ((zIn>132.*cm && zIn < 133.*cm) && 
		   (zOut>133.*cm && zOut < 134.*cm))
		 {
                  if(IDnow != IDold) 
                  {
                    IDold=IDnow ;
		   if (particleName == "proton" ||particleName == "alpha" )
		     {
		       analysis -> hadronEnergySpectrum2(particleEnergy/MeV);         
		     }
		   if (particleName == "e-" ||particleName == "e+" )
		     {
		       analysis -> leptonsEnergySpectrum2(particleEnergy/MeV);         
		     }
		   if (particleName == "gamma" )
		     {
		       analysis -> gammaEnergySpectrum2(particleEnergy/MeV);         
		     }
		  }
		 }
	     }
	 }
     }
#endif
}
    







