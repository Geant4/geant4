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
// $Id: RemSimSteppingAction.cc,v 1.1 2004-02-03 09:16:47 guatelli Exp $
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
#include "G4ParticleDefinition.hh"
#include "G4VPhysicalVolume.hh"
#include "RemSimSteppingAction.hh"
#include "RemSimRunAction.hh"
#include "RemSimPrimaryGeneratorAction.hh"
#include "RemSimDetectorConstruction.hh"

#ifdef G4ANALYSIS_USE
#include "RemSimAnalysisManager.hh"
#endif

RemSimSteppingAction::RemSimSteppingAction(RemSimPrimaryGeneratorAction* primary,RemSimRunAction* run, RemSimDetectorConstruction* det):
  primaryAction(primary),
  runAction(run),
  detector(det)
{ }

RemSimSteppingAction::~RemSimSteppingAction()
{ }

void RemSimSteppingAction::UserSteppingAction(const G4Step* aStep)
{  
  
  G4double charge =  aStep->GetTrack()->GetDynamicParticle()->
                                            GetDefinition()->GetPDGCharge(); 
 if ((aStep->GetTrack()->GetParentID()== 0)&&(charge != 0))
    {
      G4int stepNumber = aStep->GetTrack()->GetCurrentStepNumber();
      G4int runID = runAction->GetRunID();
      G4double primaryParticleEnergy = primaryAction->GetInitialEnergy();
	   
     if (stepNumber == 1)
	// Stopping Power test ...
	{ 
	  G4double initialStepXPosition = aStep->GetPreStepPoint()->GetPosition().x();
	  G4double initialStepYPosition = aStep->GetPreStepPoint()->GetPosition().y();
	  G4double initialStepZPosition = aStep->GetPreStepPoint()->GetPosition().z();
	  if (0 == initialStepXPosition && 
	      0 == initialStepYPosition &&
	      0 == initialStepZPosition)
	    { 
	      G4double energyLost = abs(aStep->GetDeltaEnergy());
	      G4double stepLength = aStep->GetTrack()-> GetStepLength();
               if (stepLength != 0) 
		{
		  // calculation of Stopping Power 
		  // reference http://physics.nist.gov/PhysRefData/Star/Text/contents.html
		  G4double totalStoppingPower = energyLost / stepLength;
                  
		  G4double targetDensity = detector->GetDensity();
           
		  G4double nistStoppingPower = totalStoppingPower/targetDensity;
#ifdef G4ANALYSIS_USE                
  RemSimAnalysisManager* analysis = RemSimAnalysisManager::getInstance();	  analysis -> StoppingPower
           (runID,primaryParticleEnergy/MeV,nistStoppingPower/(MeV*(cm2/g)));
#endif
		}
	    }
	}
    
      // CSDA range ...
      G4double  ParticleKineticEnergy = aStep->GetTrack()->GetKineticEnergy();	     
      if (0.0 == ParticleKineticEnergy) 
	{  
	  G4double xend = aStep->GetTrack()->GetPosition().x()/mm ;
	  G4double yend = aStep->GetTrack()->GetPosition().y()/mm ;
	  G4double zend = aStep->GetTrack()->GetPosition().z()/mm ;
                 
	  // calculation of CSDA range 
	  // reference http://physics.nist.gov/PhysRefData/Star/Text/contents.html
	  G4double range = (sqrt(xend*xend + yend*yend + zend*zend)); 
	  G4double nistRange = range*(detector->GetDensity());
          G4int stepNumber = aStep->GetTrack()->GetCurrentStepNumber();
          G4cout<< stepNumber <<G4endl;

#ifdef G4ANALYSIS_USE                
  RemSimAnalysisManager* analysis = RemSimAnalysisManager::getInstance();		 analysis -> CSDARange(runID,
	  			primaryParticleEnergy/MeV,
	  			nistRange/(g/cm2));
#endif
	}
    }
}






