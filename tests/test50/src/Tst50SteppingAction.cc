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
// $Id: Tst50SteppingAction.cc,v 1.43 2010-06-07 10:08:39 gcosmo Exp $
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
#include <cmath> 
#include "G4SteppingManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4VPhysicalVolume.hh"
#include "Tst50SteppingAction.hh"
#include "Tst50DetectorConstruction.hh"
#ifdef G4ANALYSIS_USE
#include "Tst50AnalysisManager.hh"
#endif
#include "Tst50PrimaryGeneratorAction.hh"
#include "Tst50RunAction.hh"

Tst50SteppingAction::Tst50SteppingAction(Tst50PrimaryGeneratorAction* primary, 
					 Tst50RunAction* run, 
					 Tst50DetectorConstruction* det):
  primaryAction(primary), 
  runAction(run), 
  detector(det) 
{ }

Tst50SteppingAction::~Tst50SteppingAction()
{ }

void Tst50SteppingAction::UserSteppingAction(const G4Step* aStep)
{
#ifdef G4ANALYSIS_USE
  Tst50AnalysisManager* analysis = Tst50AnalysisManager::getInstance();
#endif

  G4int runID = runAction->GetRunID();

  // The user can select transmission test with this flag
  G4bool trasmissionTestFlag = runAction->GetFlag();
  
  G4double  primaryParticleEnergy = primaryAction->GetInitialEnergy();
  G4double  ParticleKineticEnergy = aStep->GetTrack()->GetKineticEnergy();
  G4String primaryParticleName = primaryAction->GetParticle(); 
  G4double  particleZMomentumDirection = 
                            aStep->GetTrack()->GetMomentumDirection().z();
  
  // ---------- primary transmitted particles test -----------// 
  if (aStep ->GetTrack()->GetParentID() == 0)
    {
      G4double charge = 
      aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetPDGCharge();
      // Select primary particles for transmission test ...
      if (charge != 0 || primaryParticleName == "gamma" )
	{
          // Check if the primary particle is going outside the target ...
	  if ((aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() 
	       == "Target") && 
	      (aStep->GetTrack()->GetNextVolume()->GetName() == "World") )
	    { 
	      if (particleZMomentumDirection >= 0.)
                //transmitted primary charged particles
		{runAction -> TransmittedParticleNumber();}
	      else
                // backscattered primary charged particles
		{
                 runAction -> BackscatteredParticleNumber();
                 runAction -> BackscatteredEnergy(ParticleKineticEnergy);
                 } 
	      if ((primaryParticleEnergy == ParticleKineticEnergy) 
		  &&( particleZMomentumDirection == 1.))
		 // transmitted primary gamma 
		 {runAction -> TransmittedGammaNumber();}  
	    }
	}
    }

#ifdef G4ANALYSIS_USE

  // Stopping Power and CSDA range tests ...

  if (trasmissionTestFlag == false)
    { 
      if (aStep->GetTrack()->GetParentID() == 0)
	{
	  G4double charge = 
	    aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetPDGCharge();
	  // This test works just with charged particles ...
	  if (charge != 0) 
	    { 
	      G4int stepNumber = aStep->GetTrack()->GetCurrentStepNumber();
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
		      G4double energyLost = std::abs(aStep->GetDeltaEnergy());
		      G4double stepLength = aStep->GetTrack()-> GetStepLength();
		      if (stepLength != 0) 
			{
                          // calculation of Stopping Power 
                         // reference http://physics.nist.gov/PhysRefData/Star/Text/contents.html
			  G4double totalStoppingPower = energyLost / stepLength;
                          G4double targetDensity = detector->GetDensity();
			  G4double nistStoppingPower = totalStoppingPower/targetDensity;
			  analysis -> StoppingPower(runID,
                                                    primaryParticleEnergy/MeV,
                                                    nistStoppingPower/(MeV*(cm2/g)));
			}
		    }
		}
    
	      // CSDA range ...
	     
	      if (0.0 == ParticleKineticEnergy) 
		{  
		  G4double xend = aStep->GetTrack()->GetPosition().x()/mm ;
		  G4double yend = aStep->GetTrack()->GetPosition().y()/mm ;
		  G4double zend = aStep->GetTrack()->GetPosition().z()/mm ;
                 
                  // calculation of CSDA range 
                  // reference http://physics.nist.gov/PhysRefData/Star/Text/contents.html
		  G4double range = (std::sqrt(xend*xend + yend*yend + zend*zend)); 
		  G4double nistRange = range*(detector->GetDensity());
		  analysis -> CSDARange(runID,
                                        primaryParticleEnergy/MeV,
                                        nistRange/(g/cm2));
                  //G4cout<< nistRange/(g/cm2)<<G4endl;
		}
	    }
	}
    }

#endif

}




