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
// $Id: Tst50SteppingAction.cc,v 1.34 2003-05-18 10:42:35 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Author: Susanna Guatelli (guatelli@ge.infn.it)
//
// History:
// -----------
// 17 May  2003   S. Guatelli   1st implementation
//
// -------------------------------------------------------------------
#include "G4ios.hh"
#include <math.h> // standard c math library
#include "G4SteppingManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4VPhysicalVolume.hh"
#include "Tst50SteppingAction.hh"
#include "Tst50DetectorConstruction.hh"
#include "Tst50EventAction.hh"
#ifdef G4ANALYSIS_USE
#include "Tst50AnalysisManager.hh"
#endif
#include "Tst50PrimaryGeneratorAction.hh"
#include "Tst50RunAction.hh"

Tst50SteppingAction::Tst50SteppingAction(Tst50EventAction* event, 
					 Tst50PrimaryGeneratorAction* primary, 
					 Tst50RunAction* run, 
					 Tst50DetectorConstruction* det):
  idOld(-1),
  eventAction (event), 
  primaryAction(primary), 
  runAction(run), 
  detector(det) 
{ }

Tst50SteppingAction::~Tst50SteppingAction()
{ }

void Tst50SteppingAction::UserSteppingAction(const G4Step* Step)
{
  Tst50AnalysisManager* analysis = Tst50AnalysisManager::getInstance();

  G4int eventNumber = eventAction -> GetEventNo() ;
  G4int runID = runAction -> GetRunID();
  G4bool trasmissionTestFlag = runAction -> GetFlag();

  G4double  primaryParticleEnergy = primaryAction->GetInitialEnergy();
 
  G4int IDnow = runID + 1000 * eventNumber + 10000 * (Step->GetTrack()->GetTrackID()) +
    100000000 * (Step->GetTrack()->GetParentID()); 

  G4double  ParticleKineticEnergy = Step->GetTrack()->GetKineticEnergy();
  G4String primaryParticleName = primaryAction->GetParticle(); 
  G4double  particleZMomentumDirection = Step->GetTrack()->GetMomentumDirection().z();
    
  // ---------- energy of primary transmitted particles-----------// 

  if (primaryParticleName == "e-" || 
      primaryParticleName == "e+"|| 
      primaryParticleName == "proton" ||
      primaryParticleName == "gamma" )
    { 
      if (Step->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "Target")
	{
	  // The particle exits the Target
	  if (Step->GetTrack()->GetNextVolume()->GetName() == "World") 
	    { 
	      if (0 == Step->GetTrack()->GetParentID()) 
		{ 
		  if (IDnow != idOld) 
		    {
		      // The particle is a primary one
		      idOld = IDnow ;  
		      if (particleZMomentumDirection >= 0.)
			{
			  runAction -> TransmittedParticleNumber();
			}
		      else
			{
			  runAction -> BackscatteredParticleNumber(); 
			}
		    
		      if ( (primaryParticleEnergy == ParticleKineticEnergy) &&( particleZMomentumDirection == 1.))
			{
			  runAction -> TransmittedGammaNumber();
			}
		    }
		}
	    }
	}
    }

  // Stopping Power test

  if (trasmissionTestFlag == false)
    {
      if (primaryParticleName == "e-" || 
	  primaryParticleName == "proton" || 
	  primaryParticleName == "e+")
	{
	  if(0 == Step->GetTrack()->GetParentID() ) 
	    {
	      if (IDnow != idOld)
		{ 
		  // The particle is a primary one
		  idOld = IDnow ; 
		  // Identify the first step of the particle entirely in the Target
		  G4double initialStepXPosition = Step->GetPreStepPoint()->GetPosition().x();
		  G4double initialStepYPosition = Step->GetPreStepPoint()->GetPosition().y();
		  G4double initialStepZPosition = Step->GetPreStepPoint()->GetPosition().z();
		  if (0 == initialStepXPosition && 0 == initialStepYPosition && 0 == initialStepZPosition)
		    { 
		      G4double energyLost = abs(Step->GetDeltaEnergy());
		      G4double stepLength = Step->GetTrack()-> GetStepLength();
		      if (stepLength != 0) 
			{
			  G4double totalStoppingPower = energyLost / stepLength;
			  G4double nistStoppingPower = totalStoppingPower / (detector->GetDensity());
			  analysis -> StoppingPower(runID,primaryParticleEnergy/MeV,nistStoppingPower/(MeV*(cm2/g)));
			}
		    }
		}
       
	      // CSDA range
	     
	      if (0.0 == ParticleKineticEnergy) 
		{  
		  G4double xend = Step->GetTrack()->GetPosition().x()/mm ;
		  G4double yend = Step->GetTrack()->GetPosition().y()/mm ;
		  G4double zend = Step->GetTrack()->GetPosition().z()/mm ;
		  G4double range = (sqrt(xend*xend + yend*yend + zend*zend)); 
    
		  analysis -> CSDARange(runID,primaryParticleEnergy/MeV,(range*(detector->GetDensity()))/(g/cm2));
		}
	    }
	}
    }
}


