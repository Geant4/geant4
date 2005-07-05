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
// $Id: Tst51SteppingAction.cc,v 1.1 2005-07-05 11:06:27 guatelli Exp $
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
#include "G4VParticleChange.hh"
#include "G4ParticleDefinition.hh"
#include "G4VPhysicalVolume.hh"
#include "Tst51SteppingAction.hh"
#include "Tst51DetectorConstruction.hh"
#ifdef G4ANALYSIS_USE
#include "Tst51AnalysisManager.hh"
#endif

Tst51SteppingAction::Tst51SteppingAction(): preliminaryTest(false)
 { }

Tst51SteppingAction::~Tst51SteppingAction()
{ }

void Tst51SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  Tst51AnalysisManager* analysis = Tst51AnalysisManager::getInstance();
  
  // ---------- seond test on bremmstrahlung -----------// 
  if (preliminaryTest == false)
    {
  if ((aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() 
       == "Target")  && (aStep->GetPostStepPoint()->GetPhysicalVolume()->GetName() 
       == "World"))
    {
      G4String name = aStep -> GetTrack()-> GetDynamicParticle() 
                            -> GetDefinition() -> GetParticleName();
      if (name == "gamma")
	{

	  G4double theta =  ( aStep -> GetPreStepPoint() -> GetMomentum()).theta();
	  G4double kineticEnergy = aStep -> GetPreStepPoint() -> GetKineticEnergy();
	 
	      if (theta < 90 * deg)
		{
		  analysis -> angularDistributionTransmittedGamma(theta/deg);
		  analysis -> energyDistributionTransmittedGamma(kineticEnergy/keV);

		}
	      if (theta > 90 * deg)
		{
		  analysis -> angularDistributionBackGamma(theta/deg);
		  analysis -> energyDistributionBackGamma(kineticEnergy/keV); 
		}
	      analysis-> fillNtuple(kineticEnergy/keV, theta/deg);
	     }	
    }
    }

  // preliminary test: energy and angular distribution of photons generated in the target
  // by 70 kev electrons
  if (preliminaryTest == true)
    {
      // Electrons are generated inside the slab 
  G4String name = aStep -> GetTrack()-> GetDynamicParticle() ->
  GetDefinition() -> GetParticleName();
   if (name == "e-") 
  { 
    G4int stepNumber = aStep->GetTrack()->GetCurrentStepNumber(); 
  
    if (stepNumber == 1) 
      {                  
        // After the first step of the electrons in the slab, the information
        //  of secondary photons is retrieved
        G4SteppingManager * SM = fpSteppingManager; 
        G4TrackVector* fSecondary = SM -> GetfSecondary();
       
        for(size_t lp1=0;lp1<(*fSecondary).size(); lp1++) 
          { 
        	
          G4String name = (*fSecondary)[lp1]->GetDefinition() -> GetParticleName();
          if (name == "gamma")
	    {          
          // theta is ZMoD/MMod
          G4double theta =  ( (*fSecondary)[lp1] -> GetMomentum()).theta();
	  G4double energy = (*fSecondary)[lp1]->GetKineticEnergy(); 

          analysis -> angularDistributionPostStep(theta/deg); 
          analysis -> energyDistributionPostStep(energy/keV); 
	  }
	  }
      }
  }
    }  
}



