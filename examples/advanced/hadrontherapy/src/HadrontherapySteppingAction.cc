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
// $Id: HadrontherapyProtonSteppingAction.cc; 
//
// See more at: http://workgroup.lngs.infn.it/geant4lns
//
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, F. Di Rosa(a), S. Guatelli(b), G. Russo(a)
// 
// (a) Laboratori Nazionali del Sud 
//     of the INFN, Catania, Italy
// (b) INFN Section of Genova, Genova, Italy
// 
// * cirrone@lns.infn.it
// ----------------------------------------------------------------------------

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

#ifdef G4ANALYSIS_USE 	
#include "HadrontherapyAnalysisManager.hh"
#endif

#include "HadrontherapyRunAction.hh"

/////////////////////////////////////////////////////////////////////////////
HadrontherapySteppingAction::HadrontherapySteppingAction( HadrontherapyRunAction* run)
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
  
  // Electromagnetic and hadronic processes of primary particles in the phantom
  if ((aStep -> GetTrack() -> GetTrackID() == 1) &&
    (aStep -> GetTrack() -> GetVolume() -> GetName() == "PhantomPhys") &&
    (aStep -> GetPostStepPoint() -> GetProcessDefinedStep() != NULL))
	    {
	      G4String process = aStep -> GetPostStepPoint() -> 
		GetProcessDefinedStep() -> GetProcessName();
 
	      if ((process == "Transportation") || (process == "StepLimiter")) {;}
	      else {
		if ((process == "msc") || (process == "hLowEIoni") || (process == "hIoni")) 
		  { 
                    runAction -> AddEMProcess();
		  } 
		else  
                 {
		   runAction -> AddHadronicProcess();

                   if ( (process != "LElastic") && (process != "ProtonInelastic") && (process != "hElastic") )
		     G4cout << "Warning! Unknown proton process: "<< process << G4endl;
		 }
	      }         
	    }

 // Retrieve information about the secondaries originated in the phantom

#ifdef G4ANALYSIS_USE 	
 G4SteppingManager*  steppingManager = fpSteppingManager;
  
  // check if it is alive
  //if(theTrack-> GetTrackStatus() == fAlive) { return; }

  // Retrieve the secondary particles
  G4TrackVector* fSecondary = steppingManager -> GetfSecondary();
     
  for(size_t lp1=0;lp1<(*fSecondary).size(); lp1++)
    { 
      G4String volumeName = (*fSecondary)[lp1] -> GetVolume() -> GetName(); 
 
      if (volumeName == "PhantomPhys")
	{
	  G4String secondaryParticleName =  (*fSecondary)[lp1]->GetDefinition() -> GetParticleName();  
	  G4double secondaryParticleKineticEnergy =  (*fSecondary)[lp1] -> GetKineticEnergy();     
   
	  HadrontherapyAnalysisManager* analysis =  HadrontherapyAnalysisManager::getInstance();   
        
          if (secondaryParticleName == "e-")
	    analysis -> electronEnergyDistribution(secondaryParticleKineticEnergy/MeV);
	  	    
	  if (secondaryParticleName == "gamma")
	    analysis -> gammaEnergyDistribution(secondaryParticleKineticEnergy/MeV);

	  if (secondaryParticleName == "deuteron")
	    analysis -> deuteronEnergyDistribution(secondaryParticleKineticEnergy/MeV);
       
	  if (secondaryParticleName == "triton")
	    analysis -> tritonEnergyDistribution(secondaryParticleKineticEnergy/MeV);
           
	  if (secondaryParticleName == "alpha")
	    analysis -> alphaEnergyDistribution(secondaryParticleKineticEnergy/MeV);
	
	  G4double z = (*fSecondary)[lp1]-> GetDynamicParticle() -> GetDefinition() -> GetPDGCharge();
	  if (z > 0.)
            {      
	      G4int a = (*fSecondary)[lp1]-> GetDynamicParticle() -> GetDefinition() -> GetBaryonNumber();
	      G4int electronOccupancy = (*fSecondary)[lp1] ->  GetDynamicParticle() -> GetTotalOccupancy(); 
	      
	      // If a generic ion is originated in the detector, its baryonic number, PDG charge, 
	      // total number of electrons in the orbitals are stored in a ntuple 
	      analysis -> genericIonInformation(a, z, electronOccupancy, secondaryParticleKineticEnergy/MeV);
	    }
	}
    }
#endif
}





