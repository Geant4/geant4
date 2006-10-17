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
//    **************************************
//    *                                    *
//    *        CellSteppingAction.cc       *
//    *                                    *
//    **************************************
//
// Author: Susanna Guatelli (guatelli@ge.infn.it)
//	   Barbara Mascialino (Barbara.Mascialino@ge.infn.it)
//
// History:
// -----------
//20 September 2006   S. Guatelli, B. Mascialino   1st implementation
//
// -------------------------------------------------------------------
#include "G4ios.hh"
#include <cmath> 
#include "G4SteppingManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4VPhysicalVolume.hh"
#include "CellSteppingAction.hh"
#include "CellDetectorConstruction.hh"
#ifdef G4ANALYSIS_USE
#include "CellAnalysisManager.hh"
#endif
#include "CellPrimaryGeneratorAction.hh"
#include "CellRunAction.hh"

CellSteppingAction::CellSteppingAction(CellPrimaryGeneratorAction* primary, 
					 CellRunAction* run, 
					 CellDetectorConstruction* det):
  primaryAction(primary), 
  runAction(run), 
  detector(det) 
{ }

CellSteppingAction::~CellSteppingAction()
{ }

void CellSteppingAction::UserSteppingAction(const G4Step* aStep)
{
  // The stepping action allows to retrieve information concerning 
  // the steps of particles in the experimental set-up

  // The information are retrieve and stored in histograms contained in cell.hbk
  CellAnalysisManager* analysis = CellAnalysisManager::getInstance();
	     
  // Check if the particle is primary
  if (aStep ->GetTrack()->GetParentID() == 0)
    {
      // Retrieve the energy of the primary particles when they outgo the phantom
 
      // Check if the primary particle is going outside the target ...
      if ((aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() 
	   == "Target") && 
	  (aStep->GetTrack()->GetNextVolume()->GetName() == "World") )
	{ 
	  // Retrieve the energy
          G4double  particleKineticEnergy = aStep -> GetTrack()-> GetKineticEnergy(); 
     
	  // Store the information in the histogram
         analysis -> primary_energy_outgoing(particleKineticEnergy/MeV);
	}
    }

  // Analysis of the secondary particles generated in the target
  // Whenever a secondary particle is generated in the phantom, its information are retrieved here
  // The secondary particle can be produced by a primary particle or a secondary particle
 
 G4SteppingManager*  steppingManager = fpSteppingManager;
 G4Track* theTrack = aStep -> GetTrack();

 // check if it is alive
 if(theTrack -> GetTrackStatus() == fAlive) { return; }
 
 // Retrieve the vector with generated secondary particles
 G4TrackVector* fSecondary = steppingManager -> GetfSecondary();

for(size_t lp1=0;lp1<(*fSecondary).size(); lp1++)
     { 
       // Retrieve the name of the volume where the secondary particle is generated
       G4String volumeName = (*fSecondary)[lp1] -> GetVolume() -> GetName(); 

       // Retrieve the name of the secondary particle
       G4String secondaryParticleName =  (*fSecondary)[lp1]-> GetDefinition() -> GetParticleName();  

       // Retrieve the energy of the created secondary particle
       G4double secondaryParticleKineticEnergy =  (*fSecondary)[lp1] -> GetKineticEnergy(); 
     
       // Retrieve the process that generated the secondary particle
       G4String process = (*fSecondary)[lp1]-> GetCreatorProcess()-> GetProcessName();   
      
       // If the secondary particle is originated in the target ... 
       if (volumeName == "Target")
	 {
	   G4cout<< "A secondary " << secondaryParticleName << " was originated in "
                 << volumeName << " with energy "<<secondaryParticleKineticEnergy/MeV <<" MeV"
                 << " from " <<  process << G4endl;
	 }
     }
}

