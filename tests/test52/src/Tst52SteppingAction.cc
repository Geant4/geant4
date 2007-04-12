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
// $Id: Tst52SteppingAction.cc,v 1.1 2007-04-12 12:00:17 guatelli Exp $
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
#include "Tst52SteppingAction.hh"
#include "Tst52DetectorConstruction.hh"
#ifdef G4ANALYSIS_USE
#include "Tst52AnalysisManager.hh"
#endif
#include "Tst52PrimaryGeneratorAction.hh"
#include "Tst52RunAction.hh"

Tst52SteppingAction::Tst52SteppingAction(Tst52PrimaryGeneratorAction* primary, 
					 Tst52RunAction* run, 
					 Tst52DetectorConstruction* det):
  primaryAction(primary), 
  runAction(run), 
  detector(det) 
{ }

Tst52SteppingAction::~Tst52SteppingAction()
{ }

void Tst52SteppingAction::UserSteppingAction(const G4Step* aStep )
{
  // Analysis of the secondaries generated in the phantom
#ifdef G4ANALYSIS_USE
  Tst52AnalysisManager* analysis = Tst52AnalysisManager::getInstance();
 
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

       //  G4cout <<"process: " << process <<" to generate particle: " << secondaryParticleName << G4endl;

       if (volumeName == "Target")
	 {
	   G4cout <<"process: " << process <<" to generate particle: " << secondaryParticleName << "with energy: " <<secondaryParticleKineticEnergy/MeV<< G4endl;
            
	   if (secondaryParticleName == "e-") 
	     {
	       if ((process == "eIoni") || (process ==  "LowEnergyIoni") || (process ==  "PenelopeIoni")) 
		 analysis -> secondaryElectron(secondaryParticleKineticEnergy/MeV);
	       else
		 {
		   if ((process == "LowEnPhotoElec") || (process == "phot") || (process == "PenPhotoElec"))
		     analysis -> secondaryElectronPhoto(secondaryParticleKineticEnergy/MeV);
		   else
		     {
		       if ((process != "eIoni") && (process != "LowEnergyIoni")  && (process != "PenelopeIoni")
			   && (process != "LowEnPhotoElec") && (process != "phot") && (process != "PenPhotoElec"))
			 G4cout<< secondaryParticleName << "is generated with "<< process <<"!!!!!" << G4endl;
		     }
		 }
	       //   G4cout<<secondaryParticleName<< " " << secondaryParticleKineticEnergy/MeV << G4endl; 
	     }
	 
	   if (secondaryParticleName == "gamma") 
	     {
	       analysis -> secondaryPhoton(secondaryParticleKineticEnergy/MeV);
	       // G4cout<<secondaryParticleName<< " " << process << G4endl;
	     }
	 }
     }

#endif
}




