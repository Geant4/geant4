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
// ********************************************************************
// *                                                                  *
// * cosmicray_charging advanced example for Geant4                   *
// * (adapted simulation of test-mass charging in the LISA mission)   *
// *                                                                  *
// * Henrique Araujo (h.araujo@imperial.ac.uk) & Peter Wass           *
// * Imperial College London                                          *
// *                                                                  *
// * LISASteppingAction class                                         *
// *                                                                  *
// ********************************************************************
//
// HISTORY
// 22/02/2004: migrated from LISA-V04
//
// ********************************************************************


#include "LISASteppingAction.hh"

#include "LISASteppingActionMessenger.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4StepPoint.hh"
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleTypes.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <fstream>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
LISASteppingAction::LISASteppingAction()  {

  charge_in[0]  = 0;
  charge_in[1]  = 0;
  charge_out[0] = 0;
  charge_out[1] = 0;

  // defaults
  FlagSpectrum = false;

  // create messenger
  steppingMessenger = new LISASteppingActionMessenger(this);

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
LISASteppingAction::~LISASteppingAction() {

  delete steppingMessenger;

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void LISASteppingAction::UserSteppingAction(const G4Step* fStep) {

  if(fStep->GetTrack()->GetNextVolume()) {

    // step points
    G4TouchableHistory* thePreTouchable = 
      (G4TouchableHistory*)(fStep->GetPreStepPoint()->GetTouchable());
    G4TouchableHistory* thePostTouchable = 
      (G4TouchableHistory*)(fStep->GetPostStepPoint()->GetTouchable());
    G4String PreVol = thePreTouchable->GetVolume()->GetName();
    G4String PosVol = thePostTouchable->GetVolume()->GetName();

    // particle
    G4double charge = fStep->GetTrack()->GetDefinition()->GetPDGCharge();
    G4String pName  = fStep->GetTrack()->GetDefinition()->GetParticleName();
    G4double energy = fStep->GetTrack()->GetKineticEnergy();

    G4int tm = -1;


    // *** entering a test mass ***
    if(PreVol!="tmass_phys" && PosVol=="tmass_phys") {

      tm = thePostTouchable->GetReplicaNumber(6);
      
      charge_in[tm] += (int)charge;
      
      // energy spectra at test-mass boundary
      if(FlagSpectrum) {
	// electrons entering test mass
	if(pName=="e-") {
	  std::ofstream electrons_in("electrons_in.out",std::ios::app);
	  electrons_in << energy/eV << G4endl; 
	}
	// hadrons entering test mass
	else if(pName=="proton" || 
		pName=="pi+" || pName=="pi-" ||
		pName=="deuteron" || pName=="He3" || 
		pName=="alpha" || pName=="triton") {
	  std::ofstream hadrons_in("hadrons_in.out",std::ios::app);
	  hadrons_in << energy/eV << G4endl;
	}
      }

    } // *** entering a test mass ***
 

    // *** leaving a test mass ***
    else if(PreVol!=PosVol && PreVol=="tmass_phys") {

      tm = thePreTouchable->GetReplicaNumber(6);
      
      charge_out[tm] += (int)charge;
      
      // energy spectra at test-mass boundary
      if(FlagSpectrum) {
	// electrons leaving test mass 
	if(pName=="e-") {
	  std::ofstream electrons_out("electrons_out.out",std::ios::app);
	  electrons_out << energy/eV << G4endl; 
	}
	// hadrons leaving test mass
	else if(pName=="proton" || 
		pName=="pi+" || pName=="pi-" ||
		pName=="deuteron" || pName=="He3" || 
		pName=="alpha" || pName=="triton") {
	  std::ofstream hadrons_out("hadrons_out.out",std::ios::app);
	  hadrons_out << energy/eV << G4endl; 
	}
      } 

    } // *** leaving a test mass ***


  } // if(next volume)


  // count particles hitting probe volume
  //   if(fStep->GetTrack()->GetNextVolume()) {
  //     static int i=0;
  //     G4String PosVol=fStep->GetPostStepPoint()
  //       ->GetPhysicalVolume()->GetName();
  //     if(PosVol=="probe_phys") {
  //       fStep->GetTrack()->SetTrackStatus(fStopAndKill);
  //       G4cout << i++ << G4endl;
  //     }
  //   }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

