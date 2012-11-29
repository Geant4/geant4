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
// Code developed by:
//  S.Larsson
//
//    ***********************************
//    *                                 *
//    *    PurgMagSteppingAction.cc     *
//    *                                 *
//    ***********************************
//
// $Id$
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "PurgMagSteppingAction.hh"
#include "PurgMagRunAction.hh"
#include "PurgMagDetectorConstruction.hh"

#include "G4SteppingManager.hh"
#include "G4VTouchable.hh"
#include "G4VPhysicalVolume.hh"

#ifdef G4ANALYSIS_USE
#include"PurgMagAnalysisManager.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PurgMagSteppingAction::PurgMagSteppingAction(PurgMagRunAction* run,PurgMagDetectorConstruction* det)
:PurgMagRun(run),Detector(det)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PurgMagSteppingAction::~PurgMagSteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifdef G4ANALYSIS_USE
void PurgMagSteppingAction::UserSteppingAction(const G4Step* aStep)
  
{ 
  //Collection at SSD in N-tuples. Electrons and photons separated
  //Prestep point in World, next volume MeasureVolume, process transportation
 PurgMagAnalysisManager* analysis = PurgMagAnalysisManager::getInstance();	
  if ((aStep->GetPreStepPoint()->GetPhysicalVolume() == Detector->GetWorld())&&
      (aStep->GetTrack()->GetNextVolume() == Detector->GetMeasureVolume())&&
      //(aStep->GetTrack()->GetMomentumDirection().z()>0.)&& // only particles with positive momentum
      (aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() == "Transportation"))
    {
      G4double gx, gy, gz, ge, gpx, gpy, gpz, ex, ey, ez, ee; 
      G4double epx, epy, epz, px, py, pz, pe, ppx, ppy, ppz;
   
		 
      // Electrons 
      if(aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName()
      == "e-")
	{//Position 
         ex = (aStep->GetTrack()->GetPosition().x())/cm;
         ey = (aStep->GetTrack()->GetPosition().y())/cm; 
         ez = (aStep->GetTrack()->GetPosition().z())/cm;
	 // Energy
	 ee = (aStep->GetTrack()->GetKineticEnergy())/MeV;
	 // Momentum
	 epx = aStep->GetTrack()->GetMomentum().x();
	 epy = aStep->GetTrack()->GetMomentum().y();
	 epz = aStep->GetTrack()->GetMomentum().z();
	  
	 // Fill N-tuple electrons
	 analysis->fill_Tuple_Electrons(ex, ey, ez, ee, epx, epy, epz);
      }

      // Photons
      if (aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->
          GetParticleName() == "gamma")
	{
	  
	  // Position
          gx = (aStep->GetTrack()->GetPosition().x())/cm;
	  gy = (aStep->GetTrack()->GetPosition().y())/cm;
	  gz = (aStep->GetTrack()->GetPosition().z())/cm;
	  
	  // Energy
	  ge = (aStep->GetTrack()->GetKineticEnergy())/MeV;
	  
	  // Momentum
	  gpx = aStep->GetTrack()->GetMomentum().x();
	  gpy = aStep->GetTrack()->GetMomentum().y();
	  gpz = aStep->GetTrack()->GetMomentum().z();

	  // Fill N-tuple photons
	  analysis->fill_Tuple_Gamma(gx, gy, gz, ge, gpx, gpy, gpz);
	}


      // Positrons
      if (aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName() == "e+")
	{

	  // Position
          px = (aStep->GetTrack()->GetPosition().x())/cm;
	  py = (aStep->GetTrack()->GetPosition().y())/cm;
	  pz = (aStep->GetTrack()->GetPosition().z())/cm;
	  
	  // Energy
	  pe = (aStep->GetTrack()->GetKineticEnergy())/MeV;
	  
	  // Momentum
	  ppx = aStep->GetTrack()->GetMomentum().x();
	  ppy = aStep->GetTrack()->GetMomentum().y();
	  ppz = aStep->GetTrack()->GetMomentum().z();
	  
	  // Fill Ntuple positrons
	  analysis->fill_Tuple_Positrons(px, py, pz, pe, ppx, ppy, ppz);
	}
    }
#else
void PurgMagSteppingAction::UserSteppingAction(const G4Step* )
{ 
#endif
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
