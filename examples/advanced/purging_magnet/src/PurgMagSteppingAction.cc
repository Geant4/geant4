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
// Code developed by:
//  S.Larsson
//
//    ***********************************
//    *                                 *
//    *    PurgMagSteppingAction.cc     *
//    *                                 *
//    ***********************************
//
// $Id: PurgMagSteppingAction.cc,v 1.4 2004/12/07 13:45:02 guatelli Exp $
// GEANT4 tag $Name: geant4-08-00 $
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

void PurgMagSteppingAction::UserSteppingAction(const G4Step* aStep)
  
{ 
#ifdef G4ANALYSIS_USE
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
#endif
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
