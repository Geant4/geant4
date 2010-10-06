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
// -------------------------------------------------------------------
// $Id: SteppingAction.cc,v 1.4 2010-10-06 13:21:06 sincerti Exp $
// -------------------------------------------------------------------

#include "SteppingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SteppingAction::SteppingAction(RunAction* run,DetectorConstruction* det,PrimaryGeneratorAction* pri, HistoManager* his)
:Run(run),Detector(det),Primary(pri),Histo(his)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SteppingAction::UserSteppingAction(const G4Step* s)
  
{ 

if (Detector->GetCoef()==1) 
{

  if  ( (s->GetTrack()->GetDynamicParticle()->GetDefinition() == 
       G4Proton::ProtonDefinition())

/*
// for doublet

	 && (s->GetPostStepPoint()->GetPosition().z()/mm>-3230.2)
         && (s->GetPostStepPoint()->GetPosition().z()/mm<-3229.8) 
*/

// for triplet and whole line

	 && (s->GetPostStepPoint()->GetPosition().z()/mm>249.99999)
         && (s->GetPostStepPoint()->GetPosition().z()/mm<250.00001) 
         && (s->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "Vol")
         && (s->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "World")
     )
		
     {
     	      xIn = s->GetPostStepPoint()->GetPosition().x();
	      yIn = s->GetPostStepPoint()->GetPosition().y();
	      zIn = s->GetPostStepPoint()->GetPosition().z();
	      E   = s->GetTrack()->GetKineticEnergy();

              G4ThreeVector angleIn;
              angleIn = s->GetTrack()->GetMomentumDirection();

              thetaIn = std::asin(angleIn[0]/std::sqrt(angleIn[0]*angleIn[0]+angleIn[1]*angleIn[1]+angleIn[2]*angleIn[2]));
              phiIn = std::asin(angleIn[1]/std::sqrt(angleIn[0]*angleIn[0]+angleIn[1]*angleIn[1]+angleIn[2]*angleIn[2]));

              G4cout << "    =>IMAGE : X(microns)=" << xIn/micrometer <<" Y(microns)="<< yIn/micrometer << " THETA(mrad)=" << (thetaIn/mrad) << " PHI(mrad)=" << (phiIn/mrad) << G4endl;
	      G4cout << G4endl;

              Run->AddRow();
              Run->AddToXVector(xIn/um);
              Run->AddToYVector(yIn/um);
              Run->AddToThetaVector(thetaIn/mrad);
              Run->AddToPhiVector(phiIn/mrad);

	      Histo->FillNtuple(2, 0, xIn/um);
	      Histo->FillNtuple(2, 1, yIn/um);
	      Histo->FillNtuple(2, 2, thetaIn/mrad);
	      Histo->FillNtuple(2, 3, phiIn/mrad);
	      Histo->AddRowNtuple(2);      

     }
}

if (Detector->GetProfile()==1) 
{

	if  (
	    (s->GetTrack()->GetDynamicParticle()->GetDefinition() == G4Proton::ProtonDefinition())
         && (s->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "Vol")
         && (s->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "Vol") )
	{
     	      xIn = s->GetPostStepPoint()->GetPosition().x();
	      yIn = s->GetPostStepPoint()->GetPosition().y();
	      zIn = s->GetPostStepPoint()->GetPosition().z();
	      
	      Histo->FillNtuple(0, 0, xIn/um);
	      Histo->FillNtuple(0, 1, yIn/um);
	      Histo->FillNtuple(0, 2, zIn/mm);
	      Histo->AddRowNtuple(0);      
	}
}
	
if (Detector->GetGrid()==1) 
{

	if  (
	    (s->GetTrack()->GetDynamicParticle()->GetDefinition() == G4Proton::ProtonDefinition())
         && (s->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "ControlVol_GridShadow")
         && (s->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "World") )
	{
     	      xIn = s->GetPostStepPoint()->GetPosition().x();
	      yIn = s->GetPostStepPoint()->GetPosition().y();
              E   = s->GetTrack()->GetKineticEnergy();

	      Histo->FillNtuple(1, 0, xIn/um);
	      Histo->FillNtuple(1, 1, yIn/um);
	      Histo->FillNtuple(1, 2, E/MeV);
	      Histo->AddRowNtuple(1);
	}
}

// end
}     
