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
// $Id: Tst26SteppingVerbose.cc,v 1.2 2003-03-26 17:29:30 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Tst26SteppingVerbose.hh"

#include "G4MaterialCutsCouple.hh"
#include "G4SteppingManager.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Tst26SteppingVerbose::Tst26SteppingVerbose():G4SteppingVerbose()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Tst26SteppingVerbose::~Tst26SteppingVerbose()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst26SteppingVerbose::StepInfo()
{
  CopyState();
  G4int prec = G4cout.precision(3);
  if( verboseLevel >= 1 ) {

    if( verboseLevel >= 4 ) VerboseTrack();

    PrintStep(false);

    if(fStep->GetPostStepPoint()->GetProcessDefinedStep() != NULL){
      G4cout << fStep->GetPostStepPoint()->GetProcessDefinedStep()
	                                ->GetProcessName();
    } else {
      G4cout << "UserLimit";
    }

    G4cout << G4endl;

    if( verboseLevel >= 2 ){
      G4int tN2ndariesTot = fN2ndariesAtRestDoIt +
	                    fN2ndariesAlongStepDoIt +
	                    fN2ndariesPostStepDoIt;
      if(tN2ndariesTot>0){
	G4cout << "    :----- List of 2ndaries - "
	       << "#SpawnInStep=" << G4std::setw(3) << tN2ndariesTot
	       << "(Rest="  << G4std::setw(2) << fN2ndariesAtRestDoIt
	       << ",Along=" << G4std::setw(2) << fN2ndariesAlongStepDoIt
	       << ",Post="  << G4std::setw(2) << fN2ndariesPostStepDoIt
	       << "), "
	       << "#SpawnTotal=" << G4std::setw(3) << (*fSecondary).size()
	       << " ---------------"
	       << G4endl;

	for(size_t lp1=(*fSecondary).size()-tN2ndariesTot;
                        lp1<(*fSecondary).size(); lp1++){
	  G4cout << "    : "
		 << G4std::setw(6)
		 << G4BestUnit((*fSecondary)[lp1]->GetPosition().x(),"Length")
		 << G4std::setw(6)
		 << G4BestUnit((*fSecondary)[lp1]->GetPosition().y(),"Length")
		 << G4std::setw(6)
		 << G4BestUnit((*fSecondary)[lp1]->GetPosition().z(),"Length")
		 << G4std::setw(6)
		 << G4BestUnit((*fSecondary)[lp1]->GetKineticEnergy(),"Energy")
		 << G4std::setw(10)
		 << (*fSecondary)[lp1]->GetDefinition()->GetParticleName();
	  G4cout << G4endl;
	}

	G4cout << "    :-----------------------------"
	       << "----------------------------------"
	       << "-- EndOf2ndaries Info ---------------"
	       << G4endl;
      }
    }

    G4cout.precision(prec);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst26SteppingVerbose::TrackingStarted()
{
  CopyState();

  G4int prec = G4cout.precision(3);
  if( verboseLevel > 0 ){

    PrintStep(true);
    G4cout  << "initStep" << G4endl;
  }
  G4cout.precision(prec);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst26SteppingVerbose::PrintStep(G4bool firstStep)
{
  if( firstStep || verboseLevel >= 3 )
    {
      G4cout << G4std::setw( 5) << "#Step#"      << " "
             << G4std::setw( 8) << "X(mm)"      << " " 
	     << G4std::setw( 8) << "Y(mm)"      << " "
             << G4std::setw( 8) << "Z(mm)"      << " "
             << G4std::setw( 9) << "KinE(MeV)"  << " " 
	     << G4std::setw( 8) << "dE(keV)"    << " "
             << G4std::setw( 8) << "StepL(mm)"   << " "
	     << G4std::setw(10) << "TrackL(cm)"  << " "
             << G4std::setw( 9) << "idxCouple" << " "
             << G4std::setw( 9) << "Material" << " "
             << G4std::setw(10) << "NextVolume" << " "
	     << G4std::setw( 8) << "ProcName"   << G4endl;
    }
    G4cout << G4std::setw( 5) << fTrack->GetCurrentStepNumber() << " "
                 << G4std::setw( 8) << fTrack->GetPosition().x()/mm  << " "
                 << G4std::setw( 8) << fTrack->GetPosition().y()/mm  << " "
                 << G4std::setw( 8) << fTrack->GetPosition().z()/mm  << " "
                 << G4std::setw( 9) << fTrack->GetKineticEnergy()/MeV << " "
                 << G4std::setw( 8) << fStep->GetTotalEnergyDeposit()/keV << " "
                 << G4std::setw( 8) << fStep->GetStepLength()/mm << " "
                 << G4std::setw(10) << fTrack->GetTrackLength()/cm << "     ";

  if( firstStep ) {
    G4cout << "                    ";
  } else {
    const G4MaterialCutsCouple* couple = fTrack->GetMaterialCutsCouple();
    G4String mate = couple->GetMaterial()->GetName();
    G4int idx = couple->GetIndex();
    G4cout << G4std::setw(5) << idx << "      " << mate << "       ";
  }

  if( fTrack->GetNextVolume() != 0 ) {
      G4cout << G4std::setw(10) << fTrack->GetVolume()->GetName() << "    ";
  } else {
      G4cout << G4std::setw(10) << "OutOfWorld    ";
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

