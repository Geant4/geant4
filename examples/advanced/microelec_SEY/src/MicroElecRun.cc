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
/// \file MicroElecRun.cc
/// \brief Implementation of the MicroElecRun class
//
//

#include "MicroElecHitSey.hh"
#include "MicroElecRun.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MicroElecRun::MicroElecRun()
 : G4Run()
{
	ElecPrimScorer	= 0.0;
	ElecSecoScorer	= 0.0;
	ElecSup50Scorer = 0.0;
	ElecTotaScorer	= 0.0;
	ElecEneIncPart = 0.0;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MicroElecRun::~MicroElecRun()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//  RecordEvent is called at end of event.
//  For scoring purpose, the resultant quantity in a event,
//  is accumulated during a MicroElecRun.
void MicroElecRun::RecordEvent(const G4Event* evt)
{
	G4Run::RecordEvent(evt);
	G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  
  
  if (!HCE)
	  return;
  else {
	  G4int NbOfCollection = HCE->GetNumberOfCollections();
	  G4int NbOfHitsinCollection;
	  
	  for (G4int i = 0; i < NbOfCollection; i++) {
		  
		  NbOfHitsinCollection = HCE->GetHC(i)->GetSize();
		  for (G4int j = 0; j< NbOfHitsinCollection;j++)
		  {
			  MicroElecHitSey* localHit =  (MicroElecHitSey*) (HCE->GetHC(i)->GetHit(j));
			  
			  if (localHit->GetParentID()==0){ 
				  
				SetElecEneIncPart(localHit->GetVertexKineticEnergy());
			  }
			  
		  }

	  }

  }
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void MicroElecRun::Merge(const G4Run * aRun) {
  const MicroElecRun * localRun = static_cast<const MicroElecRun *>(aRun);
  
  ElecPrimScorer += localRun->ElecPrimScorer;
  ElecSecoScorer += localRun->ElecSecoScorer;
  ElecSup50Scorer += localRun->ElecSup50Scorer;
  ElecTotaScorer += localRun->ElecTotaScorer;
  ElecEneIncPart= localRun->ElecEneIncPart;

  G4cout << "-------------------- MicroElec Run Merge Events : begin--------------------------" << G4endl;
  G4cout << "Vertex Prim NRJ = " << GetElecEneIncPart() << G4endl;
  G4cout << "Primaries = " << GetElecPrimScorer() << G4endl;
  G4cout << "Secondaries = " << GetElecSecoScorer() << G4endl;
  G4cout << "Sec > 50 eV = " << GetElecSup50Scorer() << G4endl;
  G4cout << "Total = " << GetElecTotaScorer() << G4endl;
  G4cout << "-------------------- MicroElec Run Merge Events : end--------------------------" << G4endl;

  G4Run::Merge(aRun);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
