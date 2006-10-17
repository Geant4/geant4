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
//    **************************************
//    *                                    *
//    *          CellRunAction.cc          *
//    *                                    *
//    **************************************
//
// Author: Susanna Guatelli (guatelli@ge.infn.it)
//	   Barbara Mascialino (Barbara.Mascialino@ge.infn.it)
//
// History:
// -----------
// 20 September 2006   S. Guatelli, B. Mascialino   1st implementation
//
// -------------------------------------------------------------------
 
#include "G4ios.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "CellDetectorConstruction.hh"
#include "CellRunAction.hh"
#include "CellSurvival.hh"

CellRunAction::CellRunAction()
{
  totalEnergyDeposit = 0;
}

CellRunAction::~CellRunAction()
{

}

void CellRunAction::BeginOfRunAction(const G4Run*)
{  
  // This method allows to do "operations" at the beginning of the run
 totalEnergyDeposit = 0;
}

void CellRunAction::IntegrateEnergyDeposit(G4double energy)
{
  totalEnergyDeposit += energy;
}

void CellRunAction::EndOfRunAction(const G4Run* aRun)
{
  G4cout << "Total energy deposit of the Run: " << totalEnergyDeposit/MeV
         << " MeV" << G4endl;

  CellDetectorConstruction* detector = 
(CellDetectorConstruction*)G4RunManager::GetRunManager()->GetUserDetectorConstruction();
  
  G4double mass =detector-> GetTargetMass(); 

  G4double dose = totalEnergyDeposit/mass;

  G4cout << " The dose in the target is: " << dose/gray << " gray"<< G4endl;

  CellSurvival* survival = new CellSurvival();

  survival -> SurvivalFormula(dose);
  
  G4double result = survival -> GetSurvival();

  G4cout<< "Survival probability = " << result << G4endl;

  delete survival;

  // This method allows to do "operations" at the end of the run
  G4int numberEvents = aRun -> GetNumberOfEvent();
  G4cout << "Number of events of this run: "<< numberEvents << G4endl;
}
