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
// Authors: S. Guatelli and M. G. Pia, INFN Genova, Italy
// 
// Based on code developed by the undergraduate student G. Guerrieri 
// Note: this is a preliminary beta-version of the code; an improved 
// version will be distributed in the next Geant4 public release, compliant
// with the design in a forthcoming publication, and subject to a 
// design and code review.
//

#include "G4HumanPhantomRunAction.hh"
#include "G4HumanPhantomAnalysisManager.hh"
#include "G4ios.hh"
#include "G4Run.hh"

G4HumanPhantomRunAction::G4HumanPhantomRunAction()
{}

G4HumanPhantomRunAction::~G4HumanPhantomRunAction()
{}

void G4HumanPhantomRunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4int run_number = aRun->GetRunID();
  G4cout << "### Run " << run_number << " start." << G4endl;

#ifdef G4ANALYSIS_USE
  if (run_number == 0)
    {
  G4HumanPhantomAnalysisManager* analysis = G4HumanPhantomAnalysisManager::getInstance();
  analysis->book();
    }
#endif
}

void G4HumanPhantomRunAction::EndOfRunAction(const G4Run* aRun)
{
  G4cout << "Number of events = " << aRun->GetNumberOfEvent() << G4endl;
}



