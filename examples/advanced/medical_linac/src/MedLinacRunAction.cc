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
// $Id: MedLinacRunAction.cc,v 1.6 2006/06/29 16:04:47 gunter Exp $
//
//
// Code developed by: M. Piergentili

#include "MedLinacRunAction.hh"
#include "MedLinacEventAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "MedLinacDetectorConstruction.hh"
#include "G4SDManager.hh"
#include "G4Timer.hh"

#ifdef G4ANALYSIS_USE
#include "MedLinacAnalysisManager.hh"
#endif


MedLinacRunAction::MedLinacRunAction()
{
}


MedLinacRunAction::~MedLinacRunAction()
{ 
}

void MedLinacRunAction::BeginOfRunAction(const G4Run* aRun)
{
#ifdef G4ANALYSIS_USE
  MedLinacAnalysisManager* analysis = MedLinacAnalysisManager::getInstance();
  analysis->book();
#endif


  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

}

void MedLinacRunAction::EndOfRunAction(const G4Run* aRun)
{
  
#ifdef G4ANALYSIS_USE
  MedLinacAnalysisManager* analysis = MedLinacAnalysisManager::getInstance();
#endif
  G4cout << "number of event = " << aRun->GetNumberOfEvent() << G4endl;
  
#ifdef G4ANALYSIS_USE      
  analysis->finish();
#endif


}



