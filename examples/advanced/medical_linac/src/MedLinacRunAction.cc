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
// $Id: MedLinacRunAction.cc,v 1.5 2005/11/25 22:02:04 mpiergen Exp $
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



