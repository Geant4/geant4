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
// Rich advanced example for Geant4
// RichTbRunAction.cc for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#include "G4Timer.hh"

#include "RichTbRunAction.hh"

#include "G4ios.hh"
#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "RichTbAnalysisManager.hh"

#include "RichTbRunConfig.hh"

RichTbRunAction::RichTbRunAction(){ 

  timer = new G4Timer;

}

RichTbRunAction::~RichTbRunAction()
{
  delete timer;
}

void RichTbRunAction::BeginOfRunAction(const G4Run* aRun)
{

  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  if(G4VVisManager::GetConcreteInstance())
    {
      G4UImanager* UI = G4UImanager::GetUIpointer();
      UI->ApplyCommand("/vis/scene/notifyHandlers");
    }


  timer->Start();


#ifdef G4ANALYSIS_USE

  RichTbAnalysisManager * analysis = RichTbAnalysisManager::getInstance();

  analysis->book();

#endif

}

void RichTbRunAction::EndOfRunAction(const G4Run* aRun)
{
#ifdef G4ANALYSIS_USE

  RichTbAnalysisManager * analysis = RichTbAnalysisManager::getInstance();
#endif


  timer->Stop();
  G4cout << "number of event = " << aRun->GetNumberOfEvent() 
       << " " << *timer << G4endl;

#ifdef G4ANALYSIS_USE
  analysis->finish();
#endif

}

