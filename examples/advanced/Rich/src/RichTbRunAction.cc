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

