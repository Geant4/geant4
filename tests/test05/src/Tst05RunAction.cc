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
// $Id: Tst05RunAction.cc,v 1.7 2001-07-11 10:09:39 gunter Exp $
// ------------------------------------------------------------

// Make this appear first!
#include "G4Timer.hh"

#include "Tst05RunAction.hh"

#include "G4Run.hh"
#include "G4UImanager.hh"

Tst05RunAction::Tst05RunAction()
{
  timer = new G4Timer;
}

Tst05RunAction::~Tst05RunAction()
{
  delete timer;
}

void Tst05RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4UImanager* UI = G4UImanager::GetUIpointer();
  UI->ApplyCommand("/event/Verbose 1");
  UI->ApplyCommand("/tracking/Verbose 1");

  //G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  timer->Start();
}

void Tst05RunAction::EndOfRunAction(const G4Run* aRun)
{
  timer->Stop();
  //G4cout << "number of event = " << aRun->GetNumberOfEvent() 
  //     << " " << *timer << G4endl;
}

