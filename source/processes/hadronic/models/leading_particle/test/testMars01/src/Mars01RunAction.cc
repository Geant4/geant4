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
// $Id: Mars01RunAction.cc,v 1.1 2001-12-13 14:58:48 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4Timer.hh"
#include "Mars01RunAction.hh"

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"

Mars01RunAction::Mars01RunAction()
{
  timer = new G4Timer;
}

Mars01RunAction::~Mars01RunAction()
{
  delete timer;
}

void Mars01RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  timer->Start();
}

void Mars01RunAction::EndOfRunAction(const G4Run* aRun)
{
  timer->Stop();
  G4cout << "number of event = " << aRun->GetNumberOfEvent() << G4endl;
  G4cout << " " << *timer << G4endl;
}

