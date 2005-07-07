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
// $Id: Tst51RunAction.cc,v 1.2 2005-07-07 07:33:43 pandola Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Susanna Guatelli (guatelli@ge.infn.it)
//
// History:
// -----------
// 17 May  2003   S. Guatelli   1st implementation
// 06 Jul  2005   L. Pandola    Moved booking from main() to here in order to 
//                              set the filename via macro
//
// -------------------------------------------------------------------
 
#include "G4ios.hh"
#include <cmath>
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "Tst51RunAction.hh"
#include "Tst51AnalysisManager.hh"
#include "Tst51RunActionMessenger.hh"

Tst51RunAction::Tst51RunAction()
{
  filename = "test51.hbk";
  theMessenger = new Tst51RunActionMessenger(this);
}

Tst51RunAction::~Tst51RunAction()
{
  delete theMessenger;
}

void Tst51RunAction::BeginOfRunAction(const G4Run* aRun)
{      
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  if (G4VVisManager::GetConcreteInstance())
    {
      G4UImanager* UI = G4UImanager::GetUIpointer();
      UI->ApplyCommand("/vis/scene/notifyHandlers");
    } 
   if (aRun->GetRunID() == 0) //first run:book histos
     Tst51AnalysisManager::getInstance()->book(filename);
}

void Tst51RunAction::EndOfRunAction(const G4Run* aRun)
{
  if (G4VVisManager::GetConcreteInstance())
    {
      G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
    }

  G4double numberEvents = aRun->GetNumberOfEvent();
  G4cout << "Number of Events:"<< numberEvents << G4endl;
}

