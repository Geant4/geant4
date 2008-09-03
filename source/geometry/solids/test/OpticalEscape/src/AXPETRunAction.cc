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
// $Id: AXPETRunAction.cc,v 1.1 2008-09-03 13:34:03 gcosmo Exp $
// ------------------------------------------------------------
// Geant4 class implementation file
//
// 03/09/2008, by T.Nikitina
// ------------------------------------------------------------

#include "G4Timer.hh"

#include "G4RunManager.hh"
#include "AXPETRunAction.hh"

#include "G4Run.hh"


AXPETRunAction::AXPETRunAction(): 
energylyso(0), timelyso(-1), lyso(-1)
{
  timer = new G4Timer;
}


AXPETRunAction::~AXPETRunAction()
{
  delete timer;
}


void AXPETRunAction::BeginOfRunAction(const G4Run* aRun)
{
   G4RunManager::GetRunManager()->SetRandomNumberStore(true);
   G4cout << "Run " << aRun -> GetRunID() << " starts ..." << G4endl;
   runID = aRun->GetRunID();

   timer->Start();
}


void AXPETRunAction::EndOfRunAction(const G4Run* aRun)
{
  timer->Stop();
  G4cout << "number of event = " << aRun->GetNumberOfEvent() 
         << " " << *timer << G4endl;	 
}
