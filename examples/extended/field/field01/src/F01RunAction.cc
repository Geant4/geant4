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
// $Id: F01RunAction.cc,v 1.6 2001-11-07 16:36:31 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 


#include "F01RunAction.hh"
#include "F01RunMessenger.hh"

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "g4std/iomanip"

#include "Randomize.hh"

//////////////////////////////////////////////////////////////////////////////

F01RunAction::F01RunAction()
  : saveRndm(0)
{
  runMessenger = new F01RunMessenger(this);
}

////////////////////////////////////////////////////////////////////////////

F01RunAction::~F01RunAction()
{
  delete runMessenger;
}

/////////////////////////////////////////////////////////////////////////////

void F01RunAction::BeginOfRunAction(const G4Run* aRun)
{  
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  
  // save Rndm status
  if (saveRndm > 0)
  { 
      HepRandom::showEngineStatus();
      HepRandom::saveEngineStatus("beginOfRun.rndm");
  }  
  G4UImanager* UI = G4UImanager::GetUIpointer();
   
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

  if(pVVisManager)    UI->ApplyCommand("/vis/scene/notifyHandlers");
}

/////////////////////////////////////////////////////////////////////////////

void F01RunAction::EndOfRunAction(const G4Run* aRun)
{
  if (G4VVisManager::GetConcreteInstance())
  {
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
  }

  // save Rndm status

  if (saveRndm == 1)
  { 
    HepRandom::showEngineStatus();
    HepRandom::saveEngineStatus("endOfRun.rndm");
  }     
}

//
//
////////////////////////////////////////////////////////////////////////
