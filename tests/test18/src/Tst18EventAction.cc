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
//  File:        Tst18EventAction.cc
//  Description: Event action for radioactive decay system test 
//  Author:	 Dennis Wright (SLAC)
//                 (original by F. Lei DERA UK)
//  Date:        14 August 2013
//

#include "G4ios.hh"
#include "Tst18EventActionMessenger.hh"
#include "Tst18EventAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4UImanager.hh"
#include "G4UnitsTable.hh"


Tst18EventAction::Tst18EventAction()
 : eventMessenger(0), numberOfSecondariesPerEvent(0)
{
  eventMessenger = new Tst18EventActionMessenger(this);
}

Tst18EventAction::~Tst18EventAction()
{
  delete eventMessenger;
}

void Tst18EventAction::BeginOfEventAction(const G4Event*)
{
  numberOfSecondariesPerEvent = 0;
}

void Tst18EventAction::EndOfEventAction(const G4Event*)
{
  if (numberOfSecondariesPerEvent < 2) {
    G4cout << " Abort " << G4endl;
    G4Exception("Tst18EventAction::EndOfEventAction", "Tst18_00",
                FatalException, "No Decay Products Found");
  }
}

void Tst18EventAction::IncrementParticleNumber()
{
  numberOfSecondariesPerEvent++;
}

