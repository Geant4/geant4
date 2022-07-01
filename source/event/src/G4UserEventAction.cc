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
// G4UserEventAction class implementation
//
// Author: Makoto Asai (SLAC)
// --------------------------------------------------------------------

#include "G4UserEventAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ParticleTable.hh"
#include "globals.hh"

G4UserEventAction::G4UserEventAction()
{
 if(!(G4ParticleTable::GetParticleTable()->GetReadiness()))
 {
   G4String msg;
   msg =  " You are instantiating G4UserEventAction BEFORE your\n";
   msg += "G4VUserPhysicsList is instantiated and assigned to G4RunManager.\n";
   msg += " Such an instantiation is prohibited by Geant4 version 8.0. To fix this problem,\n";
   msg += "please make sure that your main() instantiates G4VUserPhysicsList AND\n";
   msg += "set it to G4RunManager before instantiating other user action classes\n";
   msg += "such as G4UserEventAction.";
   G4Exception("G4UserEventAction::G4UserEventAction()",
              "Event0032",FatalException,msg);
 }
}

void G4UserEventAction::SetEventManager(G4EventManager* value)
{
  fpEventManager = value;
}

void G4UserEventAction::BeginOfEventAction(const G4Event*)
{;}

void G4UserEventAction::EndOfEventAction(const G4Event*)
{;}

