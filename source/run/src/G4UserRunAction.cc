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
//
// $Id: G4UserRunAction.cc 68997 2013-04-15 09:21:24Z gcosmo $
//

#include "G4UserRunAction.hh"
#include "G4ParticleTable.hh"
#include "globals.hh"


G4UserRunAction::G4UserRunAction()
:isMaster(true)
{
 if(!(G4ParticleTable::GetParticleTable()->GetReadiness()))
 {
   G4String msg;
   msg =  " You are instantiating G4UserRunAction BEFORE your G4VUserPhysicsList is\n";
   msg += "instantiated and assigned to G4RunManager.\n";
   msg += " Such an instantiation is prohibited by Geant4 version 8.0. To fix this problem,\n";
   msg += "please make sure that your main() instantiates G4VUserPhysicsList AND\n";
   msg += "set it to G4RunManager before instantiating other user action classes\n";
   msg += "such as G4UserRunAction.";
   G4Exception("G4UserRunAction::G4UserRunAction()",
              "Run0041",FatalException,msg);
 }
}

G4UserRunAction::~G4UserRunAction()
{;}

G4Run* G4UserRunAction::GenerateRun()
{ return 0; }

void G4UserRunAction::BeginOfRunAction(const G4Run*)
{;}

void G4UserRunAction::EndOfRunAction(const G4Run*)
{;}

