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
// $Id: G4UserRunAction.cc,v 1.7 2005/11/22 19:54:53 asaim Exp $
// GEANT4 tag $Name: geant4-08-00 $
//

#include "G4UserRunAction.hh"
#include "G4ParticleTable.hh"
#include "globals.hh"


G4UserRunAction::G4UserRunAction()
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
              "Run0002",FatalException,msg);
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

