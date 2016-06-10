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
// $Id: G4UserSteppingAction.cc 66241 2012-12-13 18:34:42Z gunter $
//
// ---------------------------------------------------------------
//
// G4UserSteppingAction.cc
//
// Description:
//   This class represents actions taken place by the user at each
//   end of stepping. 
//
// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//
// ---------------------------------------------------------------

#include "G4UserSteppingAction.hh"
#include "G4ParticleTable.hh"
#include "globals.hh"

/////////////////////////////////////////////////////////
G4UserSteppingAction::G4UserSteppingAction()
/////////////////////////////////////////////////////////
   : fpSteppingManager(0)
{
 if(!(G4ParticleTable::GetParticleTable()->GetReadiness()))
 {
   G4String msg;
   msg =  " You are instantiating G4UserSteppingAction BEFORE your\n";
   msg += "G4VUserPhysicsList is instantiated and assigned to G4RunManager.\n";
   msg += " Such an instantiation is prohibited since Geant4 version 8.0. To fix this problem,\n";
   msg += "please make sure that your main() instantiates G4VUserPhysicsList AND\n";
   msg += "set it to G4RunManager before instantiating other user action classes\n";
   msg += "such as G4UserSteppingAction.";
   G4Exception("G4UserSteppingAction::G4UserSteppingAction()",
               "Tracking0002", FatalException, msg);
 }
}

/////////////////////////////////////////////////////////
G4UserSteppingAction::~G4UserSteppingAction()
/////////////////////////////////////////////////////////
{;}

/////////////////////////////////////////////////////////
void G4UserSteppingAction::
     SetSteppingManagerPointer(G4SteppingManager* pValue)
/////////////////////////////////////////////////////////
{
  fpSteppingManager = pValue;
}  
