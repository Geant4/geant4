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
// $Id: G4UserSteppingAction.cc,v 1.7 2005/11/22 21:04:07 asaim Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
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
{
 if(!(G4ParticleTable::GetParticleTable()->GetReadiness()))
 {
   G4String msg;
   msg =  " You are instantiating G4UserSteppingAction BEFORE your\n";
   msg += "G4VUserPhysicsList is instantiated and assigned to G4RunManager.\n";
   msg += " Such an instantiation is prohibited by Geant4 version 8.0. To fix this problem,\n";
   msg += "please make sure that your main() instantiates G4VUserPhysicsList AND\n";
   msg += "set it to G4RunManager before instantiating other user action classes\n";
   msg += "such as G4UserSteppingAction.";
   G4Exception("G4UserSteppingAction::G4UserSteppingAction()",
              "Tracking0002",FatalException,msg);
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




