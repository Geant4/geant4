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
// $Id: G4UserTrackingAction.cc 66241 2012-12-13 18:34:42Z gunter $
//
// ---------------------------------------------------------------
//
// G4UserTrackingAction.cc
//
// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//
// ---------------------------------------------------------------

#include "G4UserTrackingAction.hh"
#include "G4ParticleTable.hh"
#include "globals.hh"

/////////////////////////////////////////////////////////
G4UserTrackingAction::G4UserTrackingAction()
/////////////////////////////////////////////////////////
  : fpTrackingManager(0)
{
  if(!(G4ParticleTable::GetParticleTable()->GetReadiness()))
  {
   G4String msg;
   msg =  " You are instantiating G4UserTrackingAction BEFORE your\n";
   msg += "G4VUserPhysicsList is instantiated and assigned to G4RunManager.\n";
   msg += " Such an instantiation is prohibited since Geant4 version 8.0. To fix this problem,\n";
   msg += "please make sure that your main() instantiates G4VUserPhysicsList AND\n";
   msg += "set it to G4RunManager before instantiating other user action classes\n";
   msg += "such as G4UserTrackingAction.";
   G4Exception("G4UserTrackingAction::G4UserTrackingAction()",
               "Tracking0001", FatalException, msg);
  }
}

/////////////////////////////////////////////////////////
G4UserTrackingAction::~G4UserTrackingAction()
/////////////////////////////////////////////////////////
{;}

/////////////////////////////////////////////////////////
void G4UserTrackingAction::
     SetTrackingManagerPointer(G4TrackingManager* pValue)
/////////////////////////////////////////////////////////
{
  fpTrackingManager = pValue;
}
