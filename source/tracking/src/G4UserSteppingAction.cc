// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UserSteppingAction.cc,v 1.4 1999-12-15 14:54:00 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

/////////////////////////////////////////////////////////
void G4UserSteppingAction::
     SetSteppingManagerPointer(G4SteppingManager* pValue)
/////////////////////////////////////////////////////////
{
  fpSteppingManager = pValue;
}  




