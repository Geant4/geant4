// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SteppingControl.hh,v 2.1 1998/07/12 03:08:49 urbi Exp $
// GEANT4 tag $Name: geant4-00 $
//
//
//---------------------------------------------------------------
//
// G4SteppingControl  
//
// Description:
//   This enumaration specifies possible conditions to control
//   the stepping manager behavier.
//
// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//
//---------------------------------------------------------------

#ifndef G4SteppingControl_h
#define G4SteppingControl_h 1

/////////////////////
enum G4SteppingControl  
/////////////////////
{
  NormalCondition,
  AvoidHitInvocation,            
    // Hit will NOT be called 
  Debug
};

#endif


