// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SteppingControl.hh,v 1.2 1999-11-07 16:32:03 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//---------------------------------------------------------------
//
// G4SteppingControl  
//
// Class Description:
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


