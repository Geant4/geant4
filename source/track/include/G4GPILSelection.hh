// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4GPILSelection.hh,v 1.2 1999-11-07 16:31:58 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//---------------------------------------------------------------
//
// G4GPILSelection  
//
// Class Description:
// This enumaration is used to control whether a AlongStepProcess
//   can be a winner of the GPIL race or not. 
//
// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//
//---------------------------------------------------------------

#ifndef G4GPILSelection_h
#define G4GPILSelection_h 1

/////////////////////
enum G4GPILSelection  
/////////////////////
{
  CandidateForSelection,            
  // This AlongStep process partecipates in the process selection 
  // mechanism, i.e. it can be the winner of the GPIL race.
  // (this case is default)

  NotCandidateForSelection          
  // This AlongStep process does not partecipate in the
  // process selection mechanism even when it limits the Step.

};

#endif


