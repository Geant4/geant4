// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ForceCondition.hh,v 1.1 1999-01-07 16:14:20 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//---------------------------------------------------------------
//
// G4ForceCondition  
//
// Description:
//   This enumaration specifies possible conditions the three
//   types of DoIt can be assinged by physics processes.
//
// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//
//---------------------------------------------------------------

#ifndef G4ForceCondition_h
#define G4ForceCondition_h 1

/////////////////////
enum G4ForceCondition  
/////////////////////
{
  Forced,            
    // This PostStepDoIt is forced to invoke.
  NotForced,         
    // This PostStepDoIt is not forced to invoke.
  Conditionally,     
    // This PostStepDoIt is forced to invoke only when corresponding
    // AlongStepDoIt limits the Step.
  ExclusivelyForced  
    // Only this PostStepDoIt (or AtRestDoIt) is exclusively forced 
    // to invoke - all other DoIt including AlongStepDoIts are ignored.
};

#endif


