// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em1SteppingVerbose.hh,v 1.2 1999-12-15 14:48:56 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  
//---------------------------------------------------------------
//
// Em1SteppingVerbose.hh
//
// Description:
//   This class manages the vervose outputs in G4SteppingManager. 
//   
//
// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//
//---------------------------------------------------------------

class Em1SteppingVerbose;

#ifndef Em1SteppingVerbose_h
#define Em1SteppingVerbose_h 1

#include "G4SteppingVerbose.hh"

class Em1SteppingVerbose : public G4SteppingVerbose {
public:   
// Constructor/Destructor
  Em1SteppingVerbose();
 ~Em1SteppingVerbose();
//
  void StepInfo();
  void TrackingStarted();
//


};

#endif
