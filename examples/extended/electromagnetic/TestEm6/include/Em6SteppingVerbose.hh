// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em6SteppingVerbose.hh,v 1.1 1999-11-10 17:58:27 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  
//---------------------------------------------------------------
//
// Em6SteppingVerbose.hh
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

class Em6SteppingVerbose;

#ifndef Em6SteppingVerbose_h
#define Em6SteppingVerbose_h 1

#include "G4SteppingVerbose.hh"

class Em6SteppingVerbose : public G4SteppingVerbose {
public:   
// Constructor/Destructor
  Em6SteppingVerbose();
 ~Em6SteppingVerbose();
//
  void StepInfo();
  void TrackingStarted();
//


};

#endif
