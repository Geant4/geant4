// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em2SteppingVerbose.hh,v 1.1 1999-11-10 16:44:17 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  
//---------------------------------------------------------------
//
// Em2SteppingVerbose.hh
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

class Em2SteppingVerbose;

#ifndef Em2SteppingVerbose_h
#define Em2SteppingVerbose_h 1

#include "G4SteppingVerbose.hh"

class Em2SteppingVerbose : public G4SteppingVerbose {
public:   
// Constructor/Destructor
  Em2SteppingVerbose();
 ~Em2SteppingVerbose();
//
  void StepInfo();
  void TrackingStarted();
//


};

#endif
