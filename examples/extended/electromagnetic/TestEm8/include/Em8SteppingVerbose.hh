// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em8SteppingVerbose.hh,v 1.1 2000-01-07 14:50:23 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  
//---------------------------------------------------------------
//
// Em8SteppingVerbose.hh
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

class Em8SteppingVerbose;

#ifndef Em8SteppingVerbose_h
#define Em8SteppingVerbose_h 1

#include "G4SteppingVerbose.hh"

class Em8SteppingVerbose : public G4SteppingVerbose 
{
  public:   
         // Constructor/Destructor

    Em8SteppingVerbose();
   ~Em8SteppingVerbose();

    void StepInfo();
    void TrackingStarted();

};

#endif
