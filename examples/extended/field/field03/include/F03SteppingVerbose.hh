// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: F03SteppingVerbose.hh,v 1.1 2001-06-08 11:55:40 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  
//---------------------------------------------------------------
//
// F03SteppingVerbose.hh
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

class F03SteppingVerbose;

#ifndef F03SteppingVerbose_h
#define F03SteppingVerbose_h 1

#include "G4SteppingVerbose.hh"

class F03SteppingVerbose : public G4SteppingVerbose 
{
  public:   
         // Constructor/Destructor

    F03SteppingVerbose();
   ~F03SteppingVerbose();

    void StepInfo();
    void TrackingStarted();

};

#endif
