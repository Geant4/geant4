// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: F01SteppingVerbose.hh,v 1.1 2001-03-27 16:21:31 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  
//---------------------------------------------------------------
//
// F01SteppingVerbose.hh
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

class F01SteppingVerbose;

#ifndef F01SteppingVerbose_h
#define F01SteppingVerbose_h 1

#include "G4SteppingVerbose.hh"

class F01SteppingVerbose : public G4SteppingVerbose 
{
  public:   
         // Constructor/Destructor

    F01SteppingVerbose();
   ~F01SteppingVerbose();

    void StepInfo();
    void TrackingStarted();

};

#endif
