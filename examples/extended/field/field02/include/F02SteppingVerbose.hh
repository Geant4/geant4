// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: F02SteppingVerbose.hh,v 1.1 2001-03-27 16:26:21 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  
//---------------------------------------------------------------
//
// F02SteppingVerbose.hh
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

class F02SteppingVerbose;

#ifndef F02SteppingVerbose_h
#define F02SteppingVerbose_h 1

#include "G4SteppingVerbose.hh"

class F02SteppingVerbose : public G4SteppingVerbose 
{
  public:   
         // Constructor/Destructor

    F02SteppingVerbose();
   ~F02SteppingVerbose();

    void StepInfo();
    void TrackingStarted();

};

#endif
