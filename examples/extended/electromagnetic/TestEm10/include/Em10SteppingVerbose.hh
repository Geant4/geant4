// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em10SteppingVerbose.hh,v 1.1 2000-07-14 15:51:18 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  
//---------------------------------------------------------------
//
// Em10SteppingVerbose.hh
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

class Em10SteppingVerbose;

#ifndef Em10SteppingVerbose_h
#define Em10SteppingVerbose_h 1

#include "G4SteppingVerbose.hh"

class Em10SteppingVerbose : public G4SteppingVerbose 
{
  public:   
         // Constructor/Destructor

    Em10SteppingVerbose();
   ~Em10SteppingVerbose();

    void StepInfo();
    void TrackingStarted();

};

#endif
