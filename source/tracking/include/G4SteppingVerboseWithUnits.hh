// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SteppingVerboseWithUnits.hh,v 1.1 1999-07-27 09:21:50 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  
//---------------------------------------------------------------
//
// G4SteppingVerboseWithUnits.hh
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

class G4SteppingVerboseWithUnits;

#ifndef G4SteppingVerboseWithUnits_h
#define G4SteppingVerboseWithUnits_h 1

#include "G4SteppingVerbose.hh"

class G4SteppingVerboseWithUnits : public G4SteppingVerbose {
public:   
// Constructor/Destructor
//  G4SteppingVerboseWithUnits();
  G4SteppingVerboseWithUnits(G4SteppingManager* const );
 ~G4SteppingVerboseWithUnits();
//
  virtual void AtRestDoItInvoked();
  virtual void AlongStepDoItAllDone();
  virtual void PostStepDoItAllDone();
  virtual void AlongStepDoItOneByOne();
  virtual void PostStepDoItOneByOne();
  virtual void StepInfo();
  virtual void TrackingStarted();
  virtual void DPSLStarted();
  virtual void DPSLUserLimit();
  virtual void DPSLPostStep();
  virtual void DPSLAlongStep();
  virtual void VerboseTrack();
//


};


#endif
