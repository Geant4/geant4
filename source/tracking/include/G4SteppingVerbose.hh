// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SteppingVerbose.hh,v 1.8 2000-11-11 06:34:09 tsasaki Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  
//---------------------------------------------------------------
//
// G4SteppingVerbose.hh
//
// class dscription:
//   This class manages the vervose outputs in G4SteppingManager. 
//   
//
// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//
//---------------------------------------------------------------

class G4SteppingVerbose;

#ifndef G4SteppingVerose_h
#define G4SteppingVerose_h 1

#include "G4VSteppingVerbose.hh"

class G4SteppingVerbose : public G4VSteppingVerbose {
public:   // with description
// Constructor/Destructor
  G4SteppingVerbose();
 ~G4SteppingVerbose();
// methods to be invoked in the SteppingManager
  void NewStep();
  void AtRestDoItInvoked();
  void AlongStepDoItAllDone();
  void PostStepDoItAllDone();
  void AlongStepDoItOneByOne();
  void PostStepDoItOneByOne();
  void StepInfo();
  void TrackingStarted();
  void DPSLStarted();
  void DPSLUserLimit();
  void DPSLPostStep();
  void DPSLAlongStep();
//  void DPSLAlongStepDoItOneByOne();
//  void DPSLPostStepDoItOneByOne();
  void VerboseTrack();
  void VerboseParticleChange();
  void ShowStep() const;
//

};


#endif

