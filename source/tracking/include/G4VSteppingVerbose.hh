// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VSteppingVerbose.hh,v 1.1 1999-03-24 04:45:47 tsasaki Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  
//---------------------------------------------------------------
//
// G4VSteppingVerbose.hh
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

#ifndef G4VSteppingVerbose_h
#define G4VSteppingVerbose_h
class G4SteppingManager;

class G4VSteppingVerbose{

public:   
// Constructor/Destructor
  G4VSteppingVerbose();
  G4VSteppingVerbose(G4SteppingManager*);
  virtual ~G4VSteppingVerbose();
//
    virtual void NewStep() = 0;
    virtual void CopyState() = 0;
    virtual void AtRestDoItInvoked() = 0;
    virtual void AlongStepDoItAllDone() = 0;
    virtual void PostStepDoItAllDone() = 0;
    virtual void AlongStepDoItOneByOne() = 0;
    virtual void PostStepDoItOneByOne() = 0;
    virtual void StepInfo() = 0;
    virtual void TrackingStarted() = 0;
    virtual void DPSLStarted() = 0;
    virtual void DPSLUserLimit() = 0;
    virtual void DPSLPostStep() = 0;
    virtual void DPSLAlongStep() = 0;
//    virtual void DPSLAlongStepDoItOneByOne() = 0;
//    virtual void DPSLPostStepDoItOneByOne() = 0;
    virtual void VerboseTrack() = 0;
    virtual void VerboseParticleChange() = 0;
};
#endif
