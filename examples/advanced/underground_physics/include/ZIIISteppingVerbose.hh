// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ZIIISteppingVerbose.hh,v 1.1 2001-06-26 11:23:24 ahoward Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class ZIIISteppingVerbose;

#ifndef ZIIISteppingVerbose_h
#define ZIIISteppingVerbose_h 1

#include "G4SteppingVerbose.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class ZIIISteppingVerbose : public G4SteppingVerbose
{
 public:   

   ZIIISteppingVerbose();
  ~ZIIISteppingVerbose();

   void StepInfo();
   void TrackingStarted();

};

#endif
