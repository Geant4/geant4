//
// $Id: ThyroidSteppingVerboseTest.hh,v 1.1 2003-05-23 11:55:59 francy Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//   This class manages the verbose outputs in G4SteppingManager. 
//   It inherits from G4SteppingVerbose.
//   It shows how to extract informations during the tracking of a particle.
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

class ThyroidSteppingVerbose;

#ifndef ThyroidSteppingVerbose_h
#define ThyroidSteppingVerbose_h 1

#include "G4SteppingVerbose.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

class ThyroidSteppingVerbose : public G4SteppingVerbose {

public:
   
 ThyroidSteppingVerbose();
 ~ThyroidSteppingVerbose();

  void StepInfo();
  void TrackingStarted();

};

#endif










