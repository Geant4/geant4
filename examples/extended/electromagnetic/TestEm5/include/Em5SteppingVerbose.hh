// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em5SteppingVerbose.hh,v 1.5 2001-02-19 14:08:22 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//   This class manages the verbose outputs in G4SteppingManager. 
//   It inherits from G4SteppingVerbose.
//   It shows how to extract informations during the tracking of a particle.
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

#ifndef Em5SteppingVerbose_h
#define Em5SteppingVerbose_h 1

#include "G4SteppingVerbose.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

class Em5SteppingVerbose : public G4SteppingVerbose {

public:   

  Em5SteppingVerbose();
 ~Em5SteppingVerbose();

  void StepInfo();
  void TrackingStarted();

};

#endif
