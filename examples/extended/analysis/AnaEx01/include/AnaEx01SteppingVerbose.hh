// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: AnaEx01SteppingVerbose.hh,v 1.1.1.1 2000-09-14 11:37:21 barrand Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  
//---------------------------------------------------------------
//
// AnaEx01SteppingVerbose.hh
//
// Description:
//   This class manages the verbose outputs in G4SteppingManager. 
//   It inherits from G4SteppingVerbose   
//
//---------------------------------------------------------------

class AnaEx01SteppingVerbose;

#ifndef AnaEx01SteppingVerbose_h
#define AnaEx01SteppingVerbose_h 1

#include "G4SteppingVerbose.hh"

class AnaEx01SteppingVerbose : public G4SteppingVerbose {
public:   
// Constructor/Destructor
  AnaEx01SteppingVerbose();
 ~AnaEx01SteppingVerbose();
//
  void StepInfo();
  void TrackingStarted();
//


};

#endif
