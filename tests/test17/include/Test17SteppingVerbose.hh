// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Test17SteppingVerbose.hh,v 1.1 2000-05-26 06:34:28 chauvie Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  
//---------------------------------------------------------------
//
// Test17SteppingVerbose.hh
//
// Description:
//   This class manages the vervose outputs in G4SteppingManager. 
//   
//
//---------------------------------------------------------------

class Test17SteppingVerbose;

#ifndef Test17SteppingVerbose_h
#define Test17SteppingVerbose_h 1

#include "G4SteppingVerbose.hh"

class Test17SteppingVerbose : public G4SteppingVerbose {
public:   
// Constructor/Destructor
  Test17SteppingVerbose();
 ~Test17SteppingVerbose();
//
  void StepInfo();
  void TrackingStarted();
//


};

#endif
