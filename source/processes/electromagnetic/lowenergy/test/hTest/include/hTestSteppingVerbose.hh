// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: hTestSteppingVerbose.hh,v 1.1 2000-05-21 18:37:46 chauvie Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  
//---------------------------------------------------------------
//
// hTestSteppingVerbose.hh
//
// Description:
//   This class manages the vervose outputs in G4SteppingManager. 
//   
//
//---------------------------------------------------------------

class hTestSteppingVerbose;

#ifndef hTestSteppingVerbose_h
#define hTestSteppingVerbose_h 1

#include "G4SteppingVerbose.hh"

class hTestSteppingVerbose : public G4SteppingVerbose {
public:   
// Constructor/Destructor
  hTestSteppingVerbose();
 ~hTestSteppingVerbose();
//
  void StepInfo();
  void TrackingStarted();
//


};

#endif
