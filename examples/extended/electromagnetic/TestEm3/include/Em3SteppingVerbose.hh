// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em3SteppingVerbose.hh,v 1.3 2000-02-29 12:18:58 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  
//---------------------------------------------------------------
//
// Em3SteppingVerbose.hh
//
// Description:
//   This class manages the vervose outputs in G4SteppingManager. 
//   
//
//---------------------------------------------------------------

class Em3SteppingVerbose;

#ifndef Em3SteppingVerbose_h
#define Em3SteppingVerbose_h 1

#include "G4SteppingVerbose.hh"

class Em3SteppingVerbose : public G4SteppingVerbose {
public:   
// Constructor/Destructor
  Em3SteppingVerbose();
 ~Em3SteppingVerbose();
//
  void StepInfo();
  void TrackingStarted();
//


};

#endif
