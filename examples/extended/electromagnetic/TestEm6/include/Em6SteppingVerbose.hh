// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em6SteppingVerbose.hh,v 1.3 2000-02-29 12:26:36 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  
//---------------------------------------------------------------
//
// Em6SteppingVerbose.hh
//
// Description:
//   This class manages the vervose outputs in G4SteppingManager. 
//   
//
//---------------------------------------------------------------

class Em6SteppingVerbose;

#ifndef Em6SteppingVerbose_h
#define Em6SteppingVerbose_h 1

#include "G4SteppingVerbose.hh"

class Em6SteppingVerbose : public G4SteppingVerbose {
public:   
// Constructor/Destructor
  Em6SteppingVerbose();
 ~Em6SteppingVerbose();
//
  void StepInfo();
  void TrackingStarted();
//


};

#endif
