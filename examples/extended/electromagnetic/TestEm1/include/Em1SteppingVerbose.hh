// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em1SteppingVerbose.hh,v 1.3 2000-02-29 12:13:24 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  
//---------------------------------------------------------------
//
// Em1SteppingVerbose.hh
//
// Description:
//   This class manages the verbose outputs in G4SteppingManager. 
//
//---------------------------------------------------------------

class Em1SteppingVerbose;

#ifndef Em1SteppingVerbose_h
#define Em1SteppingVerbose_h 1

#include "G4SteppingVerbose.hh"

class Em1SteppingVerbose : public G4SteppingVerbose {
public:   
// Constructor/Destructor
  Em1SteppingVerbose();
 ~Em1SteppingVerbose();
//
  void StepInfo();
  void TrackingStarted();
//


};

#endif
