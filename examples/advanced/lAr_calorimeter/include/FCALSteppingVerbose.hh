// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: FCALSteppingVerbose.hh,v 1.1 2002-10-02 19:37:01 ahoward Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  
//---------------------------------------------------------------
//
// ExN03SteppingVerbose.hh
//
// Description:
//   This class manages the verbose outputs in G4SteppingManager. 
//   It inherits from G4SteppingVerbose   
//
//---------------------------------------------------------------

class ExN03SteppingVerbose;

#ifndef ExN03SteppingVerbose_h
#define ExN03SteppingVerbose_h 1

#include "G4SteppingVerbose.hh"

class ExN03SteppingVerbose : public G4SteppingVerbose {
public:   
// Constructor/Destructor
  ExN03SteppingVerbose();
 ~ExN03SteppingVerbose();
//
  void StepInfo();
  void TrackingStarted();
//


};

#endif
