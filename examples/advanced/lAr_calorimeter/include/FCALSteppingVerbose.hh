// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: FCALSteppingVerbose.hh,v 1.2 2002-10-02 19:40:09 ahoward Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  
//---------------------------------------------------------------
//
// FCALSteppingVerbose.hh
//
// Description:
//   This class manages the verbose outputs in G4SteppingManager. 
//   It inherits from G4SteppingVerbose   
//
//---------------------------------------------------------------

class FCALSteppingVerbose;

#ifndef FCALSteppingVerbose_h
#define FCALSteppingVerbose_h 1

#include "G4SteppingVerbose.hh"

class FCALSteppingVerbose : public G4SteppingVerbose {
public:   
// Constructor/Destructor
  FCALSteppingVerbose();
 ~FCALSteppingVerbose();
//
  void StepInfo();
  void TrackingStarted();
//


};

#endif
