// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN02SteppingVerbose.hh,v 1.3 2000-02-28 17:58:49 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  
//---------------------------------------------------------------
//
// ExN02SteppingVerbose.hh
//
// Description:
//   This class manages the vervose outputs in G4SteppingManager. 
//   It inherits from G4SteppingVerbose
//
//---------------------------------------------------------------

class ExN02SteppingVerbose;

#ifndef ExN02SteppingVerbose_h
#define ExN02SteppingVerbose_h 1

#include "G4SteppingVerbose.hh"

class ExN02SteppingVerbose : public G4SteppingVerbose {
public:   
// Constructor/Destructor
  ExN02SteppingVerbose();
 ~ExN02SteppingVerbose();
//
  void StepInfo();
  void TrackingStarted();
//


};

#endif
