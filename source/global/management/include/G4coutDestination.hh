// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4coutDestination.hh,v 1.2 1999-11-16 17:40:51 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ---------------------------------------------------------------
// GEANT 4 class header file
//
// G4coutDestination.hh
//
// ---------------------------------------------------------------
#ifndef G4COUTDESTINATION_HH
#define G4COUTDESTINATION_HH

#include "globals.hh"

class G4coutDestination
{
public:
  virtual int ReceiveG4cout(G4String){return 0;}
  virtual int ReceiveG4cerr(G4String){return 0;}
};

#endif

