// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPNames.hh,v 1.1 1999-01-07 16:13:07 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPNames_h
#define G4NeutronHPNames_h 1

#include "G4ios.hh"
#include <fstream.h>
#ifdef WIN32
  #include <strstrea.h>
#else
  #include <strstream.h>
#endif
#include <stdlib.h>
#include "globals.hh"
#include "G4NeutronHPDataUsed.hh"

class G4NeutronHPNames
{
  public:
  
  G4NeutronHPNames(){}
  ~G4NeutronHPNames(){}
  
  G4NeutronHPDataUsed GetName(G4int A, G4int Z, G4String base, G4String rest, G4bool & active);
  G4String GetName(G4int i) { return theString[i]; }
  
  private:
  
  static const G4String theString[99];

};
#endif
