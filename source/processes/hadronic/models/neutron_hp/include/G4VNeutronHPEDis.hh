// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VNeutronHPEDis.hh,v 1.3 1999-07-02 10:00:04 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4VNeutronHPEDis_h
#define G4VNeutronHPEDis_h 1

#include "globals.hh"
#include "G4NeutronHPVector.hh"
#include "Randomize.hh"
#include "G4ios.hh"
#include <fstream.h>

class G4VNeutronHPEDis
{
  public:
  G4VNeutronHPEDis()
  {
  }
  virtual ~G4VNeutronHPEDis()
  {
  }
  
  virtual void Init(ifstream & theData) = 0;
  
  virtual G4double GetFractionalProbability(G4double anEnergy) = 0;
  
  virtual G4double Sample(G4double anEnergy) = 0; 
  
  private:
    
};

#endif
