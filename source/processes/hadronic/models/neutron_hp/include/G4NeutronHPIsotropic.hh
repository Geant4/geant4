// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPIsotropic.hh,v 1.3 1999-07-02 09:59:17 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPIsotropic_h
#define G4NeutronHPIsotropic_h 1

#include "G4ios.hh"
#include <fstream.h>
#include "globals.hh"
#include "G4ReactionProduct.hh"
#include "G4VNeutronHPEnergyAngular.hh"

class G4NeutronHPIsotropic : public G4VNeutronHPEnergyAngular
{
  public:
  
  G4NeutronHPIsotropic(){}
  
  ~G4NeutronHPIsotropic(){}
  
  public:
  
  void Init(ifstream & aDataFile);
  G4ReactionProduct * Sample(G4double anEnergy, G4double massCode, G4double mass);
  G4double MeanEnergyOfThisInteraction()
  {
    return -1.;
  }
  private:

};
#endif
