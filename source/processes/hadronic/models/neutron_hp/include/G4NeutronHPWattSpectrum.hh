// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPWattSpectrum.hh,v 1.1 1999-01-07 16:13:10 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPWattSpectrum_h
#define G4NeutronHPWattSpectrum_h 1

#include "globals.hh"
#include "G4NeutronHPVector.hh"
#include "Randomize.hh"
#include "G4ios.hh"
#include <fstream.h>
#include "G4VNeutronHPEDis.hh"

// we will need a List of these .... one per term.

class G4NeutronHPWattSpectrum : public G4VNeutronHPEDis
{
  public:
  G4NeutronHPWattSpectrum()
  {
    expm1 = exp(-1.);
  }
  ~G4NeutronHPWattSpectrum()
  {
  }
  
  inline void Init(ifstream & aDataFile)
  {
    theFractionalProb.Init(aDataFile, eV);
    theApar.Init(aDataFile, eV);
    theBpar.Init(aDataFile, eV);
  }
  
  inline G4double GetFractionalProbability(G4double anEnergy)
  {
    return theFractionalProb.GetY(anEnergy);
  }
  
  G4double Sample(G4double anEnergy);
  
  private:
  
  inline G4double Watt(G4double anEnergy, G4double a, G4double b)
  {
    G4double energy = anEnergy/eV;
    G4double result = exp(-energy/a)*sinh(sqrt(b*energy));
    return result;
  }
  
  private:
  
  G4double expm1;
  
  G4NeutronHPVector theFractionalProb;
  
  G4NeutronHPVector theApar;
  G4NeutronHPVector theBpar;
  
};

#endif
