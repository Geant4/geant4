// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPSimpleEvapSpectrum.hh,v 1.2 1999-06-29 18:44:14 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPSimpleEvapSpectrum_h
#define G4NeutronHPSimpleEvapSpectrum_h 1

#include "globals.hh"
#include "G4NeutronHPVector.hh"
#include "Randomize.hh"
#include "G4ios.hh"
#include <fstream.h>
#include "G4VNeutronHPEDis.hh"

// we will need a List of these .... one per term.

class G4NeutronHPSimpleEvapSpectrum : public G4VNeutronHPEDis
{
  public:
  G4NeutronHPSimpleEvapSpectrum();
  ~G4NeutronHPSimpleEvapSpectrum();
  
  inline void Init(ifstream & aDataFile)
  {
    theFractionalProb.Init(aDataFile,eV);
    theThetaDist.Init(aDataFile, eV);
  }
  
  inline G4double GetFractionalProbability(G4double anEnergy)
  {
    return theFractionalProb.GetY(anEnergy);
  }
  
  inline G4double Sample(G4double anEnergy) 
  {
    G4double theta = theThetaDist.GetY(anEnergy)*eV;
    G4double random, cut, max, result;
    max = 10.*theta;
    do
    {
      random = G4UniformRand();
      result = -theta*log(random); 
      cut = G4UniformRand();
    }
    while(cut>result/max); 
    return result;
  }
  
  private:
  
  inline G4double Evapo(G4double anEnergy, G4double theta)
  {
    G4double result = (anEnergy*eV)*exp(-anEnergy*eV/theta);
    return result;
  }
  
  private:
  
  G4double expm1;
  
  G4NeutronHPVector theFractionalProb;
  
  G4NeutronHPVector theThetaDist;
  
};

#endif
