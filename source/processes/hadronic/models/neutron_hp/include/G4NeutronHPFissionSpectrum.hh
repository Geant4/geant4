// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPFissionSpectrum.hh,v 1.3 1999-07-02 09:59:07 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPFissionSpectrum_h
#define G4NeutronHPFissionSpectrum_h 1

#include "globals.hh"
#include "G4NeutronHPVector.hh"
#include "Randomize.hh"
#include "G4ios.hh"
#include <fstream.h>
#include "G4VNeutronHPEDis.hh"

// we will need a List of these .... one per term.

class G4NeutronHPFissionSpectrum : public G4VNeutronHPEDis
{
  public:
  G4NeutronHPFissionSpectrum()
  {
    expm1 = exp(-1.);
  }
  ~G4NeutronHPFissionSpectrum()
  {
  }
  
  inline void Init(ifstream & aDataFile)
  {
    theFractionalProb.Init(aDataFile, eV);
    theThetaDist.Init(aDataFile, eV);
  }
  
  inline G4double GetFractionalProbability(G4double anEnergy)
  {
    return theFractionalProb.GetY(anEnergy);
  }
  
  inline G4double Sample(G4double anEnergy) 
  {
    G4double theta = theThetaDist.GetY(anEnergy);
    // here we need to sample Maxwells distribution, if 
    // need be.
    G4double result, cut;
    G4double range =50*MeV;
    G4double max = Maxwell((theta*eV)/2., theta);
    G4double value;
    do
    {
      result = range*G4UniformRand();
      value = Maxwell(result, theta);
      cut = G4UniformRand();
    }
    while(cut > value/max);
    return result;
  }
  
  private:
 
  // this is the function to sample from. 
  inline G4double Maxwell(G4double anEnergy, G4double theta)
  {
    G4double result = sqrt(anEnergy/eV)*exp(-anEnergy/eV/theta);
    return result;
  }
  
  private:
  
  G4double expm1;
  
  G4NeutronHPVector theFractionalProb;
  
  G4NeutronHPVector theThetaDist;
  
};

#endif
