// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPEvapSpectrum.hh,v 1.4 1999-12-15 14:53:11 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPEvapSpectrum_h
#define G4NeutronHPEvapSpectrum_h 1

#include "globals.hh"
#include "G4NeutronHPVector.hh"
#include "Randomize.hh"
#include "G4ios.hh"
#include "g4std/fstream"
#include "G4VNeutronHPEDis.hh"

// we will need a List of these .... one per term.

class G4NeutronHPEvapSpectrum : public G4VNeutronHPEDis
{
  public:
  G4NeutronHPEvapSpectrum()
  {
  }
  ~G4NeutronHPEvapSpectrum()
  {
  }
  
  inline void Init(G4std::ifstream & aDataFile)
  {
    theFractionalProb.Init(aDataFile);
    theThetaDist.Init(aDataFile);
    theXDist.Init(aDataFile);
  }
  
  inline G4double GetFractionalProbability(G4double anEnergy)
  {
    return theFractionalProb.GetY(anEnergy);
  }
  
  inline G4double Sample(G4double anEnergy) 
  {
    // when this is called, theFractionalProb was used, and 'k' is sorted out already.
    G4double x = theXDist.Sample();
    G4double theta = theThetaDist.GetY(anEnergy);
    G4double result = x*theta;
    return result*eV;
  }
  
  private:
  
  G4NeutronHPVector theFractionalProb;
  
  G4NeutronHPVector theThetaDist;
  G4NeutronHPVector theXDist;
  
};

#endif
