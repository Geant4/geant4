// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPPhotonXSection.hh,v 1.2 1999-06-29 18:44:12 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Here for historic reasons only.
//
#ifndef G4NeutronHPPhotonXSection_h
#define G4NeutronHPPhotonXSection_h 1

#include "globals.hh"
#include "G4NeutronHPVector.hh"
#include "Randomize.hh"
#include "G4ios.hh"
#include <fstream.h>
#include "globals.hh"

// we will need a List of these .... one per term.

class G4NeutronHPPhotonXSection 
{
  public:
  G4NeutronHPPhotonXSection();
  ~G4NeutronHPPhotonXSection();
  
  inline void Init(ifstream & aDataFile)
  {
    aDataFile  >> nChannels >> targetMass;
    if(nChannels!=1) 
    {
      aDataFile >> theIncEnergy>>theIncShell>>theIncFlag>>theIncDisFlag;
      theInclusive.Init(aDataFile, eV);
    }
    theExclusive = new G4NeutronHPVector[nChannels];
    theExShell = new G4double[nChannels];
    theExEnergy = new G4double[nChannels];
    theExFlag = new G4int[nChannels];
    theExDisFlag = new G4int[nChannels];   
    for(G4int i=0; i<nChannels; i++)
    {
      aDataFile>>theExEnergy[i]>>theExShell[i]>>theExFlag[i]>>theExDisFlag[i];
      theExclusive[i].Init(aDataFile,eV);
    }
  }
  
  G4double Sample(G4double anEnergy)
  {
    return -1;
  }
  
  private:
   
  G4double targetMass;
  
  G4double theIncShell;
  G4double theIncEnergy;
  G4int theIncFlag;
  G4int theIncDisFlag;
  G4NeutronHPVector theInclusive;
  
  G4int nChannels;
  G4double * theExShell;
  G4double * theExEnergy;
  G4int * theExFlag;
  G4int * theExDisFlag;
  G4NeutronHPVector * theExclusive;
  
};

#endif
