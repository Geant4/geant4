// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPPhotonXSection.hh,v 1.3 1999-07-02 09:59:52 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPPhotonXSection_h
#define G4NeutronHPPhotonXSection_h 1

#include "globals.hh"
#include "G4NeutronHPVector.hh"
#include "Randomize.hh"
#include "G4ios.hh"
#include <fstream.h>
#include "globals.hh"
#include "G4VNeutronVector.hh"

// we will need a List of these .... one per term.

class G4NeutronHPPhotonXSection 
{
  public:
  G4NeutronHPPhotonXSection()
  {
    theExclusive = NULL;
    theExShell = NULL;
    theExEnergy = NULL;
    theExFlag = NULL;
    theExDisFlag = NULL;
  }
  ~G4NeutronHPPhotonXSection()
  {
    if(theExclusive!=NULL) delete [] theExclusive;
    if(theExShell != NULL) delete [] theExShell;
    if(theExEnergy != NULL) delete [] theExEnergy;
    if(theExFlag != NULL) delete [] theExFlag;
    if(theExDisFlag != NULL) delete [] theExDisFlag;
  }
  
  inline void Init(ifstream & aDataFile)
  {
    aDataFile  >> nChannels >> targetMass;
    if(nChannels!=1) 
    {
      aDataFile >> theIncEnergy>>theIncShell>>theIncFlag>>theIncDisFlag;
      theaDataFileInclusive.Init(aDataFile, eV);
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
  G4int theIncFlag
  G4int theIncDisFlag;
  G4NeutronHPVector theInclusive;
  
  G4int nChannels;
  G4double * theExShell;
  G4double * theExEnergy;
  G4int * theExFlag
  G4int * theExDisFlag;
  G4NeutronHPVector * theExclusive;
  
};

#endif
