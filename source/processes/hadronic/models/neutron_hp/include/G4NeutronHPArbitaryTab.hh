// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPArbitaryTab.hh,v 1.1 1999-01-07 16:12:54 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4NeutronHPArbitaryTab_h
#define G4NeutronHPArbitaryTab_h 1

#include "globals.hh"
#include "G4NeutronHPVector.hh"
#include "Randomize.hh"
#include "G4ios.hh"
#include <fstream.h>
#include "G4VNeutronHPEDis.hh"
#include "G4InterpolationManager.hh"

// we will need a List of these .... one per term.

class G4NeutronHPArbitaryTab : public G4VNeutronHPEDis
{
  public:
  G4NeutronHPArbitaryTab()
  {
   theDistFunc = NULL;
  }
  ~G4NeutronHPArbitaryTab()
  {
   if(theDistFunc!=NULL) delete [] theDistFunc;
  }
  
  inline void Init(ifstream & theData)
  {
    G4int i, total;
    theFractionalProb.Init(theData, eV);
    theData >> nDistFunc; // = number of incoming n energy points
    theDistFunc = new G4NeutronHPVector [nDistFunc];
    theManager.Init(theData);
    G4double currentEnergy;
    for(i=0; i<nDistFunc; i++)
    {
      theData >> currentEnergy;
      theDistFunc[i].SetLabel(currentEnergy*eV);
      theDistFunc[i].Init(theData, eV);
      theDistFunc[i].ThinOut(0.02); // @@@ optimization to be finished.
    }
  }
  
  inline G4double GetFractionalProbability(G4double anEnergy)
  {
    return theFractionalProb.GetY(anEnergy);
  }
  
  G4double Sample(G4double anEnergy) ;
  
  private:
  
  G4NeutronHPVector theFractionalProb;
  G4int nDistFunc;
  G4InterpolationManager theManager; // knows the interpolation between stores
  G4NeutronHPVector * theDistFunc; // one per incoming energy
  G4NeutronHPVector theBuffer;
  
};

#endif
