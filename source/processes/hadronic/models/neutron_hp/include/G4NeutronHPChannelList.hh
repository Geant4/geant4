// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPChannelList.hh,v 1.1 1999-01-07 16:12:55 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Process: Very Low Energy Neutron X-Sections
 // original by H.P. Wellisch, TRIUMF, 14-Feb-97
 // Builds and has the Cross-section data for one material.
 
#ifndef G4NeutronHPChannelList_h
#define G4NeutronHPChannelList_h 1

#include "globals.hh"
#include "G4NeutronHPChannel.hh"
#include "G4StableIsotopes.hh"

class G4Element;
class G4ParticleChange;
class G4Track;
class G4NeutronHPFinalState;

class G4NeutronHPChannelList
{      

  public:
  
  G4NeutronHPChannelList(G4int n);
  
  G4NeutronHPChannelList();
  
  void Init(G4int n);

  ~G4NeutronHPChannelList();
  
  G4ParticleChange * ApplyYourself(const G4Element * theElement, const G4Track & aTrack);
      
  void Init(G4Element * anElement, const G4String & dirName);
    
  void Register(G4NeutronHPFinalState *theFS, const G4String &aName);
  
  inline G4double GetXsec(G4double anEnergy)
  {
    G4double result=0;
    G4int i;
    for(i=0; i<nChannels; i++)
    {
      result+=theChannels[i]->GetXsec(anEnergy);
    }
    return result;
  }
  
  inline G4int GetNumberOfChannels() { return nChannels; }
      
  inline G4bool HasDataInAnyFinalState()
  {
    G4bool result = false;
    G4int i;
    for(i=0; i<nChannels; i++)
    {
      if(theChannels[i]->HasDataInAnyFinalState()) result = true;
    }
    return result;
  }
  inline void RestartRegistration()
  {
    allChannelsCreated = true;
    theInitCount = 0;
  }
  private:

  static G4int trycounter;
  G4NeutronHPChannel ** theChannels;
  G4int nChannels;
  G4String theDir;
  G4Element * theElement; 
  
  G4bool allChannelsCreated;
  G4int theInitCount;

  G4StableIsotopes theStableOnes;
  
};

#endif
