// This code implementation is the intellectual property of
// neutron_hp -- header file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronHPChannel.hh,v 1.4 1999-07-06 16:56:21 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Process: Very Low Energy Neutron X-Sections
 // original by H.P. Wellisch, TRIUMF, 14-Feb-97
 // Builds and has the Cross-section data for one element and channel.
//
// Bug fixes and workarounds in the destructor, F.W.Jones 06-Jul-1999
 
#ifndef G4NeutronHPChannel_h
#define G4NeutronHPChannel_h 1
#include "globals.hh"
#include "G4NeutronHPIsoData.hh"
#include "G4NeutronHPVector.hh"
#include "G4Material.hh"
#include "G4Track.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronElasticProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4StableIsotopes.hh"
#include "G4NeutronHPCaptureFS.hh"
#include "G4NeutronHPFinalState.hh"
#include "G4Element.hh"

class G4NeutronHPChannel
{
public:

  G4NeutronHPChannel()
  {
    theChannelData = new G4NeutronHPVector; 
    theBuffer = NULL;
    theIsotopeWiseData = NULL;
    theFinalStates = NULL;
    active = NULL;
    registerCount = -1;
  }
  
  ~G4NeutronHPChannel()
  {
    delete theChannelData; 
    G4int i;
    // Following statement disabled to avoid SEGV
    // theBuffer is also deleted as "theChannelData" in
    // ~G4NeutronHPIsoData.  FWJ 06-Jul-1999
    //if(theBuffer != NULL) delete theBuffer; 
    if(theIsotopeWiseData != NULL) delete [] theIsotopeWiseData;
    // Deletion of FinalStates disabled to avoid endless looping
    // in the destructor heirarchy.  FWJ 06-Jul-1999
    //if(theFinalStates != NULL)
    //{
    //  for(i=0; i<niso; i++)
    //  {
    //    delete theFinalStates[i];
    //  }
    //  delete [] theFinalStates;
    //}
    // FWJ experiment
    //if(active!=NULL) delete [] active;
    
  }
  
  G4double GetXsec(G4double energy);
  
  G4double GetWeightedXsec(G4double energy, G4int isoNumber);
  
  G4double GetFSCrossSection(G4double energy, G4int isoNumber);
  
  inline G4bool IsActive(G4int isoNumber) { return active[isoNumber]; }
  
  inline G4bool HasFSData(G4int isoNumber) { return theFinalStates[isoNumber]->HasFSData(); }
  
  inline G4bool HasAnyData(G4int isoNumber) { return theFinalStates[isoNumber]->HasAnyData(); }
  
  G4bool Register(G4NeutronHPFinalState *theFS);
  
  void Init(G4Element * theElement, const G4String dirName); 

  void Init(G4Element * theElement, const G4String dirName, const G4String fsType); 
  
  void UpdateData(G4int A, G4int Z, G4int index, G4double abundance);
  
  void Harmonise(G4NeutronHPVector *& theStore, G4NeutronHPVector * theNew);

  G4ParticleChange * ApplyYourself(const G4Track & theTrack, G4int isoNumber=-1);
    
  inline G4int GetNiso() {return niso;}
  
  inline G4bool HasDataInAnyFinalState()
  {
    G4bool result = false;
    G4int i;
    for(i=0; i<niso; i++)
    {
      if(theFinalStates[i]->HasAnyData()) result = true;
    }
    return result;
  }
  
private:
  G4NeutronHPVector * theChannelData;  // total (element) cross-section for this channel
  G4NeutronHPVector * theBuffer;
  
  G4NeutronHPIsoData * theIsotopeWiseData; // these are the isotope-wise cross-sections for each final state.
  G4NeutronHPFinalState ** theFinalStates; // also these are isotope-wise pionters, parallel to the above.
  G4bool * active;
  G4int niso;

  G4StableIsotopes theStableOnes;
  
  G4String theDir;
  G4String theFSType;
  G4Element * theElement;
  
  G4int registerCount;
    
};

#endif
