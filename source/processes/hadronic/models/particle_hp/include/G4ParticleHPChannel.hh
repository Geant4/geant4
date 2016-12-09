//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
//
 // Hadronic Process: Very Low Energy Neutron X-Sections
 // original by H.P. Wellisch, TRIUMF, 14-Feb-97
 // Builds and has the Cross-section data for one element and channel.
//
// Bug fixes and workarounds in the destructor, F.W.Jones 06-Jul-1999
// 070612 Fix memory leaking by T. Koi
//
// 080520 Delete unnecessary dependencies by T. Koi
 
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#ifndef G4ParticleHPChannel_h
#define G4ParticleHPChannel_h 1
#include "globals.hh"
#include "G4ParticleHPIsoData.hh"
#include "G4ParticleHPVector.hh"
#include "G4Material.hh"
#include "G4HadProjectile.hh"
//#include "G4NeutronInelasticProcess.hh"
//#include "G4HadronFissionProcess.hh"
//#include "G4HadronElasticProcess.hh"
//#include "G4HadronCaptureProcess.hh"
#include "G4StableIsotopes.hh"
#include "G4ParticleHPCaptureFS.hh"
#include "G4ParticleHPFinalState.hh"
#include "G4Element.hh"
#include "G4WendtFissionFragmentGenerator.hh"
class G4ParticleDefinition;

class G4ParticleHPChannel
{
public:

  G4ParticleHPChannel(G4ParticleDefinition* projectile)
  : wendtFissionGenerator(getenv("G4NEUTRON_HP_USE_WENDT_FISSION_MODEL") == NULL ? NULL : G4WendtFissionFragmentGenerator::GetInstance())
  {
    theProjectile = const_cast<G4ParticleDefinition*>(projectile);
    theChannelData = new G4ParticleHPVector; 
    theBuffer = 0;
    theIsotopeWiseData = 0;
    theFinalStates = 0;
    active = 0;
    registerCount = -1;
    niso = -1;
    theElement = NULL;
  }

  G4ParticleHPChannel()
  : wendtFissionGenerator(getenv("G4NEUTRON_HP_USE_WENDT_FISSION_MODEL") == NULL ? NULL : G4WendtFissionFragmentGenerator::GetInstance())
  {
    theProjectile = G4Neutron::Neutron();
    theChannelData = new G4ParticleHPVector; 
    theBuffer = 0;
    theIsotopeWiseData = 0;
    theFinalStates = 0;
    active = 0;
    registerCount = -1;
    niso = -1;
    theElement = NULL;
  }

  ~G4ParticleHPChannel()
  {
    delete theChannelData; 
    // Following statement disabled to avoid SEGV
    // theBuffer is also deleted as "theChannelData" in
    // ~G4ParticleHPIsoData.  FWJ 06-Jul-1999
    //if(theBuffer != 0) delete theBuffer; 
    if(theIsotopeWiseData != 0) delete [] theIsotopeWiseData;
    // Deletion of FinalStates disabled to avoid endless looping
    // in the destructor heirarchy.  FWJ 06-Jul-1999
    //if(theFinalStates != 0)
    //{
    //  for(i=0; i<niso; i++)
    //  {
    //    delete theFinalStates[i];
    //  }
    //  delete [] theFinalStates;
    //}
    // FWJ experiment
    //if(active!=0) delete [] active;
// T.K. 
   if ( theFinalStates != 0 )
   {
      for ( G4int i = 0 ; i < niso ; i++ )
      {
         delete theFinalStates[i];
      }
      delete [] theFinalStates;
   }
   if ( active != 0 ) delete [] active;
    
  }
  
  G4double GetXsec(G4double energy);
  
  G4double GetWeightedXsec(G4double energy, G4int isoNumber);
  
  G4double GetFSCrossSection(G4double energy, G4int isoNumber);
  
  inline G4bool IsActive(G4int isoNumber) { return active[isoNumber]; }
  
  inline G4bool HasFSData(G4int isoNumber) { return theFinalStates[isoNumber]->HasFSData(); }
  
  inline G4bool HasAnyData(G4int isoNumber) { return theFinalStates[isoNumber]->HasAnyData(); }
  
  G4bool Register(G4ParticleHPFinalState *theFS);
  
  void Init(G4Element * theElement, const G4String dirName); 

  void Init(G4Element * theElement, const G4String dirName, const G4String fsType); 
  
  //void UpdateData(G4int A, G4int Z, G4int index, G4double abundance);
  void UpdateData(G4int A, G4int Z, G4int index, G4double abundance, G4ParticleDefinition* projectile) { G4int M = 0; UpdateData( A, Z, M, index, abundance, projectile); };
  void UpdateData(G4int A, G4int Z, G4int M, G4int index, G4double abundance, G4ParticleDefinition* projectile);
  
  void Harmonise(G4ParticleHPVector *& theStore, G4ParticleHPVector * theNew);

  G4HadFinalState * ApplyYourself(const G4HadProjectile & theTrack, G4int isoNumber=-1);
    
  inline G4int GetNiso() {return niso;}
  
  inline G4double GetN(G4int i) {return theFinalStates[i]->GetN();}
  inline G4double GetZ(G4int i) {return theFinalStates[i]->GetZ();}
  inline G4double GetM(G4int i) {return theFinalStates[i]->GetM();}
  
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

  void DumpInfo();

  G4String GetFSType() const {
    return theFSType;
  }

  G4ParticleHPFinalState ** GetFinalStates() const {
    return theFinalStates; 
  }
  
private:
  G4ParticleDefinition* theProjectile;

  G4ParticleHPVector * theChannelData;  // total (element) cross-section for this channel
  G4ParticleHPVector * theBuffer;
  
  G4ParticleHPIsoData * theIsotopeWiseData; // these are the isotope-wise cross-sections for each final state.
  G4ParticleHPFinalState ** theFinalStates; // also these are isotope-wise pionters, parallel to the above.
  G4bool * active;
  G4int niso;

  G4StableIsotopes theStableOnes;
  
  G4String theDir;
  G4String theFSType;
  G4Element * theElement;
  
  G4int registerCount;

  G4WendtFissionFragmentGenerator* const wendtFissionGenerator;
    
};

#endif
