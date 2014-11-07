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
 // Builds and has the Cross-section data for one material.
 
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#ifndef G4ParticleHPChannelList_h
#define G4ParticleHPChannelList_h 1

#include "globals.hh"
#include "G4ParticleHPChannel.hh"
#include "G4StableIsotopes.hh"

class G4Element;
class G4HadFinalState;
class G4HadProjectile;
class G4ParticleHPFinalState;
class G4ParticleDefinition;

class G4ParticleHPChannelList
{      

  public:
  
  G4ParticleHPChannelList(G4int n, G4ParticleDefinition* projectile);
  
  G4ParticleHPChannelList();
  
  void Init(G4int n);

  ~G4ParticleHPChannelList();
  
  G4HadFinalState * ApplyYourself(const G4Element * theElement, const G4HadProjectile & aTrack);
      
  void Init(G4Element * anElement, const G4String & dirName, G4ParticleDefinition* projectile);
    
  void Register(G4ParticleHPFinalState *theFS, const G4String &aName);
  
  inline G4double GetXsec(G4double anEnergy)
  {
    G4double result=0;
    G4int i;
    for(i=0; i<nChannels; i++)
    {
      result+=std::max(0., theChannels[i]->GetXsec(anEnergy));
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

  void DumpInfo();


private:

  static G4ThreadLocal G4int trycounter;
  G4ParticleHPChannel ** theChannels;
  G4int nChannels;
  G4String theDir;
  G4Element * theElement; 
  
  G4bool allChannelsCreated;
  G4int theInitCount;

  G4StableIsotopes theStableOnes;

  G4ParticleDefinition* theProjectile;  

  G4HadFinalState unChanged;

};

#endif
