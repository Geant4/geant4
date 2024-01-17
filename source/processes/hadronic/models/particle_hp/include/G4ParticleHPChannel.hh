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
// Hadronic Process: Very Low Energy Neutron X-Sections
// original by H.P. Wellisch, TRIUMF, 14-Feb-97
// Builds and has the Cross-section data for one element and channel.
//
// Bug fixes and workarounds in the destructor, F.W.Jones 06-Jul-1999
// 070612 Fix memory leaking by T. Koi
//
// 080520 Delete unnecessary dependencies by T. Koi

// P. Arce, June-2014 Conversion neutron_hp to particle_hp
// V. Ivanchenko, July-2023 Basic revision of particle HP classes
//

#ifndef G4ParticleHPChannel_h
#define G4ParticleHPChannel_h 1

#include "G4Element.hh"
#include "G4HadProjectile.hh"
#include "G4Material.hh"
#include "G4ParticleHPCaptureFS.hh"
#include "G4ParticleHPFinalState.hh"
#include "G4ParticleHPIsoData.hh"
#include "G4ParticleHPManager.hh"
#include "G4ParticleHPVector.hh"
#include "G4StableIsotopes.hh"
#include "G4WendtFissionFragmentGenerator.hh"
#include "globals.hh"

class G4ParticleDefinition;

class G4ParticleHPChannel
{
public:
  G4ParticleHPChannel(G4ParticleDefinition* projectile = nullptr);

  ~G4ParticleHPChannel();

  G4double GetXsec(G4double energy);

  G4double GetWeightedXsec(G4double energy, G4int isoNumber);

  G4double GetFSCrossSection(G4double energy, G4int isoNumber);

  G4bool IsActive(G4int isoNumber) const
  {
    return active[isoNumber];
  }

  G4bool HasFSData(G4int isoNumber)
  { 
    return theFinalStates[isoNumber]->HasFSData();
  }

  G4bool HasAnyData(G4int isoNumber)
  {
    return theFinalStates[isoNumber]->HasAnyData();
  }

  G4bool Register(G4ParticleHPFinalState* theFS);

  void Init(G4Element* theElement, const G4String& dirName);

  void Init(G4Element* theElement, const G4String& dirName,
	    const G4String& fsType);

  void UpdateData(G4int A, G4int Z, G4int index, G4double abundance,
		  G4ParticleDefinition* projectile)
  {
    UpdateData(A, Z, 0, index, abundance, projectile);
  }

  void UpdateData(G4int A, G4int Z, G4int M, G4int index, G4double abundance,
		  G4ParticleDefinition* projectile);

  void Harmonise(G4ParticleHPVector*& theStore, G4ParticleHPVector* theNew);

  G4HadFinalState* ApplyYourself(const G4HadProjectile& theTrack,
				 G4int isoNumber = -1,
				 G4bool isElastic = false);

  G4int GetNiso() const { return niso; }

  G4double GetN(G4int i) const { return theFinalStates[i]->GetN(); }
  G4double GetZ(G4int i) const { return theFinalStates[i]->GetZ(); }
  G4double GetM(G4int i) const { return theFinalStates[i]->GetM(); }

  G4bool HasDataInAnyFinalState()
  {
    G4bool result = false;
    for (G4int i = 0; i < niso; ++i) {
      if (theFinalStates[i]->HasAnyData()) {
	result = true;
	break;
      }
    }
    return result;
  }

  void DumpInfo();

  G4String& GetFSType() { return theFSType; }

  G4ParticleHPFinalState** GetFinalStates() const { return theFinalStates; }

  G4ParticleHPChannel(G4ParticleHPChannel &) = delete;
  G4ParticleHPChannel & operator=
  (const G4ParticleHPChannel &right) = delete;

protected:
  G4ParticleHPManager* fManager;

private:
  G4ParticleDefinition* theProjectile;

  G4ParticleHPVector* theChannelData;
  G4Element* theElement{nullptr};

  // total (element) cross-section for this channel
  G4ParticleHPVector* theBuffer{nullptr};

  G4ParticleHPIsoData* theIsotopeWiseData{nullptr};
  // these are the isotope-wise cross-sections for each final state.
  G4ParticleHPFinalState** theFinalStates{nullptr};
  // also these are isotope-wise pionters, parallel to the above.

  G4WendtFissionFragmentGenerator* wendtFissionGenerator{nullptr};
  G4bool* active{nullptr};
  G4int niso{-1};
  G4int registerCount{-1};

  G4String theDir{""};
  G4String theFSType{""};
};

#endif
