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

#include "G4PionDecayMakeSpin.hh"

#include "G4Decay.hh"
#include "G4DecayProducts.hh"

#include "G4RandomDirection.hh"

// constructor

G4PionDecayMakeSpin::G4PionDecayMakeSpin(const G4String& processName)
                               : G4Decay(processName) 
{
  // set Process Sub Type
  SetProcessSubType(static_cast<int>(DECAY_PionMakeSpin));

}

G4PionDecayMakeSpin::~G4PionDecayMakeSpin() { }

void G4PionDecayMakeSpin::DaughterPolarization(const G4Track& aTrack,
                                   G4DecayProducts* products)
{
  //  This routine deals only with particles that can decay into a muon
  //                 pi+, pi-, K+, K- and K0_long

  // get particle

  const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
  const G4ParticleDefinition* aParticleDef = aParticle->GetDefinition();

  G4ParticleDefinition* aMuonPlus =
            G4ParticleTable::GetParticleTable()->FindParticle("mu+");
  G4ParticleDefinition* aMuonMinus =
            G4ParticleTable::GetParticleTable()->FindParticle("mu-");
  G4ParticleDefinition* aPionPlus = 
            G4ParticleTable::GetParticleTable()->FindParticle("pi+");
  G4ParticleDefinition* aPionMinus =
            G4ParticleTable::GetParticleTable()->FindParticle("pi-");
  G4ParticleDefinition* aKaonPlus =
            G4ParticleTable::GetParticleTable()->FindParticle("kaon+");
  G4ParticleDefinition* aKaonMinus =
            G4ParticleTable::GetParticleTable()->FindParticle("kaon-");
  G4ParticleDefinition* aKaon0Long =
            G4ParticleTable::GetParticleTable()->FindParticle("kaon0L");
  G4ParticleDefinition* aNeutrinoMu =
            G4ParticleTable::GetParticleTable()->FindParticle("nu_mu");
  G4ParticleDefinition* aAntiNeutrinoMu =
            G4ParticleTable::GetParticleTable()->FindParticle("anti_nu_mu");

  if( aParticleDef == aPionPlus   ||
      aParticleDef == aPionMinus  || 
      aParticleDef == aKaonPlus   ||
      aParticleDef == aKaonMinus  ||
      aParticleDef == aKaon0Long ) {
  } else {
      return;
  }

  G4DynamicParticle* aMuon = nullptr;

  G4double emu(0), eneutrino(0);
  G4ThreeVector p_muon, p_neutrino;

  G4int numberOfSecondaries = products->entries();

  if (numberOfSecondaries > 0) {
    for (G4int index=0; index < numberOfSecondaries; index++){
      G4DynamicParticle* aSecondary = (*products)[index];
      const G4ParticleDefinition* aSecondaryDef = aSecondary->GetDefinition();
      
      if (aSecondaryDef == aMuonPlus ||
	  aSecondaryDef == aMuonMinus ) {
	//          Muon+ or Muon-
	aMuon = aSecondary;
	emu    = aSecondary->GetTotalEnergy();
	p_muon = aSecondary->GetMomentum();
      } else if (aSecondaryDef == aNeutrinoMu ||
		 aSecondaryDef == aAntiNeutrinoMu ) {
	//          Muon-Neutrino / Muon-Anti-Neutrino
	eneutrino  = aSecondary->GetTotalEnergy();
	p_neutrino = aSecondary->GetMomentum();
      }
    }
  }

  //  This routine deals only with decays with a 
  //  muon and mu-(anti)neutrinos in the final state
  if (aMuon == nullptr) return;
  if (eneutrino==0||emu==0) return;

  G4ThreeVector spin(0,0,0);

  const G4DynamicParticle* theParentParticle = products->GetParentParticle();

  G4double amass = theParentParticle->GetMass();
  G4double emmu = aMuonPlus->GetPDGMass();

  if (numberOfSecondaries == 2 ) {
     G4double scale = - (eneutrino - ( p_muon * p_neutrino )/(emu+emmu));

     p_muon = scale * p_muon;
     p_neutrino = emmu * p_neutrino;
     spin = p_muon + p_neutrino;

     scale = 2./(amass*amass-emmu*emmu);
     spin = scale * spin;

     if (aParticle->GetCharge() < 0.0) spin = -spin;

  } else {
       spin = G4RandomDirection();

  }

  spin = spin.unit();

  aMuon->SetPolarization(spin.x(),spin.y(),spin.z());

  return; 
}

void G4PionDecayMakeSpin::ProcessDescription(std::ostream& outFile) const
{
  outFile << GetProcessName() 
	  << ": Decay of mesons that can decay into a muon \n"
	  << " i.e. pi+, pi-, K+, K- and K0_long \n"
	  << " kinematics of daughters are dertermined by DecayChannels \n"
	  << " polarization of daughter particles are take into account. \n";
}

