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
// Hadronic Process: High Precision low E neutron tracking
// original by H.P. Wellisch, TRIUMF, 14-Feb-97
// Builds and has the Cross-section data for one material.
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
// V. Ivanchenko, July-2023 Basic revision of particle HP classes
//

#ifndef G4ParticleHPInelastic_h
#define G4ParticleHPInelastic_h 1

// Final state production model for a high precision (based on evaluated data
// libraries) description of neutron inelastic scattering below 20 MeV or
// light ion inelastic interaction below 100 MeV.
// 36 exclusive final states are consideded.

#include "G4HadronicInteraction.hh"
#include "G4ParticleHPChannel.hh"
#include "G4ParticleHPChannelList.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"

class G4ParticleHPInelastic : public G4HadronicInteraction
{
public:
  G4ParticleHPInelastic(G4ParticleDefinition* p = G4Neutron::Neutron(),
			const char* name = "NeutronHPInelastic");

  ~G4ParticleHPInelastic() override;

  G4HadFinalState* ApplyYourself(const G4HadProjectile& aTrack,
				 G4Nucleus& aTargetNucleus) override;

  const std::pair<G4double, G4double>
  GetFatalEnergyCheckLevels() const override;

  void BuildPhysicsTable(const G4ParticleDefinition&) override;
  void ModelDescription(std::ostream& outFile) const override;

  G4ParticleHPInelastic(G4ParticleHPInelastic &) = delete;
  G4ParticleHPInelastic & operator=
  (const G4ParticleHPInelastic &right) = delete;

private:

  void ClearData();

  G4ParticleDefinition* theProjectile;
  G4bool isFirst{false};
  static G4bool fLock[6];

protected:
  // one List per element
  static std::vector<G4ParticleHPChannelList*>* theInelastic[6];
  G4ParticleHPManager* fManager;
  G4String dirName;
  G4int numEle{0};
  G4int indexP;
};

#endif
