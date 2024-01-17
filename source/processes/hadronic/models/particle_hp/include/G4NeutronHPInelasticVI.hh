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
// Geant4 header : G4NeutronHPInelasticVI
// Created:  15 October 2023
//
// Author  V.Ivanchenko
//
// Class Description:
// Neutron HP inelastic scattering model below 20 MeV
//

#ifndef G4NeutronHPInelasticVI_h
#define G4NeutronHPInelasticVI_h 1

// Final state production model for a high precision (based on evaluated data
// libraries) description of neutron inelastic scattering below 20 MeV or
// light ion inelastic interaction below 100 MeV.
// 36 exclusive final states are consideded.

#include "G4HadronicInteraction.hh"
#include "G4ParticleHPChannelList.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"

class G4NeutronHPInelasticVI : public G4HadronicInteraction
{
  public:
    G4NeutronHPInelasticVI();

    ~G4NeutronHPInelasticVI() override;

    G4HadFinalState* ApplyYourself(const G4HadProjectile& aTrack,
                                   G4Nucleus& aTargetNucleus) override;

    const std::pair<G4double, G4double> GetFatalEnergyCheckLevels() const override;

    void BuildPhysicsTable(const G4ParticleDefinition&) override;

    void ModelDescription(std::ostream& outFile) const override;

    G4NeutronHPInelasticVI(G4NeutronHPInelasticVI &) = delete;
    G4NeutronHPInelasticVI & operator=
    (const G4NeutronHPInelasticVI &right) = delete;

private:

    void InitialiseOnFly();

    void Initialise();

    G4ParticleHPManager* fManagerHP;
    G4bool fInitializer{false};

    static G4bool fLock;
    static constexpr G4int ZMAXHPI{101}; // 101 because Z range is 1-100
    static G4ParticleHPChannelList* theChannels[ZMAXHPI];
};

#endif
