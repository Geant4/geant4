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
// Builds and has the Cross-section data for one material.

// Class Description
// Final state production model for a high precision (based on evaluated data
// libraries) description of neutron capture below 20 MeV;
// To be used in your physics list in case you need this physics.
// In this case you want to register an object of this class with
// the corresponding process.
// Class Description - End

// P. Arce, June-2014 Conversion neutron_hp to particle_hp
// V. Ivanchenko July-2023 converted back
//
#ifndef G4NeutronHPCapture_h
#define G4NeutronHPCapture_h 1

#include "G4HadronicInteraction.hh"
#include "G4ParticleHPChannel.hh"
#include "globals.hh"

class G4NeutronHPCapture : public G4HadronicInteraction
{
  public:
    G4NeutronHPCapture();
    ~G4NeutronHPCapture() override;

    G4HadFinalState* ApplyYourself(const G4HadProjectile& aTrack,
                                   G4Nucleus& aTargetNucleus) override;

    const std::pair<G4double, G4double> GetFatalEnergyCheckLevels() const override;
    G4int GetVerboseLevel() const;
    void SetVerboseLevel(G4int);
    void BuildPhysicsTable(const G4ParticleDefinition&) override;
    void ModelDescription(std::ostream& outFile) const override;

  private:
    std::vector<G4ParticleHPChannel*>* theCapture{nullptr};
    G4String dirName;
    G4int numEle{0};

    G4HadFinalState theResult;
};

#endif
