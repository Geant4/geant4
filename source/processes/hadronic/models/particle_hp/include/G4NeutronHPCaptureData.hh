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
// 080417 Add IsZAApplicable method (return false) by T. Koi
// 080428 Add bool onFlightDB by T. Koi
// 091118 Add Ignore and Enable On Flight Doppler Broadening methods by T. Koi
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
// V. Ivanchenko July-2023 converted back
//
#ifndef G4NeutronHPCaptureData_h
#define G4NeutronHPCaptureData_h 1

// Class Description
// Cross-section data set for a high precision (based on evaluated data
// libraries) description of neutron capture below 20 MeV;
// To be used in your physics list in case you need this physics.
// In this case you want to register an object of this class with
// the corresponding process.
// Class Description - End

#include "G4DynamicParticle.hh"
#include "G4VCrossSectionDataSet.hh"

class G4ParticleHPManager;
class G4ParticleDefinition;
class G4PhysicsTable;
class G4Element;
class G4Material;

class G4NeutronHPCaptureData : public G4VCrossSectionDataSet
{
  public:
    G4NeutronHPCaptureData();
    ~G4NeutronHPCaptureData() override;

    G4bool IsIsoApplicable(const G4DynamicParticle*, G4int Z, G4int A,
                           const G4Element*, const G4Material*) override;

    G4double GetIsoCrossSection(const G4DynamicParticle*, G4int Z, G4int A,
                                const G4Isotope*, const G4Element*,
                                const G4Material*) override;

    G4double GetCrossSection(const G4DynamicParticle*, const G4Element*, G4double aT);

    void BuildPhysicsTable(const G4ParticleDefinition&) override;

    void DumpPhysicsTable(const G4ParticleDefinition&) override;

    void CrossSectionDescription(std::ostream&) const override;

  private:
    static G4PhysicsTable* theCrossSections;
    G4ParticleHPManager* fManager;
    G4double emax;

    G4double ke_cache{0.0};
    G4double xs_cache{0.0};
    const G4Element* element_cache{nullptr};
    const G4Material* material_cache{nullptr};

    G4bool isFirst{false};
    static G4bool fLock;
};

#endif
