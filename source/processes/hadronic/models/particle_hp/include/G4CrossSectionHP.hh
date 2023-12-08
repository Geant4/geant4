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
// V. Ivanchenko September 2023 
//               
// G4CrossSectionHP is a generic class implementing 
// cross sections for neutrons, protons and light ions
// It is an alternative to code developed by J.P. Wellisch & T.Koi
//

#ifndef G4CrossSectionHP_h
#define G4CrossSectionHP_h 1

#include "G4VCrossSectionDataSet.hh"
#include "globals.hh"
#include "G4ElementData.hh"
#include "G4PhysicsVector.hh"
#include "G4LorentzVector.hh"
#include "G4ThreeVector.hh"
#include <vector>

class G4DynamicParticle;
class G4ParticleDefinition;
class G4ParticleHPManager;
class G4Element;
class G4Material;

class G4CrossSectionHP : public G4VCrossSectionDataSet
{
  public:
    explicit G4CrossSectionHP(const G4ParticleDefinition*,
                              const G4String& nameData,
                              const G4String& nameDir, G4double emaxHP,
                              G4int zmin, G4int zmax);

    ~G4CrossSectionHP() override = default;

    G4bool IsIsoApplicable(const G4DynamicParticle*, G4int Z, G4int A,
                           const G4Element*, const G4Material*) override;

    G4double GetIsoCrossSection(const G4DynamicParticle*, G4int Z, G4int A,
                                const G4Isotope*, const G4Element*,
                                const G4Material*) override;

    G4double ComputeIsoCrossSection(G4double kinEnergy, G4double loge,
                                    const G4ParticleDefinition*,
                                    G4int Z, G4int A,
                                    const G4Isotope* iso,
                                    const G4Element* elm,
                                    const G4Material* mat) override;

    const G4Isotope* SelectIsotope(const G4Element*, G4double kinEnergy,
                                   G4double logE) override;

    void BuildPhysicsTable(const G4ParticleDefinition&) override;

    void DumpPhysicsTable(const G4ParticleDefinition&) override;

    G4CrossSectionHP & operator=(const G4CrossSectionHP &right) = delete;
    G4CrossSectionHP(const G4CrossSectionHP&) = delete;

  protected:

    inline G4double GetMaxHPEnergy() const;

  private:

    void Initialise(const G4int Z);

    void InitialiseOnFly(const G4int Z);

    void PrepareCache(const G4Material*);

    G4double IsoCrossSection(const G4double kinE, const G4double loge,
			     const G4int Z, const G4int A,
                             const G4double temperature);

    inline G4double GetCrossSection(const G4int Z, const G4int A);

    inline G4bool CheckCache(const G4int Z);

    const G4ParticleDefinition* fParticle;
    G4ParticleHPManager* fManagerHP;

    const G4double emax;
    const G4double emaxT;
    const G4double elimit;
    const G4double logElimit;
    G4LorentzVector fLV;
    G4ThreeVector fBoost;

    const G4int minZ;
    const G4int maxZ;
    std::size_t index{0};
    G4bool fPrinted{false};

    // cache
    const G4Material* fCurrentMat{nullptr};
    std::vector<std::pair<G4int, G4int> > fZA;
    std::vector<G4double> fIsoXS;
    std::vector<G4double> fTemp;

    const G4String fDataName;
    const G4String fDataDirectory;
    G4ElementData* fData{nullptr};
};

inline G4double
G4CrossSectionHP::GetCrossSection(const G4int Z, const G4int A)
{
  G4double res = 0.0;
  for (std::size_t i=0; i<fZA.size(); ++i) {
    if (Z == fZA[i].first && A == fZA[i].second) {
      res = fIsoXS[i];
      break;
    }
  }
  return res;
}

inline G4bool G4CrossSectionHP::CheckCache(const G4int Z)
{
  G4bool res = false;
  for (auto const & p : fZA) {
    if (Z == p.first) {
      res = true;
      break;
    }
  }
  return res;
}

inline G4double G4CrossSectionHP::GetMaxHPEnergy() const
{ 
  return emax;
}

#endif
