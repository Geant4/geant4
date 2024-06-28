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
// author: Vladimir.Grichine@cern.ch
//
// Implements data from: Barashenkov V.S., Nucleon-Nucleus Cross Section,
// Preprint JINR P2-89-770, p. 12, Dubna 1989 (scanned version from KEK)
// Based on G. Folger version of G4PiNuclearCrossSection class
//
// Modified:
// 05.03.2024 V.Ivanchenko removed obsolete methods and calls
//

#ifndef G4NucleonNuclearCrossSection_h
#define G4NucleonNuclearCrossSection_h

#include "G4VCrossSectionDataSet.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"

class G4ComponentBarNucleonNucleusXsc;
class G4ParticleDefinition;

class G4NucleonNuclearCrossSection : public G4VCrossSectionDataSet
{
public:
  
  G4NucleonNuclearCrossSection();
  ~G4NucleonNuclearCrossSection() override = default;
    
  static const char* Default_Name() { return "BarashenkovNucleonXS"; }

  G4bool IsElementApplicable(const G4DynamicParticle* aParticle,
			     G4int Z, const G4Material* mat) final;

  // return inelastic x-section
  G4double GetElementCrossSection(const G4DynamicParticle* aParticle, 
                                  G4int Z, const G4Material* mat=nullptr) final;

  void CrossSectionDescription(std::ostream&) const final;

  // return elastic x-section
  inline G4double GetElasticCrossSection(const G4DynamicParticle* aParticle,
					 G4int Z);

  // access methods should be called after ComputeCrossSection(...)
  inline G4double GetTotalXsc()     { return fTotalXsc; };
  inline G4double GetInelasticXsc() { return fInelasticXsc; };
  inline G4double GetElasticXsc()   { return fElasticXsc; };

  G4NucleonNuclearCrossSection& operator=
  (const G4NucleonNuclearCrossSection &right) = delete;
  G4NucleonNuclearCrossSection(const G4NucleonNuclearCrossSection&) = delete;

private:

  void ComputeCrossSections(const G4ParticleDefinition*, 
                            G4double kinEnergy, G4int Z);

  G4ComponentBarNucleonNucleusXsc* fBarash;

  G4double fTotalXsc{0.0};
  G4double fInelasticXsc{0.0};
  G4double fElasticXsc{0.0};
};

inline
G4double G4NucleonNuclearCrossSection::GetElasticCrossSection(
         const G4DynamicParticle* dp, G4int Z)
{
  ComputeCrossSections(dp->GetDefinition(), dp->GetKineticEnergy(), Z);
  return fElasticXsc;
}

#endif
