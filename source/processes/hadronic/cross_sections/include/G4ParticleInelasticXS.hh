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
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:    G4ParticleInelasticXS
//
// Author  Ivantchenko, Geant4, 24 May 2018
//

// Class Description:
// This is a base class for n,p,d,t,he3,he4 inelastic hadronic cross 
// section based on data files from G4PARTICLEXSDATA data set
// and Glauber-Gribov model for high energy 
//
 
#ifndef G4ParticleInelasticXS_h
#define G4ParticleInelasticXS_h 1

#include "G4VCrossSectionDataSet.hh"
#include "globals.hh"
#include "G4ElementData.hh"
#include "G4PhysicsVector.hh"
#include <vector>

class G4DynamicParticle;
class G4ParticleDefinition;
class G4Element;
class G4VComponentCrossSection;

class G4ParticleInelasticXS final : public G4VCrossSectionDataSet
{
public: 

  explicit G4ParticleInelasticXS(const G4ParticleDefinition*);

  ~G4ParticleInelasticXS() final;

  G4bool IsElementApplicable(const G4DynamicParticle*, G4int Z,
			     const G4Material* mat = nullptr) final;

  G4bool IsIsoApplicable(const G4DynamicParticle*, G4int Z, G4int A,
			 const G4Element* elm = nullptr,
                         const G4Material* mat = nullptr) final;

  G4double GetElementCrossSection(const G4DynamicParticle*, G4int Z,
                                  const G4Material* mat = nullptr) final;

  G4double ComputeCrossSectionPerElement(G4double kinEnergy, G4double loge,
                                         const G4ParticleDefinition*,
                                         const G4Element*,
                                         const G4Material*) final;
  
  G4double ComputeIsoCrossSection(G4double kinEnergy, G4double loge,
                                  const G4ParticleDefinition*,
                                  G4int Z, G4int A,
                                  const G4Isotope* iso,
                                  const G4Element* elm,
                                  const G4Material* mat) final;

  G4double GetIsoCrossSection(const G4DynamicParticle*, G4int Z, G4int A,
                              const G4Isotope* iso = nullptr,
                              const G4Element* elm = nullptr,
                              const G4Material* mat = nullptr) final;

  const G4Isotope* SelectIsotope(const G4Element*, 
                                 G4double kinEnergy, G4double logE) final;

  void BuildPhysicsTable(const G4ParticleDefinition&) final;

  void CrossSectionDescription(std::ostream&) const final;

  G4double ElementCrossSection(G4double kinEnergy, G4double loge, G4int Z);

  G4double IsoCrossSection(G4double ekin, G4double logE, G4int Z, G4int A);

  G4ParticleInelasticXS & operator=
  (const G4ParticleInelasticXS &right) = delete;
  G4ParticleInelasticXS(const G4ParticleInelasticXS&) = delete;

private: 

  void Initialise(G4int Z);

  void InitialiseOnFly(G4int Z);

  const G4String& FindDirectoryPath();

  inline const G4PhysicsVector* GetPhysicsVector(G4int Z);

  G4PhysicsVector* RetrieveVector(std::ostringstream& in, G4bool warn);

  G4VComponentCrossSection* highEnergyXsection = nullptr;
  const G4ParticleDefinition* particle;

  std::vector<G4double> temp;
  G4double elimit;

  G4int index = 0;
  G4bool isMaster = false;

  static const G4int MAXZINELP = 93;
  static G4ElementData* data[5];
  static G4double coeff[MAXZINELP][5];
  static G4String gDataDirectory[5];
};

inline
const G4PhysicsVector* G4ParticleInelasticXS::GetPhysicsVector(G4int Z)
{
  const G4PhysicsVector* pv = data[index]->GetElementData(Z);
  if(pv == nullptr) { 
    InitialiseOnFly(Z);
    pv = data[index]->GetElementData(Z);
  }
  return pv;
}

#endif
