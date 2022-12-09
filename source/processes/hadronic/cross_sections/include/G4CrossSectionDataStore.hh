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
// -------------------------------------------------------------------
// File name:     G4CrossSectionDataStore
//
// Modifications:
// 23.01.2009 V.Ivanchenko move constructor and destructor to source,
//                         use STL vector instead of C-array
//
// August 2011  Re-designed
//              by G. Folger, V. Ivantchenko, T. Koi and D.H. Wright

// Class Description
// This is the class to which cross section data sets may be registered. 
// An instance of it is contained in each hadronic process, allowing
// the use of the AddDataSet() method to tailor the cross sections to
// your application.
// Class Description - End

#ifndef G4CrossSectionDataStore_h
#define G4CrossSectionDataStore_h 1

#include "globals.hh"
#include "G4VCrossSectionDataSet.hh"
#include "G4DynamicParticle.hh"
#include "G4PhysicsVector.hh"
#include <vector>
#include <iostream>

class G4Nucleus;
class G4ParticleDefinition;
class G4Isotope;
class G4Element;
class G4Material;
class G4NistManager;

class G4CrossSectionDataStore
{
public:

  G4CrossSectionDataStore();

  ~G4CrossSectionDataStore() = default;

  // Run time cross section per unit volume
  inline G4double GetCrossSection(const G4DynamicParticle*, const G4Material*);
  G4double ComputeCrossSection(const G4DynamicParticle*, const G4Material*);

  // Cross section per element
  G4double GetCrossSection(const G4DynamicParticle*, 
			   const G4Element*, const G4Material*);

  // Cross section per isotope 
  G4double GetCrossSection(const G4DynamicParticle*, G4int Z, G4int A,
                           const G4Isotope*,
			   const G4Element*, const G4Material*);

  // Sample Z and A of a target nucleus and upload into G4Nucleus
  const G4Element* SampleZandA(const G4DynamicParticle*, const G4Material*,
			       G4Nucleus& target);

  // Initialisation before run
  void BuildPhysicsTable(const G4ParticleDefinition&);

  // Dump store to G4cout
  void DumpPhysicsTable(const G4ParticleDefinition&);

  // Dump store as html
  void DumpHtml(const G4ParticleDefinition&, std::ofstream&) const;
  void PrintCrossSectionHtml(const G4VCrossSectionDataSet *cs) const;
  
  void AddDataSet(G4VCrossSectionDataSet*);
  void AddDataSet(G4VCrossSectionDataSet*, std::size_t);
  inline const std::vector<G4VCrossSectionDataSet*>& GetDataSetList() const;

  inline void SetVerboseLevel(G4int value);

  // may be used by special processes
  inline void SetForcedElement(const G4Element*);

  G4CrossSectionDataStore & operator=
  (const G4CrossSectionDataStore &right) = delete;
  G4CrossSectionDataStore(const G4CrossSectionDataStore&) = delete;

private:

  G4double GetIsoCrossSection(const G4DynamicParticle*, G4int Z, G4int A,
			      const G4Isotope*, const G4Element*,
                              const G4Material*, const G4int index);

  G4String HtmlFileName(const G4String & in) const;

  G4NistManager* nist;
  const G4Material* currentMaterial = nullptr;
  const G4ParticleDefinition* matParticle = nullptr;
  const G4Element* forcedElement = nullptr;
  G4double matKinEnergy = 0.0;
  G4double matCrossSection = 0.0;

  G4int nDataSetList = 0;
  G4int verboseLevel = 1;

  std::vector<G4VCrossSectionDataSet*> dataSetList;
  std::vector<G4double> xsecelm;
  std::vector<G4double> xseciso;
};

inline void G4CrossSectionDataStore::SetVerboseLevel(G4int value)
{
  verboseLevel = value;
}

inline void G4CrossSectionDataStore::SetForcedElement(const G4Element* ptr)
{
  forcedElement = ptr;
}

inline const std::vector<G4VCrossSectionDataSet*>&
G4CrossSectionDataStore::GetDataSetList() const
{
  return dataSetList;
}

inline G4double 
G4CrossSectionDataStore::GetCrossSection(const G4DynamicParticle* dp,
                                         const G4Material* mat)
{
  if(dp->GetKineticEnergy() != matKinEnergy || mat != currentMaterial ||
     dp->GetDefinition() != matParticle) {
    ComputeCrossSection(dp, mat);
  }
  return matCrossSection;
}

#endif
