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
#include "G4FastPathHadronicCrossSection.hh"
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

  ~G4CrossSectionDataStore();

  // Cross section per unit volume is computed (inverse mean free path)
  inline G4double GetCrossSection(const G4DynamicParticle*, const G4Material*);
  G4double ComputeCrossSection(const G4DynamicParticle*, const G4Material*);

  // Cross section per element is computed
  G4double GetCrossSection(const G4DynamicParticle*, 
			   const G4Element*, const G4Material*);

  // Cross section per isotope is computed
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
  void AddDataSet(G4VCrossSectionDataSet*,size_t);

  inline void SetVerboseLevel(G4int value);

private:

  G4double GetIsoCrossSection(const G4DynamicParticle*, G4int Z, G4int A,
			      const G4Isotope*,
			      const G4Element*, const G4Material* aMaterial,
			      G4int index);

  G4CrossSectionDataStore & operator=(const G4CrossSectionDataStore &right);
  G4CrossSectionDataStore(const G4CrossSectionDataStore&);

  G4String HtmlFileName(const G4String & in) const;

  G4NistManager* nist;

  std::vector<G4VCrossSectionDataSet*> dataSetList;
  std::vector<G4double> xsecelm;
  std::vector<G4double> xseciso;

  const G4Material* currentMaterial;
  const G4ParticleDefinition* matParticle;
  G4double matKinEnergy;
  G4double matCrossSection;

  G4int nDataSetList;
  G4int verboseLevel;
  //Fast path: caching
public:
  inline const G4FastPathHadronicCrossSection::fastPathParameters&
  GetFastPathParameters() const { return fastPathParams; }
  inline const G4FastPathHadronicCrossSection::controlFlag&
  GetFastPathControlFlags() const { return fastPathFlags; }
  void DumpFastPath( const G4ParticleDefinition* , const G4Material* , std::ostream& os);
  void ActivateFastPath( const G4ParticleDefinition*, const G4Material* , G4double);
private:
  friend struct G4FastPathHadronicCrossSection::fastPathEntry;
  //The following method is called by the public one GetCrossSection(const G4DynamicParticle*, const G4Material*)
  //The third parameter is used to force the calculation of cross-sections skipping the fast-path mechanism
  G4double GetCrossSection(const G4DynamicParticle*, const G4Material*, G4bool requiresSlowPath);
  G4FastPathHadronicCrossSection::controlFlag fastPathFlags;
  G4FastPathHadronicCrossSection::fastPathParameters fastPathParams;
  //Counters
  G4FastPathHadronicCrossSection::getCrossSectionCount counters;
  //TODO: share this among threads
  G4FastPathHadronicCrossSection::G4CrossSectionDataStore_Cache fastPathCache;
  G4FastPathHadronicCrossSection::timing timing;
  G4FastPathHadronicCrossSection::G4CrossSectionDataStore_Requests requests;
};

inline G4double G4CrossSectionDataStore::GetCrossSection(const G4DynamicParticle* particle , const G4Material* material ) {
	//By default tries to use the fast-path mechanism
	return GetCrossSection( particle , material , false);
}

inline void G4CrossSectionDataStore::SetVerboseLevel(G4int value)
{
  verboseLevel = value;
}

#endif
