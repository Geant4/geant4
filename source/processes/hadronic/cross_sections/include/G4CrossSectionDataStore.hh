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
// $Id: G4CrossSectionDataStore.hh,v 1.15 2010-04-29 14:46:08 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4CrossSectionDataStore
//
//
// Modifications:
// 23.01.2009 V.Ivanchenko move constructor and destructor to source,
//                         use STL vector instead of C-array
//

// Class Description
// This is the class to which cross section data sets may be registered. 
// An instance of it is contained in each hadronic process, allowing
// the use of the AddDataSet() method to tailor the cross sections to
// your application.
// Class Description - End

#ifndef G4CrossSectionDataStore_h
#define G4CrossSectionDataStore_h 1

#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4VCrossSectionDataSet.hh"
#include <vector>

class G4Nucleus;

class G4CrossSectionDataStore
{
public:

  G4CrossSectionDataStore();

  ~G4CrossSectionDataStore();

  G4double GetCrossSection(const G4DynamicParticle*, 
			   const G4Element*, G4double aTemperature);

  G4double GetCrossSection(const G4DynamicParticle*, 
			   const G4Isotope*, G4double aTemperature);

  G4double GetCrossSection(const G4DynamicParticle*, 
			   G4double Z, G4double A, G4double aTemperature);

  G4double GetCrossSection(const G4DynamicParticle*, 
			   G4int Z, G4int A, G4double aTemperature);

  // to replace GetMicroscopicCrossSection
  G4double GetCrossSection(const G4DynamicParticle*, const G4Material*);

  //std::pair<G4double/*Z*/, G4double/*A*/> 
  // SelectRandomIsotope(const G4DynamicParticle*, const G4Material*);

  G4Element* SampleZandA(const G4DynamicParticle*, const G4Material*,
			 G4Nucleus& target);

  void AddDataSet(G4VCrossSectionDataSet*);

  void BuildPhysicsTable(const G4ParticleDefinition&);

  void DumpPhysicsTable(const G4ParticleDefinition&);

  inline void SetVerboseLevel(G4int value)
  {
    verboseLevel = value;
  }

  inline G4int GetVerboseLevel()
  {
    return verboseLevel;
  }

private:

  G4VCrossSectionDataSet* whichDataSetInCharge(const G4DynamicParticle*, 
					       const G4Element*);

  std::vector<G4VCrossSectionDataSet*> DataSetList;
  G4int NDataSetList;
  G4int verboseLevel;
};

#endif
