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
// $Id: G4CrossSectionDataStore.cc 78442 2013-12-18 15:04:33Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4CrossSectionDataStore
//
// Modifications:
// 23.01.2009 V.Ivanchenko add destruction of data sets
// 29.04.2010 G.Folger     modifictaions for integer A & Z
// 14.03.2011 V.Ivanchenko fixed DumpPhysicsTable
// 15.08.2011 G.Folger, V.Ivanchenko, T.Koi, D.Wright redesign the class
// 07.03.2013 M.Maire cosmetic in DumpPhysicsTable
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

#include "G4CrossSectionDataStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4HadronicException.hh"
#include "G4HadTmpUtil.hh"
#include "Randomize.hh"
#include "G4Nucleus.hh"

#include "G4DynamicParticle.hh"
#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4CrossSectionDataStore::G4CrossSectionDataStore() :
  nDataSetList(0), verboseLevel(0)
{
  nist = G4NistManager::Instance();
  currentMaterial = elmMaterial = 0;
  currentElement = 0;  //ALB 14-Aug-2012 Coverity fix.
  matParticle = elmParticle = 0;
  matKinEnergy = elmKinEnergy = matCrossSection = elmCrossSection = 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4CrossSectionDataStore::~G4CrossSectionDataStore()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4double
G4CrossSectionDataStore::GetCrossSection(const G4DynamicParticle* part,
                                         const G4Material* mat)
{
  if(mat == currentMaterial && part->GetDefinition() == matParticle
     && part->GetKineticEnergy() == matKinEnergy) 
    { return matCrossSection; }
  
  currentMaterial = mat;
  matParticle = part->GetDefinition();
  matKinEnergy = part->GetKineticEnergy();
  matCrossSection = 0;

  G4int nElements = mat->GetNumberOfElements();
  const G4double* nAtomsPerVolume = mat->GetVecNbOfAtomsPerVolume();

  if(G4int(xsecelm.size()) < nElements) { xsecelm.resize(nElements); }

  for(G4int i=0; i<nElements; ++i) {
    matCrossSection += nAtomsPerVolume[i] * 
      GetCrossSection(part, (*mat->GetElementVector())[i], mat);
    xsecelm[i] = matCrossSection;
  }
  return matCrossSection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4double
G4CrossSectionDataStore::GetCrossSection(const G4DynamicParticle* part,
                                         const G4Element* elm,
					 const G4Material* mat)
{
  if(mat == elmMaterial && elm == currentElement &&
     part->GetDefinition() == elmParticle &&
     part->GetKineticEnergy() == elmKinEnergy) 
    { return elmCrossSection; }

  elmMaterial = mat;
  currentElement = elm;
  elmParticle = part->GetDefinition();
  elmKinEnergy = part->GetKineticEnergy();
  elmCrossSection = 0.0;

  G4int i = nDataSetList-1;  
  G4int Z = G4lrint(elm->GetZ());
  if (elm->GetNaturalAbundanceFlag() &&
      dataSetList[i]->IsElementApplicable(part, Z, mat)) {

    // element wise cross section
    elmCrossSection = dataSetList[i]->GetElementCrossSection(part, Z, mat);

    //G4cout << "Element wise " << elmParticle->GetParticleName() 
    //	   << " xsec(barn)= " <<  elmCrossSection/barn 
    //	   << "  E(MeV)= " << elmKinEnergy/MeV 
    //	   << " Z= " << Z << " AbundFlag= " << elm->GetNaturalAbandancesFlag()
    //	   <<G4endl;

  } else {
    // isotope wise cross section
    G4int nIso = elm->GetNumberOfIsotopes();    
    G4Isotope* iso = 0;

    // user-defined isotope abundances        
    G4IsotopeVector* isoVector = elm->GetIsotopeVector();
    G4double* abundVector = elm->GetRelativeAbundanceVector();

    for (G4int j = 0; j<nIso; ++j) {
      if(abundVector[j] > 0.0) {
	iso = (*isoVector)[j];
	elmCrossSection += abundVector[j]*
	  GetIsoCrossSection(part, Z, iso->GetN(), iso, elm, mat, i);
	//G4cout << "Isotope wise " << elmParticle->GetParticleName() 
	//       << " xsec(barn)= " <<  elmCrossSection/barn 
	//       << "  E(MeV)= " << elmKinEnergy/MeV 
	//       << " Z= " << Z << "  A= " << iso->GetN() << "  j= " << j << G4endl;
      }
    }
  }
  //G4cout << "  E(MeV)= " << elmKinEnergy/MeV 
  //	 << "xsec(barn)= " <<  elmCrossSection/barn <<G4endl;
  return elmCrossSection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4double
G4CrossSectionDataStore::GetIsoCrossSection(const G4DynamicParticle* part,
					    G4int Z, G4int A, 
					    const G4Isotope* iso,
					    const G4Element* elm,
					    const G4Material* mat, 
					    G4int idx)
{
  // this methods is called after the check that dataSetList[idx] 
  // depend on isotopes, so for this DataSet only isotopes are checked

  // isotope-wise cross section does exist
  if(dataSetList[idx]->IsIsoApplicable(part, Z, A, elm, mat) ) {
    return dataSetList[idx]->GetIsoCrossSection(part, Z, A, iso, elm, mat);

  } else {
    // seach for other dataSet
    for (G4int j = nDataSetList-1; j >= 0; --j) { 
      if (dataSetList[j]->IsElementApplicable(part, Z, mat)) {
	return dataSetList[j]->GetElementCrossSection(part, Z, mat);
      } else if (dataSetList[j]->IsIsoApplicable(part, Z, A, elm, mat)) {
	return dataSetList[j]->GetIsoCrossSection(part, Z, A, iso, elm, mat);
      }
    }
  }
  G4cout << "G4CrossSectionDataStore::GetCrossSection ERROR: "
	 << " no isotope cross section found"
	 << G4endl;
  G4cout << "  for " << part->GetDefinition()->GetParticleName() 
	 << " off Element " << elm->GetName()
         << "  in " << mat->GetName() 
	 << " Z= " << Z << " A= " << A
	 << " E(MeV)= " << part->GetKineticEnergy()/MeV << G4endl; 
  throw G4HadronicException(__FILE__, __LINE__, 
                      " no applicable data set found for the isotope");
  return 0.0;
  //return dataSetList[idx]->ComputeCrossSection(part, elm, mat);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4double
G4CrossSectionDataStore::GetCrossSection(const G4DynamicParticle* part,
                                         G4int Z, G4int A,
					 const G4Isotope* iso,
                                         const G4Element* elm,
					 const G4Material* mat)
{
  for (G4int i = nDataSetList-1; i >= 0; --i) {
    if (dataSetList[i]->IsIsoApplicable(part, Z, A, elm, mat) ) {
      return dataSetList[i]->GetIsoCrossSection(part, Z, A, iso, elm, mat);
    }
  }
  G4cout << "G4CrossSectionDataStore::GetCrossSection ERROR: "
	 << " no isotope cross section found"
	 << G4endl;
  G4cout << "  for " << part->GetDefinition()->GetParticleName() 
	 << " off Element " << elm->GetName()
         << "  in " << mat->GetName() 
	 << " Z= " << Z << " A= " << A
	 << " E(MeV)= " << part->GetKineticEnergy()/MeV << G4endl; 
  throw G4HadronicException(__FILE__, __LINE__, 
                      " no applicable data set found for the isotope");
  return 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4Element*
G4CrossSectionDataStore::SampleZandA(const G4DynamicParticle* part, 
                                     const G4Material* mat,
				     G4Nucleus& target)
{
  G4int nElements = mat->GetNumberOfElements();
  const G4ElementVector* theElementVector = mat->GetElementVector();
  G4Element* anElement = (*theElementVector)[0];

  G4double cross = GetCrossSection(part, mat);

  // select element from a compound 
  if(1 < nElements) {
    cross *= G4UniformRand();
    for(G4int i=0; i<nElements; ++i) {
      if(cross <= xsecelm[i]) {
	anElement = (*theElementVector)[i];
        break;
      }
    }
  }

  G4int Z = G4lrint(anElement->GetZ());
  G4Isotope* iso = 0;

  G4int i = nDataSetList-1; 
  if (dataSetList[i]->IsElementApplicable(part, Z, mat)) {

    //----------------------------------------------------------------
    // element-wise cross section
    // isotope cross section is not computed
    //----------------------------------------------------------------
    G4int nIso = anElement->GetNumberOfIsotopes();
    if (0 >= nIso) { 
      G4cout << " Element " << anElement->GetName() << " Z= " << Z 
	     << " has no isotopes " << G4endl; 
      throw G4HadronicException(__FILE__, __LINE__, 
                      " Isotope vector is not defined");
      return anElement;
    }
    // isotope abundances        
    G4IsotopeVector* isoVector = anElement->GetIsotopeVector();
    iso = (*isoVector)[0];

    // more than 1 isotope
    if(1 < nIso) { 
      iso = dataSetList[i]->SelectIsotope(anElement, part->GetKineticEnergy());
    }

  } else {

    //----------------------------------------------------------------
    // isotope-wise cross section
    // isotope cross section is computed
    //----------------------------------------------------------------
    G4int nIso = anElement->GetNumberOfIsotopes();
    cross = 0.0;

    if (0 >= nIso) { 
      G4cout << " Element " << anElement->GetName() << " Z= " << Z 
	     << " has no isotopes " << G4endl; 
      throw G4HadronicException(__FILE__, __LINE__, 
                      " Isotope vector is not defined");
      return anElement;
    }

    // user-defined isotope abundances        
    G4IsotopeVector* isoVector = anElement->GetIsotopeVector();
    iso = (*isoVector)[0];

    // more than 1 isotope
    if(1 < nIso) {
      G4double* abundVector = anElement->GetRelativeAbundanceVector();
      if(G4int(xseciso.size()) < nIso) { xseciso.resize(nIso); }

      for (G4int j = 0; j<nIso; ++j) {
	G4double xsec = 0.0;
	if(abundVector[j] > 0.0) {
	  iso = (*isoVector)[j];
	  xsec = abundVector[j]*
	    GetIsoCrossSection(part, Z, iso->GetN(), iso, anElement, mat, i);
	}
	cross += xsec;
	xseciso[j] = cross;
      }
      cross *= G4UniformRand();
      for (G4int j = 0; j<nIso; ++j) {
	if(cross <= xseciso[j]) {
	  iso = (*isoVector)[j];
	  break;
	}
      }
    }
  }
  target.SetIsotope(iso);
  return anElement;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void
G4CrossSectionDataStore::BuildPhysicsTable(const G4ParticleDefinition& aParticleType)
{
  if (nDataSetList == 0) 
    {
      throw G4HadronicException(__FILE__, __LINE__, 
				"G4CrossSectionDataStore: no data sets registered");
      return;
    }
  for (G4int i=0; i<nDataSetList; ++i) {
    dataSetList[i]->BuildPhysicsTable(aParticleType);
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void 
G4CrossSectionDataStore::DumpPhysicsTable(const G4ParticleDefinition& aParticleType)
{
  // Print out all cross section data sets used and the energies at
  // which they apply

  if (nDataSetList == 0) {
    G4cout << "WARNING - G4CrossSectionDataStore::DumpPhysicsTable: "
	   << " no data sets registered" << G4endl;
    return;
  }

  for (G4int i = nDataSetList-1; i >= 0; --i) {
    G4double e1 = dataSetList[i]->GetMinKinEnergy();
    G4double e2 = dataSetList[i]->GetMaxKinEnergy();
     G4cout 
      << "     Cr_sctns: " << std::setw(25) << dataSetList[i]->GetName() << ": "
      <<  G4BestUnit(e1, "Energy")
      << " ---> "
      <<  G4BestUnit(e2, "Energy") << "\n";
      if (dataSetList[i]->GetName() == "G4CrossSectionPairGG") {
        dataSetList[i]->DumpPhysicsTable(aParticleType);
      }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4CrossSectionDataStore::DumpHtml(const G4ParticleDefinition&,
                                       std::ofstream& outFile)
{
  // Write cross section data set info to html physics list
  // documentation page

  G4double ehi = 0;
  G4double elo = 0;
  for (G4int i = nDataSetList-1; i > 0; i--) {
    elo = dataSetList[i]->GetMinKinEnergy()/GeV;
    ehi = dataSetList[i]->GetMaxKinEnergy()/GeV;
    outFile << "      <li><b><a href=\"" << dataSetList[i]->GetName() << ".html\"> "
            << dataSetList[i]->GetName() << "</a> from "
            << elo << " GeV to " << ehi << " GeV </b></li>\n";
  }

  G4double defaultHi = dataSetList[0]->GetMaxKinEnergy()/GeV;
  if (ehi < defaultHi) {
    outFile << "      <li><b><a href=\"" << dataSetList[0]->GetName() << ".html\"> "
            << dataSetList[0]->GetName() << "</a> from "
            << ehi << " GeV to " << defaultHi << " GeV </b></li>\n";
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
