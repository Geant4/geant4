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
#include "Randomize.hh"
#include "G4Nucleus.hh"

#include "G4DynamicParticle.hh"
#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4NistManager.hh"
#include <algorithm>
#include <typeinfo>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4CrossSectionDataStore::G4CrossSectionDataStore()
  : nist(G4NistManager::Instance())
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4double 
G4CrossSectionDataStore::ComputeCrossSection(const G4DynamicParticle* dp,
					     const G4Material* mat)
{
  currentMaterial = mat;
  matParticle = dp->GetDefinition();
  matKinEnergy = dp->GetKineticEnergy();
  matCrossSection = 0.0;

  std::size_t nElements = mat->GetNumberOfElements();
  const G4double* nAtomsPerVolume = mat->GetVecNbOfAtomsPerVolume();

  if(xsecelm.size() < nElements) { xsecelm.resize(nElements); }

  for(G4int i=0; i<(G4int)nElements; ++i) {
    G4double xs = 
      nAtomsPerVolume[i]*GetCrossSection(dp, mat->GetElement(i), mat);
    matCrossSection += std::max(xs, 0.0); 
    xsecelm[i] = matCrossSection;
  }
  return matCrossSection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4double G4CrossSectionDataStore::GetCrossSection(const G4DynamicParticle* dp,
                                                  const G4Element* elm,
                                                  const G4Material* mat)
{
  // first check the most last cross section
  G4int i = nDataSetList-1;
  G4int Z = elm->GetZasInt();

  if(elm->GetNaturalAbundanceFlag() &&
     dataSetList[i]->IsElementApplicable(dp, Z, mat))
  {
    // element wise cross section
    return dataSetList[i]->GetElementCrossSection(dp, Z, mat);
  }

  // isotope wise cross section
  G4int nIso = (G4int)elm->GetNumberOfIsotopes();

  // user-defined isotope abundances
  const G4double* abundVector = elm->GetRelativeAbundanceVector();

  G4double sigma = 0.0;

  // isotope and element wise cross sections
  for(G4int j = 0; j < nIso; ++j) 
  {
    const G4Isotope* iso = elm->GetIsotope(j);
    sigma += abundVector[j] *
      GetIsoCrossSection(dp, Z, iso->GetN(), iso, elm, mat, i);
  }
  return sigma;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4double
G4CrossSectionDataStore::GetIsoCrossSection(const G4DynamicParticle* dp,
					    G4int Z, G4int A, 
					    const G4Isotope* iso,
					    const G4Element* elm,
					    const G4Material* mat, 
					    G4int idx)
{
  // this methods is called after the check that dataSetList[idx] 
  // depend on isotopes, so first isotopes are checked
  if(dataSetList[idx]->IsIsoApplicable(dp, Z, A, elm, mat) ) {
    return dataSetList[idx]->GetIsoCrossSection(dp, Z, A, iso, elm, mat);
  }

  // no isotope wise cross section - check other datasets
  for (G4int j = nDataSetList-1; j >= 0; --j) {
    if(dataSetList[j]->IsElementApplicable(dp, Z, mat)) {
      return dataSetList[j]->GetElementCrossSection(dp, Z, mat);
    } else if (dataSetList[j]->IsIsoApplicable(dp, Z, A, elm, mat)) {
      return dataSetList[j]->GetIsoCrossSection(dp, Z, A, iso, elm, mat);
    }
  }
  G4ExceptionDescription ed;
  ed << "No isotope cross section found for " 
     << dp->GetDefinition()->GetParticleName()
     << " off target Element " << elm->GetName() 
     << " Z= " << Z << " A= " << A;
  if(nullptr != mat) ed << " from " << mat->GetName();
  ed << " E(MeV)=" << dp->GetKineticEnergy()/MeV << G4endl; 
  G4Exception("G4CrossSectionDataStore::GetIsoCrossSection", "had001", 
              FatalException, ed);
  return 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4double
G4CrossSectionDataStore::GetCrossSection(const G4DynamicParticle* dp,
                                         G4int Z, G4int A,
					 const G4Isotope* iso,
                                         const G4Element* elm,
					 const G4Material* mat)
{
  for (G4int i = nDataSetList-1; i >= 0; --i) {
    if (dataSetList[i]->IsIsoApplicable(dp, Z, A, elm, mat) ) {
      return dataSetList[i]->GetIsoCrossSection(dp, Z, A, iso, elm, mat);
    } else if(dataSetList[i]->IsElementApplicable(dp, Z, mat)) {
      return dataSetList[i]->GetElementCrossSection(dp, Z, mat);
    }
  }
  G4ExceptionDescription ed;
  ed << "No isotope cross section found for " 
     << dp->GetDefinition()->GetParticleName()
     << " off target Element " << elm->GetName() 
     << " Z= " << Z << " A= " << A;
  if(nullptr != mat) ed << " from " << mat->GetName();
  ed << " E(MeV)=" << dp->GetKineticEnergy()/MeV << G4endl; 
  G4Exception("G4CrossSectionDataStore::GetCrossSection", "had001", 
              FatalException, ed);
  return 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

const G4Element*
G4CrossSectionDataStore::SampleZandA(const G4DynamicParticle* dp,
                                     const G4Material* mat,
				     G4Nucleus& target)
{
  if(nullptr != forcedElement) { return forcedElement; }
  std::size_t nElements = mat->GetNumberOfElements();
  const G4Element* anElement = mat->GetElement(0);

  // select element from a compound 
  if(1 < nElements) {
    G4double cross = matCrossSection*G4UniformRand();
    for(G4int i=0; i<(G4int)nElements; ++i) {
      if(cross <= xsecelm[i]) {
        anElement = mat->GetElement(i);
        break;
      }
    }
  }

  G4int Z = anElement->GetZasInt();
  const G4Isotope* iso = nullptr;

  G4int i = nDataSetList-1;
  if (dataSetList[i]->IsElementApplicable(dp, Z, mat)) {

    //----------------------------------------------------------------
    // element-wise cross section
    // isotope cross section is not computed
    //----------------------------------------------------------------
    std::size_t nIso = anElement->GetNumberOfIsotopes();
    iso = anElement->GetIsotope(0);

    // more than 1 isotope
    if(1 < nIso) { 
      iso = dataSetList[i]->SelectIsotope(anElement, 
                                          dp->GetKineticEnergy(),
					  dp->GetLogKineticEnergy());
    }
  } else {

    //----------------------------------------------------------------
    // isotope-wise cross section
    // isotope cross section is computed
    //----------------------------------------------------------------
    std::size_t nIso = anElement->GetNumberOfIsotopes();
    iso = anElement->GetIsotope(0);

    // more than 1 isotope
    if(1 < nIso) {
      const G4double* abundVector = anElement->GetRelativeAbundanceVector();
      if(xseciso.size() < nIso) { xseciso.resize(nIso); }

      G4double cross = 0.0;
      G4int j;
      for (j = 0; j<(G4int)nIso; ++j) {
	G4double xsec = 0.0;
	if(abundVector[j] > 0.0) {
	  iso = anElement->GetIsotope(j);
	  xsec = abundVector[j]*
	    GetIsoCrossSection(dp, Z, iso->GetN(), iso, anElement, mat, i);
	}
	cross += xsec;
	xseciso[j] = cross;
      }
      cross *= G4UniformRand();
      for (j = 0; j<(G4int)nIso; ++j) {
	if(cross <= xseciso[j]) {
	  iso = anElement->GetIsotope(j);
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
G4CrossSectionDataStore::BuildPhysicsTable(const G4ParticleDefinition& part)
{
  if (nDataSetList == 0) {
    G4ExceptionDescription ed;
    ed << "No cross section is registered for " 
       << part.GetParticleName() << G4endl;
    G4Exception("G4CrossSectionDataStore::BuildPhysicsTable", "had001", 
                FatalException, ed);
    return;
  }
  matParticle = &part;
  for (G4int i=0; i<nDataSetList; ++i) {
    dataSetList[i]->BuildPhysicsTable(part);
  }
  const G4MaterialTable* theMatTable = G4Material::GetMaterialTable();
  std::size_t nelm = 0;
  std::size_t niso = 0;
  for(auto mat : *theMatTable) {
    std::size_t nElements = mat->GetNumberOfElements();
    nelm = std::max(nelm, nElements);
    for(G4int j=0; j<(G4int)nElements; ++j) {
      niso = std::max(niso, mat->GetElement(j)->GetNumberOfIsotopes());
    }
  }
  // define vectors for a run
  xsecelm.resize(nelm, 0.0);
  xseciso.resize(niso, 0.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void 
G4CrossSectionDataStore::DumpPhysicsTable(const G4ParticleDefinition& part)
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
      << G4BestUnit(e1, "Energy") << " ---> "
      <<  G4BestUnit(e2, "Energy") << "\n";
    if (dataSetList[i]->GetName() == "G4CrossSectionPairGG") {
      dataSetList[i]->DumpPhysicsTable(part);
    }
    G4cout << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4CrossSectionDataStore::DumpHtml(const G4ParticleDefinition& /* pD */,
                                       std::ofstream& outFile) const
{
  // Write cross section data set info to html physics list
  // documentation page

  G4double ehi = 0;
  G4double elo = 0;
  G4String physListName(std::getenv("G4PhysListName"));
  for (G4int i = nDataSetList-1; i > 0; i--) {
    elo = dataSetList[i]->GetMinKinEnergy()/GeV;
    ehi = dataSetList[i]->GetMaxKinEnergy()/GeV;
    outFile << "      <li><b><a href=\"" << physListName << "_"
	    << dataSetList[i]->GetName() << ".html\"> "
            << dataSetList[i]->GetName() << "</a> from "
            << elo << " GeV to " << ehi << " GeV </b></li>\n";
    PrintCrossSectionHtml(dataSetList[i]);			
  }

  G4double defaultHi = dataSetList[0]->GetMaxKinEnergy()/GeV;
  if (ehi < defaultHi) {
    outFile << "      <li><b><a href=\"" << dataSetList[0]->GetName() 
	    << ".html\"> "
            << dataSetList[0]->GetName() << "</a> from "
            << ehi << " GeV to " << defaultHi << " GeV </b></li>\n";
    PrintCrossSectionHtml(dataSetList[0]);			
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4CrossSectionDataStore::PrintCrossSectionHtml(const G4VCrossSectionDataSet *cs) const
{
  G4String dirName(std::getenv("G4PhysListDocDir"));
  G4String physListName(std::getenv("G4PhysListName"));

  G4String pathName = dirName + "/" + physListName + "_" + HtmlFileName(cs->GetName());
  std::ofstream outCS;
  outCS.open(pathName);
  outCS << "<html>\n";
  outCS << "<head>\n";
  outCS << "<title>Description of " << cs->GetName() 
	<< "</title>\n";
  outCS << "</head>\n";
  outCS << "<body>\n";
  
  cs->CrossSectionDescription(outCS);
  
  outCS << "</body>\n";
  outCS << "</html>\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4String G4CrossSectionDataStore::HtmlFileName(const G4String & in) const
{
   G4String str(in);
   // replace blanks by _  C++11 version:
   std::transform(str.begin(), str.end(), str.begin(), [](char ch) {
       return ch == ' ' ? '_' : ch;
   });
   str=str + ".html";		
   return str;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4CrossSectionDataStore::AddDataSet(G4VCrossSectionDataSet* p)
{
  if(p->ForAllAtomsAndEnergies()) { 
    dataSetList.clear();
    nDataSetList = 0;
  }
  dataSetList.push_back(p);
  ++nDataSetList;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4CrossSectionDataStore::AddDataSet(G4VCrossSectionDataSet* p, std::size_t i)
{
  if(p->ForAllAtomsAndEnergies()) {
    dataSetList.clear();
    dataSetList.push_back(p);
    nDataSetList = 1;
  } else if ( i >= dataSetList.size() ) {
    dataSetList.push_back(p);
    ++nDataSetList;
  } else {
    std::vector< G4VCrossSectionDataSet* >::iterator it = dataSetList.end() - i;
    dataSetList.insert(it , p);
    ++nDataSetList;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
