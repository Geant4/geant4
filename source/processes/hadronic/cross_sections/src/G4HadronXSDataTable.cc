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
//---------------------------------------------------------------------------
//
// GEANT4 Class file
//
// Description: Data structure for cross sections per materials
//
// Author:      V.Ivanchenko 31.05.2018
//
// Modifications:
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4HadronXSDataTable.hh"
#include "G4PhysicsLogVector.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4DynamicParticle.hh"
#include "G4CrossSectionDataStore.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4HadElementSelector::G4HadElementSelector(G4DynamicParticle* dp, 
					   G4CrossSectionDataStore* xs, 
					   const G4Material* mat, 
					   G4int bins, G4double emin, 
					   G4double emax, G4bool)
{
  std::size_t n = mat->GetNumberOfElements();
  nElmMinusOne = G4int(n - 1);
  theElementVector = mat->GetElementVector();
  if(nElmMinusOne > 0) {
    G4PhysicsVector* first = nullptr;
    xSections.resize(n, first);
    first = new G4PhysicsLogVector(emin,emax,bins,false);
    xSections[0] = first;
    for(std::size_t i=1; i<n; ++i) {
      xSections[i] = new G4PhysicsVector(*first);
    }
    std::vector<G4double> temp;
    temp.resize(n, 0.0);
    for(G4int j=0; j<=bins; ++j) {
      G4double cross = 0.0;
      G4double e = first->Energy(j);
      dp->SetKineticEnergy(e);
      for(std::size_t i=0; i<n; ++i) {
	cross += xs->GetCrossSection(dp, (*theElementVector)[i], mat);
        temp[i] = cross;
      }
      G4double fact = (cross > 0.0) ? 1.0/cross : 0.0;
      for(std::size_t i=0; i<n; ++i) {
	G4double y = (i<n-1) ? temp[i]*fact : 1.0; 
        xSections[i]->PutValue(j, y);
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4HadElementSelector::~G4HadElementSelector()
{
  if(nElmMinusOne > 0) {
    for(G4int i=0; i<=nElmMinusOne; ++i) { delete xSections[i]; }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4HadElementSelector::Dump()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4HadronXSDataTable::G4HadronXSDataTable() : nMaterials(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4HadronXSDataTable::Initialise(G4DynamicParticle* dp, 
				     G4CrossSectionDataStore* xs, 
				     G4int bins, G4double emin, G4double emax, 
				     G4bool spline)
{
  std::size_t nn = G4Material::GetNumberOfMaterials();
  if(nn > nMaterials) {
    if(0 == nMaterials) {
      xsData.reserve(nn);
      elmSelectors.reserve(nn);
    }
    G4PhysicsLogVector* first = nullptr;
    G4int sbins = std::max(10, bins/5);
    const G4MaterialTable* mtable = G4Material::GetMaterialTable();
    for(std::size_t i=nMaterials; i<nn; ++i) {
      const G4Material* mat = (*mtable)[i];
      G4PhysicsVector* v = nullptr;
      G4HadElementSelector* es = nullptr;
      // create real vector only for complex materials
      if(mat->GetNumberOfElements() > 1) {
	if(nullptr == first) {
          first = new G4PhysicsLogVector(emin, emax, bins, spline);
	  v = first;
	} else {
	  v = new G4PhysicsVector(*first);
	}
        for(G4int j=0; j<=bins; ++j) {
          G4double e = first->Energy(j);
          dp->SetKineticEnergy(e);
          G4double cros = xs->ComputeCrossSection(dp, mat);
          v->PutValue(j, cros);
	}
        if(spline) v->FillSecondDerivatives();
	elmSelectors[i] = new G4HadElementSelector(dp, xs, mat, sbins, emin, emax, spline);
      }
      xsData.push_back(v);
      elmSelectors.push_back(es);
    }
    nMaterials = nn;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4HadronXSDataTable::~G4HadronXSDataTable()
{
  for(std::size_t i=0; i<nMaterials; ++i) {
    delete xsData[i];
    delete elmSelectors[i];
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4HadronXSDataTable::Dump()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
