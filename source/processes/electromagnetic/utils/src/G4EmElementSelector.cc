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
// File name:     G4EmElementSelector
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 29.05.2008
//
// Modifications:
//
// Class Description:
//
// Generic helper class for the random selection of an element

// -------------------------------------------------------------------
//

#include "G4EmElementSelector.hh"
#include "G4VEmModel.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmElementSelector::G4EmElementSelector(G4VEmModel* mod, 
                                         const G4Material* mat, 
                                         G4int bins, 
                                         G4double emin, 
                                         G4double emax,
                                         G4bool):
  model(mod), material(mat), nbins(bins), cutEnergy(-1.0), 
  lowEnergy(emin), highEnergy(emax)
{
  G4int n = (G4int)material->GetNumberOfElements();
  nElmMinusOne = n - 1;
  theElementVector = material->GetElementVector();
  if(nElmMinusOne > 0) {
    xSections.reserve(n);
    auto v0 = new G4PhysicsLogVector(lowEnergy,highEnergy,nbins,false);
    xSections.push_back(v0);
    for(G4int i=1; i<n; ++i) {
      auto v = new G4PhysicsLogVector(*v0);
      xSections.push_back(v);
    }
  }
  /*  
  G4cout << "G4EmElementSelector for " << mat->GetName() << " n= " << n
         << " nbins= " << nbins << "  Emin= " << lowEnergy 
         << " Emax= " << highEnergy << G4endl;
  */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmElementSelector::~G4EmElementSelector()
{
  if(nElmMinusOne > 0) {
    for(G4int i=0; i<=nElmMinusOne; ++i) { delete xSections[i]; }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmElementSelector::Initialise(const G4ParticleDefinition* part, 
                                     G4double cut)
{
  //G4cout << "G4EmElementSelector initialise for " << material->GetName()
  //  << G4endl;
  if(0 == nElmMinusOne || cut == cutEnergy) { return; }

  cutEnergy = cut;
  //G4cout << "cut(keV)= " << cut/keV << G4endl;
  G4double cross;

  const G4double* theAtomNumDensityVector = 
    material->GetVecNbOfAtomsPerVolume();

  // loop over bins
  for(G4int j=0; j<=nbins; ++j) {
    G4double e = (xSections[0])->Energy(j);
    model->SetupForMaterial(part, material, e);
    cross = 0.0;
    //G4cout << "j= " << j << " e(MeV)= " << e/MeV << G4endl;
    for (G4int i=0; i<=nElmMinusOne; ++i) {
      cross += theAtomNumDensityVector[i]*      
        model->ComputeCrossSectionPerAtom(part, (*theElementVector)[i], e, 
                                          cutEnergy, e);
      xSections[i]->PutValue(j, cross);
    }
  }
  // xSections start from null, so use probabilities from the next bin
  if(0.0 == (*xSections[nElmMinusOne])[0]) {
    for (G4int i=0; i<=nElmMinusOne; ++i) {
      xSections[i]->PutValue(0, (*xSections[i])[1]);
    }
  }
  // xSections ends with null, so use probabilities from the previous bin
  if(0.0 == (*xSections[nElmMinusOne])[nbins]) {
    for (G4int i=0; i<=nElmMinusOne; ++i) {
      xSections[i]->PutValue(nbins, (*xSections[i])[nbins-1]);
    }
  }
  // perform normalization
  for(G4int j=0; j<=nbins; ++j) {
    cross = (*xSections[nElmMinusOne])[j];
    // only for positive X-section 
    if(cross > 0.0) {
      for (G4int i=0; i<nElmMinusOne; ++i) {
        G4double x = (*xSections[i])[j]/cross;
        xSections[i]->PutValue(j, x);
      }
    }
  }
  /*
  G4cout << "======== G4EmElementSelector for the " << model->GetName() 
         << G4endl;
  for (G4int i=0; i<=nElmMinusOne; ++i) {
    G4cout << "#### i=" << i << G4endl;
    G4cout << (*xSections[i]) << G4endl;
  }
  */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

const G4Element* 
G4EmElementSelector::SelectRandomAtom(const G4double e, 
                                      const G4double loge) const
{
  const G4Element* element = (*theElementVector)[nElmMinusOne];
  if (nElmMinusOne > 0) {
    // 1. Determine energy index (only once)
    // handle cases below/above the enrgy grid (by ekin, idx that gives a=0/1)
    // ekin = x[0]   if e<=x[0]   and idx will be   0 ^ a=0 => so y=y0
    // ekin = x[N-1] if e>=x[N-1] and idx will be N-2 ^ a=1 => so y=y_{N-1}
    G4double ekin = e;
    std::size_t idx = 0;
    if(e <= (xSections[0])->Energy(0)) {
      ekin = (xSections[0])->Energy(0);
    } else if(e < (xSections[0])->GetMaxEnergy()) {
      idx = (xSections[0])->ComputeLogVectorBin(loge);
    } else {
      ekin = (xSections[0])->GetMaxEnergy();
      idx = (xSections[0])->GetVectorLength() - 2;
    }
    // 2. Do the linear interp.(corner cases are already excluded)
    const G4double x1 = (xSections[0])->Energy(idx);
    const G4double  a = (ekin - x1)/((xSections[0])->Energy(idx+1) - x1);
    const G4double urnd = G4UniformRand();
    for (G4int i = 0; i < nElmMinusOne; ++i) {
      const G4double  y1 = (*xSections[i])[idx];
      if (urnd <= y1 + a*((*xSections[i])[idx+1] - y1)) {
        element = (*theElementVector)[i];
        break;
      }
    }
  }
  return element;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmElementSelector::Dump(const G4ParticleDefinition* part)
{
  G4cout << "======== G4EmElementSelector for the " << model->GetName();
  if(part) G4cout << " and " << part->GetParticleName();
  G4cout << " for " << material->GetName() << " ========" << G4endl;
  if(0 < nElmMinusOne) {
    for(G4int i=0; i<nElmMinusOne; i++) {
      G4cout << "      " << (*theElementVector)[i]->GetName() << " : " << G4endl;
      G4cout << *(xSections[i]) << G4endl;
    }
  }  
  G4cout << "Last Element in element vector " 
         << (*theElementVector)[nElmMinusOne]->GetName() 
         << G4endl;
  G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


