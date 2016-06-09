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
// $Id: G4EmElementSelector.cc,v 1.4 2008/08/21 18:53:32 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-02 $
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmElementSelector::G4EmElementSelector(G4VEmModel* mod, 
					 const G4Material* mat, 
					 G4int bins, 
					 G4double emin, 
					 G4double emax,
					 G4bool spline):
  model(mod), material(mat), nbins(bins), cutEnergy(-1.0), 
  lowEnergy(emin), highEnergy(emax)
{
  G4int n = material->GetNumberOfElements();
  nElmMinusOne = n - 1;
  theElementVector = material->GetElementVector();
  if(nElmMinusOne > 0) {
    for(G4int i=0; i<nElmMinusOne; i++) {
      G4PhysicsLogVector* v = new G4PhysicsLogVector(lowEnergy,highEnergy,nbins);
      v->SetSpline(spline);
      xSections.push_back(v);
    }
  }
  //G4cout << "G4EmElementSelector for " << mat->GetName() << " n= " << n << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmElementSelector::~G4EmElementSelector()
{
  if(nElmMinusOne > 0) {
    for(G4int i=0; i<nElmMinusOne; i++) {
      delete xSections[i];
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmElementSelector::Initialise(const G4ParticleDefinition* part, 
				     G4double cut)
{
  //G4cout << "G4EmElementSelector initialise for " << material->GetName() << G4endl;
  if(0 == nElmMinusOne || cut == cutEnergy) return;

  cutEnergy = cut;
  //G4cout << "cut(keV)= " << cut/keV << G4endl;
  G4double cross;

  const G4double* theAtomNumDensityVector = material->GetVecNbOfAtomsPerVolume();

  G4int i;
  G4int n = nElmMinusOne + 1;
  G4double* xsec = new G4double[n];  

  // loop over bins
  for(G4int j=0; j<nbins; j++) {
    G4double e = (xSections[0])->GetLowEdgeEnergy(j);
    model->SetupForMaterial(part, material, e);
    cross = 0.0;
    //G4cout << "j= " << j << " e(MeV)= " << e/MeV << G4endl;
    for (i=0; i<n; i++) {
      cross += theAtomNumDensityVector[i]*      
	model->ComputeCrossSectionPerAtom(part, (*theElementVector)[i], e, 
					  cutEnergy, e);
      xsec[i] = cross;
    }
    if(DBL_MIN >= cross) cross = 1.0;
    // normalise cross section sum 
    for (i=0; i<nElmMinusOne; i++) {
      xSections[i]->PutValue(j, xsec[i]/cross);
      //G4cout << "i= " << i << " xs= " << xsec[i]/cross << G4endl;
    }
  }
  delete [] xsec;
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
  G4cout << "Last Element in element vector" 
	 << (*theElementVector)[nElmMinusOne]->GetName() 
	 << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


