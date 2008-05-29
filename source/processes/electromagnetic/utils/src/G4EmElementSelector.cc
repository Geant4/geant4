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
// $Id: G4EmElementSelector.cc,v 1.1 2008-05-29 13:38:05 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
#include "G4LossTableManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmElementSelector::G4EmElementSelector(G4VEmModel* mod, 
  const G4Material* mat, G4int bins, G4double emin, G4double emax):
  model(mod), material(mat), nbins(bins), cutEnergy(-1.0), 
  lowEnergy(emin), highEnergy(emax)
{
  nElements = material->GetNumberOfElements();
  theElementVector = material->GetElementVector();
  lastElement = (*theElementVector)[0];
  if(nElements > 1) {
    xSections.resize(nElements);
    for(G4int i=0; i<nElements; i++) {
      G4PhysicsLogVector* v = new G4PhysicsLogVector(lowEnergy,highEnergy,nbins);
      v->SetSpline((G4LossTableManager::Instance())->SplineFlag());
      xSections.push_back(v);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmElementSelector::~G4EmElementSelector()
{
  if(nElements > 1) {
    for(G4int i=0; i<nElements; i++) {
      delete xSections[i];
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmElementSelector::Initialise(const G4ParticleDefinition* part, 
				     G4double cut)
{
  if(1 == nElements || cut == cutEnergy) return;

  cutEnergy = cut;
  model->SetupForMaterial(part, material);
  G4double cross;
  const G4Element* elm;

  const G4double* theAtomNumDensityVector = material->GetVecNbOfAtomsPerVolume();

  for(G4int j=0; j<nbins; j++) {
    G4double e = (xSections[0])->GetLowEdgeEnergy(j);
    cross = 0.0;
    for (G4int i=0; i<nElements; i++) {
      elm = (*theElementVector)[i];
      cross += theAtomNumDensityVector[i]*      
	model->ComputeCrossSectionPerAtom(part, e, elm->GetZ(), elm->GetN(),
					  cutEnergy, e);
      xSections[i]->PutValue(j, cross);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmElementSelector::Dump(const G4ParticleDefinition* part)
{
  G4cout << "======== G4EmElementSelector for the " << model->GetName();
  if(part) G4cout << " and " << part->GetParticleName();
  G4cout << " for " << material->GetName() << " ========" << G4endl;
  if(nElements > 1) {
    for(G4int i=0; i<nElements; i++) {
      G4cout << "      " << (*theElementVector)[i]->GetName() << " : " << G4endl;
      G4cout << &(xSections[i]) << G4endl;
    }
  }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


