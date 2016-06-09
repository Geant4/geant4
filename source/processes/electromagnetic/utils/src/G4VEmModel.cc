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
// $Id: G4VEmModel.cc,v 1.20 2008/11/13 23:13:18 schaelic Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4VEmModel
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 25.07.2005
//
// Modifications:
// 25.10.2005 Set default highLimit=100.TeV (V.Ivanchenko)
// 06.02.2006 add method ComputeMeanFreePath() (mma)
//
//
// Class Description:
//
// Abstract interface to energy loss models

// -------------------------------------------------------------------
//

#include "G4VEmModel.hh"
#include "G4LossTableManager.hh"
#include "G4ProductionCutsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VEmModel::G4VEmModel(const G4String& nam):
  fluc(0), name(nam), lowLimit(0.1*keV), highLimit(100.0*TeV), 
  polarAngleLimit(0.0),secondaryThreshold(DBL_MAX),theLPMflag(false),
  pParticleChange(0),nuclearStopping(false),nsec(5) 
{
  xsec.resize(nsec);
  nSelectors = 0;
  G4LossTableManager::Instance()->Register(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VEmModel::~G4VEmModel()
{
  G4LossTableManager::Instance()->DeRegister(this);
  G4int n = elmSelectors.size();
  if(n > 0) {
    for(G4int i=0; i<n; i++) { 
      delete elmSelectors[i]; 
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4VEmModel::CrossSectionPerVolume(const G4Material* material,
					   const G4ParticleDefinition* p,
					   G4double ekin,
					   G4double emin,
					   G4double emax)
{
  SetupForMaterial(p, material, ekin);
  G4double cross = 0.0;
  const G4ElementVector* theElementVector = material->GetElementVector();
  const G4double* theAtomNumDensityVector = material->GetVecNbOfAtomsPerVolume();
  G4int nelm = material->GetNumberOfElements(); 
  if(nelm > nsec) {
    xsec.resize(nelm);
    nsec = nelm;
  }
  for (G4int i=0; i<nelm; i++) {
    cross += theAtomNumDensityVector[i]*
      ComputeCrossSectionPerAtom(p,(*theElementVector)[i],ekin,emin,emax);
    xsec[i] = cross;
  }
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4VEmModel::ComputeMeanFreePath(const G4ParticleDefinition* p,
					 G4double ekin,
					 const G4Material* material,     
					 G4double emin,
					 G4double emax)
{
  G4double mfp = DBL_MAX;
  G4double cross = CrossSectionPerVolume(material,p,ekin,emin,emax);
  if (cross > DBL_MIN) mfp = 1./cross;
  return mfp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4VEmModel::InitialiseElementSelectors(const G4ParticleDefinition* p, 
					    const G4DataVector& cuts)
{
  G4int nbins = G4int(std::log10(highLimit/lowLimit) + 0.5);
  if(nbins < 3) nbins = 3;
  G4bool spline = G4LossTableManager::Instance()->SplineFlag();

  G4ProductionCutsTable* theCoupleTable=
    G4ProductionCutsTable::GetProductionCutsTable();
  G4int numOfCouples = theCoupleTable->GetTableSize();

  // prepare vector
  if(numOfCouples > nSelectors) {
    elmSelectors.resize(numOfCouples);
    nSelectors = numOfCouples;
  }

  // initialise vector
  for(G4int i=0; i<numOfCouples; i++) {
    const G4MaterialCutsCouple* couple =
      theCoupleTable->GetMaterialCutsCouple(i);
    const G4Material* material = couple->GetMaterial();
    G4int idx = couple->GetIndex();

    // selector already exist check if should be deleted
    G4bool create = true;
    if(elmSelectors[i]) {
      if(material == elmSelectors[i]->GetMaterial()) create = false;
      else delete elmSelectors[i];
    }
    if(create) {
      elmSelectors[i] = new G4EmElementSelector(this,material,nbins,
						lowLimit,highLimit,spline);
    }
    elmSelectors[i]->Initialise(p, cuts[idx]);
    //elmSelectors[i]->Dump(p);
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


