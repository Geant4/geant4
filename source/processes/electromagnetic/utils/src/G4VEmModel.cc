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
// $Id: G4VEmModel.cc,v 1.37 2010-10-14 16:27:35 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
// 16.02.2009 Move implementations of virtual methods to source (VI)
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
#include "G4ParticleChangeForLoss.hh"
#include "G4ParticleChangeForGamma.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VEmModel::G4VEmModel(const G4String& nam):
  flucModel(0),anglModel(0), name(nam), lowLimit(0.1*keV), highLimit(100.0*TeV), 
  eMinActive(0.0),eMaxActive(DBL_MAX),
  polarAngleLimit(0.0),secondaryThreshold(DBL_MAX),theLPMflag(false),
  pParticleChange(0),
  currentCouple(0),currentElement(0),
  nsec(5),flagDeexcitation(false) 
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
    for(G4int i=0; i<n; ++i) { 
      delete elmSelectors[i]; 
    }
  }
  delete anglModel;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ParticleChangeForLoss* G4VEmModel::GetParticleChangeForLoss()
{
  G4ParticleChangeForLoss* p = 0;
  if (pParticleChange) {
    p = static_cast<G4ParticleChangeForLoss*>(pParticleChange);
  } else {
    p = new G4ParticleChangeForLoss();
    pParticleChange = p;
  }
  return p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ParticleChangeForGamma* G4VEmModel::GetParticleChangeForGamma()
{
  G4ParticleChangeForGamma* p = 0;
  if (pParticleChange) {
    p = static_cast<G4ParticleChangeForGamma*>(pParticleChange);
  } else {
    p = new G4ParticleChangeForGamma();
    pParticleChange = p;
  }
  return p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4VEmModel::InitialiseElementSelectors(const G4ParticleDefinition* p, 
					    const G4DataVector& cuts)
{
  // initialise before run
  G4LossTableManager* man = G4LossTableManager::Instance();
  // G4bool spline = man->SplineFlag();
  G4bool spline = false;

  // two times less bins because probability functon is normalized 
  // so correspondingly is more smooth
  G4int nbins = G4int(man->GetNumberOfBinsPerDecade()
		      * std::log10(highLimit/lowLimit) / 6.0);
  if(nbins < 5) { nbins = 5; }

  G4ProductionCutsTable* theCoupleTable=
    G4ProductionCutsTable::GetProductionCutsTable();
  G4int numOfCouples = theCoupleTable->GetTableSize();

  // prepare vector
  if(numOfCouples > nSelectors) { 
    elmSelectors.reserve(numOfCouples); 
    for(G4int i=nSelectors; i<numOfCouples; ++i) { elmSelectors.push_back(0); }
    nSelectors = numOfCouples;
  }

  // initialise vector
  for(G4int i=0; i<numOfCouples; ++i) {
    currentCouple = theCoupleTable->GetMaterialCutsCouple(i);
    const G4Material* material = currentCouple->GetMaterial();
    G4int idx = currentCouple->GetIndex();

    // selector already exist check if should be deleted
    G4bool create = true;
    if(elmSelectors[i]) {
      if(material == elmSelectors[i]->GetMaterial()) { create = false; }
      else { delete elmSelectors[i]; }
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

G4double G4VEmModel::ComputeDEDXPerVolume(const G4Material*,
					  const G4ParticleDefinition*,
					  G4double,G4double)
{
  return 0.0;
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

const G4Element* G4VEmModel::SelectRandomAtom(const G4Material* material,
					      const G4ParticleDefinition* pd,
					      G4double kinEnergy,
					      G4double tcut,
					      G4double tmax)
{
  const G4ElementVector* theElementVector = material->GetElementVector();
  G4int n = material->GetNumberOfElements() - 1;
  currentElement = (*theElementVector)[n];
  if (n > 0) {
    G4double x = G4UniformRand()*
                 G4VEmModel::CrossSectionPerVolume(material,pd,kinEnergy,tcut,tmax);
    for(G4int i=0; i<n; ++i) {
      if (x <= xsec[i]) {
	currentElement = (*theElementVector)[i];
	break;
      }
    }
  }
  return currentElement;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4VEmModel::ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
						G4double, G4double, G4double,
						G4double, G4double)
{
  return 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4VEmModel::DefineForRegion(const G4Region*) 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4VEmModel::MinEnergyCut(const G4ParticleDefinition*,
				  const G4MaterialCutsCouple*)
{
  return 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4VEmModel::ChargeSquareRatio(const G4Track& track)
{
  return GetChargeSquareRatio(track.GetParticleDefinition(), 
			      track.GetMaterial(), track.GetKineticEnergy());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4VEmModel::GetChargeSquareRatio(const G4ParticleDefinition* p,
					  const G4Material*, G4double)
{
  G4double q = p->GetPDGCharge()/CLHEP::eplus;
  return q*q;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4VEmModel::GetParticleCharge(const G4ParticleDefinition* p,
				       const G4Material*, G4double)
{
  return p->GetPDGCharge();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4VEmModel::CorrectionsAlongStep(const G4MaterialCutsCouple*,
				      const G4DynamicParticle*,
				      G4double&,G4double&,G4double)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4VEmModel::SampleDeexcitationAlongStep(const G4Material*,
					     const G4Track&,
					     G4double& )
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4VEmModel::MaxSecondaryEnergy(const G4ParticleDefinition*,
					G4double kineticEnergy)
{
  return kineticEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4VEmModel::SetupForMaterial(const G4ParticleDefinition*,
				  const G4Material*, G4double)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void 
G4VEmModel::SetParticleChange(G4VParticleChange* p, G4VEmFluctuationModel* f)
{
  if(p && pParticleChange != p) { pParticleChange = p; }
  flucModel = f;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
