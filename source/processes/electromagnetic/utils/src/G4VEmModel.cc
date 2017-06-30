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
// $Id: G4VEmModel.cc 104457 2017-05-31 15:52:37Z gcosmo $
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
#include "G4ElementData.hh"
#include "G4LossTableManager.hh"
#include "G4ProductionCutsTable.hh"
#include "G4ParticleChangeForLoss.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4SystemOfUnits.hh"
#include "G4Log.hh"
#include "Randomize.hh"
#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const G4double G4VEmModel::inveplus = 1.0/CLHEP::eplus;

G4VEmModel::G4VEmModel(const G4String& nam):
  flucModel(nullptr),anglModel(nullptr), name(nam), lowLimit(0.1*CLHEP::keV), 
  highLimit(100.0*CLHEP::TeV),eMinActive(0.0),eMaxActive(DBL_MAX),
  polarAngleLimit(CLHEP::pi),secondaryThreshold(DBL_MAX),
  theLPMflag(false),flagDeexcitation(false),flagForceBuildTable(false),
  isMaster(true),fElementData(nullptr),pParticleChange(nullptr),
  xSectionTable(nullptr),theDensityFactor(nullptr),theDensityIdx(nullptr),
  lossFlucFlag(true),fCurrentCouple(nullptr),
  fCurrentElement(nullptr),fCurrentIsotope(nullptr),
  fTripletModel(nullptr),nsec(5) 
{
  xsec.resize(nsec);
  nSelectors = 0;
  elmSelectors = nullptr;
  localElmSelectors = true;
  localTable = true;
  useAngularGenerator = false;
  isLocked = false;
  idxTable = 0;

  fEmManager = G4LossTableManager::Instance();
  fEmManager->Register(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VEmModel::~G4VEmModel()
{
  if(localElmSelectors) { 
    if(nSelectors > 0) {
      for(G4int i=0; i<nSelectors; ++i) { 
        delete (*elmSelectors)[i]; 
      }
    }
    delete elmSelectors; 
  }
  delete anglModel;
  
  if(localTable && xSectionTable) { 
    xSectionTable->clearAndDestroy();
    delete xSectionTable;
    xSectionTable = nullptr; 
  }
  if(isMaster && fElementData) {
    delete fElementData;
    fElementData = nullptr;
  }
  fEmManager->DeRegister(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ParticleChangeForLoss* G4VEmModel::GetParticleChangeForLoss()
{
  G4ParticleChangeForLoss* p = nullptr;
  if (pParticleChange) {
    p = static_cast<G4ParticleChangeForLoss*>(pParticleChange);
  } else {
    p = new G4ParticleChangeForLoss();
    pParticleChange = p;
  }
  if(fTripletModel) { fTripletModel->SetParticleChange(p); }
  return p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ParticleChangeForGamma* G4VEmModel::GetParticleChangeForGamma()
{
  G4ParticleChangeForGamma* p = nullptr;
  if (pParticleChange) {
    p = static_cast<G4ParticleChangeForGamma*>(pParticleChange);
  } else {
    p = new G4ParticleChangeForGamma();
    pParticleChange = p;
  }
  if(fTripletModel) { fTripletModel->SetParticleChange(p); }
  return p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4VEmModel::InitialiseElementSelectors(const G4ParticleDefinition* part, 
                                            const G4DataVector& cuts)
{
  // using spline for element selectors should be investigated in details
  // because small number of points may provide biased results
  // large number of points requires significant increase of memory
  G4bool spline = false;
  
  //G4cout << "IES: for " << GetName() << " Emin(MeV)= " << lowLimit/MeV 
  //         << " Emax(MeV)= " << highLimit/MeV << G4endl;
  
  // two times less bins because probability functon is normalized 
  // so correspondingly is more smooth
  if(highLimit <= lowLimit) { return; }

  G4int nbinsPerDec = G4EmParameters::Instance()->NumberOfBinsPerDecade();

  G4ProductionCutsTable* theCoupleTable=
    G4ProductionCutsTable::GetProductionCutsTable();
  G4int numOfCouples = theCoupleTable->GetTableSize();

  // prepare vector
  if(!elmSelectors) {
    elmSelectors = new std::vector<G4EmElementSelector*>;
  }
  if(numOfCouples > nSelectors) { 
    for(G4int i=nSelectors; i<numOfCouples; ++i) { 
      elmSelectors->push_back(nullptr); 
    }
    nSelectors = numOfCouples;
  }

  // initialise vector
  for(G4int i=0; i<numOfCouples; ++i) {

    // no need in element selectors for infionite cuts
    if(cuts[i] == DBL_MAX) { continue; }
   
    fCurrentCouple = theCoupleTable->GetMaterialCutsCouple(i);
    const G4Material* material = fCurrentCouple->GetMaterial();

    // selector already exist check if should be deleted
    G4bool create = true;
    if((*elmSelectors)[i]) {
      if(material == ((*elmSelectors)[i])->GetMaterial()) { create = false; }
      else { delete (*elmSelectors)[i]; }
    }
    if(create) {
      G4double emin = std::max(lowLimit, 
                               MinPrimaryEnergy(material, part, cuts[i]));
      G4double emax = std::max(highLimit, 10*emin);
      static const G4double invlog106 = 1.0/(6*G4Log(10.));
      G4int nbins = (G4int)(nbinsPerDec*G4Log(emax/emin)*invlog106);
      nbins = std::max(nbins, 3);

      (*elmSelectors)[i] = new G4EmElementSelector(this,material,nbins,
                                                   emin,emax,spline);
    }
    ((*elmSelectors)[i])->Initialise(part, cuts[i]);
    /*      
      G4cout << "G4VEmModel::InitialiseElmSelectors i= " << i
             << " idx= " << fCurrentCouple->GetIndex() 
             << "  "  << part->GetParticleName() 
             << " for " << GetName() << "  cut= " << cuts[i] 
             << "  " << (*elmSelectors)[i] << G4endl;      
      ((*elmSelectors)[i])->Dump(part);
    */
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4VEmModel::InitialiseLocal(const G4ParticleDefinition*, 
                                 G4VEmModel*)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4VEmModel::InitialiseForMaterial(const G4ParticleDefinition* part,
                                       const G4Material* material)
{
  if(material) {
    const G4ElementVector* theElementVector = material->GetElementVector();
    G4int n = material->GetNumberOfElements();
    for(G4int i=0; i<n; ++i) {
      G4int Z = ((*theElementVector)[i])->GetZasInt();
      InitialiseForElement(part, Z);
    }
  } else {
    //G4cout << "G4VEmModel::InitialiseForMaterial for " << GetName();
    //if(part) { G4cout << " and  " << part->GetParticleName(); }
    //G4cout << " with no material" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4VEmModel::InitialiseForElement(const G4ParticleDefinition*, G4int)
{}

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
  const G4double* theAtomNumDensityVector = 
    material->GetVecNbOfAtomsPerVolume();
  G4int nelm = material->GetNumberOfElements(); 
  if(nelm > nsec) {
    xsec.resize(nelm);
    nsec = nelm;
  }
  for (G4int i=0; i<nelm; ++i) {
    cross += theAtomNumDensityVector[i]*
      ComputeCrossSectionPerAtom(p,(*theElementVector)[i],ekin,emin,emax);
    xsec[i] = cross;
  }
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4VEmModel::GetPartialCrossSection(const G4Material*, G4int,
                                            const G4ParticleDefinition*,
                                            G4double)
{ 
  return 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4VEmModel::StartTracking(G4Track*)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const G4Element* G4VEmModel::SelectRandomAtom(const G4Material* material,
                                              const G4ParticleDefinition* pd,
                                              G4double kinEnergy,
                                              G4double tcut,
                                              G4double tmax)
{
  const G4ElementVector* theElementVector = material->GetElementVector();
  G4int n = material->GetNumberOfElements() - 1;
  fCurrentElement = (*theElementVector)[n];
  if (n > 0) {
    G4double x = G4UniformRand()*
      G4VEmModel::CrossSectionPerVolume(material,pd,kinEnergy,tcut,tmax);
    for(G4int i=0; i<n; ++i) {
      if (x <= xsec[i]) {
        fCurrentElement = (*theElementVector)[i];
        break;
      }
    }
  }
  return fCurrentElement;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4VEmModel::ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
                                                G4double, G4double, G4double,
                                                G4double, G4double)
{
  return 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double 
G4VEmModel::ComputeCrossSectionPerShell(const G4ParticleDefinition*,
                                        G4int, G4int, 
                                        G4double, G4double, G4double)
{
  return 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4VEmModel::DefineForRegion(const G4Region*) 
{}

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
  G4double q = p->GetPDGCharge()*inveplus;
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

G4double G4VEmModel::Value(const G4MaterialCutsCouple* couple,
                           const G4ParticleDefinition* p, G4double e)
{
  SetCurrentCouple(couple);
  return e*e*CrossSectionPerVolume(couple->GetMaterial(),p,e,0.0,DBL_MAX);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4VEmModel::MinPrimaryEnergy(const G4Material*,
                                      const G4ParticleDefinition*,
                                      G4double)
{
  return 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4VEmModel::MinEnergyCut(const G4ParticleDefinition*, 
                                  const G4MaterialCutsCouple*)
{
  return 0.0;
}

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
  if(flucModel != f) { flucModel = f; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4VEmModel::SetCrossSectionTable(G4PhysicsTable* p, G4bool isLocal)
{
  if(p != xSectionTable) {
    if(xSectionTable && localTable) { 
      xSectionTable->clearAndDestroy(); 
      delete xSectionTable;
    }
    xSectionTable = p;
  }
  localTable = isLocal;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void  G4VEmModel::ModelDescription(std::ostream& outFile) const
{
  outFile << "The description for this model has not been written yet.\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
