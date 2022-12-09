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
//---------------------------------------------------------------------------
//
// ClassName:    G4TablesForExtrapolator
//  
// Description:  This class provide calculation of energy loss, fluctuation, 
//               and msc angle
//
// Author:       09.12.04 V.Ivanchenko 
//
// Modification: 
// 08-04-05 Rename Propogator -> Extrapolator (V.Ivanchenko)
// 16-03-06 Add muon tables and fix bug in units (V.Ivanchenko)
// 21-03-06 Add verbosity defined in the constructor and Initialisation
//          start only when first public method is called (V.Ivanchenko)
// 03-05-06 Remove unused pointer G4Material* from number of methods (VI)
// 12-05-06 SEt linLossLimit=0.001 (VI)
//
//----------------------------------------------------------------------------
//
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4TablesForExtrapolator.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4LossTableManager.hh"
#include "G4PhysicsLogVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4EmParameters.hh"
#include "G4MollerBhabhaModel.hh"
#include "G4BetheBlochModel.hh"
#include "G4eBremsstrahlungRelModel.hh"
#include "G4MuPairProductionModel.hh"
#include "G4MuBremsstrahlungModel.hh"
#include "G4ProductionCuts.hh"
#include "G4LossTableBuilder.hh"
#include "G4WentzelVIModel.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4TablesForExtrapolator::G4TablesForExtrapolator(G4int verb, G4int bins, 
                                                 G4double e1, G4double e2)
  : emin(e1), emax(e2), verbose(verb), nbins(bins)
{
  electron = G4Electron::Electron();
  positron = G4Positron::Positron();
  proton   = G4Proton::Proton();
  muonPlus = G4MuonPlus::MuonPlus();
  muonMinus= G4MuonMinus::MuonMinus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4TablesForExtrapolator::~G4TablesForExtrapolator()
{
  if(nullptr != dedxElectron) {
    dedxElectron->clearAndDestroy();
    delete dedxElectron;
  }
  if(nullptr != dedxPositron) {
    dedxPositron->clearAndDestroy();
    delete dedxPositron;
  }
  if(nullptr != dedxProton) {
    dedxProton->clearAndDestroy();
    delete dedxProton;
  }
  if(nullptr != dedxMuon) {
    dedxMuon->clearAndDestroy();
    delete dedxMuon;
  }
  if(nullptr != rangeElectron) {
    rangeElectron->clearAndDestroy();
    delete rangeElectron;
  }
  if(nullptr != rangePositron) {
    rangePositron->clearAndDestroy();
    delete rangePositron;
  }
  if(nullptr != rangeProton) {
    rangeProton->clearAndDestroy();
    delete rangeProton;
  }
  if(nullptr != rangeMuon) {
    rangeMuon->clearAndDestroy();
    delete rangeMuon;
  }
  if(nullptr != invRangeElectron) {
    invRangeElectron->clearAndDestroy();
    delete invRangeElectron;
  }
  if(nullptr != invRangePositron) {
    invRangePositron->clearAndDestroy();
    delete invRangePositron;
  }
  if(nullptr != invRangeProton) {
    invRangeProton->clearAndDestroy();
    delete invRangeProton;
  }
  if(nullptr != invRangeMuon) {
    invRangeMuon->clearAndDestroy();
    delete invRangeMuon;
  }
  if(nullptr != mscElectron) {
    mscElectron->clearAndDestroy();
    delete mscElectron;
  }
  delete pcuts;
  delete builder;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4PhysicsTable* 
G4TablesForExtrapolator::GetPhysicsTable(ExtTableType type) const
{
  const G4PhysicsTable* table = nullptr;
  switch (type) 
    {
    case fDedxElectron:
      table = dedxElectron;
      break;
    case fDedxPositron:
      table = dedxPositron;
      break;
    case fDedxProton:
      table = dedxProton;
      break;
    case fDedxMuon:
      table = dedxMuon;
      break;
    case fRangeElectron:
      table = rangeElectron;
      break;
    case fRangePositron:
      table = rangePositron;
      break;
    case fRangeProton:
      table = rangeProton;
      break;
    case fRangeMuon:
      table = rangeMuon;
      break;
    case fInvRangeElectron:
      table = invRangeElectron;
      break;
    case fInvRangePositron:
      table = invRangePositron;
      break;
    case fInvRangeProton:
      table = invRangeProton;
      break;
    case fInvRangeMuon:
      table = invRangeMuon;
      break;
    case fMscElectron:
      table = mscElectron;
    }
  return table;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4TablesForExtrapolator::Initialisation()
{
  if(verbose>1) {
    G4cout << "### G4TablesForExtrapolator::Initialisation" << G4endl;
  }
  G4int num = (G4int)G4Material::GetNumberOfMaterials();
  if(nmat == num) { return; }
  nmat = num;
  cuts.resize(nmat, DBL_MAX);
  couples.resize(nmat, nullptr);

  const G4MaterialTable* mtable = G4Material::GetMaterialTable();
  if(!pcuts) { pcuts = new G4ProductionCuts(); }

  for(G4int i=0; i<nmat; ++i) {
    couples[i] = new G4MaterialCutsCouple((*mtable)[i],pcuts);
  }

  dedxElectron     = PrepareTable(dedxElectron);
  dedxPositron     = PrepareTable(dedxPositron);
  dedxMuon         = PrepareTable(dedxMuon);
  dedxProton       = PrepareTable(dedxProton);
  rangeElectron    = PrepareTable(rangeElectron);
  rangePositron    = PrepareTable(rangePositron);
  rangeMuon        = PrepareTable(rangeMuon);
  rangeProton      = PrepareTable(rangeProton);
  invRangeElectron = PrepareTable(invRangeElectron);
  invRangePositron = PrepareTable(invRangePositron);
  invRangeMuon     = PrepareTable(invRangeMuon);
  invRangeProton   = PrepareTable(invRangeProton);
  mscElectron      = PrepareTable(mscElectron);

  builder = new G4LossTableBuilder(true);
  builder->SetBaseMaterialActive(false);

  if(verbose>1) {
    G4cout << "### G4TablesForExtrapolator Builds electron tables" 
	   << G4endl;
  }
  ComputeElectronDEDX(electron, dedxElectron);
  builder->BuildRangeTable(dedxElectron,rangeElectron);  
  builder->BuildInverseRangeTable(rangeElectron, invRangeElectron);  

  if(verbose>1) {
    G4cout << "### G4TablesForExtrapolator Builds positron tables" 
	   << G4endl;
  }
  ComputeElectronDEDX(positron, dedxPositron);
  builder->BuildRangeTable(dedxPositron, rangePositron);  
  builder->BuildInverseRangeTable(rangePositron, invRangePositron);  

  if(verbose>1) {
    G4cout << "### G4TablesForExtrapolator Builds muon tables" << G4endl;
  }
  ComputeMuonDEDX(muonPlus, dedxMuon);
  builder->BuildRangeTable(dedxMuon, rangeMuon);  
  builder->BuildInverseRangeTable(rangeMuon, invRangeMuon);  
  
  if(verbose>2) {
    G4cout << "DEDX MUON" << G4endl;
    G4cout << *dedxMuon << G4endl;
    G4cout << "RANGE MUON" << G4endl;
    G4cout << *rangeMuon << G4endl;
    G4cout << "INVRANGE MUON" << G4endl;
    G4cout << *invRangeMuon << G4endl;
  }
  if(verbose>1) {
    G4cout << "### G4TablesForExtrapolator Builds proton tables" 
	   << G4endl;
  }
  ComputeProtonDEDX(proton, dedxProton);
  builder->BuildRangeTable(dedxProton, rangeProton);  
  builder->BuildInverseRangeTable(rangeProton, invRangeProton);  

  ComputeTrasportXS(electron, mscElectron);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsTable* G4TablesForExtrapolator::PrepareTable(G4PhysicsTable* ptr)
{
  G4PhysicsTable* table = ptr;
  if(nullptr == ptr) { table = new G4PhysicsTable(); }
  G4int n = (G4int)table->length();
  for(G4int i=n; i<nmat; ++i) {  
    G4PhysicsVector* v = new G4PhysicsLogVector(emin, emax, nbins, splineFlag);
    table->push_back(v);
  }
  return table;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4TablesForExtrapolator::ComputeElectronDEDX(
                                  const G4ParticleDefinition* part,
				  G4PhysicsTable* table) 
{
  G4MollerBhabhaModel* ioni = new G4MollerBhabhaModel();
  G4eBremsstrahlungRelModel* brem = new G4eBremsstrahlungRelModel();
  ioni->Initialise(part, cuts);
  brem->Initialise(part, cuts);
  ioni->SetUseBaseMaterials(false);
  brem->SetUseBaseMaterials(false);

  mass    = electron_mass_c2;
  charge2 = 1.0;
  currentParticle = part;

  const G4MaterialTable* mtable = G4Material::GetMaterialTable();
  if(0<verbose) {
    G4cout << "G4TablesForExtrapolator::ComputeElectronDEDX for " 
           << part->GetParticleName() 
           << G4endl;
  }
  for(G4int i=0; i<nmat; ++i) {  

    const G4Material* mat = (*mtable)[i];
    if(1<verbose) {
      G4cout << "i= " << i << "  mat= " << mat->GetName() << G4endl;
    }
    G4PhysicsVector* aVector = (*table)[i];

    for(G4int j=0; j<=nbins; ++j) {
        
       G4double e = aVector->Energy(j);
       G4double dedx = ioni->ComputeDEDXPerVolume(mat,part,e,e) 
	 + brem->ComputeDEDXPerVolume(mat,part,e,e);
       if(1<verbose) {
         G4cout << "j= " << j
                << "  e(MeV)= " << e/MeV  
                << " dedx(Mev/cm)= " << dedx*cm/MeV
                << " dedx(Mev.cm2/g)= " 
		<< dedx/((MeV*mat->GetDensity())/(g/cm2)) 
		<< G4endl;
       }
       aVector->PutValue(j,dedx);
    }
    if(splineFlag) { aVector->FillSecondDerivatives(); }
  }
  delete ioni;
  delete brem;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4TablesForExtrapolator::ComputeMuonDEDX(const G4ParticleDefinition* part, 
					 G4PhysicsTable* table)
{
  G4BetheBlochModel* ioni = new G4BetheBlochModel();
  G4MuPairProductionModel* pair = new G4MuPairProductionModel();
  G4MuBremsstrahlungModel* brem = new G4MuBremsstrahlungModel();
  ioni->Initialise(part, cuts);
  pair->Initialise(part, cuts);
  brem->Initialise(part, cuts);
  ioni->SetUseBaseMaterials(false);
  pair->SetUseBaseMaterials(false);
  brem->SetUseBaseMaterials(false);

  mass    = part->GetPDGMass();
  charge2 = 1.0;
  currentParticle = part;

  const G4MaterialTable* mtable = G4Material::GetMaterialTable();

  if(0<verbose) {
    G4cout << "G4TablesForExtrapolator::ComputeMuonDEDX for " 
	   << part->GetParticleName() 
           << G4endl;
  }
 
  for(G4int i=0; i<nmat; ++i) {  

    const G4Material* mat = (*mtable)[i];
    if(1<verbose) {
      G4cout << "i= " << i << "  mat= " << mat->GetName() << G4endl;
    }
    G4PhysicsVector* aVector = (*table)[i];
    for(G4int j=0; j<=nbins; ++j) {
        
       G4double e = aVector->Energy(j);
       G4double dedx = ioni->ComputeDEDXPerVolume(mat,part,e,e) +
	               pair->ComputeDEDXPerVolume(mat,part,e,e) +
	               brem->ComputeDEDXPerVolume(mat,part,e,e);
       aVector->PutValue(j,dedx);
       if(1<verbose) {
         G4cout << "j= " << j
                << "  e(MeV)= " << e/MeV  
                << " dedx(Mev/cm)= " << dedx*cm/MeV
                << " dedx(Mev/(g/cm2)= " 
                << dedx/((MeV*mat->GetDensity())/(g/cm2))
		<< G4endl;
       }
    }
    if(splineFlag) { aVector->FillSecondDerivatives(); }
  }
  delete ioni;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4TablesForExtrapolator::ComputeProtonDEDX(const G4ParticleDefinition* part, 
                        		   G4PhysicsTable* table)
{
  G4BetheBlochModel* ioni = new G4BetheBlochModel();
  ioni->Initialise(part, cuts);
  ioni->SetUseBaseMaterials(false);

  mass    = part->GetPDGMass();
  charge2 = 1.0;
  currentParticle = part;

  const G4MaterialTable* mtable = G4Material::GetMaterialTable();

  if(0<verbose) {
    G4cout << "G4TablesForExtrapolator::ComputeProtonDEDX for " 
	   << part->GetParticleName() 
           << G4endl;
  }
 
  for(G4int i=0; i<nmat; ++i) {  

    const G4Material* mat = (*mtable)[i];
    if(1<verbose)
      G4cout << "i= " << i << "  mat= " << mat->GetName() << G4endl;
    G4PhysicsVector* aVector = (*table)[i];
    for(G4int j=0; j<=nbins; ++j) {
        
       G4double e = aVector->Energy(j);
       G4double dedx = ioni->ComputeDEDXPerVolume(mat,part,e,e);
       aVector->PutValue(j,dedx);
       if(1<verbose) {
         G4cout << "j= " << j
                << "  e(MeV)= " << e/MeV  
                << " dedx(Mev/cm)= " << dedx*cm/MeV
                << " dedx(Mev.cm2/g)= " << dedx/((mat->GetDensity())/(g/cm2)) 
		<< G4endl;
       }
    }
    if(splineFlag) { aVector->FillSecondDerivatives(); }
  }
  delete ioni;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4TablesForExtrapolator::ComputeTrasportXS(const G4ParticleDefinition* part, 
					   G4PhysicsTable* table)
{
  G4WentzelVIModel* msc = new G4WentzelVIModel();
  msc->SetPolarAngleLimit(CLHEP::pi);
  msc->Initialise(part, cuts);
  msc->SetUseBaseMaterials(false);

  mass    = part->GetPDGMass();
  charge2 = 1.0;
  currentParticle = part;

  const G4MaterialTable* mtable = G4Material::GetMaterialTable();

  if(0<verbose) {
    G4cout << "G4TablesForExtrapolator::ComputeTransportXS for " 
	   << part->GetParticleName() 
           << G4endl;
  }
 
  for(G4int i=0; i<nmat; ++i) {  

    const G4Material* mat = (*mtable)[i];
    msc->SetCurrentCouple(couples[i]);
    if(1<verbose)
      G4cout << "i= " << i << "  mat= " << mat->GetName() << G4endl;
    G4PhysicsVector* aVector = (*table)[i];
    for(G4int j=0; j<=nbins; ++j) {
        
       G4double e = aVector->Energy(j);
       G4double xs = msc->CrossSectionPerVolume(mat,part,e);
       aVector->PutValue(j,xs);
       if(1<verbose) {
         G4cout << "j= " << j << "  e(MeV)= " << e/MeV  
                << " xs(1/mm)= " << xs*mm << G4endl;
       }
    }
    if(splineFlag) { aVector->FillSecondDerivatives(); }
  }
  delete msc;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

