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
// $Id: G4EnergyLossForExtrapolator.cc,v 1.21 2010-11-04 12:40:29 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:    G4EnergyLossForExtrapolator
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

#include "G4EnergyLossForExtrapolator.hh"
#include "G4PhysicsLogVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4ParticleTable.hh"
#include "G4LossTableBuilder.hh"
#include "G4MollerBhabhaModel.hh"
#include "G4BetheBlochModel.hh"
#include "G4eBremsstrahlungModel.hh"
#include "G4MuPairProductionModel.hh"
#include "G4MuBremsstrahlungModel.hh"
#include "G4ProductionCuts.hh"
#include "G4LossTableManager.hh"
#include "G4WentzelVIModel.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EnergyLossForExtrapolator::G4EnergyLossForExtrapolator(G4int verb)
  :maxEnergyTransfer(DBL_MAX),verbose(verb),isInitialised(false)
{
  currentParticle = 0;
  currentMaterial = 0;

  linLossLimit = 0.01;
  emin         = 1.*MeV;
  emax         = 10.*TeV;
  nbins        = 70;

  nmat = index = 0;
  cuts = 0;

  mass = charge2 = electronDensity = radLength = bg2 = beta2 
    = kineticEnergy = tmax = 0;
  gam = 1.0;

  dedxElectron = dedxPositron = dedxProton = rangeElectron 
    = rangePositron = rangeProton = invRangeElectron = invRangePositron 
    = invRangeProton = mscElectron = dedxMuon = rangeMuon = invRangeMuon = 0;
  cuts = 0;
  electron = positron = proton = muonPlus = muonMinus = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EnergyLossForExtrapolator:: ~G4EnergyLossForExtrapolator()
{
  for(G4int i=0; i<nmat; i++) {delete couples[i];}
  delete dedxElectron;
  delete dedxPositron;
  delete dedxProton;
  delete dedxMuon;
  delete rangeElectron;
  delete rangePositron;
  delete rangeProton;
  delete rangeMuon;
  delete invRangeElectron;
  delete invRangePositron;
  delete invRangeProton;
  delete invRangeMuon;
  delete mscElectron;
  delete cuts;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EnergyLossForExtrapolator::EnergyAfterStep(G4double kinEnergy, 
						      G4double stepLength, 
						      const G4Material* mat, 
						      const G4ParticleDefinition* part)
{
  if(!isInitialised) Initialisation();
  G4double kinEnergyFinal = kinEnergy;
  if(SetupKinematics(part, mat, kinEnergy)) {
    G4double step = TrueStepLength(kinEnergy,stepLength,mat,part);
    G4double r  = ComputeRange(kinEnergy,part);
    if(r <= step) {
      kinEnergyFinal = 0.0;
    } else if(step < linLossLimit*r) {
      kinEnergyFinal -= step*ComputeDEDX(kinEnergy,part);
    } else {  
      G4double r1 = r - step;
      kinEnergyFinal = ComputeEnergy(r1,part);
    }
  }
  return kinEnergyFinal;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EnergyLossForExtrapolator::EnergyBeforeStep(G4double kinEnergy, 
						       G4double stepLength, 
						       const G4Material* mat, 
						       const G4ParticleDefinition* part)
{
  if(!isInitialised) Initialisation();
  G4double kinEnergyFinal = kinEnergy;

  if(SetupKinematics(part, mat, kinEnergy)) {
    G4double step = TrueStepLength(kinEnergy,stepLength,mat,part);
    G4double r  = ComputeRange(kinEnergy,part);

    if(step < linLossLimit*r) {
      kinEnergyFinal += step*ComputeDEDX(kinEnergy,part);
    } else {  
      G4double r1 = r + step;
      kinEnergyFinal = ComputeEnergy(r1,part);
    }
  }
  return kinEnergyFinal;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EnergyLossForExtrapolator::TrueStepLength(G4double kinEnergy, 
						     G4double stepLength,
						     const G4Material* mat, 
						     const G4ParticleDefinition* part)
{
  G4double res = stepLength;
  if(!isInitialised) Initialisation();
  if(SetupKinematics(part, mat, kinEnergy)) {
    if(part == electron || part == positron) {
      G4double x = stepLength*ComputeValue(kinEnergy, mscElectron);
      if(x < 0.2) res *= (1.0 + 0.5*x + x*x/3.0);
      else if(x < 0.9999) res = -std::log(1.0 - x)*stepLength/x;
      else res = ComputeRange(kinEnergy,part);
    
    } else {
      res = ComputeTrueStep(mat,part,kinEnergy,stepLength);
    }
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4EnergyLossForExtrapolator::SetupKinematics(const G4ParticleDefinition* part, 
						    const G4Material* mat, 
						    G4double kinEnergy)
{
  if(!part || !mat || kinEnergy < keV) return false;
  if(!isInitialised) Initialisation();
  G4bool flag = false;
  if(part != currentParticle) {
    flag = true;
    currentParticle = part;
    mass = part->GetPDGMass();
    G4double q = part->GetPDGCharge()/eplus;
    charge2 = q*q;
  }
  if(mat != currentMaterial) {
    G4int i = mat->GetIndex();
    if(i >= nmat) {
      G4cout << "### G4EnergyLossForExtrapolator WARNING:index i= " 
	     << i << " is out of table - NO extrapolation" << G4endl;
    } else {
      flag = true;
      currentMaterial = mat;
      electronDensity = mat->GetElectronDensity();
      radLength       = mat->GetRadlen();
      index           = i;
    }
  }
  if(flag || kinEnergy != kineticEnergy) {
    kineticEnergy = kinEnergy;
    G4double tau  = kinEnergy/mass;

    gam   = tau + 1.0;
    bg2   = tau * (tau + 2.0);
    beta2 = bg2/(gam*gam);
    tmax  = kinEnergy;
    if(part == electron) tmax *= 0.5;
    else if(part != positron) {
      G4double r = electron_mass_c2/mass;
      tmax = 2.0*bg2*electron_mass_c2/(1.0 + 2.0*gam*r + r*r);
    }
    if(tmax > maxEnergyTransfer) tmax = maxEnergyTransfer;
  }
  return true;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EnergyLossForExtrapolator::Initialisation()
{
  isInitialised = true;
  if(verbose>1) {
    G4cout << "### G4EnergyLossForExtrapolator::Initialisation" << G4endl;
  }
  currentParticle = 0;
  currentMaterial = 0;
  kineticEnergy   = 0.0;

  electron = G4Electron::Electron();
  positron = G4Positron::Positron();
  proton   = G4Proton::Proton();
  muonPlus = G4MuonPlus::MuonPlus();
  muonMinus= G4MuonMinus::MuonMinus();

  currentParticleName = "";

  nmat = G4Material::GetNumberOfMaterials();
  const G4MaterialTable* mtable = G4Material::GetMaterialTable();
  cuts = new G4ProductionCuts();

  const G4MaterialCutsCouple* couple;
  for(G4int i=0; i<nmat; i++) {
    couple = new G4MaterialCutsCouple((*mtable)[i],cuts);  
    couples.push_back(couple);
  }

  dedxElectron     = PrepareTable();
  dedxPositron     = PrepareTable();
  dedxMuon         = PrepareTable();
  dedxProton       = PrepareTable();
  rangeElectron    = PrepareTable();
  rangePositron    = PrepareTable();
  rangeMuon        = PrepareTable();
  rangeProton      = PrepareTable();
  invRangeElectron = PrepareTable();
  invRangePositron = PrepareTable();
  invRangeMuon     = PrepareTable();
  invRangeProton   = PrepareTable();
  mscElectron      = PrepareTable();

  G4LossTableBuilder builder; 

  if(verbose>1) 
    G4cout << "### G4EnergyLossForExtrapolator Builds electron tables" << G4endl;

  ComputeElectronDEDX(electron, dedxElectron);
  builder.BuildRangeTable(dedxElectron,rangeElectron);  
  builder.BuildInverseRangeTable(rangeElectron, invRangeElectron);  

  if(verbose>1) 
    G4cout << "### G4EnergyLossForExtrapolator Builds positron tables" << G4endl;

  ComputeElectronDEDX(positron, dedxPositron);
  builder.BuildRangeTable(dedxPositron, rangePositron);  
  builder.BuildInverseRangeTable(rangePositron, invRangePositron);  

  if(verbose>1) 
    G4cout << "### G4EnergyLossForExtrapolator Builds muon tables" << G4endl;

  ComputeMuonDEDX(muonPlus, dedxMuon);
  builder.BuildRangeTable(dedxMuon, rangeMuon);  
  builder.BuildInverseRangeTable(rangeMuon, invRangeMuon);  

  if(verbose>1) 
    G4cout << "### G4EnergyLossForExtrapolator Builds proton tables" << G4endl;

  ComputeProtonDEDX(proton, dedxProton);
  builder.BuildRangeTable(dedxProton, rangeProton);  
  builder.BuildInverseRangeTable(rangeProton, invRangeProton);  

  ComputeTrasportXS(electron, mscElectron);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsTable* G4EnergyLossForExtrapolator::PrepareTable()
{
  G4PhysicsTable* table = new G4PhysicsTable();

  for(G4int i=0; i<nmat; i++) {  

    G4PhysicsVector* v = new G4PhysicsLogVector(emin, emax, nbins);
    v->SetSpline(G4LossTableManager::Instance()->SplineFlag());
    table->push_back(v);
  }
  return table;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4ParticleDefinition* G4EnergyLossForExtrapolator::FindParticle(const G4String& name)
{
  const G4ParticleDefinition* p = 0;
  if(name != currentParticleName) {
    p = G4ParticleTable::GetParticleTable()->FindParticle(name);
    if(!p) {
      G4cout << "### G4EnergyLossForExtrapolator WARNING:FindParticle fails to find " 
	     << name << G4endl;
    }
  } else {
    p = currentParticle;
  }
  return p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EnergyLossForExtrapolator::ComputeDEDX(G4double kinEnergy, 
						  const G4ParticleDefinition* part)
{
  G4double x = 0.0;
  if(part == electron)      x = ComputeValue(kinEnergy, dedxElectron);
  else if(part == positron) x = ComputeValue(kinEnergy, dedxPositron);
  else if(part == muonPlus || part == muonMinus) {
    x = ComputeValue(kinEnergy, dedxMuon);
  } else {
    G4double e = kinEnergy*proton_mass_c2/mass;
    x = ComputeValue(e, dedxProton)*charge2;
  }
  return x;
}
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EnergyLossForExtrapolator::ComputeRange(G4double kinEnergy, 
						   const G4ParticleDefinition* part)
{
  G4double x = 0.0;
  if(part == electron)      x = ComputeValue(kinEnergy, rangeElectron);
  else if(part == positron) x = ComputeValue(kinEnergy, rangePositron);
  else if(part == muonPlus || part == muonMinus) 
    x = ComputeValue(kinEnergy, rangeMuon);
  else {
    G4double massratio = proton_mass_c2/mass;
    G4double e = kinEnergy*massratio;
    x = ComputeValue(e, rangeProton)/(charge2*massratio);
  }
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EnergyLossForExtrapolator::ComputeEnergy(G4double range, 
						    const G4ParticleDefinition* part)
{
  G4double x = 0.0;
  if(part == electron)      x = ComputeValue(range, invRangeElectron);
  else if(part == positron) x = ComputeValue(range, invRangePositron);
  else if(part == muonPlus || part == muonMinus) 
    x = ComputeValue(range, invRangeMuon);
  else {
    G4double massratio = proton_mass_c2/mass;
    G4double r = range*massratio*charge2;
    x = ComputeValue(r, invRangeProton)/massratio;
  }
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EnergyLossForExtrapolator::ComputeElectronDEDX(const G4ParticleDefinition* part, 
						      G4PhysicsTable* table) 
{
  G4DataVector v;
  G4MollerBhabhaModel* ioni = new G4MollerBhabhaModel();
  G4eBremsstrahlungModel* brem = new G4eBremsstrahlungModel();
  ioni->Initialise(part, v);
  brem->Initialise(part, v);

  mass    = electron_mass_c2;
  charge2 = 1.0;
  currentParticle = part;

  const G4MaterialTable* mtable = G4Material::GetMaterialTable();
  if(0<verbose) {
    G4cout << "G4EnergyLossForExtrapolator::ComputeElectronDEDX for " 
           << part->GetParticleName() 
           << G4endl;
  }
  for(G4int i=0; i<nmat; i++) {  

    const G4Material* mat = (*mtable)[i];
    if(1<verbose)
      G4cout << "i= " << i << "  mat= " << mat->GetName() << G4endl;
    const G4MaterialCutsCouple* couple = couples[i];
    G4PhysicsVector* aVector = (*table)[i];

    for(G4int j=0; j<=nbins; j++) {
        
       G4double e = aVector->Energy(j);
       G4double dedx = ioni->ComputeDEDX(couple,part,e,e) + brem->ComputeDEDX(couple,part,e,e);
       if(1<verbose) {
         G4cout << "j= " << j
                << "  e(MeV)= " << e/MeV  
                << " dedx(Mev/cm)= " << dedx*cm/MeV
                << " dedx(Mev.cm2/g)= " << dedx/((MeV*mat->GetDensity())/(g/cm2)) << G4endl;
       }
       aVector->PutValue(j,dedx);
    }
  }
  delete ioni;
  delete brem;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EnergyLossForExtrapolator::ComputeMuonDEDX(const G4ParticleDefinition* part, 
						  G4PhysicsTable* table)
{
  G4DataVector v;
  G4BetheBlochModel* ioni = new G4BetheBlochModel();
  G4MuPairProductionModel* pair = new G4MuPairProductionModel();
  G4MuBremsstrahlungModel* brem = new G4MuBremsstrahlungModel();
  ioni->Initialise(part, v);
  pair->Initialise(part, v);
  brem->Initialise(part, v);

  mass    = part->GetPDGMass();
  charge2 = 1.0;
  currentParticle = part;

  const G4MaterialTable* mtable = G4Material::GetMaterialTable();

  if(0<verbose) {
    G4cout << "G4EnergyLossForExtrapolator::ComputeMuonDEDX for " << part->GetParticleName() 
           << G4endl;
  }
 
  for(G4int i=0; i<nmat; i++) {  

    const G4Material* mat = (*mtable)[i];
    if(1<verbose)
      G4cout << "i= " << i << "  mat= " << mat->GetName() << G4endl;
    const G4MaterialCutsCouple* couple = couples[i];
    G4PhysicsVector* aVector = (*table)[i];
    for(G4int j=0; j<=nbins; j++) {
        
       G4double e = aVector->Energy(j);
       G4double dedx = ioni->ComputeDEDX(couple,part,e,e) +
	               pair->ComputeDEDX(couple,part,e,e) +
	               brem->ComputeDEDX(couple,part,e,e);
       aVector->PutValue(j,dedx);
       if(1<verbose) {
         G4cout << "j= " << j
                << "  e(MeV)= " << e/MeV  
                << " dedx(Mev/cm)= " << dedx*cm/MeV
                << " dedx(Mev/(g/cm2)= " << dedx/((MeV*mat->GetDensity())/(g/cm2)) << G4endl;
       }
    }
  }
  delete ioni;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EnergyLossForExtrapolator::ComputeProtonDEDX(const G4ParticleDefinition* part, 
						    G4PhysicsTable* table)
{
  G4DataVector v;
  G4BetheBlochModel* ioni = new G4BetheBlochModel();
  ioni->Initialise(part, v);

  mass    = part->GetPDGMass();
  charge2 = 1.0;
  currentParticle = part;

  const G4MaterialTable* mtable = G4Material::GetMaterialTable();

  if(0<verbose) {
    G4cout << "G4EnergyLossForExtrapolator::ComputeProtonDEDX for " << part->GetParticleName() 
           << G4endl;
  }
 
  for(G4int i=0; i<nmat; i++) {  

    const G4Material* mat = (*mtable)[i];
    if(1<verbose)
      G4cout << "i= " << i << "  mat= " << mat->GetName() << G4endl;
    const G4MaterialCutsCouple* couple = couples[i];
    G4PhysicsVector* aVector = (*table)[i];
    for(G4int j=0; j<=nbins; j++) {
        
       G4double e = aVector->Energy(j);
       G4double dedx = ioni->ComputeDEDX(couple,part,e,e);
       aVector->PutValue(j,dedx);
       if(1<verbose) {
         G4cout << "j= " << j
                << "  e(MeV)= " << e/MeV  
                << " dedx(Mev/cm)= " << dedx*cm/MeV
                << " dedx(Mev.cm2/g)= " << dedx/((mat->GetDensity())/(g/cm2)) << G4endl;
       }
    }
  }
  delete ioni;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EnergyLossForExtrapolator::ComputeTrasportXS(const G4ParticleDefinition* part, 
						    G4PhysicsTable* table)
{
  G4DataVector v;
  G4WentzelVIModel* msc = new G4WentzelVIModel();
  msc->SetPolarAngleLimit(CLHEP::pi);
  msc->Initialise(part, v);

  mass    = part->GetPDGMass();
  charge2 = 1.0;
  currentParticle = part;

  const G4MaterialTable* mtable = G4Material::GetMaterialTable();

  if(0<verbose) {
    G4cout << "G4EnergyLossForExtrapolator::ComputeProtonDEDX for " << part->GetParticleName() 
           << G4endl;
  }
 
  for(G4int i=0; i<nmat; i++) {  

    const G4Material* mat = (*mtable)[i];
    if(1<verbose)
      G4cout << "i= " << i << "  mat= " << mat->GetName() << G4endl;
    G4PhysicsVector* aVector = (*table)[i];
    for(G4int j=0; j<=nbins; j++) {
        
       G4double e = aVector->Energy(j);
       G4double xs = msc->CrossSectionPerVolume(mat,part,e);
       aVector->PutValue(j,xs);
       if(1<verbose) {
         G4cout << "j= " << j << "  e(MeV)= " << e/MeV  
                << " xs(1/mm)= " << xs*mm << G4endl;
       }
    }
  }
  delete msc;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

