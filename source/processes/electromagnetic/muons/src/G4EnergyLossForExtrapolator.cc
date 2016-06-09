//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: G4EnergyLossForExtrapolator.cc,v 1.5 2005/06/27 15:29:29 gunter Exp $
// GEANT4 tag $Name: geant4-07-01 $
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
// 08-04-05 Rename Propogator -> Extrapolator
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4EnergyLossForExtrapolator.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"
#include "G4ParticleTable.hh"
#include "G4LossTableBuilder.hh"
#include "G4MollerBhabhaModel.hh"
#include "G4BetheBlochModel.hh"
#include "G4eBremsstrahlungModel.hh"
#include "G4ProductionCuts.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EnergyLossForExtrapolator::G4EnergyLossForExtrapolator()
{
  Initialisation();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EnergyLossForExtrapolator:: ~G4EnergyLossForExtrapolator()
{
  for(G4int i=0; i<nmat; i++) {delete couples[i];}
  delete dedxElectron;
  delete dedxPositron;
  delete dedxProton;
  delete rangeElectron;
  delete rangePositron;
  delete rangeProton;
  delete invRangeElectron;
  delete invRangePositron;
  delete invRangeProton;
  delete cuts;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EnergyLossForExtrapolator::EnergyAfterStep(G4double kinEnergy, 
						      G4double step, 
						      const G4Material* mat, 
						      const G4ParticleDefinition* part)
{
  G4double kinEnergyFinal = kinEnergy;
  if(mat && part) {
    G4double r  = ComputeRange(kinEnergy,mat,part);
    if(r <= step) {
      kinEnergyFinal = 0.0;
    } else if(step < linLossLimit*r) {
      kinEnergyFinal -= step*ComputeDEDX(kinEnergy,mat,part);
    } else {  
      G4double r1 = r - step;
      kinEnergyFinal = ComputeEnergy(r1,mat,part);
    }
  }
  return kinEnergyFinal;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EnergyLossForExtrapolator::EnergyBeforeStep(G4double kinEnergy, 
						       G4double step, 
						       const G4Material* mat, 
						       const G4ParticleDefinition* part)
{
  G4double kinEnergyFinal = kinEnergy;
  if(mat && part) {
    G4double r  = ComputeRange(kinEnergy,mat,part);
    if(step < linLossLimit*r) {
      kinEnergyFinal += step*ComputeDEDX(kinEnergy,mat,part);
    } else {  
      G4double r1 = r + step;
      kinEnergyFinal = ComputeEnergy(r1,mat,part);
    }
  }
  return kinEnergyFinal;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EnergyLossForExtrapolator::AverageScatteringAngle(G4double kinEnergy, 
							     G4double step, 
							     const G4Material* mat, 
							     const G4ParticleDefinition* part)
{
  G4double theta = 0.0;
  if(mat && part && kinEnergy > 0.0) {
    G4double t = step/mat->GetRadlen();
    G4double M = part->GetPDGMass();
    G4double z = std::abs(part->GetPDGCharge()/eplus);
    G4double bpc = kinEnergy*(kinEnergy + 2.0*M)/(kinEnergy+M);
    theta = 13.6*MeV*z*std::sqrt(t)/bpc;
  }
  return theta;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
G4double G4EnergyLossForExtrapolator::EnergyDispersion(G4double kinEnergy, 
						       G4double step, 
						       const G4Material* mat, 
						       const G4ParticleDefinition* part)
{
  G4double sig2 = 0.0;
  if(part) {
    G4double tmax  = kinEnergy;
    G4double mass  = part->GetPDGMass();
    G4double tau   = kinEnergy/mass;
    G4double gam   = tau + 1.0;
    G4double bg2   = tau * (tau+2.0);
    G4double beta2 = bg2/(gam*gam);
    if(part == electron) tmax *= 0.5;
    else if(part != positron) {
      G4double r = electron_mass_c2/mass;
      tmax = 2.0*bg2*electron_mass_c2/(1.0 + 2.0*gam*r + r*r);
    }
    G4double electronDensity = mat->GetElectronDensity();
    G4double q = part->GetPDGCharge()/eplus; 
    sig2 = (1.0/beta2 - 0.5)* twopi_mc2_rcl2 * tmax*step * electronDensity * q *q;
  }
  return sig2;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EnergyLossForExtrapolator::Initialisation()
{
  currentParticle = 0;
  electron = G4Electron::Electron();
  positron = G4Positron::Positron();
  proton   = G4Proton::Proton();

  currentParticleName = "";

  linLossLimit = 0.05;
  emin         = 1.*MeV;
  emax         = 100.*GeV;
  nbins        = 50;
  verbose      = 1;

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
  dedxProton       = PrepareTable();
  rangeElectron    = PrepareTable();
  rangePositron    = PrepareTable();
  rangeProton      = PrepareTable();
  invRangeElectron = PrepareTable();
  invRangePositron = PrepareTable();
  invRangeProton   = PrepareTable();

  G4LossTableBuilder builder; 

  ComputeElectronDEDX(electron, dedxElectron);
  builder.BuildRangeTable(dedxElectron,rangeElectron);  
  builder.BuildInverseRangeTable(rangeElectron, invRangeElectron);  

  ComputeElectronDEDX(positron, dedxPositron);
  builder.BuildRangeTable(dedxPositron, rangePositron);  
  builder.BuildInverseRangeTable(rangePositron, invRangePositron);  

  ComputeProtonDEDX(proton, dedxProton);
  builder.BuildRangeTable(dedxProton, rangeProton);  
  builder.BuildInverseRangeTable(rangeProton, invRangeProton);  

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsTable* G4EnergyLossForExtrapolator::PrepareTable()
{
  G4PhysicsTable* table = new G4PhysicsTable();

  for(G4int i=0; i<nmat; i++) {  

    G4PhysicsVector* v = new G4PhysicsLogVector(emin, emax, nbins);
    table->push_back(v);
  }
  return table;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4ParticleDefinition* G4EnergyLossForExtrapolator::FindParticle(const G4String& name)
{
  if(name != currentParticleName) {
    currentParticle = G4ParticleTable::GetParticleTable()->FindParticle(name);
    if(!currentParticle) {
      G4cout << "### WARNING: G4EmCalculator::FindParticle fails to find " << name << G4endl;
    }
    currentParticleName = name;
  }
  return currentParticle;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EnergyLossForExtrapolator::ComputeDEDX(G4double kinEnergy, 
						  const G4Material* mat, 
						  const G4ParticleDefinition* part)
{
  G4double x = 0.0;
  if(part == electron)      x = ComputeValue(kinEnergy, mat, dedxElectron);
  else if(part == positron) x = ComputeValue(kinEnergy, mat, dedxPositron);
  else {
    G4double e = kinEnergy*proton_mass_c2/part->GetPDGMass();
    G4double z = part->GetPDGCharge()/eplus;
    x = ComputeValue(e, mat, dedxProton)*z*z;
  }
  return x;
}
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EnergyLossForExtrapolator::ComputeRange(G4double kinEnergy, 
						   const G4Material* mat, 
						   const G4ParticleDefinition* part)
{
  G4double x = 0.0;
  if(part == electron)      x = ComputeValue(kinEnergy, mat, rangeElectron);
  else if(part == positron) x = ComputeValue(kinEnergy, mat, rangePositron);
  else {
    G4double massratio = proton_mass_c2/part->GetPDGMass();
    G4double e = kinEnergy*massratio;
    G4double z = part->GetPDGCharge()/eplus;
    x = ComputeValue(e, mat, rangeProton)/(z*z*massratio);
  }
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EnergyLossForExtrapolator::ComputeEnergy(G4double range, 
						    const G4Material* mat, 
						    const G4ParticleDefinition* part)
{
  G4double x = 0.0;
  if(part == electron)      x = ComputeValue(range, mat, invRangeElectron);
  else if(part == positron) x = ComputeValue(range, mat, invRangePositron);
  else {
    G4double massratio = proton_mass_c2/part->GetPDGMass();
    G4double z = part->GetPDGCharge()/eplus;
    G4double r = range*massratio*z*z;
    x = ComputeValue(r, mat, invRangeProton)/massratio;
  }
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EnergyLossForExtrapolator::ComputeValue(G4double x, 
						   const G4Material* mat, 
						   const G4PhysicsTable* table)
{
  G4double res = 0.0;
  bool b;
  G4int i = mat->GetIndex();
  if(i<nmat) res = ((*table)[i])->GetValue(x, b);

  return res;
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

    for(G4int j=0; j<nbins; j++) {
        
       G4double e = aVector->GetLowEdgeEnergy(j);
       G4double dedx = ioni->ComputeDEDX(couple,part,e,e) + brem->ComputeDEDX(couple,part,e,e);
       dedx /= (MeV/cm);
       if(1<verbose) {
         G4cout << "j= " << j
                << "  e(MeV)= " << e/MeV  
                << " dedx(Mev/cm)= " << dedx
                << " dedx(Mev.cm2/g)= " << dedx/((mat->GetDensity())/(g/cm3)) << G4endl;
       }
       aVector->PutValue(j,dedx);
    }
  }
  delete ioni;
  delete brem;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EnergyLossForExtrapolator::ComputeProtonDEDX(const G4ParticleDefinition* part, 
						    G4PhysicsTable* table)
{
  G4DataVector v;
  G4BetheBlochModel* ioni = new G4BetheBlochModel();
  ioni->Initialise(part, v);

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
    for(G4int j=0; j<nbins; j++) {
        
       G4double e = aVector->GetLowEdgeEnergy(j);
       G4double dedx = ioni->ComputeDEDX(couple,part,e,e);
       dedx /= (MeV/cm);
       aVector->PutValue(j,dedx);
       if(1<verbose) {
         G4cout << "j= " << j
                << "  e(MeV)= " << e/MeV  
                << " dedx(Mev/cm)= " << dedx
                << " dedx(Mev.cm2/g)= " << dedx/((mat->GetDensity())/(g/cm3)) << G4endl;
       }
    }
  }
  delete ioni;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

