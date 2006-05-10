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
// $Id: G4ionIonisation.cc,v 1.35 2006-05-10 19:38:43 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4ionIonisation
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 07.05.2002
//
// Modifications:
//
// 23-12-02 Change interface in order to move to cut per region (V.Ivanchenko)
// 26-12-02 Secondary production moved to derived classes (V.Ivanchenko)
// 13-02-03 SubCutoff regime is assigned to a region (V.Ivanchenko)
// 18-04-03 Use IonFluctuations (V.Ivanchenko)
// 03-08-03 Add effective charge (V.Ivanchenko)
// 12-11-03 G4EnergyLossSTD -> G4EnergyLossProcess (V.Ivanchenko)
// 27-05-04 Set integral to be a default regime (V.Ivanchenko)
// 08-11-04 Migration to new interface of Store/Retrieve tables (V.Ivantchenko)
// 08-04-05 Major optimisation of internal interfaces (V.Ivantchenko)
// 10-01-06 SetStepLimits -> SetStepFunction (V.Ivantchenko)
// 10-05-06 Add a possibility to download user data (V.Ivantchenko)
//
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4ionIonisation.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4GenericIon.hh"
#include "G4BraggModel.hh"
#include "G4BraggIonModel.hh"
#include "G4BetheBlochModel.hh"
#include "G4IonFluctuations.hh"
#include "G4UniversalFluctuation.hh"
#include "G4UnitsTable.hh"
#include "G4LossTableManager.hh"
#include "G4NistManager.hh"
#include "G4LPhysicsFreeVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4ionIonisation::G4ionIonisation(const G4String& name)
  : G4VEnergyLossProcess(name),
    theParticle(0),
    theBaseParticle(0),
    isInitialised(false),
    stopDataActive(false)
{
  SetDEDXBinning(120);
  SetLambdaBinning(120);
  SetMinKinEnergy(0.1*keV);
  SetMaxKinEnergy(100.0*TeV);
  SetLinearLossLimit(0.15);
  SetStepFunction(0.1, 0.1*mm);
  SetIntegral(true);
  SetVerboseLevel(1);
  curMaterial = 0;
  curParticle = 0;
  theBraggModel = 0;
  corr = G4LossTableManager::Instance()->EmCorrections();  
  nist = G4NistManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ionIonisation::~G4ionIonisation()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ionIonisation::InitialiseEnergyLossProcess(
		      const G4ParticleDefinition* part,
		      const G4ParticleDefinition* bpart)
{
  if(isInitialised) return;

  theParticle = part;

  if(part == bpart || part == G4GenericIon::GenericIon()) theBaseParticle = 0;
  else if(bpart == 0) theBaseParticle = G4GenericIon::GenericIon();
  else                theBaseParticle = bpart;

  SetBaseParticle(theBaseParticle);
  SetSecondaryParticle(G4Electron::Electron());

  if(theBaseParticle) baseMass = theBaseParticle->GetPDGMass();
  else                baseMass = theParticle->GetPDGMass();

  flucModel = new G4IonFluctuations();

  eth = 2.0*MeV;

  theBraggModel = new G4BraggIonModel();
  theBraggModel->SetLowEnergyLimit(0.1*keV);
  theBraggModel->SetHighEnergyLimit(eth);
  AddEmModel(1, theBraggModel, flucModel);
  G4VEmModel* em1 = new G4BetheBlochModel();
  em1->SetLowEnergyLimit(eth);
  em1->SetHighEnergyLimit(100.0*TeV);
  AddEmModel(2, em1, flucModel);
  isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ionIonisation::PrintInfo()
{
  G4cout << "      Scaling relation is used to proton dE/dx and range"
         << G4endl
         << "      Bether-Bloch model for Escaled > " << eth << " MeV, ICRU49 "
         << "parametrisation for alpha particles below.";
  if(stopDataActive) 
    G4cout << G4endl << "      Stopping Power data for " << nIons
	   << " ion/material pairs are used.";
  G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ionIonisation::GetMeanFreePath(const G4Track& track,
					  G4double, 
					  G4ForceCondition* cond)
{
  DefineMassCharge(track.GetDefinition(), 
		   track.GetMaterial(),
		   track.GetDynamicParticle()->GetMass(),
		   track.GetKineticEnergy());
  return G4VEnergyLossProcess::GetMeanFreePath(track, 0.0, cond);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ionIonisation::EffectiveChargeCorrection(const G4ParticleDefinition* p,
						    const G4Material* mat,
						    G4double kinEnergy)
{
  G4double factor = 1.0;
  if(p->GetPDGCharge() <= 2.5*eplus) return factor;
  G4PhysicsVector* v = 0;
  if(p != curParticle || mat != curMaterial) {
    curParticle = p;
    curMaterial = 0;
    G4int Z = p->GetAtomicNumber();
    G4int A = p->GetAtomicMass();
    massFactor = proton_mass_c2/nist->GetIsotopeMass(Z,A);
    idx = 0;
    for(; idx<nIons; idx++) {
      if(Z == Zion[idx] && A == Aion[idx]) {
        if(materialList[idx] == mat) {
	  curMaterial = mat;
          v = stopData[idx];
	  break;
        } else if(materialList[idx] == 0) {
	  if(materialName[idx] == mat->GetName()) 
	    v = InitialiseMaterial(mat);
	}
      }
    }
  }
  if(v) {
    G4bool b;
    factor = v->GetValue(kinEnergy*massFactor,b);
  }
  return factor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ionIonisation::AddStoppingData(G4int Z, G4int A, 
				      const G4String& mname,
				      G4PhysicsVector& dVector)
{
  idx = 0;
  for(; idx<nIons; idx++) {
    if(Z == Zion[idx] && A == Aion[idx] && mname == materialName[idx]) 
      break;
  }
  if(idx == nIons) {
    Zion.push_back(Z);
    Aion.push_back(A);
    materialName.push_back(mname);
    materialList.push_back(0);
    stopData.push_back(0);
    nIons++;
  } else {
    delete stopData[idx];
  }
  massFactor = proton_mass_c2/nist->GetIsotopeMass(Z,A);
  size_t nbins = dVector.GetVectorLength();
  size_t n = 0;
  for(; n<nbins; n++) {
    if(dVector.GetLowEdgeEnergy(n)*massFactor > eth) break;
  } 
  if(n < nbins) nbins = n + 1;
  G4LPhysicsFreeVector* v = 
    new G4LPhysicsFreeVector(nbins, 
			     dVector.GetLowEdgeEnergy(0)*massFactor,
			     dVector.GetLowEdgeEnergy(nbins-1)*massFactor);
  G4bool b;
  for(size_t i=0; i<nbins; i++) {
    G4double e = dVector.GetLowEdgeEnergy(n);
    G4double dedx = dVector.GetValue(e, b);
    v->PutValues(i, e*massFactor, dedx);
  }
  stopData[idx] = v;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsVector* G4ionIonisation::InitialiseMaterial(const G4Material* mat)
{
  G4PhysicsVector* v = 0;
  const G4Material* m = nist->FindOrBuildMaterial(materialName[idx],false);
  if(m) {
    materialList[idx] = m;
    curMaterial = mat;
    v = stopData[idx];
    size_t nbins = v->GetVectorLength();
    const G4ParticleDefinition* p = G4Proton::Proton();
    G4bool b;
    for(size_t i=0; i<nbins; i++) {
      G4double e = v->GetLowEdgeEnergy(i);
      G4double dedx = v->GetValue(e, b);
      G4double dedx1= theBraggModel->ComputeDEDXPerVolume(mat, p, e, e)*
	effCharge.EffectiveChargeSquareRatio(curParticle,mat,e/massFactor);
      v->PutValue(i, dedx/dedx1);
    }
  }
  return v;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
