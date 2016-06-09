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
// $Id: G4EmCalculator.cc,v 1.19 2005/05/30 08:55:51 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4EmCalculator
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 28.06.2004
//
// Modifications:
// 12.09.2004 Add verbosity (V.Ivanchenko)
// 17.11.2004 Change signature of methods, add new methods (V.Ivanchenko)
// 08.04.2005 Major optimisation of internal interfaces (V.Ivantchenko)
// 08.05.2005 Use updated interfaces (V.Ivantchenko)
//
// Class Description:
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4EmCalculator.hh"
#include "G4LossTableManager.hh"
#include "G4VEmProcess.hh"
#include "G4VEnergyLossProcess.hh"
#include "G4VMultipleScattering.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicsTable.hh"
#include "G4ProductionCutsTable.hh"
#include "G4ProcessManager.hh"
#include "G4ionEffectiveCharge.hh"
#include "G4RegionStore.hh"
#include "G4Element.hh"
#include "G4EmCorrections.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmCalculator::G4EmCalculator()
{
  manager = G4LossTableManager::Instance();
  corr    = manager->EmCorrections();
  nLocalMaterials    = 0;
  verbose            = 0;
  currentCoupleIndex = 0;
  currentCouple      = 0;
  currentMaterial    = 0;
  currentParticle    = 0;
  currentLambda      = 0;
  chargeSquare       = 1.0;
  massRatio          = 1.0;
  currentParticleName= "";
  currentMaterialName= "";
  ionEffCharge       = new G4ionEffectiveCharge();
  isIon              = false;
  isApplicable       = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmCalculator::~G4EmCalculator()
{
  delete ionEffCharge;
  for (G4int i=0; i<nLocalMaterials; i++) {
    delete localCouples[i];
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::GetDEDX(G4double kinEnergy, const G4ParticleDefinition* p,
                                 const G4Material* mat, const G4Region* region)
{
  G4double res = 0.0;
  const G4MaterialCutsCouple* couple = FindCouple(mat, region);
  if(couple && UpdateParticle(p, kinEnergy) ) {
    res = manager->GetDEDX(p, kinEnergy, couple);
    if(verbose>0) {
      G4cout << "E(MeV)= " << kinEnergy/MeV
	     << " DEDX(MeV/mm)= " << res*mm/MeV
	     << " DEDX(MeV*cm^2/g)= " << res*gram/(MeV*cm2*mat->GetDensity())
	     << "  " <<  p->GetParticleName()
	     << " in " <<  mat->GetName()
	     << G4endl;
    }
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::GetDEDX(G4double kinEnergy, const G4String& particle,
                                 const G4String& material, const G4String& reg)
{
  return GetDEDX(kinEnergy,FindParticle(particle),FindMaterial(material),FindRegion(reg));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::GetRange(G4double kinEnergy, const G4ParticleDefinition* p,
                                  const G4Material* mat, const G4Region* region)
{
  G4double res = 0.0;
  const G4MaterialCutsCouple* couple = FindCouple(mat,region);
  if(couple && UpdateParticle(p, kinEnergy)) {
    res = manager->GetRange(p, kinEnergy, couple);
    if(verbose>0) {
      G4cout << "E(MeV)= " << kinEnergy/MeV
	     << " range(mm)= " << res/mm
	     << "  " <<  p->GetParticleName()
	     << " in " <<  mat->GetName()
	     << G4endl;
    }
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::GetRange(G4double kinEnergy, const G4String& particle,
                                  const G4String& material, const G4String& reg)
{
  return GetRange(kinEnergy,FindParticle(particle),FindMaterial(material),FindRegion(reg));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::GetKinEnergy(G4double range, const G4ParticleDefinition* p,
                                      const G4Material* mat, const G4Region* region)
{
  G4double res = 0.0;
  const G4MaterialCutsCouple* couple = FindCouple(mat,region);
  if(couple && UpdateParticle(p, 1.0*GeV)) {
    res = manager->GetEnergy(p, range, couple);
    if(verbose>0) {
      G4cout << "Range(mm)= " << range/mm
	     << " KinE(MeV)= " << res/MeV
	     << "  " <<  p->GetParticleName()
	     << " in " <<  mat->GetName()
	     << G4endl;
    }
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::GetKinEnergy(G4double range, const G4String& particle,
                                      const G4String& material, const G4String& reg)
{
  return GetKinEnergy(range,FindParticle(particle),FindMaterial(material),FindRegion(reg));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::GetCrossSectionPerVolume(G4double kinEnergy,
                                            const G4ParticleDefinition* p,
                                            const G4String& processName,
					    const G4Material* mat,
					    const G4Region* region)
{
  G4double res = 0.0;
  const G4MaterialCutsCouple* couple = FindCouple(mat,region);

  if(couple) {
    G4int idx = couple->GetIndex();
    FindLambdaTable(p, processName);
    if(currentLambda && UpdateParticle(p, kinEnergy)) {
      G4bool b;
      G4double e = kinEnergy*massRatio;
      res = (((*currentLambda)[idx])->GetValue(e,b))*chargeSquare;
      res /= mat->GetTotNbOfAtomsPerVolume();
      if(verbose>0) {
	G4cout << "E(MeV)= " << kinEnergy/MeV
	       << " cross(cm^2/g)= " << res*gram/cm2
	       << "  " <<  p->GetParticleName()
	       << " in " <<  mat->GetName()
	       << G4endl;
      }
    }
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::GetCrossSectionPerVolume(G4double kinEnergy,
                                            const G4String& particle,
					    const G4String& processName,
                                            const G4String& material,
					    const G4String& reg)
{
  return GetCrossSectionPerVolume(kinEnergy,FindParticle(particle),processName,
                                  FindMaterial(material),FindRegion(reg));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::GetCrossSectionPerAtom(G4double kinEnergy,
                                          const G4ParticleDefinition* p,
                                          const G4String& processName,
					  const G4Material* mat,
					  const G4Region* region)
{
  G4double res = GetCrossSectionPerVolume(kinEnergy,p,processName,mat,region);
  if(mat) res /= mat->GetTotNbOfAtomsPerVolume();
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::GetCrossSectionPerAtom(  G4double kinEnergy,
                                            const G4String& particle,
					    const G4String& processName,
                                            const G4String& material,
					    const G4String& reg)
{
  return GetCrossSectionPerAtom(kinEnergy,FindParticle(particle),processName,
                                FindMaterial(material),FindRegion(reg));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::GetMeanFreePath(G4double kinEnergy,
                                         const G4ParticleDefinition* p,
                                         const G4String& processName,
					 const G4Material* mat,
                                         const G4Region* region)
{
  G4double res = DBL_MAX;
  G4double x = GetCrossSectionPerVolume(kinEnergy,p, processName, mat,region);
  if(x > 0.0) res = 1.0/x;
  if(verbose>1) {
    G4cout << "E(MeV)= " << kinEnergy/MeV
	   << " MFP(mm)= " << res/mm
	   << "  " <<  p->GetParticleName()
	   << " in " <<  mat->GetName()
	   << G4endl;
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::GetMeanFreePath(G4double kinEnergy,
                                         const G4String& particle,
					 const G4String& processName,
                                         const G4String& material,
					 const G4String& reg)
{
  return GetMeanFreePath(kinEnergy,FindParticle(particle),processName,
                         FindMaterial(material),FindRegion(reg));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmCalculator::PrintDEDXTable(const G4ParticleDefinition* p)
{
  const G4VEnergyLossProcess* elp = FindEnergyLossProcess(p);
  G4cout << "##### DEDX Table for " << p->GetParticleName() << G4endl;
  if(elp) G4cout << *(elp->DEDXTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmCalculator::PrintRangeTable(const G4ParticleDefinition* p)
{
  const G4VEnergyLossProcess* elp = FindEnergyLossProcess(p);
  G4cout << "##### Range Table for " << p->GetParticleName() << G4endl;
  if(elp) G4cout << *(elp->RangeTableForLoss()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmCalculator::PrintInverseRangeTable(const G4ParticleDefinition* p)
{
  const G4VEnergyLossProcess* elp = FindEnergyLossProcess(p);
  G4cout << "### G4EmCalculator: Inverse Range Table for " << p->GetParticleName() << G4endl;
  if(elp) G4cout << *(elp->InverseRangeTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::ComputeDEDX(G4double kinEnergy,
                                     const G4ParticleDefinition* p,
                                     const G4String& processName,
				     const G4Material* mat,
                                           G4double cut)
{
  G4double res = 0.0;
  if(verbose > 1) {
    G4cout << "ComputeDEDX: " << p->GetParticleName()
           << " in " << mat->GetName()
           << " e(MeV)= " << kinEnergy/MeV << "  cut(MeV)= " << cut/MeV << G4endl;
  }
  if(FindEmModel(p, processName, kinEnergy)) {
    //    G4cout << "currentModel= " << currentModel << G4endl;
    if(UpdateParticle(p, kinEnergy)) {
      G4double escaled = kinEnergy*massRatio; 
      if(baseParticle) {
        res = currentModel->ComputeDEDXPerVolume(mat, baseParticle, escaled, cut)
            * chargeSquare;
        if(verbose > 1) {
          G4cout <<  baseParticle->GetParticleName() << " E(MeV)= " << escaled
                 << " DEDX(MeV/mm)= " << res*mm/MeV
                 << " DEDX(MeV*cm^2/g)= " << res*gram/(MeV*cm2*mat->GetDensity())
                 << G4endl;
	}          
      } else {
        res = currentModel->ComputeDEDXPerVolume(mat, p, kinEnergy, cut);
      }
      if(isIon && currentModel->HighEnergyLimit() > 100.*MeV) 
        res += corr->HighOrderCorrections(p,mat,kinEnergy);

      // emulate boundary region for different parameterisations
      G4double eth = currentModel->LowEnergyLimit();
      if(eth > 0.05*MeV && eth < 10.*MeV && escaled > eth) {
        G4double res1 = 0.0;
        if(baseParticle) {
          res1 = currentModel->ComputeDEDXPerVolume(mat, baseParticle, eth, cut)
               * chargeSquare;
	} else {
	  res1 = currentModel->ComputeDEDXPerVolume(mat, p, eth, cut);
	}
        if(isIon) res1 += corr->HighOrderCorrections(p,mat,eth/massRatio);  
        G4double res0 = res1;
        if(FindEmModel(p, processName, eth-1.0*keV)) {
	  if(baseParticle) {
	    res0 = currentModel->ComputeDEDXPerVolume(mat, baseParticle, eth, cut)
	         * chargeSquare;
	  } else {
	    res0 = currentModel->ComputeDEDXPerVolume(mat, p, eth, cut);
	  }
	}
        //G4cout << "eth= " << eth << " escaled= " << escaled << " res0= " << res0 << " res1= " 
        //       << res1 <<  "  q2= " << chargeSquare << G4endl; 
        res *= (1.0 + (res0/res1 - 1.0)*eth/escaled);
      }
      if(verbose > 0) {
        G4cout << "E(MeV)= " << kinEnergy/MeV
               << " DEDX(MeV/mm)= " << res*mm/MeV
               << " DEDX(MeV*cm^2/g)= " << res*gram/(MeV*cm2*mat->GetDensity())
               << " cut(MeV)= " << cut/MeV
               << "  " <<  p->GetParticleName()
               << " in " <<  mat->GetName()
               << " Zi^2= " << chargeSquare
               << G4endl;
      }
    }
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::ComputeDEDX(G4double kinEnergy,
                                     const G4String& particle,
				     const G4String& processName,
                                     const G4String& material,
                                           G4double cut)
{
  return ComputeDEDX(kinEnergy,FindParticle(particle),processName,
                     FindMaterial(material),cut);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::ComputeNuclearDEDX(G4double kinEnergy,
                                      const G4ParticleDefinition* p,
				      const G4Material* mat)
{

  G4double res = corr->NuclearDEDX(p, mat, kinEnergy, false);

  if(verbose > 1) {
    G4cout <<  p->GetParticleName() << " E(MeV)= " << kinEnergy/MeV
	   << " NuclearDEDX(MeV/mm)= " << res*mm/MeV
	   << " NuclearDEDX(MeV*cm^2/g)= " << res*gram/(MeV*cm2*mat->GetDensity())
	   << G4endl;
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::ComputeNuclearDEDX(G4double kinEnergy,
                                      const G4String& particle,
				      const G4String& material)
{
  return ComputeNuclearDEDX(kinEnergy,FindParticle(particle),FindMaterial(material));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::ComputeCrossSectionPerVolume(
                                                   G4double kinEnergy,
                                             const G4ParticleDefinition* p,
                                             const G4String& processName,
					     const G4Material* mat,
                                                   G4double cut)
{
  G4double res = 0.0;
  if(FindEmModel(p, processName, kinEnergy)) {
    UpdateParticle(p, kinEnergy);
    G4double e = kinEnergy;
    if(baseParticle) {
      e *= kinEnergy*massRatio;
      res = currentModel->CrossSectionPerVolume(mat, baseParticle, e, cut, e)*chargeSquare;
    } else {
      res = currentModel->CrossSectionPerVolume(mat, p, e, cut, e);
    }
    if(verbose>0) {
      G4cout << "E(MeV)= " << kinEnergy/MeV
	     << " cross(cm^2/g)= " << res*gram/cm2
	     << "  " <<  p->GetParticleName()
	     << " in " <<  mat->GetName()
	     << G4endl;
    }
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::ComputeCrossSectionPerVolume(
                                                   G4double kinEnergy,
                                             const G4String& particle,
					     const G4String& processName,
                                             const G4String& material,
                                                   G4double cut)
{
  return ComputeCrossSectionPerVolume(kinEnergy,FindParticle(particle),processName,
                                      FindMaterial(material),cut);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::ComputeCrossSectionPerAtom(
                                                   G4double kinEnergy,
					     const G4ParticleDefinition* p,
                                             const G4String& processName,
					           G4double Z, G4double A,
		                                   G4double cut)
{
  G4double res = 0.0;
  if(FindEmModel(p, processName, kinEnergy)) {
    UpdateParticle(p, kinEnergy);
    G4double e = kinEnergy;
    if(baseParticle) {
      e *= kinEnergy*massRatio;
      res = currentModel->ComputeCrossSectionPerAtom(baseParticle, e, Z, A, cut)*chargeSquare;
    } else {
      res = currentModel->ComputeCrossSectionPerAtom(p, e, Z, A, cut);
    }
    if(verbose>0) {
      G4cout << "E(MeV)= " << kinEnergy/MeV
	     << " cross(barn)= " << res/barn
	     << "  " <<  p->GetParticleName()
	     << " Z= " <<  Z << " A= " << A
	     << G4endl;
    }
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::ComputeCrossSectionPerAtom(G4double kinEnergy,
                                              const G4String& particle,
                                              const G4String& processName,
 					      const G4Element* elm,
		                                    G4double cut)
{
  return ComputeCrossSectionPerAtom(kinEnergy,FindParticle(particle),processName,
                                    elm->GetZ(),elm->GetN(),cut);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::ComputeMeanFreePath(G4double kinEnergy,
                                             const G4ParticleDefinition* p,
                                             const G4String& processName,
					     const G4Material* mat,
                                                   G4double cut)
{
  G4double mfp = DBL_MAX;
  G4double x = ComputeCrossSectionPerVolume(kinEnergy, p, processName, mat, cut);
  if(x > 0.0) mfp = 1.0/x;
  if(verbose>1) {
    G4cout << "E(MeV)= " << kinEnergy/MeV
	   << " MFP(mm)= " << mfp/mm
	   << "  " <<  p->GetParticleName()
	   << " in " <<  mat->GetName()
	   << G4endl;
  }
  return mfp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::ComputeMeanFreePath(G4double kinEnergy,
                                             const G4String& particle,
                                             const G4String& processName,
                                             const G4String& material,
                                                   G4double cut)
{
  return ComputeMeanFreePath(kinEnergy,FindParticle(particle),processName,
                             FindMaterial(material),cut);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4EmCalculator::UpdateParticle(const G4ParticleDefinition* p, G4double kinEnergy)
{
  if(p != currentParticle) {
    currentParticle = p;
    baseParticle    = 0;
    currentParticleName = p->GetParticleName();
  }
  isApplicable    = false;
  if(p->GetPDGCharge() != 0.0) {
    const G4VEnergyLossProcess* proc = FindEnergyLossProcess(p);
    if(proc) {
      isApplicable    = true;
      baseParticle    = proc->BaseParticle();
      massRatio       = 1.0;
      G4double qbase  = 1.0;
      G4double q      = p->GetPDGCharge()/eplus;

      if(baseParticle) {
        massRatio = baseParticle->GetPDGMass()/p->GetPDGMass();
        qbase = baseParticle->GetPDGCharge()/eplus;
      }

      if(p->GetParticleType() == "nucleus") {
	isIon = true;
        chargeSquare =
          ionEffCharge->EffectiveChargeSquareRatio(p, currentMaterial, kinEnergy);
	//G4double emax = currentModel->HighEnergyLimit(p);
	//if(emax < 100.*MeV) chargeSquare *=  q*q /
        //   ionEffCharge->EffectiveChargeSquareRatio(p, currentMaterial, emax/massRatio);

      } else {
	isIon = false;
	chargeSquare = q*q;
      }

      chargeSquare /= (qbase*qbase);

    }
  }
  return isApplicable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4ParticleDefinition* G4EmCalculator::FindParticle(const G4String& name)
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

const G4Material* G4EmCalculator::FindMaterial(const G4String& name)
{
  if(name != currentMaterialName) {
    currentMaterial = G4Material::GetMaterial(name);
    currentMaterialName = name;
    if(!currentMaterial)
      G4cout << "### WARNING: G4EmCalculator::FindMaterial fails to find " << name << G4endl;
  }
  return currentMaterial;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4Region* G4EmCalculator::FindRegion(const G4String& reg)
{
  return G4RegionStore::GetInstance()->GetRegion(reg);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4MaterialCutsCouple* G4EmCalculator::FindCouple(const G4Material* material,
                                                       const G4Region* region)
{
  if(!material) return 0;
  currentMaterial = material;
  currentMaterialName = material->GetName();
  // Access to materials
  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  const G4Region* r = region;
  if(!r) r = G4RegionStore::GetInstance()->GetRegion("DefaultRegionForTheWorld");

  return theCoupleTable->GetMaterialCutsCouple(material,r->GetProductionCuts());

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4EmCalculator::UpdateCouple(const G4Material* material, G4double cut)
{
  if(!material) return false;
  currentMaterial = material;
  currentMaterialName = material->GetName();
  for (G4int i=0; i<nLocalMaterials; i++) {
    if(material == localMaterials[i] && cut == localCuts[i]) {
      currentCouple = localCouples[i];
      currentCoupleIndex = currentCouple->GetIndex();
      currentCut = cut;
      return true;
    }
  }
  const G4MaterialCutsCouple* cc = new G4MaterialCutsCouple(material);
  localMaterials.push_back(material);
  localCouples.push_back(cc);
  localCuts.push_back(cut);
  nLocalMaterials++;
  currentCouple = cc;
  currentCoupleIndex = currentCouple->GetIndex();
  currentCut = cut;
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmCalculator::FindLambdaTable(const G4ParticleDefinition* p,
                                     const G4String& processName)
{
  // Search for the process
  if (p != currentParticle || processName != currentName) {
    currentName     = processName;
    currentLambda   = 0;
    G4LossTableManager* lManager = G4LossTableManager::Instance();
    const std::vector<G4VEnergyLossProcess*> vel = lManager->GetEnergyLossProcessVector();
    G4int n = vel.size();
    for(G4int i=0; i<n; i++) {
      if((vel[i])->GetProcessName() == currentName && (vel[i])->Particle() == p) {
        currentLambda = (vel[i])->LambdaTable();
	isApplicable    = true;
	break;
      }
    }
    if(!currentLambda) {
      const std::vector<G4VEmProcess*> vem = lManager->GetEmProcessVector();
      G4int n = vem.size();
      for(G4int i=0; i<n; i++) {
        if((vem[i])->GetProcessName() == currentName && (vem[i])->Particle() == p) {
          currentLambda = (vem[i])->LambdaTable();
          isApplicable    = true;
	  break;
        }
      }
    }
    if(!currentLambda) {
      const std::vector<G4VMultipleScattering*> vmsc = lManager->GetMultipleScatteringVector();
      G4int n = vmsc.size();
      for(G4int i=0; i<n; i++) {
        if((vmsc[i])->GetProcessName() == currentName && (vmsc[i])->Particle() == p) {
          currentLambda = (vmsc[i])->LambdaTable();
	  isApplicable    = true;
	  break;
        }
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4EmCalculator::FindEmModel(const G4ParticleDefinition* p,
                                   const G4String& processName,
	  	 		         G4double kinEnergy)
{
  G4bool res = false;
  if(verbose > 1) {
    G4cout << "G4EmCalculator::FindEmModel for " << p->GetParticleName()
           << " and " << processName << " at e(MeV)= " << kinEnergy
           << G4endl;
  }
  // Search for the process
  currentName = processName;
  currentModel = 0;
  size_t idx   = 0;
  G4LossTableManager* lManager = G4LossTableManager::Instance();
  const std::vector<G4VEnergyLossProcess*> vel = lManager->GetEnergyLossProcessVector();
  G4int n = vel.size();
  for(G4int i=0; i<n; i++) {
    if((vel[i])->GetProcessName() == currentName && (vel[i])->Particle() == p) {
      const G4ParticleDefinition* bp = (vel[i])->BaseParticle();
      //      G4cout << "i= " << i << " bp= " << bp << G4endl;
      if(!bp) {
        currentModel = (vel[i])->SelectModelForMaterial(kinEnergy, idx);
      } else {
        G4double e = kinEnergy*bp->GetPDGMass()/p->GetPDGMass();
        for(G4int j=0; j<n; j++) {
          if((vel[j])->Particle() == bp) {
            currentModel = (vel[j])->SelectModelForMaterial(e, idx);
            break;
	  }
	}
      }
      break;
    }
  }
  if(!currentModel) {
    const std::vector<G4VEmProcess*> vem = lManager->GetEmProcessVector();
    G4int n = vem.size();
    for(G4int i=0; i<n; i++) {
      if((vem[i])->GetProcessName() == currentName && (vem[i])->Particle() == p) {
        currentModel = (vem[i])->SelectModelForMaterial(kinEnergy, idx);
	isApplicable    = true;
        break;
      }
    }
  }
  if(!currentModel) {
    const std::vector<G4VMultipleScattering*> vmsc = lManager->GetMultipleScatteringVector();
    G4int n = vmsc.size();
    for(G4int i=0; i<n; i++) {
      if((vmsc[i])->GetProcessName() == currentName && (vmsc[i])->Particle() == p) {
        currentModel = (vmsc[i])->SelectModelForMaterial(kinEnergy, idx);
	isApplicable    = true;
        break;
      }
    }
  }
  if(currentModel) res = true;
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4VEnergyLossProcess* G4EmCalculator::FindEnergyLossProcess(
                      const G4ParticleDefinition* p)
{
  const G4VEnergyLossProcess* elp = 0;
  
  G4LossTableManager* lManager = G4LossTableManager::Instance();
  const std::vector<G4VEnergyLossProcess*> vel = lManager->GetEnergyLossProcessVector();
  G4int n = vel.size();
  for(G4int i=0; i<n; i++) {
    if((vel[i])->Particle() == p) {
      elp = vel[i];
      break;
    }
  }
  return elp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmCalculator::SetVerbose(G4int verb)
{
  verbose = verb;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

