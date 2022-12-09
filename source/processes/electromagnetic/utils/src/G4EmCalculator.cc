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
//
// Class Description: V.Ivanchenko & M.Novak
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4EmCalculator.hh"
#include "G4SystemOfUnits.hh"
#include "G4LossTableManager.hh"
#include "G4EmParameters.hh"
#include "G4NistManager.hh"
#include "G4DynamicParticle.hh"
#include "G4VEmProcess.hh"
#include "G4VEnergyLossProcess.hh"
#include "G4VMultipleScattering.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4PhysicsTable.hh"
#include "G4ProductionCutsTable.hh"
#include "G4ProcessManager.hh"
#include "G4ionEffectiveCharge.hh"
#include "G4RegionStore.hh"
#include "G4Element.hh"
#include "G4EmCorrections.hh"
#include "G4GenericIon.hh"
#include "G4ProcessVector.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4EmUtility.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmCalculator::G4EmCalculator()
{
  manager = G4LossTableManager::Instance();
  nist    = G4NistManager::Instance();
  theParameters = G4EmParameters::Instance();
  corr    = manager->EmCorrections();
  cutenergy[0] = cutenergy[1] = cutenergy[2] = DBL_MAX;
  theGenericIon = G4GenericIon::GenericIon();
  ionEffCharge  = new G4ionEffectiveCharge();
  dynParticle   = new G4DynamicParticle();
  ionTable      = G4ParticleTable::GetParticleTable()->GetIonTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmCalculator::~G4EmCalculator()
{
  delete ionEffCharge;
  delete dynParticle;
  for (G4int i=0; i<nLocalMaterials; ++i) {
    delete localCouples[i];
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::GetDEDX(G4double kinEnergy, 
                                 const G4ParticleDefinition* p,
                                 const G4Material* mat, 
                                 const G4Region* region)
{
  G4double res = 0.0;
  const G4MaterialCutsCouple* couple = FindCouple(mat, region);
  if(nullptr != couple && UpdateParticle(p, kinEnergy) ) {
    res = manager->GetDEDX(p, kinEnergy, couple);
    
    if(isIon) {
      if(FindEmModel(p, currentProcessName, kinEnergy)) {
        G4double length = CLHEP::nm;
        G4double eloss = res*length;
        //G4cout << "### GetDEDX: E= " << kinEnergy << " dedx0= " << res 
        //       << " de= " << eloss << G4endl;; 
        dynParticle->SetKineticEnergy(kinEnergy);
        currentModel->GetChargeSquareRatio(p, mat, kinEnergy);
        currentModel->CorrectionsAlongStep(couple,dynParticle,length,eloss);
        res = eloss/length; 
             //G4cout << " de1= " << eloss << " res1= " << res 
        //       << " " << p->GetParticleName() <<G4endl;;
      }
    } 
    
    if(verbose>0) {
      G4cout << "G4EmCalculator::GetDEDX: E(MeV)= " << kinEnergy/MeV
             << " DEDX(MeV/mm)= " << res*mm/MeV
             << " DEDX(MeV*cm^2/g)= " << res*gram/(MeV*cm2*mat->GetDensity())
             << "  " <<  p->GetParticleName()
             << " in " <<  mat->GetName()
             << " isIon= " << isIon
             << G4endl;
    }
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::GetRangeFromRestricteDEDX(G4double kinEnergy, 
                                                   const G4ParticleDefinition* p,
                                                   const G4Material* mat,
                                                   const G4Region* region)
{
  G4double res = 0.0;
  const G4MaterialCutsCouple* couple = FindCouple(mat,region);
  if(couple && UpdateParticle(p, kinEnergy)) {
    res = manager->GetRangeFromRestricteDEDX(p, kinEnergy, couple);
    if(verbose>1) {
      G4cout << " G4EmCalculator::GetRangeFromRestrictedDEDX: E(MeV)= " 
	     << kinEnergy/MeV
             << " range(mm)= " << res/mm
             << "  " <<  p->GetParticleName()
             << " in " <<  mat->GetName()
             << G4endl;
    }
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::GetCSDARange(G4double kinEnergy, 
                                      const G4ParticleDefinition* p,
                                      const G4Material* mat, 
                                      const G4Region* region)
{
  G4double res = 0.0;
  if(!theParameters->BuildCSDARange()) {
    G4ExceptionDescription ed;
    ed << "G4EmCalculator::GetCSDARange: CSDA table is not built; " 
       << " use UI command: /process/eLoss/CSDARange true";
    G4Exception("G4EmCalculator::GetCSDARange", "em0077",
                JustWarning, ed);
    return res;
  }

  const G4MaterialCutsCouple* couple = FindCouple(mat,region);
  if(nullptr != couple && UpdateParticle(p, kinEnergy)) {
    res = manager->GetCSDARange(p, kinEnergy, couple);
    if(verbose>1) {
      G4cout << " G4EmCalculator::GetCSDARange: E(MeV)= " << kinEnergy/MeV
             << " range(mm)= " << res/mm
             << "  " <<  p->GetParticleName()
             << " in " <<  mat->GetName()
             << G4endl;
    }
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::GetRange(G4double kinEnergy, 
                                  const G4ParticleDefinition* p,
                                  const G4Material* mat, 
                                  const G4Region* region)
{
  G4double res = 0.0;
  if(theParameters->BuildCSDARange()) {
    res = GetCSDARange(kinEnergy, p, mat, region);
  } else {
    res = GetRangeFromRestricteDEDX(kinEnergy, p, mat, region);
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::GetKinEnergy(G4double range, 
                                      const G4ParticleDefinition* p,
                                      const G4Material* mat,
                                      const G4Region* region)
{
  G4double res = 0.0;
  const G4MaterialCutsCouple* couple = FindCouple(mat,region);
  if(nullptr != couple && UpdateParticle(p, 1.0*GeV)) {
    res = manager->GetEnergy(p, range, couple);
    if(verbose>0) {
      G4cout << "G4EmCalculator::GetKinEnergy: Range(mm)= " << range/mm
             << " KinE(MeV)= " << res/MeV
             << "  " <<  p->GetParticleName()
             << " in " <<  mat->GetName()
             << G4endl;
    }
  }
  return res;
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

  if(nullptr != couple && UpdateParticle(p, kinEnergy)) {
    if(FindEmModel(p, processName, kinEnergy)) {
      G4int idx      = couple->GetIndex();
      G4int procType = -1;
      FindLambdaTable(p, processName, kinEnergy, procType);

      G4VEmProcess* emproc = FindDiscreteProcess(p, processName);
      if(nullptr != emproc) {
	res = emproc->GetCrossSection(kinEnergy, couple);
      } else if(currentLambda) {
        // special tables are built for Msc models
	// procType is set in FindLambdaTable
        if(procType==2) {
          auto mscM = static_cast<G4VMscModel*>(currentModel);
          mscM->SetCurrentCouple(couple);
          G4double tr1Mfp = mscM->GetTransportMeanFreePath(p, kinEnergy);
          if (tr1Mfp<DBL_MAX) {
            res = 1./tr1Mfp;
          }
        } else {  
          G4double e = kinEnergy*massRatio;
          res = (((*currentLambda)[idx])->Value(e))*chargeSquare;
        } 
      } else {
      	res = ComputeCrossSectionPerVolume(kinEnergy, p, processName, mat, kinEnergy);
      }
      if(verbose>0) {
      	G4cout << "G4EmCalculator::GetXSPerVolume: E(MeV)= " << kinEnergy/MeV
      	       << " cross(cm-1)= " << res*cm
      	       << "  " <<  p->GetParticleName()
      	       << " in " <<  mat->GetName();
      	if(verbose>1) 
      	  G4cout << "  idx= " << idx << "  Escaled((MeV)= " 
      		 << kinEnergy*massRatio 
      		 << "  q2= " << chargeSquare; 
      	G4cout << G4endl;
      } 
    }
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::GetShellIonisationCrossSectionPerAtom(
                                         const G4String& particle, 
                                         G4int Z, 
                                         G4AtomicShellEnumerator shell,
                                         G4double kinEnergy)
{
  G4double res = 0.0;
  const G4ParticleDefinition* p = FindParticle(particle);
  G4VAtomDeexcitation* ad = manager->AtomDeexcitation();
  if(nullptr != p && nullptr != ad) { 
    res = ad->GetShellIonisationCrossSectionPerAtom(p, Z, shell, kinEnergy); 
  }
  return res;
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
  if(x > 0.0) { res = 1.0/x; }
  if(verbose>1) {
    G4cout << "G4EmCalculator::GetMeanFreePath: E(MeV)= " << kinEnergy/MeV
           << " MFP(mm)= " << res/mm
           << "  " <<  p->GetParticleName()
           << " in " <<  mat->GetName()
           << G4endl;
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmCalculator::PrintDEDXTable(const G4ParticleDefinition* p)
{
  const G4VEnergyLossProcess* elp = manager->GetEnergyLossProcess(p);
  G4cout << "##### DEDX Table for " << p->GetParticleName() << G4endl;
  if(nullptr != elp) G4cout << *(elp->DEDXTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmCalculator::PrintRangeTable(const G4ParticleDefinition* p)
{
  const G4VEnergyLossProcess* elp = manager->GetEnergyLossProcess(p);
  G4cout << "##### Range Table for " << p->GetParticleName() << G4endl;
  if(nullptr != elp) G4cout << *(elp->RangeTableForLoss()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmCalculator::PrintInverseRangeTable(const G4ParticleDefinition* p)
{
  const G4VEnergyLossProcess* elp = manager->GetEnergyLossProcess(p);
  G4cout << "### G4EmCalculator: Inverse Range Table for " 
         << p->GetParticleName() << G4endl;
  if(nullptr != elp) G4cout << *(elp->InverseRangeTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::ComputeDEDX(G4double kinEnergy,
                                     const G4ParticleDefinition* p,
                                     const G4String& processName,
                                     const G4Material* mat,
                                           G4double cut)
{
  SetupMaterial(mat);
  G4double res = 0.0;
  if(verbose > 1) {
    G4cout << "### G4EmCalculator::ComputeDEDX: " << p->GetParticleName()
           << " in " << currentMaterialName
           << " e(MeV)= " << kinEnergy/MeV << "  cut(MeV)= " << cut/MeV
           << G4endl;
  }
  if(UpdateParticle(p, kinEnergy)) {
    if(FindEmModel(p, processName, kinEnergy)) {
      G4double escaled = kinEnergy*massRatio;
      if(nullptr != baseParticle) {
	res = currentModel->ComputeDEDXPerVolume(mat, baseParticle,
                                                 escaled, cut) * chargeSquare;
	if(verbose > 1) {
	  G4cout << "Particle: " << p->GetParticleName()
		 << " E(MeV)=" << kinEnergy
		 << " Base particle: " << baseParticle->GetParticleName()
		 << " Escaled(MeV)= " << escaled 
		 << " q2=" << chargeSquare << G4endl;
	} 
      } else {
	res = currentModel->ComputeDEDXPerVolume(mat, p, kinEnergy, cut);
	if(verbose > 1) { 
	  G4cout << "Particle: " << p->GetParticleName()
		 << " E(MeV)=" << kinEnergy << G4endl;
	}
      }
      if(verbose > 1) {
	G4cout << currentModel->GetName() << ": DEDX(MeV/mm)= " << res*mm/MeV
	       << " DEDX(MeV*cm^2/g)= "
	       << res*gram/(MeV*cm2*mat->GetDensity())
	       << G4endl;
      }
      // emulate smoothing procedure
      if(applySmoothing && nullptr != loweModel) {
	G4double eth = currentModel->LowEnergyLimit();
	G4double res0 = 0.0;
	G4double res1 = 0.0;
	if(nullptr != baseParticle) {
	  res1 = chargeSquare*
	    currentModel->ComputeDEDXPerVolume(mat, baseParticle, eth, cut);
	  res0 = chargeSquare*
	    loweModel->ComputeDEDXPerVolume(mat, baseParticle, eth, cut);
	} else {
	  res1 = currentModel->ComputeDEDXPerVolume(mat, p, eth, cut);
	  res0 = loweModel->ComputeDEDXPerVolume(mat, p, eth, cut);
	}
	if(res1 > 0.0 && escaled > 0.0) {
	  res *= (1.0 + (res0/res1 - 1.0)*eth/escaled);
	}  
	if(verbose > 1) {
	  G4cout << "At boundary energy(MeV)= " << eth/MeV
		 << " DEDX(MeV/mm)= " << res0*mm/MeV << "  " << res1*mm/MeV
		 << " after correction DEDX(MeV/mm)=" << res*mm/MeV << G4endl;
	}
      } 
      // correction for ions
      if(isIon) {
	const G4double length = CLHEP::nm;
	if(UpdateCouple(mat, cut)) {
	  G4double eloss = res*length;
	  dynParticle->SetKineticEnergy(kinEnergy);
	  currentModel->CorrectionsAlongStep(currentCouple,dynParticle,
                                             length,eloss);
	  res = eloss/length; 
        
	  if(verbose > 1) {
	    G4cout << "After Corrections: DEDX(MeV/mm)= " << res*mm/MeV
		   << " DEDX(MeV*cm^2/g)= "
		   << res*gram/(MeV*cm2*mat->GetDensity()) << G4endl;
	  } 
	}
      }
      if(verbose > 0) {
	G4cout << "## E(MeV)= " << kinEnergy/MeV
	       << " DEDX(MeV/mm)= " << res*mm/MeV
	       << " DEDX(MeV*cm^2/g)= " << res*gram/(MeV*cm2*mat->GetDensity())
	       << " cut(MeV)= " << cut/MeV
	       << "  " <<  p->GetParticleName()
	       << " in " <<  currentMaterialName
	       << " Zi^2= " << chargeSquare
	       << " isIon=" << isIon
	       << G4endl;
      }
    }
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::ComputeElectronicDEDX(G4double kinEnergy,
                                               const G4ParticleDefinition* part,
                                               const G4Material* mat,
                                               G4double cut)
{
  SetupMaterial(mat);
  G4double dedx = 0.0;
  if(UpdateParticle(part, kinEnergy)) {

    G4LossTableManager* lManager = G4LossTableManager::Instance();
    const std::vector<G4VEnergyLossProcess*> vel =
      lManager->GetEnergyLossProcessVector();
    std::size_t n = vel.size();

    //G4cout << "ComputeElectronicDEDX for " << part->GetParticleName() 
    //           << " n= " << n << G4endl;
 
    for(std::size_t i=0; i<n; ++i) {
      if(vel[i]) {
        auto p = static_cast<G4VProcess*>(vel[i]);
        if(ActiveForParticle(part, p)) {
          //G4cout << "idx= " << i << " " << (vel[i])->GetProcessName()
          //  << "  " << (vel[i])->Particle()->GetParticleName() << G4endl; 
          dedx += ComputeDEDX(kinEnergy,part,(vel[i])->GetProcessName(),mat,cut);
        }
      }
    }
  }
  return dedx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4EmCalculator::ComputeDEDXForCutInRange(G4double kinEnergy,
                                         const G4ParticleDefinition* part,
                                         const G4Material* mat,
                                         G4double rangecut)
{
  SetupMaterial(mat);
  G4double dedx = 0.0;
  if(UpdateParticle(part, kinEnergy)) {

    G4LossTableManager* lManager = G4LossTableManager::Instance();
    const std::vector<G4VEnergyLossProcess*> vel =
      lManager->GetEnergyLossProcessVector();
    std::size_t n = vel.size();

    if(mat != cutMaterial) {
      cutMaterial = mat;
      cutenergy[0] =
        ComputeEnergyCutFromRangeCut(rangecut, G4Gamma::Gamma(), mat);
      cutenergy[1] =
        ComputeEnergyCutFromRangeCut(rangecut, G4Electron::Electron(), mat);
      cutenergy[2] =
        ComputeEnergyCutFromRangeCut(rangecut, G4Positron::Positron(), mat);
    }

    //G4cout << "ComputeElectronicDEDX for " << part->GetParticleName() 
    //           << " n= " << n << G4endl;
 
    for(std::size_t i=0; i<n; ++i) {
      if(vel[i]) {
        auto p = static_cast<G4VProcess*>(vel[i]);
        if(ActiveForParticle(part, p)) {
          //G4cout << "idx= " << i << " " << (vel[i])->GetProcessName()
          // << "  " << (vel[i])->Particle()->GetParticleName() << G4endl; 
          const G4ParticleDefinition* sec = (vel[i])->SecondaryParticle();
          std::size_t idx = 0;
          if(sec == G4Electron::Electron()) { idx = 1; }
          else if(sec == G4Positron::Positron()) { idx = 2; }

          dedx += ComputeDEDX(kinEnergy,part,(vel[i])->GetProcessName(),
                              mat,cutenergy[idx]);
        }
      }
    }
  }
  return dedx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::ComputeTotalDEDX(G4double kinEnergy, 
                                          const G4ParticleDefinition* part,
                                          const G4Material* mat, 
                                          G4double cut)
{
  G4double dedx = ComputeElectronicDEDX(kinEnergy,part,mat,cut);
  if(mass > 700.*MeV) { dedx += ComputeNuclearDEDX(kinEnergy,part,mat); }
  return dedx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::ComputeNuclearDEDX(G4double kinEnergy,
                                      const G4ParticleDefinition* p,
                                      const G4Material* mat)
{
  G4double res = 0.0;
  G4VEmProcess* nucst = FindDiscreteProcess(p, "nuclearStopping");
  if(nucst) {
    G4VEmModel* mod = nucst->EmModel();
    if(mod) {
      mod->SetFluctuationFlag(false);
      res = mod->ComputeDEDXPerVolume(mat, p, kinEnergy);
    }
  }  

  if(verbose > 1) {
    G4cout <<  p->GetParticleName() << " E(MeV)= " << kinEnergy/MeV
           << " NuclearDEDX(MeV/mm)= " << res*mm/MeV
           << " NuclearDEDX(MeV*cm^2/g)= "
           << res*gram/(MeV*cm2*mat->GetDensity())
           << G4endl;
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::ComputeCrossSectionPerVolume(
                                                   G4double kinEnergy,
                                             const G4ParticleDefinition* p,
                                             const G4String& processName,
                                             const G4Material* mat,
                                                   G4double cut)
{
  SetupMaterial(mat);
  G4double res = 0.0;
  if(UpdateParticle(p, kinEnergy)) {
    if(FindEmModel(p, processName, kinEnergy)) {
      G4double e = kinEnergy;
      G4double aCut = std::max(cut, theParameters->LowestElectronEnergy()); 
      if(baseParticle) {
        e *= kinEnergy*massRatio;
        res = currentModel->CrossSectionPerVolume(
              mat, baseParticle, e, aCut, e) * chargeSquare;
      } else {
        res = currentModel->CrossSectionPerVolume(mat, p, e, aCut, e);
      }
      if(verbose>0) {
        G4cout << "G4EmCalculator::ComputeXSPerVolume: E(MeV)= "
               << kinEnergy/MeV
               << " cross(cm-1)= " << res*cm
               << " cut(keV)= " << aCut/keV
               << "  " <<  p->GetParticleName()
               << " in " <<  mat->GetName()
               << G4endl;
      }
    }
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4EmCalculator::ComputeCrossSectionPerAtom(G4double kinEnergy,
                                           const G4ParticleDefinition* p,
                                           const G4String& processName,
                                           G4double Z, G4double A,
                                           G4double cut)
{
  G4double res = 0.0;
  if(UpdateParticle(p, kinEnergy)) {
    G4int iz = G4lrint(Z);
    CheckMaterial(iz);
    if(FindEmModel(p, processName, kinEnergy)) {
      G4double e = kinEnergy;
      G4double aCut = std::max(cut, theParameters->LowestElectronEnergy()); 
      if(baseParticle) {
        e *= kinEnergy*massRatio;
        currentModel->InitialiseForElement(baseParticle, iz);
        res = currentModel->ComputeCrossSectionPerAtom(
              baseParticle, e, Z, A, aCut) * chargeSquare;
      } else {
        currentModel->InitialiseForElement(p, iz);
        res = currentModel->ComputeCrossSectionPerAtom(p, e, Z, A, aCut);
      }
      if(verbose>0) {
        G4cout << "E(MeV)= " << kinEnergy/MeV
               << " cross(barn)= " << res/barn
               << "  " <<  p->GetParticleName()
               << " Z= " <<  Z << " A= " << A/(g/mole) << " g/mole"
               << " cut(keV)= " << aCut/keV
               << G4endl;
      }
    }
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double
G4EmCalculator::ComputeCrossSectionPerShell(G4double kinEnergy,
                                            const G4ParticleDefinition* p,
                                            const G4String& processName,
                                            G4int Z, G4int shellIdx,
                                            G4double cut)
{
  G4double res = 0.0;
  if(UpdateParticle(p, kinEnergy)) {
    CheckMaterial(Z);
    if(FindEmModel(p, processName, kinEnergy)) {
      G4double e = kinEnergy;
      G4double aCut = std::max(cut, theParameters->LowestElectronEnergy()); 
      if(nullptr != baseParticle) {
        e *= kinEnergy*massRatio;
        currentModel->InitialiseForElement(baseParticle, Z);
        res =
          currentModel->ComputeCrossSectionPerShell(baseParticle, Z, shellIdx, 
                                                    e, aCut) * chargeSquare;
      } else {
        currentModel->InitialiseForElement(p, Z);
        res = currentModel->ComputeCrossSectionPerAtom(p, Z, shellIdx, e, aCut);
      }
      if(verbose>0) {
        G4cout << "E(MeV)= " << kinEnergy/MeV
               << " cross(barn)= " << res/barn
               << "  " <<  p->GetParticleName()
               << " Z= " <<  Z << " shellIdx= " << shellIdx 
               << " cut(keV)= " << aCut/keV
	       << G4endl;
      }
    }
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4EmCalculator::ComputeGammaAttenuationLength(G4double kinEnergy, 
                                              const G4Material* mat)
{
  G4double res = 0.0;
  const G4ParticleDefinition* gamma = G4Gamma::Gamma();
  res += ComputeCrossSectionPerVolume(kinEnergy, gamma, "conv", mat, 0.0);
  res += ComputeCrossSectionPerVolume(kinEnergy, gamma, "compt", mat, 0.0);
  res += ComputeCrossSectionPerVolume(kinEnergy, gamma, "phot", mat, 0.0);
  res += ComputeCrossSectionPerVolume(kinEnergy, gamma, "Rayl", mat, 0.0);
  if(res > 0.0) { res = 1.0/res; }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::ComputeShellIonisationCrossSectionPerAtom(
                                         const G4String& particle, 
                                         G4int Z, 
                                         G4AtomicShellEnumerator shell,
                                         G4double kinEnergy,
                                         const G4Material* mat)
{
  G4double res = 0.0;
  const G4ParticleDefinition* p = FindParticle(particle);
  G4VAtomDeexcitation* ad = manager->AtomDeexcitation();
  if(p && ad) { 
    res = ad->ComputeShellIonisationCrossSectionPerAtom(p, Z, shell, 
                                                        kinEnergy, mat); 
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::ComputeMeanFreePath(G4double kinEnergy,
                                             const G4ParticleDefinition* p,
                                             const G4String& processName,
                                             const G4Material* mat,
                                             G4double cut)
{
  G4double mfp = DBL_MAX;
  G4double x =
    ComputeCrossSectionPerVolume(kinEnergy, p, processName, mat, cut);
  if(x > 0.0) { mfp = 1.0/x; }
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

G4double G4EmCalculator::ComputeEnergyCutFromRangeCut(
                         G4double range, 
                         const G4ParticleDefinition* part,
                         const G4Material* mat)
{
  return G4ProductionCutsTable::GetProductionCutsTable()->
    ConvertRangeToEnergy(part, mat, range);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4EmCalculator::UpdateParticle(const G4ParticleDefinition* p,
                                      G4double kinEnergy)
{
  if(p != currentParticle) {

    // new particle
    currentParticle = p;
    dynParticle->SetDefinition(const_cast<G4ParticleDefinition*>(p));
    dynParticle->SetKineticEnergy(kinEnergy);
    baseParticle    = nullptr;
    currentParticleName = p->GetParticleName();
    massRatio       = 1.0;
    mass            = p->GetPDGMass();
    chargeSquare    = 1.0;
    currentProcess  = manager->GetEnergyLossProcess(p);
    currentProcessName = "";
    isIon = false;

    // ionisation process exist
    if(nullptr != currentProcess) {
      currentProcessName = currentProcess->GetProcessName();
      baseParticle = currentProcess->BaseParticle();
      if(currentProcessName == "ionIoni" && p->GetParticleName() != "alpha") {
        baseParticle = theGenericIon;
        isIon = true;
      }

      // base particle is used
      if(nullptr != baseParticle) {
        massRatio = baseParticle->GetPDGMass()/p->GetPDGMass();
        G4double q = p->GetPDGCharge()/baseParticle->GetPDGCharge();
        chargeSquare = q*q;
      } 
    }
  }
  // Effective charge for ions
  if(isIon && nullptr != currentProcess) {
    chargeSquare =
      corr->EffectiveChargeSquareRatio(p, currentMaterial, kinEnergy);
    currentProcess->SetDynamicMassCharge(massRatio,chargeSquare);
    if(verbose>1) {
      G4cout <<"\n NewIon: massR= "<< massRatio << "   q2= " 
	     << chargeSquare << "  " << currentProcess << G4endl;
    }
  }
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4ParticleDefinition* G4EmCalculator::FindParticle(const G4String& name)
{
  const G4ParticleDefinition* p = nullptr;
  if(name != currentParticleName) {
    p = G4ParticleTable::GetParticleTable()->FindParticle(name);
    if(nullptr == p) {
      G4cout << "### WARNING: G4EmCalculator::FindParticle fails to find " 
             << name << G4endl;
    }
  } else {
    p = currentParticle;
  }
  return p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4ParticleDefinition* G4EmCalculator::FindIon(G4int Z, G4int A)
{
  const G4ParticleDefinition* p = ionTable->GetIon(Z,A,0);
  return p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4Material* G4EmCalculator::FindMaterial(const G4String& name)
{
  if(name != currentMaterialName) {
    SetupMaterial(G4Material::GetMaterial(name, false));
    if(nullptr == currentMaterial) {
      G4cout << "### WARNING: G4EmCalculator::FindMaterial fails to find " 
             << name << G4endl;
    }
  }
  return currentMaterial;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4Region* G4EmCalculator::FindRegion(const G4String& reg)
{
  return G4EmUtility::FindRegion(reg);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4MaterialCutsCouple* G4EmCalculator::FindCouple(
                            const G4Material* material,
                            const G4Region* region)
{
  const G4MaterialCutsCouple* couple = nullptr;
  SetupMaterial(material);
  if(nullptr != currentMaterial) {
    // Access to materials
    const G4ProductionCutsTable* theCoupleTable=
      G4ProductionCutsTable::GetProductionCutsTable();
    const G4Region* r = region;
    if(nullptr != r) {
      couple = theCoupleTable->GetMaterialCutsCouple(material,
                                                     r->GetProductionCuts());
    } else {
      G4RegionStore* store = G4RegionStore::GetInstance();
      std::size_t nr = store->size();
      if(0 < nr) {
        for(std::size_t i=0; i<nr; ++i) {
          couple = theCoupleTable->GetMaterialCutsCouple(
            material, ((*store)[i])->GetProductionCuts());
          if(nullptr != couple) { break; }
        }
      }
    }
  }
  if(nullptr == couple) {
    G4ExceptionDescription ed;
    ed << "G4EmCalculator::FindCouple: fail for material <" 
       << currentMaterialName << ">";
    if(region) { ed << " and region " << region->GetName(); }
    G4Exception("G4EmCalculator::FindCouple", "em0078",
                FatalException, ed);
  }
  return couple;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4EmCalculator::UpdateCouple(const G4Material* material, G4double cut)
{
  SetupMaterial(material);
  if(!currentMaterial) { return false; }
  for (G4int i=0; i<nLocalMaterials; ++i) {
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
  ++nLocalMaterials;
  currentCouple = cc;
  currentCoupleIndex = currentCouple->GetIndex();
  currentCut = cut;
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmCalculator::FindLambdaTable(const G4ParticleDefinition* p,
                                     const G4String& processName,
                                     G4double kinEnergy, G4int& proctype)
{
  // Search for the process
  if (!currentLambda || p != lambdaParticle || processName != lambdaName) {
    lambdaName     = processName;
    currentLambda  = nullptr;
    lambdaParticle = p;
    isApplicable   = false;

    const G4ParticleDefinition* part = (isIon) ? theGenericIon : p;

    // Search for energy loss process
    currentName = processName;
    currentModel = nullptr;
    loweModel = nullptr;

    G4VEnergyLossProcess* elproc = FindEnLossProcess(part, processName);
    if(nullptr != elproc) {
      currentLambda = elproc->LambdaTable();
      proctype      = 0;
      if(nullptr != currentLambda) {
        isApplicable = true;
        if(verbose>1) { 
          G4cout << "G4VEnergyLossProcess is found out: " << currentName 
                 << G4endl;
        }
      }
      curProcess = elproc;
      return;
    }

    // Search for discrete process 
    G4VEmProcess* proc = FindDiscreteProcess(part, processName);
    if(nullptr != proc) {
      currentLambda = proc->LambdaTable();
      proctype      = 1;
      if(nullptr != currentLambda) {
        isApplicable = true;
        if(verbose>1) { 
          G4cout << "G4VEmProcess is found out: " << currentName << G4endl;
        }
      }
      curProcess = proc;
      return;
    }

    // Search for msc process
    G4VMultipleScattering* msc = FindMscProcess(part, processName);
    if(nullptr != msc) {
      currentModel = msc->SelectModel(kinEnergy,0);
      proctype     = 2;
      if(nullptr != currentModel) {
        currentLambda = currentModel->GetCrossSectionTable();
        if(nullptr != currentLambda) {
          isApplicable = true;
          if(verbose>1) { 
            G4cout << "G4VMultipleScattering is found out: " << currentName 
                   << G4endl;
          }
        }
      }
      curProcess = msc;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4EmCalculator::FindEmModel(const G4ParticleDefinition* p,
                                   const G4String& processName,
                                            G4double kinEnergy)
{
  isApplicable = false;
  if(nullptr == p || nullptr == currentMaterial) {
    G4cout << "G4EmCalculator::FindEmModel WARNING: no particle" 
           << " or materail defined; particle: " << p << G4endl;
    return isApplicable;
  }
  G4String partname =  p->GetParticleName();
  G4double scaledEnergy = kinEnergy*massRatio;
  const G4ParticleDefinition* part = (isIon) ? theGenericIon : p; 

  if(verbose > 1) {
    G4cout << "## G4EmCalculator::FindEmModel for " << partname
           << " (type= " << p->GetParticleType()
           << ") and " << processName << " at E(MeV)= " << scaledEnergy 
           << G4endl;
    if(p != part) { G4cout << "  GenericIon is the base particle" << G4endl; }
  }

  // Search for energy loss process
  currentName = processName;
  currentModel = nullptr;
  loweModel = nullptr;
  std::size_t idx = 0;

  G4VEnergyLossProcess* elproc = FindEnLossProcess(part, processName);
  if(nullptr != elproc) {
    currentModel = elproc->SelectModelForMaterial(scaledEnergy, idx);
    currentModel->InitialiseForMaterial(part, currentMaterial);
    currentModel->SetupForMaterial(part, currentMaterial, scaledEnergy);
    G4double eth = currentModel->LowEnergyLimit();
    if(eth > 0.0) {
      loweModel = elproc->SelectModelForMaterial(eth - CLHEP::eV, idx);
      if(loweModel == currentModel) { loweModel = nullptr; }
      else { 
        loweModel->InitialiseForMaterial(part, currentMaterial);
        loweModel->SetupForMaterial(part, currentMaterial, eth - CLHEP::eV); 
      }
    }
  }

  // Search for discrete process 
  if(nullptr == currentModel) {
    G4VEmProcess* proc = FindDiscreteProcess(part, processName);
    if(nullptr != proc) {
      currentModel = proc->SelectModelForMaterial(kinEnergy, idx);
      currentModel->InitialiseForMaterial(part, currentMaterial);
      currentModel->SetupForMaterial(part, currentMaterial, kinEnergy);
      G4double eth = currentModel->LowEnergyLimit();
      if(eth > 0.0) {
        loweModel = proc->SelectModelForMaterial(eth - CLHEP::eV, idx);
        if(loweModel == currentModel) { loweModel = nullptr; }
        else { 
          loweModel->InitialiseForMaterial(part, currentMaterial);
          loweModel->SetupForMaterial(part, currentMaterial, eth - CLHEP::eV); 
        }
      }
    }
  }

  // Search for msc process
  if(nullptr == currentModel) {
    G4VMultipleScattering* proc = FindMscProcess(part, processName);
    if(nullptr != proc) {
      currentModel = proc->SelectModel(kinEnergy, idx);
      loweModel = nullptr;
    }
  }
  if(nullptr != currentModel) {
    if(loweModel == currentModel) { loweModel = nullptr; }
    isApplicable = true;
    currentModel->InitialiseForMaterial(part, currentMaterial);
    if(loweModel) {
      loweModel->InitialiseForMaterial(part, currentMaterial);
    }
    if(verbose > 1) {
      G4cout << "   Model <" << currentModel->GetName() 
             << "> Emin(MeV)= " << currentModel->LowEnergyLimit()/MeV
             << " for " << part->GetParticleName();
      if(nullptr != elproc) { 
        G4cout << " and " << elproc->GetProcessName() << "  " << elproc 
               << G4endl;
      }
      if(nullptr != loweModel) { 
        G4cout << " LowEnergy model <" << loweModel->GetName() << ">"; 
      }
      G4cout << G4endl;
    } 
  }
  return isApplicable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEnergyLossProcess* 
G4EmCalculator::FindEnLossProcess(const G4ParticleDefinition* part,
                                  const G4String& processName)
{
  G4VEnergyLossProcess* proc = nullptr;
  const std::vector<G4VEnergyLossProcess*> v = 
    manager->GetEnergyLossProcessVector();
  std::size_t n = v.size();
  for(std::size_t i=0; i<n; ++i) {
    if((v[i])->GetProcessName() == processName) {
      auto p = static_cast<G4VProcess*>(v[i]);
      if(ActiveForParticle(part, p)) {
        proc = v[i];
        break;
      }
    }
  }
  return proc;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEmProcess* 
G4EmCalculator::FindDiscreteProcess(const G4ParticleDefinition* part,
                                    const G4String& processName)
{
  G4VEmProcess* proc = nullptr;
  auto v = manager->GetEmProcessVector();
  std::size_t n = v.size();
  for(std::size_t i=0; i<n; ++i) {
    auto pName = v[i]->GetProcessName();
    if(pName == "GammaGeneralProc") {
      proc = v[i]->GetEmProcess(processName);
      break;
    } else if(pName == processName) {
      auto p = static_cast<G4VProcess*>(v[i]);
      if(ActiveForParticle(part, p)) {
        proc = v[i];
        break;
      }
    }
  }
  return proc;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VMultipleScattering* 
G4EmCalculator::FindMscProcess(const G4ParticleDefinition* part,
                               const G4String& processName)
{
  G4VMultipleScattering* proc = nullptr;
  const std::vector<G4VMultipleScattering*> v = 
    manager->GetMultipleScatteringVector();
  std::size_t n = v.size();
  for(std::size_t i=0; i<n; ++i) {
    if((v[i])->GetProcessName() == processName) {
      auto p = static_cast<G4VProcess*>(v[i]);
      if(ActiveForParticle(part, p)) {
        proc = v[i];
        break;
      }
    }
  }
  return proc;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VProcess* G4EmCalculator::FindProcess(const G4ParticleDefinition* part,
                                        const G4String& processName)
{
  G4VProcess* proc = nullptr;
  const G4ProcessManager* procman = part->GetProcessManager();
  G4ProcessVector* pv = procman->GetProcessList();
  G4int nproc = (G4int)pv->size();
  for(G4int i=0; i<nproc; ++i) {
    if(processName == (*pv)[i]->GetProcessName()) {
      proc = (*pv)[i];
      break;
    }
  }
  return proc;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4EmCalculator::ActiveForParticle(const G4ParticleDefinition* part,
                                         G4VProcess* proc)
{
  G4ProcessManager* pm = part->GetProcessManager();
  G4ProcessVector* pv = pm->GetProcessList();
  G4int n = (G4int)pv->size();
  G4bool res = false;
  for(G4int i=0; i<n; ++i) {
    if((*pv)[i] == proc) {
      if(pm->GetProcessActivation(i)) { res = true; }
      break;
    }
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmCalculator::SetupMaterial(const G4Material* mat)
{
  if(mat) {
    currentMaterial = mat;
    currentMaterialName = mat->GetName();
  } else {
    currentMaterial = nullptr;
    currentMaterialName = "";
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmCalculator::SetupMaterial(const G4String& mname)
{
  SetupMaterial(nist->FindOrBuildMaterial(mname));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmCalculator::CheckMaterial(G4int Z)
{
  G4bool isFound = false;
  if(nullptr != currentMaterial) {
    G4int nn = (G4int)currentMaterial->GetNumberOfElements();
    for(G4int i=0; i<nn; ++i) { 
      if(Z == currentMaterial->GetElement(i)->GetZasInt()) {
        isFound = true;
        break;
      }
    }
  }
  if(!isFound) {
    SetupMaterial(nist->FindOrBuildSimpleMaterial(Z));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmCalculator::SetVerbose(G4int verb)
{
  verbose = verb;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

