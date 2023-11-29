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

#include "G4AdjointCSManager.hh"

#include "G4AdjointCSMatrix.hh"
#include "G4AdjointElectron.hh"
#include "G4AdjointGamma.hh"
#include "G4AdjointInterpolator.hh"
#include "G4AdjointProton.hh"
#include "G4Electron.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4Gamma.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsTable.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Proton.hh"
#include "G4SystemOfUnits.hh"
#include "G4VEmAdjointModel.hh"
#include "G4VEmProcess.hh"
#include "G4VEnergyLossProcess.hh"

G4ThreadLocal G4AdjointCSManager* G4AdjointCSManager::fInstance = nullptr;

constexpr G4double G4AdjointCSManager::fTmin;
constexpr G4double G4AdjointCSManager::fTmax;
constexpr G4int G4AdjointCSManager::fNbins;

///////////////////////////////////////////////////////
G4AdjointCSManager* G4AdjointCSManager::GetAdjointCSManager()
{
  if(fInstance == nullptr)
  {
    static G4ThreadLocalSingleton<G4AdjointCSManager> inst;
    fInstance = inst.Instance();
  }
  return fInstance;
}

///////////////////////////////////////////////////////
G4AdjointCSManager::G4AdjointCSManager()
{
  RegisterAdjointParticle(G4AdjointElectron::AdjointElectron());
  RegisterAdjointParticle(G4AdjointGamma::AdjointGamma());
  RegisterAdjointParticle(G4AdjointProton::AdjointProton());
}

///////////////////////////////////////////////////////
G4AdjointCSManager::~G4AdjointCSManager()
{
  for (auto& p : fAdjointCSMatricesForProdToProj) {
    for (auto p1 : p) {
      if (p1) {
        delete p1;
        p1 = nullptr;
      }
    }
    p.clear();
  }
  fAdjointCSMatricesForProdToProj.clear();

  for (auto& p : fAdjointCSMatricesForScatProjToProj) {
    for (auto p1 : p) {
      if (p1) {
        delete p1;
        p1 = nullptr;
      }
    }
    p.clear();
  }
  fAdjointCSMatricesForScatProjToProj.clear();

  for (auto p : fAdjointModels) {
    if (p) {
      delete p;
      p = nullptr;
    }
  }
  fAdjointModels.clear();

  for (auto p : fTotalAdjSigmaTable) {
    p->clearAndDestroy();
    delete p;
    p = nullptr;
  }
  fTotalAdjSigmaTable.clear();

  for (auto p : fSigmaTableForAdjointModelScatProjToProj) {
    p->clearAndDestroy();
    delete p;
    p = nullptr;
  }
  fSigmaTableForAdjointModelScatProjToProj.clear();

  for (auto p : fSigmaTableForAdjointModelProdToProj) {
    p->clearAndDestroy();
    delete p;
    p = nullptr;
  }
  fSigmaTableForAdjointModelProdToProj.clear();

  for (auto p : fTotalFwdSigmaTable) {
    p->clearAndDestroy();
    delete p;
    p = nullptr;
  }
  fTotalFwdSigmaTable.clear();

  for (auto p : fForwardProcesses) {
    delete p;
    p = nullptr;
  }
  fForwardProcesses.clear();

  for (auto p : fForwardLossProcesses) {
    delete p;
    p = nullptr;
  }
  fForwardLossProcesses.clear();
}

///////////////////////////////////////////////////////
std::size_t G4AdjointCSManager::RegisterEmAdjointModel(G4VEmAdjointModel* aModel)
{
  fAdjointModels.push_back(aModel);
  fSigmaTableForAdjointModelScatProjToProj.push_back(new G4PhysicsTable);
  fSigmaTableForAdjointModelProdToProj.push_back(new G4PhysicsTable);
  return fAdjointModels.size() - 1;
}

///////////////////////////////////////////////////////
void G4AdjointCSManager::RegisterEmProcess(G4VEmProcess* aProcess,
                                           G4ParticleDefinition* aFwdPartDef)
{
  G4ParticleDefinition* anAdjPartDef =
    GetAdjointParticleEquivalent(aFwdPartDef);
  if(anAdjPartDef && aProcess)
  {
    RegisterAdjointParticle(anAdjPartDef);

    for(std::size_t i = 0; i < fAdjointParticlesInAction.size(); ++i)
    {
      if(anAdjPartDef->GetParticleName() ==
         fAdjointParticlesInAction[i]->GetParticleName())
    	  fForwardProcesses[i]->push_back(aProcess);
    }
  }
}

///////////////////////////////////////////////////////
void G4AdjointCSManager::RegisterEnergyLossProcess(
  G4VEnergyLossProcess* aProcess, G4ParticleDefinition* aFwdPartDef)
{
  G4ParticleDefinition* anAdjPartDef =
    GetAdjointParticleEquivalent(aFwdPartDef);
  if(anAdjPartDef && aProcess)
  {
    RegisterAdjointParticle(anAdjPartDef);
    for(std::size_t i = 0; i < fAdjointParticlesInAction.size(); ++i)
    {
      if(anAdjPartDef->GetParticleName() ==
         fAdjointParticlesInAction[i]->GetParticleName())
    	                        fForwardLossProcesses[i]->push_back(aProcess);

    }
  }
}

///////////////////////////////////////////////////////
void G4AdjointCSManager::RegisterAdjointParticle(G4ParticleDefinition* aPartDef)
{
  G4bool found = false;
  for(auto p : fAdjointParticlesInAction)
  {
    if(p->GetParticleName() == aPartDef->GetParticleName())
    {
      found = true;
    }
  }
  if(!found)
  {
    fForwardLossProcesses.push_back(new std::vector<G4VEnergyLossProcess*>());
    fTotalFwdSigmaTable.push_back(new G4PhysicsTable);
    fTotalAdjSigmaTable.push_back(new G4PhysicsTable);
    fForwardProcesses.push_back(new std::vector<G4VEmProcess*>());
    fAdjointParticlesInAction.push_back(aPartDef);
    fEminForFwdSigmaTables.push_back(std::vector<G4double>());
    fEminForAdjSigmaTables.push_back(std::vector<G4double>());
    fEkinofFwdSigmaMax.push_back(std::vector<G4double>());
    fEkinofAdjSigmaMax.push_back(std::vector<G4double>());
  }
}

///////////////////////////////////////////////////////
void G4AdjointCSManager::BuildCrossSectionMatrices()
{
  if(fCSMatricesBuilt)
    return;
  // The Tcut and Tmax matrices will be computed probably just once.
  // When Tcut changes, some PhysicsTable will be recomputed
  // for each MaterialCutCouple but not all the matrices.
  // The Tcut defines a lower limit in the energy of the projectile before
  // scattering. In the Projectile to Scattered Projectile case we have
  // 			E_ScatProj<E_Proj-Tcut
  // Therefore in the adjoint case we have
  //			Eproj> E_ScatProj+Tcut
  // This implies that when computing the adjoint CS we should integrate over
  // Epro from E_ScatProj+Tcut to Emax
  // In the Projectile to Secondary case Tcut plays a role only in the fact that
  // Esecond should be greater than Tcut to have the possibility to have any
  // adjoint process.
  // To avoid recomputing matrices for all changes of MaterialCutCouple,
  // we propose to compute the matrices only once for the minimum possible Tcut
  // and then to interpolate the probability for a new Tcut (implemented in
  // G4VAdjointEmModel)

  fAdjointCSMatricesForScatProjToProj.clear();
  fAdjointCSMatricesForProdToProj.clear();
  const G4ElementTable* theElementTable   = G4Element::GetElementTable();
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

  G4cout << "========== Computation of cross section matrices for adjoint "
            "models =========="
         << G4endl;
  for(const auto& aModel : fAdjointModels)
  {
    G4cout << "Build adjoint cross section matrices for " << aModel->GetName()
           << G4endl;
    if(aModel->GetUseMatrix())
    {
      std::vector<G4AdjointCSMatrix*>* aListOfMat1 =
        new std::vector<G4AdjointCSMatrix*>();
      std::vector<G4AdjointCSMatrix*>* aListOfMat2 =
        new std::vector<G4AdjointCSMatrix*>();
      if(aModel->GetUseMatrixPerElement())
      {
        if(aModel->GetUseOnlyOneMatrixForAllElements())
        {
          std::vector<G4AdjointCSMatrix*> two_matrices =
            BuildCrossSectionsModelAndElement(aModel, 1, 1, 80);
          aListOfMat1->push_back(two_matrices[0]);
          aListOfMat2->push_back(two_matrices[1]);
        }
        else
        {
          for(const auto& anElement : *theElementTable)
          {
            G4int Z = G4lrint(anElement->GetZ());
            G4int A = G4lrint(anElement->GetN());
            std::vector<G4AdjointCSMatrix*> two_matrices =
              BuildCrossSectionsModelAndElement(aModel, Z, A, 40);
            aListOfMat1->push_back(two_matrices[0]);
            aListOfMat2->push_back(two_matrices[1]);
          }
        }
      }
      else
      {  // Per material case
        for(const auto& aMaterial : *theMaterialTable)
        {
          std::vector<G4AdjointCSMatrix*> two_matrices =
            BuildCrossSectionsModelAndMaterial(aModel, aMaterial, 40);
          aListOfMat1->push_back(two_matrices[0]);
          aListOfMat2->push_back(two_matrices[1]);
        }
      }
      fAdjointCSMatricesForProdToProj.push_back(*aListOfMat1);
      fAdjointCSMatricesForScatProjToProj.push_back(*aListOfMat2);
      aModel->SetCSMatrices(aListOfMat1, aListOfMat2);
    }
    else
    {
      G4cout << "The model " << aModel->GetName()
             << " does not use cross section matrices" << G4endl;
      std::vector<G4AdjointCSMatrix*> two_empty_matrices;
      fAdjointCSMatricesForProdToProj.push_back(two_empty_matrices);
      fAdjointCSMatricesForScatProjToProj.push_back(two_empty_matrices);
    }
  }
  G4cout << "              All adjoint cross section matrices are computed!"
         << G4endl;
  G4cout
    << "======================================================================"
    << G4endl;

  fCSMatricesBuilt = true;
}

///////////////////////////////////////////////////////
void G4AdjointCSManager::BuildTotalSigmaTables()
{
  if(fSigmaTableBuilt)
    return;

  const G4ProductionCutsTable* theCoupleTable =
    G4ProductionCutsTable::GetProductionCutsTable();

  // Prepare the Sigma table for all AdjointEMModel, will be filled later on
  for(std::size_t i = 0; i < fAdjointModels.size(); ++i)
  {
    fSigmaTableForAdjointModelScatProjToProj[i]->clearAndDestroy();
    fSigmaTableForAdjointModelProdToProj[i]->clearAndDestroy();
    for(std::size_t j = 0; j < theCoupleTable->GetTableSize(); ++j)
    {
      fSigmaTableForAdjointModelScatProjToProj[i]->push_back(
        new G4PhysicsLogVector(fTmin, fTmax, fNbins));
      fSigmaTableForAdjointModelProdToProj[i]->push_back(
        new G4PhysicsLogVector(fTmin, fTmax, fNbins));
    }
  }

  for(std::size_t i = 0; i < fAdjointParticlesInAction.size(); ++i)
  {
    G4ParticleDefinition* thePartDef = fAdjointParticlesInAction[i];
    DefineCurrentParticle(thePartDef);
    fTotalFwdSigmaTable[i]->clearAndDestroy();
    fTotalAdjSigmaTable[i]->clearAndDestroy();
    fEminForFwdSigmaTables[i].clear();
    fEminForAdjSigmaTables[i].clear();
    fEkinofFwdSigmaMax[i].clear();
    fEkinofAdjSigmaMax[i].clear();

    for(std::size_t j = 0; j < theCoupleTable->GetTableSize(); ++j)
    {
      const G4MaterialCutsCouple* couple =
        theCoupleTable->GetMaterialCutsCouple((G4int)j);

      // make first the total fwd CS table for FwdProcess
      G4PhysicsVector* aVector = new G4PhysicsLogVector(fTmin, fTmax, fNbins);
      G4bool Emin_found        = false;
      G4double sigma_max       = 0.;
      G4double e_sigma_max     = 0.;
      for(std::size_t l = 0; l < fNbins; ++l)
      {
        G4double totCS = 0.;
        G4double e     = aVector->Energy(l);
        for(std::size_t k = 0; k < fForwardProcesses[i]->size(); ++k)
        {
          totCS += (*fForwardProcesses[i])[k]->GetCrossSection(e, couple);
        }
        for(std::size_t k = 0; k < fForwardLossProcesses[i]->size(); ++k)
        {
          if(thePartDef == fAdjIon)
          {  // e is considered already as the scaled energy
            std::size_t mat_index = couple->GetIndex();
            G4VEmModel* currentModel =
              (*fForwardLossProcesses[i])[k]->SelectModelForMaterial(e,
                                                                     mat_index);
            G4double chargeSqRatio = currentModel->GetChargeSquareRatio(
              fFwdIon, couple->GetMaterial(), e / fMassRatio);
            (*fForwardLossProcesses[i])[k]->SetDynamicMassCharge(fMassRatio,
                                                                 chargeSqRatio);
          }
          G4double e1 = e / fMassRatio;
          totCS += (*fForwardLossProcesses[i])[k]->GetLambda(e1, couple);
        }
        aVector->PutValue(l, totCS);
        if(totCS > sigma_max)
        {
          sigma_max   = totCS;
          e_sigma_max = e;
        }
        if(totCS > 0 && !Emin_found)
        {
          fEminForFwdSigmaTables[i].push_back(e);
          Emin_found = true;
        }
      }

      fEkinofFwdSigmaMax[i].push_back(e_sigma_max);

      if(!Emin_found)
        fEminForFwdSigmaTables[i].push_back(fTmax);

      fTotalFwdSigmaTable[i]->push_back(aVector);

      Emin_found                = false;
      sigma_max                 = 0;
      e_sigma_max               = 0.;
      G4PhysicsVector* aVector1 = new G4PhysicsLogVector(fTmin, fTmax, fNbins);
      for(std::size_t eindex = 0; eindex < fNbins; ++eindex)
      {
        G4double e     = aVector1->Energy(eindex);
        G4double totCS = ComputeTotalAdjointCS(
          couple, thePartDef,
          e * 0.9999999 / fMassRatio);  // fMassRatio needed for ions
        aVector1->PutValue(eindex, totCS);
        if(totCS > sigma_max)
        {
          sigma_max   = totCS;
          e_sigma_max = e;
        }
        if(totCS > 0 && !Emin_found)
        {
          fEminForAdjSigmaTables[i].push_back(e);
          Emin_found = true;
        }
      }
      fEkinofAdjSigmaMax[i].push_back(e_sigma_max);
      if(!Emin_found)
        fEminForAdjSigmaTables[i].push_back(fTmax);

      fTotalAdjSigmaTable[i]->push_back(aVector1);
    }
  }
  fSigmaTableBuilt = true;
}

///////////////////////////////////////////////////////
G4double G4AdjointCSManager::GetTotalAdjointCS(
  G4ParticleDefinition* aPartDef, G4double Ekin,
  const G4MaterialCutsCouple* aCouple)
{
  DefineCurrentMaterial(aCouple);
  DefineCurrentParticle(aPartDef);
  return (((*fTotalAdjSigmaTable[fCurrentParticleIndex])[fCurrentMatIndex])
            ->Value(Ekin * fMassRatio));
}

///////////////////////////////////////////////////////
G4double G4AdjointCSManager::GetTotalForwardCS(
  G4ParticleDefinition* aPartDef, G4double Ekin,
  const G4MaterialCutsCouple* aCouple)
{
  DefineCurrentMaterial(aCouple);
  DefineCurrentParticle(aPartDef);
  return (((*fTotalFwdSigmaTable[fCurrentParticleIndex])[fCurrentMatIndex])
            ->Value(Ekin * fMassRatio));
}

///////////////////////////////////////////////////////
G4double G4AdjointCSManager::GetAdjointSigma(
  G4double Ekin_nuc, std::size_t index_model, G4bool is_scat_proj_to_proj,
  const G4MaterialCutsCouple* aCouple)
{
  DefineCurrentMaterial(aCouple);
  if(is_scat_proj_to_proj)
    return (((*fSigmaTableForAdjointModelScatProjToProj[index_model])
               [fCurrentMatIndex])->Value(Ekin_nuc));
  else
    return (
      ((*fSigmaTableForAdjointModelProdToProj[index_model])[fCurrentMatIndex])
        ->Value(Ekin_nuc));
}

///////////////////////////////////////////////////////
void G4AdjointCSManager::GetEminForTotalCS(G4ParticleDefinition* aPartDef,
                                           const G4MaterialCutsCouple* aCouple,
                                           G4double& emin_adj,
                                           G4double& emin_fwd)
{
  DefineCurrentMaterial(aCouple);
  DefineCurrentParticle(aPartDef);
  emin_adj = fEminForAdjSigmaTables[fCurrentParticleIndex][fCurrentMatIndex] /
             fMassRatio;
  emin_fwd = fEminForFwdSigmaTables[fCurrentParticleIndex][fCurrentMatIndex] /
             fMassRatio;
}

///////////////////////////////////////////////////////
void G4AdjointCSManager::GetMaxFwdTotalCS(G4ParticleDefinition* aPartDef,
                                          const G4MaterialCutsCouple* aCouple,
                                          G4double& e_sigma_max,
                                          G4double& sigma_max)
{
  DefineCurrentMaterial(aCouple);
  DefineCurrentParticle(aPartDef);
  e_sigma_max = fEkinofFwdSigmaMax[fCurrentParticleIndex][fCurrentMatIndex];
  sigma_max = ((*fTotalFwdSigmaTable[fCurrentParticleIndex])[fCurrentMatIndex])
                ->Value(e_sigma_max);
  e_sigma_max /= fMassRatio;
}

///////////////////////////////////////////////////////
void G4AdjointCSManager::GetMaxAdjTotalCS(G4ParticleDefinition* aPartDef,
                                          const G4MaterialCutsCouple* aCouple,
                                          G4double& e_sigma_max,
                                          G4double& sigma_max)
{
  DefineCurrentMaterial(aCouple);
  DefineCurrentParticle(aPartDef);
  e_sigma_max = fEkinofAdjSigmaMax[fCurrentParticleIndex][fCurrentMatIndex];
  sigma_max = ((*fTotalAdjSigmaTable[fCurrentParticleIndex])[fCurrentMatIndex])
                ->Value(e_sigma_max);
  e_sigma_max /= fMassRatio;
}

///////////////////////////////////////////////////////
G4double G4AdjointCSManager::GetCrossSectionCorrection(
  G4ParticleDefinition* aPartDef, G4double PreStepEkin,
  const G4MaterialCutsCouple* aCouple, G4bool& fwd_is_used)
{
  static G4double lastEkin = 0.;
  static G4ParticleDefinition* lastPartDef;

  G4double corr_fac = 1.;
  if(fForwardCSMode && aPartDef)
  {
    if(lastEkin != PreStepEkin || aPartDef != lastPartDef ||
       aCouple != fCurrentCouple)
    {
      DefineCurrentMaterial(aCouple);
      G4double preadjCS = GetTotalAdjointCS(aPartDef, PreStepEkin, aCouple);
      G4double prefwdCS = GetTotalForwardCS(aPartDef, PreStepEkin, aCouple);
      lastEkin          = PreStepEkin;
      lastPartDef       = aPartDef;
      if(prefwdCS > 0. && preadjCS > 0.)
      {
        fForwardCSUsed          = true;
        fLastCSCorrectionFactor = prefwdCS / preadjCS;
      }
      else
      {
        fForwardCSUsed          = false;
        fLastCSCorrectionFactor = 1.;
      }
    }
    corr_fac = fLastCSCorrectionFactor;
  }
  else
  {
    fForwardCSUsed          = false;
    fLastCSCorrectionFactor = 1.;
  }
  fwd_is_used = fForwardCSUsed;
  return corr_fac;
}

///////////////////////////////////////////////////////
G4double G4AdjointCSManager::GetContinuousWeightCorrection(
  G4ParticleDefinition* aPartDef, G4double PreStepEkin, G4double AfterStepEkin,
  const G4MaterialCutsCouple* aCouple, G4double step_length)
{
  G4double corr_fac    = 1.;
  G4double after_fwdCS = GetTotalForwardCS(aPartDef, AfterStepEkin, aCouple);
  G4double pre_adjCS   = GetTotalAdjointCS(aPartDef, PreStepEkin, aCouple);
  if(!fForwardCSUsed || pre_adjCS == 0. || after_fwdCS == 0.)
  {
    G4double pre_fwdCS = GetTotalForwardCS(aPartDef, PreStepEkin, aCouple);
    corr_fac *= std::exp((pre_adjCS - pre_fwdCS) * step_length);
    fLastCSCorrectionFactor = 1.;
  }
  else
  {
    fLastCSCorrectionFactor = after_fwdCS / pre_adjCS;
  }
  return corr_fac;
}

///////////////////////////////////////////////////////
G4double G4AdjointCSManager::GetPostStepWeightCorrection()
{
  return 1. / fLastCSCorrectionFactor;
}

///////////////////////////////////////////////////////
G4double G4AdjointCSManager::ComputeAdjointCS(
  G4Material* aMaterial, G4VEmAdjointModel* aModel, G4double PrimEnergy,
  G4double Tcut, G4bool isScatProjToProj, std::vector<G4double>& CS_Vs_Element)
{
  G4double EminSec = 0.;
  G4double EmaxSec = 0.;

  static G4double lastPrimaryEnergy = 0.;
  static G4double lastTcut          = 0.;
  static G4Material* lastMaterial   = nullptr;

  if(isScatProjToProj)
  {
    EminSec = aModel->GetSecondAdjEnergyMinForScatProjToProj(PrimEnergy, Tcut);
    EmaxSec = aModel->GetSecondAdjEnergyMaxForScatProjToProj(PrimEnergy);
  }
  else if(PrimEnergy > Tcut || !aModel->GetApplyCutInRange())
  {
    EminSec = aModel->GetSecondAdjEnergyMinForProdToProj(PrimEnergy);
    EmaxSec = aModel->GetSecondAdjEnergyMaxForProdToProj(PrimEnergy);
  }
  if(EminSec >= EmaxSec)
    return 0.;

  G4bool need_to_compute = false;
  if(aMaterial != lastMaterial || PrimEnergy != lastPrimaryEnergy ||
     Tcut != lastTcut)
  {
    lastMaterial      = aMaterial;
    lastPrimaryEnergy = PrimEnergy;
    lastTcut          = Tcut;
    fIndexOfAdjointEMModelInAction.clear();
    fIsScatProjToProj.clear();
    fLastAdjointCSVsModelsAndElements.clear();
    need_to_compute = true;
  }

  std::size_t ind = 0;
  if(!need_to_compute)
  {
    need_to_compute = true;
    for(std::size_t i = 0; i < fIndexOfAdjointEMModelInAction.size(); ++i)
    {
      std::size_t ind1 = fIndexOfAdjointEMModelInAction[i];
      if(aModel == fAdjointModels[ind1] &&
         isScatProjToProj == fIsScatProjToProj[i])
      {
        need_to_compute = false;
        CS_Vs_Element   = fLastAdjointCSVsModelsAndElements[ind];
      }
      ++ind;
    }
  }

  if(need_to_compute)
  {
    std::size_t ind_model = 0;
    for(std::size_t i = 0; i < fAdjointModels.size(); ++i)
    {
      if(aModel == fAdjointModels[i])
      {
        ind_model = i;
        break;
      }
    }
    G4double Tlow = Tcut;
    if(!fAdjointModels[ind_model]->GetApplyCutInRange())
      Tlow = fAdjointModels[ind_model]->GetLowEnergyLimit();
    fIndexOfAdjointEMModelInAction.push_back(ind_model);
    fIsScatProjToProj.push_back(isScatProjToProj);
    CS_Vs_Element.clear();
    if(!aModel->GetUseMatrix())
    {
      CS_Vs_Element.push_back(aModel->AdjointCrossSection(
        fCurrentCouple, PrimEnergy, isScatProjToProj));
    }
    else if(aModel->GetUseMatrixPerElement())
    {
      std::size_t n_el = aMaterial->GetNumberOfElements();
      if(aModel->GetUseOnlyOneMatrixForAllElements())
      {
        G4AdjointCSMatrix* theCSMatrix;
        if(isScatProjToProj)
        {
          theCSMatrix = fAdjointCSMatricesForScatProjToProj[ind_model][0];
        }
        else
          theCSMatrix = fAdjointCSMatricesForProdToProj[ind_model][0];
        G4double CS = 0.;
        if(PrimEnergy > Tlow)
          CS = ComputeAdjointCS(PrimEnergy, theCSMatrix, Tlow);
        G4double factor = 0.;
        for(G4int i = 0; i < (G4int)n_el; ++i)
        {  // this could be computed only once
          factor += aMaterial->GetElement(i)->GetZ() *
                    aMaterial->GetVecNbOfAtomsPerVolume()[i];
        }
        CS *= factor;
        CS_Vs_Element.push_back(CS);
      }
      else
      {
        for(G4int i = 0; i < (G4int)n_el; ++i)
        {
          std::size_t ind_el = aMaterial->GetElement(i)->GetIndex();
          G4AdjointCSMatrix* theCSMatrix;
          if(isScatProjToProj)
          {
            theCSMatrix =
              fAdjointCSMatricesForScatProjToProj[ind_model][ind_el];
          }
          else
            theCSMatrix = fAdjointCSMatricesForProdToProj[ind_model][ind_el];
          G4double CS = 0.;
          if(PrimEnergy > Tlow)
            CS = ComputeAdjointCS(PrimEnergy, theCSMatrix, Tlow);
          CS_Vs_Element.push_back(CS *
                                  (aMaterial->GetVecNbOfAtomsPerVolume()[i]));
        }
      }
    }
    else
    {
      std::size_t ind_mat = aMaterial->GetIndex();
      G4AdjointCSMatrix* theCSMatrix;
      if(isScatProjToProj)
      {
        theCSMatrix = fAdjointCSMatricesForScatProjToProj[ind_model][ind_mat];
      }
      else
        theCSMatrix = fAdjointCSMatricesForProdToProj[ind_model][ind_mat];
      G4double CS = 0.;
      if(PrimEnergy > Tlow)
        CS = ComputeAdjointCS(PrimEnergy, theCSMatrix, Tlow);
      CS_Vs_Element.push_back(CS);
    }
    fLastAdjointCSVsModelsAndElements.push_back(CS_Vs_Element);
  }

  G4double CS = 0.;
  for(const auto& cs_vs_el : CS_Vs_Element)
  {
    // We could put the progressive sum of the CS instead of the CS of an
    // element itself
    CS += cs_vs_el;
  }
  return CS;
}

///////////////////////////////////////////////////////
G4Element* G4AdjointCSManager::SampleElementFromCSMatrices(
  G4Material* aMaterial, G4VEmAdjointModel* aModel, G4double PrimEnergy,
  G4double Tcut, G4bool isScatProjToProj)
{
  std::vector<G4double> CS_Vs_Element;
  G4double CS    = ComputeAdjointCS(aMaterial, aModel, PrimEnergy, Tcut,
                                 isScatProjToProj, CS_Vs_Element);
  G4double SumCS = 0.;
  std::size_t ind     = 0;
  for(std::size_t i = 0; i < CS_Vs_Element.size(); ++i)
  {
    SumCS += CS_Vs_Element[i];
    if(G4UniformRand() <= SumCS / CS)
    {
      ind = i;
      break;
    }
  }

  return const_cast<G4Element*>(aMaterial->GetElement((G4int)ind));
}

///////////////////////////////////////////////////////
G4double G4AdjointCSManager::ComputeTotalAdjointCS(
  const G4MaterialCutsCouple* aCouple, G4ParticleDefinition* aPartDef,
  G4double Ekin)
{
  G4double TotalCS = 0.;

  DefineCurrentMaterial(aCouple);

  std::vector<G4double> CS_Vs_Element;
  G4double CS;
  G4VEmAdjointModel* adjModel = nullptr;
  for(std::size_t i = 0; i < fAdjointModels.size(); ++i)
  {
    G4double Tlow = 0.;
    adjModel      = fAdjointModels[i];
    if(!adjModel->GetApplyCutInRange())
      Tlow = adjModel->GetLowEnergyLimit();
    else
    {
      G4ParticleDefinition* theDirSecondPartDef = GetForwardParticleEquivalent(
        adjModel->GetAdjointEquivalentOfDirectSecondaryParticleDefinition());
      std::size_t idx = 56;
      if(theDirSecondPartDef->GetParticleName() == "gamma")
        idx = 0;
      else if(theDirSecondPartDef->GetParticleName() == "e-")
        idx = 1;
      else if(theDirSecondPartDef->GetParticleName() == "e+")
        idx = 2;
      if(idx < 56)
      {
        const std::vector<G4double>* aVec =
          G4ProductionCutsTable::GetProductionCutsTable()->GetEnergyCutsVector(
            idx);
        Tlow = (*aVec)[aCouple->GetIndex()];
      }
    }
    if(Ekin <= adjModel->GetHighEnergyLimit() &&
       Ekin >= adjModel->GetLowEnergyLimit())
    {
      if(aPartDef ==
         adjModel->GetAdjointEquivalentOfDirectPrimaryParticleDefinition())
      {
        CS = ComputeAdjointCS(fCurrentMaterial, adjModel, Ekin, Tlow, true,
                              CS_Vs_Element);
        TotalCS += CS;
        (*fSigmaTableForAdjointModelScatProjToProj[i])[fCurrentMatIndex]
          ->PutValue(fNbins, CS);
      }
      if(aPartDef ==
         adjModel->GetAdjointEquivalentOfDirectSecondaryParticleDefinition())
      {
        CS = ComputeAdjointCS(fCurrentMaterial, adjModel, Ekin, Tlow, false,
                              CS_Vs_Element);
        TotalCS += CS;
        (*fSigmaTableForAdjointModelProdToProj[i])[fCurrentMatIndex]->PutValue(
          fNbins, CS);
      }
    }
    else
    {
      (*fSigmaTableForAdjointModelScatProjToProj[i])[fCurrentMatIndex]
        ->PutValue(fNbins, 0.);
      (*fSigmaTableForAdjointModelProdToProj[i])[fCurrentMatIndex]->PutValue(
        fNbins, 0.);
    }
  }
  return TotalCS;
}

///////////////////////////////////////////////////////
std::vector<G4AdjointCSMatrix*>
G4AdjointCSManager::BuildCrossSectionsModelAndElement(G4VEmAdjointModel* aModel,
                                                      G4int Z, G4int A,
                                                      G4int nbin_pro_decade)
{
  G4AdjointCSMatrix* theCSMatForProdToProjBackwardScattering =
    new G4AdjointCSMatrix(false);
  G4AdjointCSMatrix* theCSMatForScatProjToProjBackwardScattering =
    new G4AdjointCSMatrix(true);

  // make the vector of primary energy of the adjoint particle. 
  G4double EkinMin        = aModel->GetLowEnergyLimit();
  G4double EkinMaxForScat = aModel->GetHighEnergyLimit() * 0.999;
  G4double EkinMaxForProd = aModel->GetHighEnergyLimit() * 0.999;
  if(aModel->GetSecondPartOfSameType())
    EkinMaxForProd = EkinMaxForProd / 2.;

  // Product to projectile backward scattering
  G4double dE = std::pow(10., 1. / nbin_pro_decade);
  G4double E2 =
    std::pow(10., double(int(std::log10(EkinMin) * nbin_pro_decade) + 1) /
                    nbin_pro_decade) / dE;
  G4double E1 = EkinMin;
  // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
  while(E1 < EkinMaxForProd)
  {
    E1 = std::max(EkinMin, E2);
    E1 = std::min(EkinMaxForProd, E1);
    std::vector<std::vector<double>*> aMat =
      aModel->ComputeAdjointCrossSectionVectorPerAtomForSecond(E1, Z, A,
                                                               nbin_pro_decade);
    if(aMat.size() >= 2)
    {
      std::vector<double>* log_ESecVec = aMat[0];
      std::vector<double>* log_CSVec   = aMat[1];
      G4double log_adjointCS           = log_CSVec->back();
      // normalise CSVec such that it becomes a probability vector
      for(std::size_t j = 0; j < log_CSVec->size(); ++j)
      {
        if(j == 0)
          (*log_CSVec)[j] = 0.;
        else
          (*log_CSVec)[j] =
            std::log(1. - std::exp((*log_CSVec)[j] - log_adjointCS) + 1e-50);
      }
      (*log_CSVec)[log_CSVec->size() - 1] =
        (*log_CSVec)[log_CSVec->size() - 2] - std::log(1000.);
      theCSMatForProdToProjBackwardScattering->AddData(
        std::log(E1), log_adjointCS, log_ESecVec, log_CSVec, 0);
    }
    E1 = E2;
    E2 *= dE;
  }

  // Scattered projectile to projectile backward scattering
  E2 = std::pow(10., G4double(int(std::log10(EkinMin) * nbin_pro_decade) + 1) /
                       nbin_pro_decade) / dE;
  E1 = EkinMin;
  // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
  while(E1 < EkinMaxForScat)
  {
    E1 = std::max(EkinMin, E2);
    E1 = std::min(EkinMaxForScat, E1);
    std::vector<std::vector<G4double>*> aMat =
      aModel->ComputeAdjointCrossSectionVectorPerAtomForScatProj(
        E1, Z, A, nbin_pro_decade);
    if(aMat.size() >= 2)
    {
      std::vector<G4double>* log_ESecVec = aMat[0];
      std::vector<G4double>* log_CSVec   = aMat[1];
      G4double log_adjointCS             = log_CSVec->back();
      // normalise CSVec such that it becomes a probability vector
      for(std::size_t j = 0; j < log_CSVec->size(); ++j)
      {
        if(j == 0)
          (*log_CSVec)[j] = 0.;
        else
          (*log_CSVec)[j] =
            std::log(1. - std::exp((*log_CSVec)[j] - log_adjointCS) + 1e-50);
      }
      (*log_CSVec)[log_CSVec->size() - 1] =
        (*log_CSVec)[log_CSVec->size() - 2] - std::log(1000.);
      theCSMatForScatProjToProjBackwardScattering->AddData(
        std::log(E1), log_adjointCS, log_ESecVec, log_CSVec, 0);
    }
    E1 = E2;
    E2 *= dE;
  }

  std::vector<G4AdjointCSMatrix*> res;
  res.push_back(theCSMatForProdToProjBackwardScattering);
  res.push_back(theCSMatForScatProjToProjBackwardScattering);

  return res;
}

///////////////////////////////////////////////////////
std::vector<G4AdjointCSMatrix*>
G4AdjointCSManager::BuildCrossSectionsModelAndMaterial(
  G4VEmAdjointModel* aModel, G4Material* aMaterial, G4int nbin_pro_decade)
{
  G4AdjointCSMatrix* theCSMatForProdToProjBackwardScattering =
    new G4AdjointCSMatrix(false);
  G4AdjointCSMatrix* theCSMatForScatProjToProjBackwardScattering =
    new G4AdjointCSMatrix(true);

  G4double EkinMin        = aModel->GetLowEnergyLimit();
  G4double EkinMaxForScat = aModel->GetHighEnergyLimit() * 0.999;
  G4double EkinMaxForProd = aModel->GetHighEnergyLimit() * 0.999;
  if(aModel->GetSecondPartOfSameType())
    EkinMaxForProd /= 2.;

  // Product to projectile backward scattering
  G4double dE = std::pow(10., 1. / nbin_pro_decade);
  G4double E2 =
    std::pow(10., double(int(std::log10(EkinMin) * nbin_pro_decade) + 1) /
                    nbin_pro_decade) / dE;
  G4double E1 = EkinMin;
  // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
  while(E1 < EkinMaxForProd)
  {
    E1 = std::max(EkinMin, E2);
    E1 = std::min(EkinMaxForProd, E1);
    std::vector<std::vector<G4double>*> aMat =
      aModel->ComputeAdjointCrossSectionVectorPerVolumeForSecond(
        aMaterial, E1, nbin_pro_decade);
    if(aMat.size() >= 2)
    {
      std::vector<G4double>* log_ESecVec = aMat[0];
      std::vector<G4double>* log_CSVec   = aMat[1];
      G4double log_adjointCS             = log_CSVec->back();

      // normalise CSVec such that it becomes a probability vector
      for(std::size_t j = 0; j < log_CSVec->size(); ++j)
      {
        if(j == 0)
          (*log_CSVec)[j] = 0.;
        else
          (*log_CSVec)[j] =
            std::log(1. - std::exp((*log_CSVec)[j] - log_adjointCS));
      }
      (*log_CSVec)[log_CSVec->size() - 1] =
        (*log_CSVec)[log_CSVec->size() - 2] - std::log(1000.);
      theCSMatForProdToProjBackwardScattering->AddData(
        std::log(E1), log_adjointCS, log_ESecVec, log_CSVec, 0);
    }

    E1 = E2;
    E2 *= dE;
  }

  // Scattered projectile to projectile backward scattering
  E2 =
    std::pow(10., G4double(G4int(std::log10(EkinMin) * nbin_pro_decade) + 1) /
                    nbin_pro_decade) /
    dE;
  E1 = EkinMin;
  while(E1 < EkinMaxForScat)
  {
    E1 = std::max(EkinMin, E2);
    E1 = std::min(EkinMaxForScat, E1);
    std::vector<std::vector<G4double>*> aMat =
      aModel->ComputeAdjointCrossSectionVectorPerVolumeForScatProj(
        aMaterial, E1, nbin_pro_decade);
    if(aMat.size() >= 2)
    {
      std::vector<G4double>* log_ESecVec = aMat[0];
      std::vector<G4double>* log_CSVec   = aMat[1];
      G4double log_adjointCS             = log_CSVec->back();

      for(std::size_t j = 0; j < log_CSVec->size(); ++j)
      {
        if(j == 0)
          (*log_CSVec)[j] = 0.;
        else
          (*log_CSVec)[j] =
            std::log(1. - std::exp((*log_CSVec)[j] - log_adjointCS));
      }
      (*log_CSVec)[log_CSVec->size() - 1] =
        (*log_CSVec)[log_CSVec->size() - 2] - std::log(1000.);

      theCSMatForScatProjToProjBackwardScattering->AddData(
        std::log(E1), log_adjointCS, log_ESecVec, log_CSVec, 0);
    }
    E1 = E2;
    E2 *= dE;
  }

  std::vector<G4AdjointCSMatrix*> res;
  res.push_back(theCSMatForProdToProjBackwardScattering);
  res.push_back(theCSMatForScatProjToProjBackwardScattering);

  return res;
}

///////////////////////////////////////////////////////
G4ParticleDefinition* G4AdjointCSManager::GetAdjointParticleEquivalent(
  G4ParticleDefinition* theFwdPartDef)
{
  if(theFwdPartDef->GetParticleName() == "e-")
    return G4AdjointElectron::AdjointElectron();
  else if(theFwdPartDef->GetParticleName() == "gamma")
    return G4AdjointGamma::AdjointGamma();
  else if(theFwdPartDef->GetParticleName() == "proton")
    return G4AdjointProton::AdjointProton();
  else if(theFwdPartDef == fFwdIon)
    return fAdjIon;
  return nullptr;
}

///////////////////////////////////////////////////////
G4ParticleDefinition* G4AdjointCSManager::GetForwardParticleEquivalent(
  G4ParticleDefinition* theAdjPartDef)
{
  if(theAdjPartDef->GetParticleName() == "adj_e-")
    return G4Electron::Electron();
  else if(theAdjPartDef->GetParticleName() == "adj_gamma")
    return G4Gamma::Gamma();
  else if(theAdjPartDef->GetParticleName() == "adj_proton")
    return G4Proton::Proton();
  else if(theAdjPartDef == fAdjIon)
    return fFwdIon;
  return nullptr;
}

///////////////////////////////////////////////////////
void G4AdjointCSManager::DefineCurrentMaterial(
  const G4MaterialCutsCouple* couple)
{
  if(couple != fCurrentCouple)
  {
    fCurrentCouple          = const_cast<G4MaterialCutsCouple*>(couple);
    fCurrentMaterial        = const_cast<G4Material*>(couple->GetMaterial());
    fCurrentMatIndex        = couple->GetIndex();
    fLastCSCorrectionFactor = 1.;
  }
}

///////////////////////////////////////////////////////
void G4AdjointCSManager::DefineCurrentParticle(
  const G4ParticleDefinition* aPartDef)
{
  static G4ParticleDefinition* currentParticleDef = nullptr;

  if(aPartDef != currentParticleDef)
  {
    currentParticleDef = const_cast<G4ParticleDefinition*>(aPartDef);
    fMassRatio         = 1.;
    if(aPartDef == fAdjIon)
      fMassRatio = proton_mass_c2 / aPartDef->GetPDGMass();
    fCurrentParticleIndex = 1000000;
    for(std::size_t i = 0; i < fAdjointParticlesInAction.size(); ++i)
    {
      if(aPartDef == fAdjointParticlesInAction[i])
        fCurrentParticleIndex = i;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
G4double G4AdjointCSManager::ComputeAdjointCS(
  G4double aPrimEnergy, G4AdjointCSMatrix* anAdjointCSMatrix, G4double Tcut)
{
  std::vector<double>* theLogPrimEnergyVector =
    anAdjointCSMatrix->GetLogPrimEnergyVector();
  if(theLogPrimEnergyVector->empty())
  {
    G4cout << "No data are contained in the given AdjointCSMatrix!" << G4endl;
    return 0.;
  }
  G4double log_Tcut = std::log(Tcut);
  G4double log_E    = std::log(aPrimEnergy);

  if(aPrimEnergy <= Tcut || log_E > theLogPrimEnergyVector->back())
    return 0.;

  G4AdjointInterpolator* theInterpolator = G4AdjointInterpolator::GetInstance();

  std::size_t ind =
    theInterpolator->FindPositionForLogVector(log_E, *theLogPrimEnergyVector);
  G4double aLogPrimEnergy1, aLogPrimEnergy2;
  G4double aLogCS1, aLogCS2;
  G4double log01, log02;
  std::vector<G4double>* aLogSecondEnergyVector1 = nullptr;
  std::vector<G4double>* aLogSecondEnergyVector2 = nullptr;
  std::vector<G4double>* aLogProbVector1         = nullptr;
  std::vector<G4double>* aLogProbVector2         = nullptr;
  std::vector<std::size_t>* aLogProbVectorIndex1 = nullptr;
  std::vector<std::size_t>* aLogProbVectorIndex2 = nullptr;

  anAdjointCSMatrix->GetData((G4int)ind, aLogPrimEnergy1, aLogCS1, log01,
                             aLogSecondEnergyVector1, aLogProbVector1,
                             aLogProbVectorIndex1);
  anAdjointCSMatrix->GetData(G4int(ind + 1), aLogPrimEnergy2, aLogCS2, log02,
                             aLogSecondEnergyVector2, aLogProbVector2,
                             aLogProbVectorIndex2);
  if (! (aLogProbVector1 && aLogProbVector2 &&
  		       aLogSecondEnergyVector1 && aLogSecondEnergyVector2)){
  	 return  0.;
  }

  if(anAdjointCSMatrix->IsScatProjToProj())
  {  // case where the Tcut plays a role
    G4double log_minimum_prob1, log_minimum_prob2;
    log_minimum_prob1 = theInterpolator->InterpolateForLogVector(
      log_Tcut, *aLogSecondEnergyVector1, *aLogProbVector1);
    log_minimum_prob2 = theInterpolator->InterpolateForLogVector(
      log_Tcut, *aLogSecondEnergyVector2, *aLogProbVector2);
    aLogCS1 += log_minimum_prob1;
    aLogCS2 += log_minimum_prob2;
  }

  G4double log_adjointCS = theInterpolator->LinearInterpolation(
    log_E, aLogPrimEnergy1, aLogPrimEnergy2, aLogCS1, aLogCS2);
  return std::exp(log_adjointCS);
}
