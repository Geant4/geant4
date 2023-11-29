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
#include "G4DNAScavengerProcess.hh"
#include <G4VScheduler.hh>
#include <memory>
#include "G4Molecule.hh"
#include "G4DNAMolecularMaterial.hh"
#include "G4MolecularConfiguration.hh"
#include "G4UnitsTable.hh"
#include "G4TrackingInformation.hh"
#include "G4DNAMolecularReactionTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4DNABoundingBox.hh"
#include "G4DNAScavengerMaterial.hh"
#include "G4MoleculeFinder.hh"
#include "G4Scheduler.hh"

#ifndef State
#  define State(theXInfo) (GetState<G4DNAScavengerProcessState>()->theXInfo)
#endif
G4DNAScavengerProcess::G4DNAScavengerProcess(const G4String& aName,
                                             const G4DNABoundingBox& box,
                                             G4ProcessType type)
  : G4VITProcess(aName, type)
  , fpBoundingBox(&box)
  , fpScavengerMaterial(nullptr)
{
  pParticleChange     = &fParticleChange;
  enableAtRestDoIt    = false;
  enableAlongStepDoIt = false;
  enablePostStepDoIt  = true;
  SetProcessSubType(66);
  G4VITProcess::SetInstantiateProcessState(false);
  fIsInitialized           = false;
  fpMolecularConfiguration = nullptr;
  fpMaterialConf           = nullptr;
  fProposesTimeStep        = true;
  fReturnedValue           = DBL_MAX;
  verboseLevel             = 0;
}

G4DNAScavengerProcess::~G4DNAScavengerProcess()
{
  for(auto& iter : fConfMap)
  {
    for(auto& iter2 : iter.second)
    {
      delete iter2.second;
    }
  }
}

G4DNAScavengerProcess::G4DNAScavengerProcessState::G4DNAScavengerProcessState()
  : G4ProcessState()
{
  fPreviousTimeAtPreStepPoint = -1;
  fIsInGoodMaterial           = false;
}
void G4DNAScavengerProcess::BuildPhysicsTable(const G4ParticleDefinition&)
{
  fpScavengerMaterial = dynamic_cast<G4DNAScavengerMaterial*>(
    G4Scheduler::Instance()->GetScavengerMaterial());
  if(fpScavengerMaterial == nullptr)
  {
    return;
  }
  fIsInitialized = true;
}

void G4DNAScavengerProcess::StartTracking(G4Track* track)
{
  G4VProcess::StartTracking(track);
  G4VITProcess::fpState = std::make_shared<G4DNAScavengerProcessState>();
  G4VITProcess::StartTracking(track);
}

void G4DNAScavengerProcess::SetReaction(MolType molConf, Data* pData)
{
  if(fIsInitialized)
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription
      << "G4DNASecondOrderReaction was already initialised. ";
    exceptionDescription << "You cannot set a reaction after initialisation.";
    G4Exception("G4DNASecondOrderReaction::SetReaction",
                "G4DNASecondOrderReaction001", FatalErrorInArgument,
                exceptionDescription);
  }
  auto materialConf = pData->GetReactant1() == molConf ? pData->GetReactant2()
                                                       : pData->GetReactant1();
  if(verboseLevel > 0)
  {
    G4cout << "G4DNAScavengerProcess::SetReaction : " << molConf->GetName()
           << "   materialConf : " << materialConf->GetName() << G4endl;
  }
  fConfMap[molConf][materialConf] = pData;
}

G4double G4DNAScavengerProcess::PostStepGetPhysicalInteractionLength(
  const G4Track& track, G4double /*previousStepSize*/,
  G4ForceCondition* pForceCond)
{
  G4Molecule* molecule = GetMolecule(track);
  auto molConf         = molecule->GetMolecularConfiguration();
  // reset
  fpMolecularConfiguration = nullptr;
  fpMaterialConf           = nullptr;

  // this because process for moleculeDifinition not for configuration
  // TODO: need change this
  auto it = fConfMap.find(molConf);
  if(it == fConfMap.end())
  {
    return DBL_MAX;
  }

  fpMolecularConfiguration = molConf;
  auto MaterialMap         = it->second;

  G4double r1 = G4UniformRand();
  std::map<G4double /*Propensity*/, std::pair<MolType, G4double>>
    ReactionDataMap;
  G4double alpha0 = 0;

  for(const auto& mat_it : MaterialMap)
  {
    auto matConf = mat_it.first;
    G4double numMol =
      fpScavengerMaterial->GetNumberMoleculePerVolumeUnitForMaterialConf(
        matConf);
    if(numMol == 0.0)  // ie : not found
    {
      continue;
    }
    if(verboseLevel > 1)
    {
      G4cout << " Material of " << matConf->GetName() << " : " << numMol
             << G4endl;
    }
    // auto data = fReactionMap[mat_it];
    auto data         = mat_it.second;
    auto reactionRate = data->GetObservedReactionRateConstant();  //_const
    G4double propensity =
      numMol * reactionRate / (fpBoundingBox->Volume() * Avogadro);
    auto reactionData = std::make_pair(mat_it.first, propensity);
    if(propensity == 0)
    {
      continue;
    }
    alpha0 += propensity;
    ReactionDataMap[alpha0] = reactionData;
  }
  if(alpha0 == 0)
  {
    if(State(fIsInGoodMaterial))
    {
      ResetNumberOfInteractionLengthLeft();
      State(fIsInGoodMaterial) = false;
    }
    return DBL_MAX;
  }
  auto rSelectedIter = ReactionDataMap.upper_bound(r1 * alpha0);

  fpMaterialConf = rSelectedIter->second.first;

  State(fIsInGoodMaterial) = true;
  G4double previousTimeStep(-1.);

  if(State(fPreviousTimeAtPreStepPoint) != -1)
  {
    previousTimeStep =
      track.GetGlobalTime() - State(fPreviousTimeAtPreStepPoint);
  }

  State(fPreviousTimeAtPreStepPoint) = track.GetGlobalTime();

  // condition is set to "Not Forced"
  *pForceCond = NotForced;

  if((previousTimeStep <= 0.0) ||
     (fpState->theNumberOfInteractionLengthLeft <= 0.0))
  {
    ResetNumberOfInteractionLengthLeft();
  }
  else if(previousTimeStep > 0.0)
  {
    SubtractNumberOfInteractionLengthLeft(previousTimeStep);
  }

  fpState->currentInteractionLength = 1 / (rSelectedIter->second.second);
  G4double value                    = DBL_MAX;
  if(fpState->currentInteractionLength < DBL_MAX)
  {
    value = fpState->theNumberOfInteractionLengthLeft *
            (fpState->currentInteractionLength);
  }
#ifdef G4VERBOSE
  if(verboseLevel > 2)
  {
    G4cout << "G4DNAScavengerProcess::PostStepGetPhysicalInteractionLength:: "
           << molConf->GetName() << G4endl;
    G4cout << "theNumberOfInteractionLengthLeft : "
           << fpState->theNumberOfInteractionLengthLeft << G4endl;
    G4cout << "currentInteractionLength : " << fpState->currentInteractionLength
           << G4endl;
    G4cout << "Material : " << fpMaterialConf->GetName()
           << " ID: " << track.GetTrackID()
           << " Track Time : " << track.GetGlobalTime()
           << " name : " << molConf->GetName()
           << " Track Position : " << track.GetPosition()
           << " Returned time : " << G4BestUnit(value, "Time") << G4endl;
  }
#endif

  if(value < fReturnedValue)
  {
    fReturnedValue = value;
  }

  return value * -1;
  // multiple by -1 to indicate to the tracking system that we are returning a
  // time
}

G4VParticleChange* G4DNAScavengerProcess::PostStepDoIt(const G4Track& track,
                                                       const G4Step& /*step*/)
{
  G4Molecule* molecule = GetMolecule(track);
  auto molConf         = molecule->GetMolecularConfiguration();
  if(fpMolecularConfiguration != molConf)
  {
    fReturnedValue = DBL_MAX;
    fParticleChange.Initialize(track);
    State(fPreviousTimeAtPreStepPoint) = -1;
    return &fParticleChange;
  }
  std::vector<G4Track*> products;
#ifdef G4VERBOSE
  if(verboseLevel > 1)
  {
    G4cout << "___________" << G4endl;
    G4cout << ">>> Beginning of G4DNAScavengerProcess verbose" << G4endl;
    G4cout << ">>> Returned value : " << G4BestUnit(fReturnedValue, "Time")
           << G4endl;
    G4cout << ">>> Time Step : "
           << G4BestUnit(G4VScheduler::Instance()->GetTimeStep(), "Time")
           << G4endl;
    G4cout << ">>> Global Time : "
           << G4BestUnit(G4VScheduler::Instance()->GetGlobalTime(), "Time")
           << G4endl;
    G4cout << ">>> Global Time Track : "
           << G4BestUnit(track.GetGlobalTime(), "Time") << G4endl;
    G4cout << ">>> Track Position : " << track.GetPosition() << G4endl;
    G4cout << ">>> Reaction : " << molecule->GetName() << "("
           << track.GetTrackID() << ") + " << fpMaterialConf->GetName()
           << G4endl;
    G4cout << ">>> End of G4DNAScavengerProcess verbose <<<" << G4endl;
  }
#endif

  G4double reactionTime = track.GetGlobalTime();
  auto data = fConfMap[fpMolecularConfiguration][fpMaterialConf];

  auto nbSecondaries = data->GetNbProducts();

  for(G4int j = 0; j < nbSecondaries; ++j)
  {
    auto pProduct = new G4Molecule(data->GetProduct(j));
    auto pProductTrack =
      pProduct->BuildTrack(reactionTime, track.GetPosition());
    pProductTrack->SetTrackStatus(fAlive);
    G4ITTrackHolder::Instance()->Push(pProductTrack);
    G4MoleculeFinder::Instance()->Push(pProductTrack);
    products.push_back(pProductTrack);
  }

#ifdef G4VERBOSE
  if(verboseLevel != 0)
  {
    G4cout << "At time : " << std::setw(7) << G4BestUnit(reactionTime, "Time")
           << " Reaction : " << GetIT(track)->GetName() << " ("
           << track.GetTrackID() << ") + " << fpMaterialConf->GetName() << " ("
           << "B"
           << ") -> ";
  }
#endif
  if(nbSecondaries > 0)
  {
    for(G4int i = 0; i < nbSecondaries; ++i)
    {
#ifdef G4VERBOSE
      if((verboseLevel != 0) && i != 0)
      {
        G4cout << " + ";
      }

      if(verboseLevel != 0)
      {
        G4cout << GetIT(products.at(i))->GetName() << " ("
               << products.at(i)->GetTrackID() << ")";
      }
#endif
    }
  }
  else
  {
#ifdef G4VERBOSE
    if(verboseLevel != 0)
    {
      G4cout << "No product";
    }
#endif
  }
#ifdef G4VERBOSE
  if(verboseLevel != 0)
  {
    G4cout << G4endl;
  }
#endif

  fReturnedValue = DBL_MAX;
  fParticleChange.Initialize(track);
  fParticleChange.ProposeTrackStatus(fStopAndKill);
  fpScavengerMaterial->ReduceNumberMoleculePerVolumeUnitForMaterialConf(
    fpMaterialConf, reactionTime);
  State(fPreviousTimeAtPreStepPoint) = -1;
  return &fParticleChange;
}
