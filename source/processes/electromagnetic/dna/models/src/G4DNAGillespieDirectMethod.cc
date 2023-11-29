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

#include "G4DNAGillespieDirectMethod.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include <memory>
#include <tuple>
#include "G4DNAEventSet.hh"
#include "G4UnitsTable.hh"
#include "G4DNAScavengerMaterial.hh"
#include "G4Scheduler.hh"
#include "G4DNAMolecularReactionTable.hh"

G4DNAGillespieDirectMethod::G4DNAGillespieDirectMethod()
  : fMolecularReactions(G4DNAMolecularReactionTable::Instance())
{}

G4DNAGillespieDirectMethod::~G4DNAGillespieDirectMethod() = default;

void G4DNAGillespieDirectMethod::SetEventSet(G4DNAEventSet* pEventSet)
{
  fpEventSet = pEventSet;
}

//#define DEBUG 1

G4double G4DNAGillespieDirectMethod::VolumeOfNode(const Voxel& voxel)
{
  auto box     = std::get<1>(voxel);
  auto LengthY = box.Getyhi() - box.Getylo();
  auto LengthX = box.Getxhi() - box.Getxlo();
  auto LengthZ = box.Getzhi() - box.Getzlo();
  G4double V   = LengthY * LengthX * LengthZ;
  if(V <= 0)
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "V > 0 !! ";
    G4Exception("G4DNAGillespieDirectMethod::VolumeOfNode",
                "G4DNAGillespieDirectMethod03", FatalErrorInArgument,
                exceptionDescription);
  }
  return V;
}
G4double G4DNAGillespieDirectMethod::PropensityFunction(const Voxel& voxel,
                                                        MolType moleType)
{
  if(moleType->GetDiffusionCoefficient() == 0)
  {
    return 0.;
  }
  const auto& node = std::get<2>(voxel);
  const auto& box  = std::get<1>(voxel);

  G4double alpha = 0;
  auto it        = node.find(moleType);
  if(it != node.end())
  {
    auto LengthY = box.Getyhi() - box.Getylo();
    G4double d   = it->first->GetDiffusionCoefficient() / std::pow(LengthY, 2);
    alpha        = d * it->second;

#ifdef DEBUG
    G4cout << it->first->GetName() << " " << it->second
           << " D : " << it->first->GetDiffusionCoefficient()
           << " LengthY : " << LengthY << " PropensityFunction : " << alpha
           << G4endl;
#endif
  }
  return alpha;
}

G4double G4DNAGillespieDirectMethod::PropensityFunction(const Voxel& voxel,
                                                        ReactionData* data)
{
  G4double value;
  auto ConfA               = data->GetReactant1();
  auto ConfB               = data->GetReactant2();
  G4double scavengerNumber = 0;
  auto typeANumber         = FindScavenging(voxel, ConfA, scavengerNumber)
                               ? scavengerNumber
                               : ComputeNumberInNode(voxel, ConfA);

  auto typeBNumber = FindScavenging(voxel, ConfB, scavengerNumber)
                       ? scavengerNumber
                       : ComputeNumberInNode(voxel, ConfB);

  if(typeANumber == 0 || typeBNumber == 0)
  {
    return 0;
  }

  auto k =
    data->GetObservedReactionRateConstant() / (Avogadro * VolumeOfNode(voxel));
  if(ConfA == ConfB)
  {
    value = typeANumber * (typeBNumber - 1) * k;
  }
  else
  {
    value = typeANumber * typeBNumber * k;
  }

  if(value < 0)
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription
      << "G4DNAGillespieDirectMethod::PropensityFunction for : "
      << ConfA->GetName() << "(" << typeANumber << ") + " << ConfB->GetName()
      << "(" << typeBNumber << ") : propensity : " << value
      << " GetObservedReactionRateConstant : "
      << data->GetObservedReactionRateConstant()
      << " GetEffectiveReactionRadius : "
      << G4BestUnit(data->GetEffectiveReactionRadius(), "Length")
      << " k : " << k << " volume : " << VolumeOfNode(voxel) << G4endl;
    G4Exception("G4DNAGillespieDirectMethod::PropensityFunction",
                "G4DNAGillespieDirectMethod013", FatalErrorInArgument,
                exceptionDescription);
  }

#ifdef DEBUG
  if(value > 0)
    G4cout << "G4DNAGillespieDirectMethod::PropensityFunction for : "
           << ConfA->GetName() << "(" << typeANumber << ") + "
           << ConfB->GetName() << "(" << typeBNumber
           << ") : propensity : " << value
           << "  Time to Reaction : " << G4BestUnit(timeToReaction, "Time")
           << " k : " << k << " Index : " << index << G4endl;
#endif

  return value;
}

void G4DNAGillespieDirectMethod::Initialize()
{
  // for Scavenger
  fpScavengerMaterial = dynamic_cast<G4DNAScavengerMaterial*>(
    G4Scheduler::Instance()->GetScavengerMaterial());

  auto begin = fpMesh->begin();
  auto end   = fpMesh->end();
  for(; begin != end; begin++)
  {
    auto index = std::get<0>(*begin);
#ifdef DEBUG
    fpMesh->PrintVoxel(index);
#endif
    CreateEvent(index);
  }
}

void G4DNAGillespieDirectMethod::SetTimeStep(const G4double& stepTime)
{
  fTimeStep = stepTime;
}
void G4DNAGillespieDirectMethod::CreateEvent(const Index& index)
{
  const auto& voxel = fpMesh->GetVoxel(index);
  if(std::get<2>(voxel).empty())
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "This voxel : " << index
                         << " is not ready to make event" << G4endl;
    G4Exception("G4DNAGillespieDirectMethod::CreateEvent",
                "G4DNAGillespieDirectMethod05", FatalErrorInArgument,
                exceptionDescription);
  }
  G4double r1         = G4UniformRand();
  G4double r2         = G4UniformRand();
  G4double dAlpha0    = DiffusiveJumping(voxel);
  G4double rAlpha0    = Reaction(voxel);
  G4double alphaTotal = dAlpha0 + rAlpha0;

  if(alphaTotal == 0)
  {
    return;
  }
  auto timeStep = ((1.0 / (alphaTotal)) * std::log(1.0 / r1)) + fTimeStep;

#ifdef DEBUG
  G4cout << "r2 : " << r2 << " rAlpha0 : " << rAlpha0
         << " dAlpha0 : " << dAlpha0 << "    rAlpha0 / (dAlpha0 + rAlpha0) : "
         << rAlpha0 / (dAlpha0 + rAlpha0) << G4endl;
#endif
  if(r2 < rAlpha0 / alphaTotal)
  {
    if(fVerbose > 1)
    {
      G4cout << "=>>>>reaction at : " << timeStep << " timeStep : "
             << G4BestUnit(((1.0 / alphaTotal) * std::log(1.0 / r1)), "Time")
             << G4endl;
    }
    auto rSelectedIter = fReactionDataMap.upper_bound(r2 * alphaTotal);
    fpEventSet->CreateEvent(timeStep, index, rSelectedIter->second);
  }
  else if(dAlpha0 > 0)
  {
    if(fVerbose > 1)
    {
      G4cout << "=>>>>jumping at : " << timeStep << " timeStep : "
             << G4BestUnit(((1.0 / alphaTotal) * std::log(1.0 / r1)), "Time")
             << G4endl;
    }

    auto dSelectedIter = fJumpingDataMap.upper_bound(r2 * alphaTotal - rAlpha0);
    auto pDSelected =
      std::make_unique<std::pair<MolType, Index>>(dSelectedIter->second);
    fpEventSet->CreateEvent(timeStep, index, std::move(pDSelected));
  }
#ifdef DEBUG
  G4cout << G4endl;
#endif
}

G4double G4DNAGillespieDirectMethod::Reaction(const Voxel& voxel)
{
  fReactionDataMap.clear();
  G4double alpha0 = 0;
  const auto& dataList =
    fMolecularReactions->GetVectorOfReactionData();  // shoud make a member
  if(dataList.empty())
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "MolecularReactionTable empty" << G4endl;
    G4Exception("G4DNAGillespieDirectMethod::Reaction",
                "G4DNAGillespieDirectMethod01", FatalErrorInArgument,
                exceptionDescription);
  }

  for(const auto& it : dataList)
  {
    auto propensity = PropensityFunction(voxel, it);
    if(propensity == 0)
    {
      continue;
    }
    alpha0 += propensity;
    fReactionDataMap[alpha0] = it;
  }
#ifdef DEBUG
  G4cout << "Reaction :alpha0 :  " << alpha0 << G4endl;
#endif
  return alpha0;
}

G4double G4DNAGillespieDirectMethod::DiffusiveJumping(const Voxel& voxel)
{
  fJumpingDataMap.clear();
  G4double alpha0        = 0;
  auto index             = std::get<0>(voxel);
  auto NeighboringVoxels = fpMesh->FindNeighboringVoxels(index);
  if(NeighboringVoxels.empty())
  {
    return 0;
  }
  auto iter = G4MoleculeTable::Instance()->GetConfigurationIterator();
  while(iter())
  {
    const auto* conf = iter.value();
    auto propensity  = PropensityFunction(voxel, conf);
    if(propensity == 0)
    {
      continue;
    }
    for(const auto& it_Neighbor : NeighboringVoxels)
    {
      alpha0 += propensity;
      fJumpingDataMap[alpha0] = std::make_pair(conf, it_Neighbor);
#ifdef DEBUG
      G4cout << "mole : " << conf->GetName()
             << " number : " << ComputeNumberInNode(index, conf)
             << " propensity : " << propensity << "  alpha0 : " << alpha0
             << G4endl;
#endif
    }
  }
#ifdef DEBUG
  G4cout << "DiffusiveJumping :alpha0 :  " << alpha0 << G4endl;
#endif
  return alpha0;
}

G4double G4DNAGillespieDirectMethod::ComputeNumberInNode(
  const Voxel& voxel, MolType type)  // depend node ?
{
  if(type->GetDiffusionCoefficient() != 0)
  {
    const auto& node = std::get<2>(voxel);
    const auto& it   = node.find(type);
    return (it != node.end()) ? (it->second) : 0;
  }
  else
  {
    return 0;
  }
}

G4bool G4DNAGillespieDirectMethod::FindScavenging(const Voxel& voxel,
                                                  MolType moletype,
                                                  G4double& numberOfScavenger)
{
  numberOfScavenger = 0;
  if(fpScavengerMaterial == nullptr)
  {
    return false;
  }
  auto volumeOfNode = VolumeOfNode(voxel);
  if(G4MoleculeTable::Instance()->GetConfiguration("H2O") == moletype)
  {
    auto factor       = Avogadro * volumeOfNode;
    numberOfScavenger = factor;
    return true;
  }

  G4double totalNumber =
    fpScavengerMaterial->GetNumberMoleculePerVolumeUnitForMaterialConf(
      moletype);
  if(totalNumber == 0)
  {
    return false;
  }
  else
  {
    G4double numberInDouble = volumeOfNode * std::floor(totalNumber) /
                              fpMesh->GetBoundingBox().Volume();
    auto numberInInterg = (int64_t) (std::floor(numberInDouble));
    G4double change     = numberInDouble - numberInInterg;
    G4UniformRand() > change ? numberOfScavenger = numberInInterg
                                 : numberOfScavenger = numberInInterg + 1;
    return true;
  }
}
