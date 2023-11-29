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
#include "G4DNAUpdateSystemModel.hh"
#include "G4Molecule.hh"
#include "G4DNAMolecularReactionTable.hh"
#include "G4UnitsTable.hh"
#include "G4MoleculeCounter.hh"
#include "G4DNAScavengerMaterial.hh"
#include "G4Scheduler.hh"
G4DNAUpdateSystemModel::G4DNAUpdateSystemModel()
  : G4VUpdateSystemModel()
{}

void G4DNAUpdateSystemModel::SetMesh(G4DNAMesh* pMesh) { fpMesh = pMesh; }
void G4DNAUpdateSystemModel::KillMolecule(const Index& index, MolType type)
{
  // kill normal molecule
  auto& node = fpMesh->GetVoxelMapList(index);
  auto iter  = node.find(type);
  if(iter != node.end())
  {
    if(iter->second <= 0)
    {
      G4ExceptionDescription exceptionDescription;
      exceptionDescription
        << "G4DNAUpdateSystemModel::KillMolecule::molecule : "
        << type->GetName() << " index : " << index
        << " number : " << iter->second << G4endl;
      G4Exception("G4DNAEventScheduler::Stepping", "G4DNAEventScheduler002",
                  FatalErrorInArgument, exceptionDescription);
    }
    iter->second--;
    if(G4VMoleculeCounter::Instance()->InUse())
    {
      G4VMoleculeCounter::Instance()->RemoveAMoleculeAtTime(type, fGlobalTime);
    }
  }
  else
  {
    auto pScavengerMaterial = dynamic_cast<G4DNAScavengerMaterial*>(
      G4Scheduler::Instance()->GetScavengerMaterial());
    if(pScavengerMaterial != nullptr)
    {
      pScavengerMaterial->ReduceNumberMoleculePerVolumeUnitForMaterialConf(
        type, fGlobalTime);
    }
    else
    {
      G4ExceptionDescription exceptionDescription;
      exceptionDescription
        << "index : " << index << " " << type->GetName()
        << "  This molecule is not belong scavengers or particle-base"
        << G4endl;
      G4Exception("G4DNAEventScheduler::Stepping", "G4DNAEventScheduler002",
                  FatalErrorInArgument, exceptionDescription);
    }
  }
}

void G4DNAUpdateSystemModel::JumpTo(const Index& index, MolType type)
{
  auto& node = fpMesh->GetVoxelMapList(index);
  auto iter  = node.find(type);
  if(iter != node.end())
  {
    if(iter->second <= 0)
    {
      G4ExceptionDescription exceptionDescription;
      exceptionDescription << "G4DNAUpdateSystemModel::JumpTo::molecule : "
                           << type->GetName() << " index : " << index
                           << " number : " << iter->second;
      G4Exception("G4DNAUpdateSystemModel::JumpTo", "G4DNAUpdateSystemModel001",
                  FatalErrorInArgument, exceptionDescription);
    }
    iter->second--;
  }
  else
  {
    fpMesh->PrintVoxel(index);
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "index : " << index << " " << type->GetName()
                         << " There is no this type";
    G4Exception("G4DNAUpdateSystemModel::JumpTo", "G4DNAUpdateSystemModel002",
                FatalErrorInArgument, exceptionDescription);
  }
}

void G4DNAUpdateSystemModel::CreateMolecule(const Index& index, MolType type)
{
  // for scavenger
  auto pScavengerMaterial = dynamic_cast<G4DNAScavengerMaterial*>(
    G4Scheduler::Instance()->GetScavengerMaterial());
  if(pScavengerMaterial != nullptr && pScavengerMaterial->find(type))
  {
    pScavengerMaterial->AddNumberMoleculePerVolumeUnitForMaterialConf(
      type, fGlobalTime);
    return;
  }
  // for molecule
  auto& node = fpMesh->GetVoxelMapList(index);
  auto iter  = node.find(type);
  if(iter != node.end())
  {
    iter->second++;
  }
  else
  {
    node[type] = 1;
  }

  if(G4VMoleculeCounter::Instance()->InUse())
  {
    G4VMoleculeCounter::Instance()->AddAMoleculeAtTime(type, fGlobalTime);
  }
}

void G4DNAUpdateSystemModel::JumpIn(const Index& index, MolType type)
{
  // for molecule
  auto& node = fpMesh->GetVoxelMapList(index);
  auto iter  = node.find(type);
  if(iter != node.end())
  {
    iter->second++;
  }
  else
  {
    node[type] = 1;
  }
}

void G4DNAUpdateSystemModel::UpdateSystem(const Index& index,
                                          const ReactionData& data)
{
  auto reactant1 = data.GetReactant1();
  auto reactant2 = data.GetReactant2();
#ifdef G4VERBOSE
  if(fVerbose != 0)
  {
    G4cout << "At time : " << std::setw(7) << G4BestUnit(fGlobalTime, "Time")
           << " Reaction : " << reactant1->GetName() << " + "
           << reactant2->GetName() << " -> ";
  }
#endif
  const G4int nbProducts = data.GetNbProducts();
  if(nbProducts != 0)
  {
    for(G4int j = 0; j < nbProducts; ++j)
    {
#ifdef G4VERBOSE
      if((fVerbose != 0) && j != 0)
      {
        G4cout << " + ";
      }
      if(fVerbose != 0)
      {
        G4cout << data.GetProduct(j)->GetName();
      }
#endif
      CreateMolecule(index, data.GetProduct(j));
    }
  }
  else
  {
#ifdef G4VERBOSE
    if(fVerbose != 0)
    {
      G4cout << "No product";
    }
#endif
  }
#ifdef G4VERBOSE
  if(fVerbose != 0)
  {
    G4cout << G4endl;
  }
#endif
  KillMolecule(index, reactant1);
  KillMolecule(index, reactant2);
}

void G4DNAUpdateSystemModel::UpdateSystem(const Index& index,
                                          const JumpingData& data)
{
  auto reactant    = std::get<0>(data);
  auto JunpToIndex = std::get<1>(data);
#ifdef G4VERBOSE
  if(fVerbose > 1)
  {
    G4cout << "At time : " << std::setw(7) << G4BestUnit(fGlobalTime, "Time")
           << " Jumping : " << reactant->GetName() << " from " << index
           << " -> " << JunpToIndex << G4endl;
  }
#endif
  JumpTo(index, reactant);
  JumpIn(JunpToIndex, reactant);
}
