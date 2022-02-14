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
  , fpMesh(nullptr)
  , fVerbose(0)
  , fGlobalTime(DBL_MAX)
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
      G4cout << "G4DNAUpdateSystemModel::KillMolecule::molecule : "
             << type->GetName() << " index : " << index
             << " number : " << iter->second << G4endl;
      assert(false);
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
      G4cout << "index : " << index << " " << type->GetName() << G4endl;
      G4cout << "This molecule is not belong scavengers or particle-base"
             << G4endl;
      assert(false);
    }
  }
}

void G4DNAUpdateSystemModel::JumpTo(const Index& index, MolType type)
{
  auto& node = fpMesh->GetVoxelMapList(index);

  auto iter = node.find(type);
  if(iter != node.end())
  {
    if(iter->second <= 0)
    {
      G4cout << "G4DNAUpdateSystemModel::KillMolecule::molecule : "
             << type->GetName() << " index : " << index
             << " number : " << iter->second << G4endl;
      assert(false);
    }
    iter->second--;
  }
  else
  {
    G4cout << "index : " << index << " " << type->GetName() << G4endl;
    G4cout << "This molecule is not belong  particle-base" << G4endl;
    assert(false);
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
    for(size_t j = 0; j < (size_t) nbProducts; ++j)
    {
#ifdef G4VERBOSE
      if((fVerbose != 0) && j != 0)
      {
        G4cout << " + ";
      }
      if(fVerbose != 0)
      {
        G4cout << data.GetProduct(j)->GetName();
        // for test
        // G4cout<<"  fGlobalTime : "<<fGlobalTime;
        // end fortest
      }
#endif
      CreateMolecule(index, data.GetProduct(j));
      //#define DEBUG 1

#ifdef DEBUG
      if(G4MoleculeCounter::Instance()->InUse())
        if(fpMesh->GetNumberOfType(data.GetProduct(j)) !=
           G4MoleculeCounter::Instance()->GetCurrentNumberOf(
             data.GetProduct(j)))
        {
          G4cout << "*********G4DNAUpdateSystemModel::DEBUG::GetNumberOfType("
                 << data.GetProduct(j)->GetName()
                 << ") : " << fpMesh->GetNumberOfType(data.GetProduct(j))
                 << G4endl;
          G4cout << "G4MoleculeCounter::GetCurrentNumberOf ("
                 << data.GetProduct(j)->GetName() << ") : "
                 << G4MoleculeCounter::Instance()->GetCurrentNumberOf(
                      data.GetProduct(j))
                 << G4endl;
          G4MoleculeCounter::Instance()->Dump();
          throw;
        }

#endif
    }
  }
  else
  {
#ifdef G4VERBOSE
    if(fVerbose != 0)
    {
      G4cout << "No product";
      // for test
      // G4cout<<"  fGlobalTime : "<<fGlobalTime;
      // end fortest
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
#ifdef DEBUG
  if(G4MoleculeCounter::Instance()->InUse())
    if(fpMesh->GetNumberOfType(reactant1) !=
       G4MoleculeCounter::Instance()->GetCurrentNumberOf(reactant1))
    {
      G4cout << "*********G4DNAUpdateSystemModel::DEBUG::GetNumberOfType("
             << reactant1->GetName()
             << ") : " << fpMesh->GetNumberOfType(reactant1) << G4endl;
      G4cout << "G4MoleculeCounter::GetCurrentNumberOf ("
             << reactant1->GetName() << ") : "
             << G4MoleculeCounter::Instance()->GetCurrentNumberOf(reactant1)
             << G4endl;
      G4MoleculeCounter::Instance()->Dump();
      throw;
    }
#endif
  KillMolecule(index, reactant2);
#ifdef DEBUG

  if(G4MoleculeCounter::Instance()->InUse())
    if(fpMesh->GetNumberOfType(reactant2) !=
       G4MoleculeCounter::Instance()->GetCurrentNumberOf(reactant2))
    {
      G4cout << "*********G4DNAUpdateSystemModel::DEBUG::GetNumberOfType("
             << reactant2->GetName()
             << ") : " << fpMesh->GetNumberOfType(reactant2) << G4endl;
      G4cout << "G4MoleculeCounter::GetCurrentNumberOf ("
             << reactant2->GetName() << ") : "
             << G4MoleculeCounter::Instance()->GetCurrentNumberOf(reactant2)
             << G4endl;
      G4MoleculeCounter::Instance()->Dump();
      throw;
    }
#endif
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