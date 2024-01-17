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
/// \file radiobiology/src/Manager.cc
/// \brief Implementation of the RadioBio::Manager class

#include "Manager.hh"

#include "G4AccumulableManager.hh"

#include <map>

#include <mutex>

std::mutex init_mutex;

namespace RadioBio
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#define width 15L

// To create the instance the DetectorConstruction must
// be passed (to account for correct voxelization)
Manager* Manager::CreateInstance()
{
  if (fInstance) {
    G4Exception("RadioBioManager::createInstance", "RecreatingRadioBioManager", FatalException,
                "Creating another, new, instance of RadioBioManager");
    delete fInstance;
  }

  fInstance = new Manager();
  return fInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Manager* Manager::GetInstance()
{
  return fInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Destructor deletes all the quantities
Manager::~Manager()
{
  for (auto q : fQuantities)
    delete q.second;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Manager::InitializeAll()
{

  std::lock_guard<std::mutex> lock(init_mutex);
  for (auto const& q : fQuantities)
    (q.second)->Initialize();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Manager::ComputeAll()
{
  for (auto const& q : fQuantities)
    (q.second)->Compute();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Manager::ResetAll()
{
  for (auto const& q : fQuantities)
    (q.second)->Reset();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Manager::StoreAll()
{
  for (auto const& q : fQuantities)
    (q.second)->Store();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

VRadiobiologicalQuantity* Manager::GetQuantity(G4String str)
{
  return fQuantities.find(str)->second;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Manager::PrintParameters()
{
  G4cout << "*******************************************" << G4endl
         << "*** right now registered quantities are ***" << G4endl;
  for (auto const& q : fQuantities)
    G4cout << "*** " << q.first;
  G4cout << "*** but their calculation might be not ****" << G4endl
         << "*** active. Ask for parameters of each ****" << G4endl
         << "*******************************************" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Manager::DigestAccumulables()
{
  for (auto q : fQuantities) {
    // Hook in the accumulable manager the one named as the quantity (eg "Dose")
    G4VAccumulable* GenAcc = G4AccumulableManager::Instance()->GetAccumulable(q.first);

    if (!GenAcc) {
      G4Exception("RadioBioManager::AddFromAccumulable", "NoAccumulable", FatalException, q.first);
    }

    // If calculation is not set enabled, exit
    if (!q.second->IsCalculationEnabled()) continue;

    // Add from the accumulable.
    q.second->AddFromAccumulable(GenAcc);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool Manager::Register(VRadiobiologicalQuantity* q, G4String name)
{
  if (q == nullptr) {
    G4Exception("RadioBioManager::Register", "RegisteringNullptr", JustWarning,
                "Asking to register a quantity with null pointer!");
    return false;
  }

  if (fQuantities.find(name) != fQuantities.end()) {
    G4Exception("RadioBioManager::Register", "RegisteringSameQuantity", FatalException,
                "Registering two radiobiological quantities with the same name!");
    return false;
  }
  fQuantities[name] = q;
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace RadioBio