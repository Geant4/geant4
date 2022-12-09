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
// G4VUserDetectorConstruction implementation
//
// Original author: M.Asai, 1999
// --------------------------------------------------------------------

#include <sstream>
#include <map>
#include <assert.h>

#include "G4VUserDetectorConstruction.hh"
#include "G4FieldManager.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4MultiSensitiveDetector.hh"
#include "G4SDManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSensitiveDetector.hh"
#include "G4VUserParallelWorld.hh"
#include "G4RunManager.hh"

// --------------------------------------------------------------------
G4VUserDetectorConstruction::G4VUserDetectorConstruction()
{
}

// --------------------------------------------------------------------
G4VUserDetectorConstruction::~G4VUserDetectorConstruction()
{
}

// --------------------------------------------------------------------
void G4VUserDetectorConstruction::
RegisterParallelWorld(G4VUserParallelWorld* aPW)
{
  auto pwItr = std::find(parallelWorld.cbegin(), parallelWorld.cend(), aPW);
  if (pwItr != parallelWorld.cend())
  {
    G4String eM = "A parallel world <";
    eM += aPW->GetName();
    eM += "> is already registered to the user detector construction.";
    G4Exception("G4VUserDetectorConstruction::RegisterParallelWorld",
                "Run0051", FatalErrorInArgument, eM);
  }
  parallelWorld.push_back(aPW);
}

// --------------------------------------------------------------------
G4int G4VUserDetectorConstruction::ConstructParallelGeometries()
{
  G4int nP = 0;
  for(auto pwItr = parallelWorld.cbegin();
           pwItr != parallelWorld.cend(); ++pwItr)
  {
    (*pwItr)->Construct();
    ++nP;
  }
  return nP;
}

// --------------------------------------------------------------------
void G4VUserDetectorConstruction::ConstructParallelSD()
{
  for(auto pwItr = parallelWorld.cbegin();
           pwItr != parallelWorld.cend(); ++pwItr)
  {
    (*pwItr)->ConstructSD();
  }
}

// --------------------------------------------------------------------
G4int G4VUserDetectorConstruction::GetNumberOfParallelWorld() const
{
  return (G4int)parallelWorld.size();
}

// --------------------------------------------------------------------
G4VUserParallelWorld*
G4VUserDetectorConstruction::GetParallelWorld(G4int i) const
{
  if(i < 0 || i >= GetNumberOfParallelWorld())
    return nullptr;
  return parallelWorld[i];
}

// --------------------------------------------------------------------
void G4VUserDetectorConstruction::ConstructSDandField()
{
  // G4RunManager::RMType rmtype =
  // G4RunManager::GetRunManager()->GetRunManagerType(); if(rmtype !=
  // G4RunManager::sequentialRM)
  // {
  //  G4cout
  //  << "User-derived detector construction class does not implement \n"
  //  << "ConstructSDandFiled method: workers will not have SD and fields!\n"
  //  << "The user can safely ignore this message if (s)he has no sensitive\n"
  //  << "detector or field in her/his application." << G4endl;
  // }
}

// --------------------------------------------------------------------
void G4VUserDetectorConstruction::CloneF()
{
  using FMtoFMmap = std::map<G4FieldManager*, G4FieldManager*>;
  using FMpair = std::pair<G4FieldManager*, G4FieldManager*>;

  FMtoFMmap masterToWorker;
  G4LogicalVolumeStore* const logVolStore = G4LogicalVolumeStore::GetInstance();
  for(auto it = logVolStore->cbegin(); it != logVolStore->cend(); ++it)
  {
    G4LogicalVolume* g4LogicalVolume = *it;
    // Use shadow of master to get instance of FM
    G4FieldManager* masterFM = nullptr;  // g4LogicalVolume->fFieldManager;
    G4FieldManager* clonedFM = nullptr;
    if(masterFM != nullptr)
    {
      auto fmFound = masterToWorker.find(masterFM);
      if(fmFound == masterToWorker.cend())
      {
        // First time we see this SD, let's clone and remember...
        try
        {
          auto insertedEl =
               masterToWorker.insert(FMpair(masterFM, masterFM->Clone()));
          clonedFM = (insertedEl.first)->second;
        } catch(...)
        {
          G4ExceptionDescription msg;
          msg << "Cloning of G4FieldManager failed."
              << " But derived class does not implement cloning. Cannot "
                 "continue.";
          G4Exception("G4VUserDetectorConstruction::CloneSD()", "Run0053",
                      FatalException, msg);
        }
      }
      else
      {
        // We have already seen this SD attached to a different LogicalVolume,
        // let's re-use previous clone
        clonedFM = (*fmFound).second;
      }
    }  // masterFM != 0
    // Note that we do not push FM to daughters (false argument), however, since
    // we area looping on all logical volumes and we implemented the "trick" of
    // the map master<->cloned the final effect is the same as using here the
    // correct Boolean flag: log-volumes that originally were sharing the same
    // FM they will have cloned ones
    g4LogicalVolume->SetFieldManager(clonedFM, false);
  }
}

// --------------------------------------------------------------------
void G4VUserDetectorConstruction::CloneSD()
{
  // Loop on ALL logial volumes to search for attached SD
  G4LogicalVolumeStore* const logVolStore = G4LogicalVolumeStore::GetInstance();

  using SDtoSDmap = std::map<G4VSensitiveDetector*, G4VSensitiveDetector*>;
  using SDpair = std::pair<G4VSensitiveDetector*, G4VSensitiveDetector*>;
  SDtoSDmap masterToWorker;

  for(auto it = logVolStore->cbegin(); it != logVolStore->cend(); ++it)
  {
    G4LogicalVolume* g4LogicalVolume = *it;
    // Use shadow of master to get the instance of SD
    G4VSensitiveDetector* masterSD = nullptr;
    G4VSensitiveDetector* clonedSD = nullptr;
    if(masterSD != nullptr)
    {
      auto sdFound = masterToWorker.find(masterSD);
      if(sdFound == masterToWorker.cend())
      {
        // First time we see this SD, let's clone and remember...
        try
        {
          auto insertedEl =
               masterToWorker.insert(SDpair(masterSD, masterSD->Clone()));
          clonedSD = (insertedEl.first)->second;
        } catch(...)
        {
          G4ExceptionDescription msg;
          msg << "Cloning of G4VSensitiveDetector requested for:"
              << masterSD->GetName() << "\n"
#ifndef WIN32
              << " (full path name: " << masterSD->GetFullPathName() << ").\n"
#endif
              << " But derived class does not implement cloning. Cannot "
                 "continue.";
          G4Exception("G4VUserDetectorConstruction::CloneSD()", "Run0053",
                      FatalException, msg);
        }
      }
      else
      {
        // We have already seen this SD attached to a different LogicalVolume,
        // let's re-use previous clone
        clonedSD = (*sdFound).second;
      }
    }  // masterSD!=0
    g4LogicalVolume->SetSensitiveDetector(clonedSD);
  }
}

// --------------------------------------------------------------------
void G4VUserDetectorConstruction::
SetSensitiveDetector(const G4String& logVolName,
                     G4VSensitiveDetector* aSD, G4bool multi)
{
  G4bool found                = false;
  G4LogicalVolumeStore* store = G4LogicalVolumeStore::GetInstance();
  auto volmap = store->GetMap();
  auto pos = volmap.find(logVolName);
  if(pos != volmap.cend())
  {
    if ((pos->second.size()>1) && !multi)
    {
      G4String eM = "More than one logical volumes of name <";
      eM += pos->first;
      eM += "> are found and thus the sensitive detector <";
      eM += aSD->GetName();
      eM += "> cannot be uniquely assigned.";
      G4Exception("G4VUserDetectorConstruction::SetSensitiveDetector()",
                  "Run0052", FatalErrorInArgument, eM);
    }
    found = true;
    for (std::size_t i = 0; i < pos->second.size(); ++i)
    {
      SetSensitiveDetector(pos->second[i], aSD);
    }
  }
  if(!found)
  {
    G4String eM2 = "No logical volume of name <";
    eM2 += logVolName;
    eM2 += "> is found. The specified sensitive detector <";
    eM2 += aSD->GetName();
    eM2 += "> couldn't be assigned to any volume.";
    G4Exception("G4VUserDetectorConstruction::SetSensitiveDetector()",
                "Run0053", FatalErrorInArgument, eM2);
  }
}

// --------------------------------------------------------------------
void G4VUserDetectorConstruction::
SetSensitiveDetector(G4LogicalVolume* logVol, G4VSensitiveDetector* aSD)
{
  assert(logVol != nullptr && aSD != nullptr);

  // The aSD has already been added by user to the manager if needed
  // G4SDManager::GetSDMpointer()->AddNewDetector(aSD);

  // New Logic: allow for "multiple" SDs being attached to a single LV.
  // To do that we use a special proxy SD called G4MultiSensitiveDetector

  // Get existing SD if already set and check if it is of the special type
  G4VSensitiveDetector* originalSD = logVol->GetSensitiveDetector();
  if(originalSD == aSD)
  {
    G4ExceptionDescription msg;
    msg << "Attempting to add multiple times the same sensitive detector (\"";
    msg << originalSD->GetName() << "\") is not allowed, skipping.";
    G4Exception("G4VUserDetectorConstruction::SetSensitiveDetector", "Run0054",
                JustWarning, msg);
    return;
  }
  if(originalSD == nullptr)
  {
    logVol->SetSensitiveDetector(aSD);
  }
  else
  {
    G4MultiSensitiveDetector* msd =
      dynamic_cast<G4MultiSensitiveDetector*>(originalSD);
    if(msd != nullptr)
    {
      msd->AddSD(aSD);
    }
    else
    {
      std::ostringstream mn;
      mn << "/MultiSD_" << logVol->GetName() << "_" << logVol;
      const G4String msdname = mn.str();
      msd                    = new G4MultiSensitiveDetector(msdname);
      // We need to register the proxy to have correct handling of IDs
      G4SDManager::GetSDMpointer()->AddNewDetector(msd);
      msd->AddSD(originalSD);
      msd->AddSD(aSD);
      logVol->SetSensitiveDetector(msd);
    }
  }
}
