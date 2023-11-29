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
// G4IStore implementation
//
// Author: Michael Dressel (CERN), 2002
// Modified: Alex Howard (CERN), 2013 - Changed class to a 'singleton'
// ----------------------------------------------------------------------

#include "G4IStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4GeometryCell.hh"
#include "G4GeometryCellStepStream.hh"
#include "G4LogicalVolume.hh"
#include "G4TransportationManager.hh"

#include "G4AutoLock.hh"

namespace
{
  G4Mutex IStoreMutex = G4MUTEX_INITIALIZER;
}

// ***************************************************************************
// Static class variable: ptr to single instance of class
// ***************************************************************************
G4ThreadLocal G4IStore* G4IStore::fInstance = nullptr;

G4IStore::G4IStore()
  : fWorldVolume(G4TransportationManager::GetTransportationManager()
                 ->GetNavigatorForTracking()->GetWorldVolume())
{
}

G4IStore::G4IStore(const G4String& ParallelWorldName)
  : fWorldVolume(G4TransportationManager::GetTransportationManager()
                 ->GetParallelWorld(ParallelWorldName))
{
#ifdef G4VERBOSE
  G4cout << " G4IStore:: ParallelWorldName = "
         << ParallelWorldName << G4endl;
  G4cout << " G4IStore:: fParallelWorldVolume = "
         << fWorldVolume->GetName() << G4endl;
#endif
}

G4IStore::~G4IStore() = default;

void G4IStore::Clear()
{
  fGeometryCelli.clear();
}

void G4IStore::SetWorldVolume()
{
  G4cout << " G4IStore:: SetWorldVolume " << G4endl;
  fWorldVolume = G4TransportationManager::GetTransportationManager()
                 ->GetNavigatorForTracking()->GetWorldVolume();
  G4cout << " World volume is: " << fWorldVolume->GetName() << G4endl;
  // fGeometryCelli = new G4GeometryCellImportance;
}

void G4IStore::SetParallelWorldVolume(const G4String& paraName)
{
  G4cout << " G4IStore:: SetParallelWorldVolume " << G4endl;
  fWorldVolume = G4TransportationManager::GetTransportationManager()
                 ->GetParallelWorld(paraName);
  G4cout << " ParallelWorld volume is: " << fWorldVolume->GetName() << G4endl;
  // fGeometryCelli = new G4GeometryCellImportance;
}

const G4VPhysicalVolume& G4IStore::GetWorldVolume() const
{
  return *fWorldVolume;
}

const G4VPhysicalVolume* G4IStore::GetParallelWorldVolumePointer() const
{
  return fWorldVolume;
}

void G4IStore::SetInternalIterator(const G4GeometryCell& gCell) const
{
  fCurrentIterator = fGeometryCelli.find(gCell);
}

void G4IStore::AddImportanceGeometryCell(G4double importance,
                                         const G4GeometryCell& gCell)
{
  if (importance < 0 )
  {
    Error("AddImportanceGeometryCell() - Invalid importance value given.");
  }  
  if (!IsInWorld(gCell.GetPhysicalVolume()) )
  {
    Error("AddImportanceGeometryCell() - Physical volume not found!");
  }
  SetInternalIterator(gCell);
  if (fCurrentIterator != fGeometryCelli.cend())
  {
    Error("AddImportanceGeometryCell() - Region already existing!");
  }
  fGeometryCelli[gCell] = importance;
}

void G4IStore::AddImportanceGeometryCell(G4double importance,
                                         const G4VPhysicalVolume& aVolume,
                                         G4int aRepNum)
{
  AddImportanceGeometryCell(importance, G4GeometryCell(aVolume, aRepNum));
}

void G4IStore::ChangeImportance(G4double importance,
                                const G4GeometryCell& gCell)
{
  if (importance < 0 )
  {
    Error("ChangeImportance() - Invalid importance value given.");
  }
  if (!IsInWorld(gCell.GetPhysicalVolume()))
  {
    Error("ChangeImportance() - Physical volume not found!");
  }
  SetInternalIterator(gCell);
  if (fCurrentIterator == fGeometryCelli.cend())
  {
    Error("ChangeImportance() - Region does not exist!");
  }
  fGeometryCelli[gCell] = importance;

}

void G4IStore::ChangeImportance(G4double importance,
                                const G4VPhysicalVolume& aVolume,
                                G4int aRepNum)
{
  ChangeImportance(importance, G4GeometryCell(aVolume, aRepNum));
}

G4double G4IStore::GetImportance(const G4VPhysicalVolume& aVolume,
                                 G4int aRepNum) const
{  
  G4AutoLock l(&IStoreMutex);
  SetInternalIterator(G4GeometryCell(aVolume, aRepNum));
  auto gCellIterator = fCurrentIterator;
  if (gCellIterator == fGeometryCelli.cend())
  {
    Error("GetImportance() - Region does not exist!");
    return 0.;
  }
  G4double importance_value = (*fCurrentIterator).second;
  l.unlock();

  return importance_value;
}

G4double G4IStore::GetImportance(const G4GeometryCell& gCell) const
{
  G4AutoLock l(&IStoreMutex);
  SetInternalIterator(gCell);
  auto gCellIterator = fCurrentIterator;
  if (gCellIterator == fGeometryCelli.cend())
  {
    std::ostringstream err_mess;
    err_mess << "GetImportance() - Region does not exist!" << G4endl
             << "Geometry cell, " << gCell
             << ", not found in: " << fGeometryCelli << ".";
    Error(err_mess.str());
    return 0.;
  }
  G4double importance_value = (*fCurrentIterator).second;
  l.unlock();

  return importance_value;
}

G4bool G4IStore::IsKnown(const G4GeometryCell& gCell) const
{
  G4AutoLock l(&IStoreMutex);
  G4bool inWorldKnown(IsInWorld(gCell.GetPhysicalVolume()));
                      
  if ( inWorldKnown )
  {
    SetInternalIterator(gCell);
    inWorldKnown = (fCurrentIterator != fGeometryCelli.cend());
  }
  l.unlock();

  return inWorldKnown;
}

G4bool G4IStore::IsInWorld(const G4VPhysicalVolume& aVolume) const
{
  G4bool isIn(true);
  if (!(aVolume == *fWorldVolume))
  {
    isIn = fWorldVolume->GetLogicalVolume()->IsAncestor(&aVolume);
  }
  return isIn;
}

void G4IStore::Error(const G4String& msg) const
{
  G4Exception("G4IStore::Error()", "GeomBias0002", FatalException, msg);
}

// ***************************************************************************
// Returns the instance of the singleton.
// Creates it in case it's called for the first time.
// ***************************************************************************
//
G4IStore* G4IStore::GetInstance()
{
  if (fInstance == nullptr)
  {
#ifdef G4VERBOSE
    G4cout << "G4IStore:: Creating new MASS IStore " << G4endl;
#endif
    fInstance = new G4IStore();
  }
  return fInstance;    
}

// ***************************************************************************
// Returns the instance of the singleton.
// Creates it in case it's called for the first time.
// ***************************************************************************
//
G4IStore* G4IStore::GetInstance(const G4String& ParallelWorldName)
{
  if (fInstance == nullptr)
  {
#ifdef G4VERBOSE
    G4cout << "G4IStore:: Creating new Parallel IStore "
           << ParallelWorldName << G4endl;
#endif
    fInstance = new G4IStore(ParallelWorldName);
  }
  return fInstance;    
}
