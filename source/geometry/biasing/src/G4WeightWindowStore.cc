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
// $Id: G4WeightWindowStore.cc 88801 2015-03-10 14:38:27Z gcosmo $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4WeightWindowStore
//
// ----------------------------------------------------------------------


#include "G4WeightWindowStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4GeometryCellStepStream.hh"
#include "G4TransportationManager.hh"

// ***************************************************************************
// Static class variable: ptr to single instance of class
// ***************************************************************************
G4ThreadLocal G4WeightWindowStore* G4WeightWindowStore::fInstance = 0;

G4WeightWindowStore::
G4WeightWindowStore() :
fWorldVolume(G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->GetWorldVolume()),
  fGeneralUpperEnergyBounds(),
  fCellToUpEnBoundLoWePairsMap(),
  fCurrentIterator(fCellToUpEnBoundLoWePairsMap.end())
{}

G4WeightWindowStore::
G4WeightWindowStore(const G4String& ParallelWorldName) :
fWorldVolume(G4TransportationManager::GetTransportationManager()->GetParallelWorld(ParallelWorldName)),
  fGeneralUpperEnergyBounds(),
  fCellToUpEnBoundLoWePairsMap(),
  fCurrentIterator(fCellToUpEnBoundLoWePairsMap.end())
{}

G4WeightWindowStore::~G4WeightWindowStore()
{}


G4double G4WeightWindowStore::GetLowerWeight(const G4GeometryCell &gCell, 
					     G4double partEnergy) const
{
  SetInternalIterator(gCell);
  G4GeometryCellWeight::const_iterator gCellIterator = fCurrentIterator;
  if (gCellIterator ==  fCellToUpEnBoundLoWePairsMap.end()) {
    Error("GetLowerWitgh() - Cell does not exist!");
    return 0.;
  }
  G4UpperEnergyToLowerWeightMap upEnLoWeiPairs =
    fCurrentIterator->second;
  G4double lowerWeight = -1;
  G4bool found = false;
  for (G4UpperEnergyToLowerWeightMap::iterator it = 
	 upEnLoWeiPairs.begin(); it != upEnLoWeiPairs.end(); it++) {
    if (partEnergy < it->first) {
      lowerWeight = it->second;
      found = true;
      break;
    }
  }
  if (!found) {
    std::ostringstream err_mess;
    err_mess << "GetLowerWitgh() - Couldn't find lower weight bound." << G4endl
             << "Energy: " << partEnergy << ".";
    Error(err_mess.str());
  }
  return lowerWeight;


}

void G4WeightWindowStore::
SetInternalIterator(const G4GeometryCell &gCell) const
{
  fCurrentIterator = fCellToUpEnBoundLoWePairsMap.find(gCell);
}

G4bool G4WeightWindowStore::
IsInWorld(const G4VPhysicalVolume &aVolume) const
{
  G4bool isIn(true);
  if (!(aVolume == *fWorldVolume)) {
    isIn = fWorldVolume->GetLogicalVolume()->IsAncestor(&aVolume);
  }
  return isIn;
}


G4bool G4WeightWindowStore::IsKnown(const G4GeometryCell &gCell) const
{
  G4bool inWorldKnown(IsInWorld(gCell.GetPhysicalVolume()));
		      
  if ( inWorldKnown ) {
    SetInternalIterator(gCell);
    inWorldKnown = (fCurrentIterator!=fCellToUpEnBoundLoWePairsMap.end());
  }
  return inWorldKnown;
}

void G4WeightWindowStore::Clear()
{
  fCellToUpEnBoundLoWePairsMap.clear();
}

void G4WeightWindowStore::SetWorldVolume()
{
  G4cout << " G4IStore:: SetWorldVolume " << G4endl;
  fWorldVolume = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->GetWorldVolume();
  G4cout << " World volume is: " << fWorldVolume->GetName() << G4endl;
  //  fGeometryCelli = new G4GeometryCellImportance;
}

void G4WeightWindowStore::SetParallelWorldVolume(G4String paraName)
{
  fWorldVolume = G4TransportationManager::GetTransportationManager()->GetParallelWorld(paraName);
    //  fGeometryCelli = new G4GeometryCellImportance;
}


const G4VPhysicalVolume &G4WeightWindowStore::GetWorldVolume() const
{
  return *fWorldVolume;
}

const G4VPhysicalVolume* G4WeightWindowStore::GetParallelWorldVolumePointer() const
{
  return fWorldVolume;
}



void G4WeightWindowStore::
AddLowerWeights(const G4GeometryCell & gCell,
                const std::vector<G4double> &lowerWeights)
{
  if (fGeneralUpperEnergyBounds.empty()) {
    Error("AddLowerWeights() - No general upper energy limits set!");
  }
  if (IsKnown(gCell)) {
    Error("AddLowerWeights() - Cell already in the store.");
  }
  if (lowerWeights.size() != fGeneralUpperEnergyBounds.size()) {
    std::ostringstream err_mess;
    err_mess << "AddLowerWeights() - Mismatch between "
             << "number of lower weights (" << lowerWeights.size()
             << ") and energy bounds (" << fGeneralUpperEnergyBounds.size()
             << ")!";
    Error(err_mess.str());
  }
  G4UpperEnergyToLowerWeightMap map;
  G4int i = 0;
  for (std::set<G4double, std::less<G4double> >::iterator it = 
	 fGeneralUpperEnergyBounds.begin(); 
       it != fGeneralUpperEnergyBounds.end();
       it++) {
    map[*it] = lowerWeights[i];
    i++;
  }
  fCellToUpEnBoundLoWePairsMap[gCell] = map;
}

 
void G4WeightWindowStore::
AddUpperEboundLowerWeightPairs(const G4GeometryCell &gCell,
                               const G4UpperEnergyToLowerWeightMap& enWeMap)
{
  if (IsKnown(gCell)) {
    Error("AddUpperEboundLowerWeightPairs() - Cell already in the store.");
  }
  if (IsKnown(gCell)) {
    Error("AddUpperEboundLowerWeightPairs() - Cell already in the store.");
  }
  fCellToUpEnBoundLoWePairsMap[gCell] = enWeMap;

}


void G4WeightWindowStore::
SetGeneralUpperEnergyBounds(const std::set<G4double,
                            std::less<G4double> > &enBounds)
{
  if (!fGeneralUpperEnergyBounds.empty()) {
    Error("SetGeneralUpperEnergyBounds() - Energy bounds already set.");
  }
  fGeneralUpperEnergyBounds = enBounds;
}

  
void G4WeightWindowStore::Error(const G4String &msg) const
{
  G4Exception("G4WeightWindowStore::Error()",
              "GeomBias0002", FatalException, msg);
}


// ***************************************************************************
// Returns the instance of the singleton.
// Creates it in case it's called for the first time.
// ***************************************************************************
//
G4WeightWindowStore* G4WeightWindowStore::GetInstance()
{
  if (!fInstance)
  {
    fInstance = new G4WeightWindowStore();
  }
  return fInstance;    
}

// ***************************************************************************
// Returns the instance of the singleton.
// Creates it in case it's called for the first time.
// ***************************************************************************
//
G4WeightWindowStore* G4WeightWindowStore::GetInstance(const G4String& ParallelWorldName)
{
  if (!fInstance)
  {
    G4cout << "G4IStore:: Creating new Parallel IStore " << ParallelWorldName << G4endl;
    fInstance = new G4WeightWindowStore(ParallelWorldName);
  }
  return fInstance;    
}

