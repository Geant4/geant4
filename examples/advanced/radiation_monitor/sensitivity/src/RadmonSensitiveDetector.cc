//
// File name:     RadmonSensitiveDetector.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonSensitiveDetector.cc,v 1.3 2006-01-06 12:52:32 guatelli Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonSensitiveDetector.hh"
#include "RadmonSensitiveDetectorDataStorer.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"

#define RADMONSENSITIVEDETECTOR_COLLECTIONNAME "radmon_default_"

RadmonSensitiveDetector :: RadmonSensitiveDetector(const G4String & name)
  :
  G4VSensitiveDetector(name),
  hitsCollection(0)
{
  collName=name;
  collName+=RADMONSENSITIVEDETECTOR_COLLECTIONNAME;
  collectionName.insert(collName);
}
 
RadmonSensitiveDetector :: ~RadmonSensitiveDetector()
{
}
 
void RadmonSensitiveDetector :: ClearDataStorersList()
{
  dataStorersSet.clear();
}
 
void RadmonSensitiveDetector :: AttachDataStorer(RadmonSensitiveDetectorDataStorer * observer)
{
  dataStorersSet.insert(observer);
}
 
void RadmonSensitiveDetector :: Initialize(G4HCofThisEvent * hitsCollections)
{
  hitsCollection = new RadmonHitsCollection(GetName(), collName);
  hitsCollections -> AddHitsCollection(GetCollectionID(0), hitsCollection);
}
 
G4bool RadmonSensitiveDetector :: ProcessHits(G4Step * step, G4TouchableHistory * /* touchableHistory */)
{
  if (dataStorersSet.empty())
    return false;
  
  RadmonHit * hit(new RadmonHit);
 
  hitsCollection->insert(hit);

  DataStorersSet::iterator i(dataStorersSet.begin());
  const DataStorersSet::iterator end(dataStorersSet.end());
 
  while (i!=end)
    {
      (*i)->StoreIntoHit(step, hit);
      i++;
    }

  return true;
}

RadmonHitsCollection* RadmonSensitiveDetector::GetDetectorCollection(void) const
{
  G4SDManager * sdManager(G4SDManager::GetSDMpointer());
 
  if (!sdManager)
    return 0;
 
  G4RunManager * rManager(G4RunManager::GetRunManager());
 
  if (!rManager)
    return 0;
 
  const G4Event * currentEvent(rManager->GetCurrentEvent());
  G4HCofThisEvent * hitsCollections(currentEvent->GetHCofThisEvent());
  return reinterpret_cast<RadmonHitsCollection *>(hitsCollections->GetHC(sdManager->GetCollectionID(collName)));
}

