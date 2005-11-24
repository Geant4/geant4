//
// File name:     RadmonSensitiveDetector.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonSensitiveDetector.cc,v 1.1 2005-11-24 02:31:47 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonSensitiveDetector.hh"
#include "RadmonSensitiveDetectorDataStorer.hh"

#define RADMONSENSITIVEDETECTOR_COLLECTIONNAME "radmon_default"

                                                RadmonSensitiveDetector :: RadmonSensitiveDetector(const G4String & name)
:
 G4VSensitiveDetector(name),
 hitCollection(0)
{
 collectionName.insert(RADMONSENSITIVEDETECTOR_COLLECTIONNAME);
}
 
 
 
                                                RadmonSensitiveDetector :: ~RadmonSensitiveDetector()
{
}
 
 
 


void                                            RadmonSensitiveDetector :: ClearDataStorersList()
{
 dataStorersSet.clear();
}
 
 
 
void                                            RadmonSensitiveDetector :: AttachDataStorer(RadmonSensitiveDetectorDataStorer * observer)
{
 dataStorersSet.insert(observer);
}
 
 
 


void                                            RadmonSensitiveDetector :: Initialize(G4HCofThisEvent * hitsCollections)
{
 hitCollection=new RadmonHitsCollection(GetName(), RADMONSENSITIVEDETECTOR_COLLECTIONNAME);
 hitsCollections->AddHitsCollection(GetCollectionID(0), hitCollection);
}
 
 
 
G4bool                                          RadmonSensitiveDetector :: ProcessHits(G4Step * step, G4TouchableHistory * /* touchableHistory */)
{
 if (dataStorersSet.empty())
  return false;
  
 RadmonHit * hit(new RadmonHit);
 
 hitCollection->insert(hit);

 DataStorersSet::iterator i(dataStorersSet.begin());
 const DataStorersSet::iterator end(dataStorersSet.end());
 
 while (i!=end)
 {
  (*i)->StoreIntoHit(step, hit);
  i++;
 }

 return true;
}
