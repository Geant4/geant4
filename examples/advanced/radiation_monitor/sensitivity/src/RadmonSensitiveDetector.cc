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
// File name:     RadmonSensitiveDetector.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonSensitiveDetector.cc,v 1.2.2.2 2006/06/29 16:20:34 gunter Exp $
// Tag:           $Name: geant4-09-01 $
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
 hitsCollection=new RadmonHitsCollection(GetName(), collName);
 hitsCollections->AddHitsCollection(GetCollectionID(0), hitsCollection);
}
 
 
 
G4bool                                          RadmonSensitiveDetector :: ProcessHits(G4Step * step, G4TouchableHistory * /* touchableHistory */)
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





RadmonHitsCollection *                          RadmonSensitiveDetector :: GetDetectorCollection(void) const
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

