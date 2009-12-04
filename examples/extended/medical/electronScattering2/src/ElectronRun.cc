//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//

#include "ElectronRun.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4SDManager.hh"
#include "G4VPrimitiveScorer.hh"
#include <assert.h>

ElectronRun::ElectronRun(const G4String& detectorName)
{
  // Get the sensitive detector manager
  G4SDManager* manager = G4SDManager::GetSDMpointer();

  // Get the sensitive detector
  G4MultiFunctionalDetector* detector =
    dynamic_cast<G4MultiFunctionalDetector*>(manager->FindSensitiveDetector(detectorName));

  // Abort if detector doesn't exist
  assert (0 != detector);
  
  G4int i(0);

  // Loop over primitive scorers registered with the detector
  for (i=0; i<detector->GetNumberOfPrimitives(); i++) {
    
    // Get scorer
    G4VPrimitiveScorer* scorer = detector->GetPrimitive(i);

    // Need to form the full collection name = detector name + "/"+ scorer name 
    // to get the collection id number
    G4String fullCollectionName = detectorName+"/"+scorer->GetName();

    G4int id = manager->GetCollectionID(fullCollectionName);

    // Abort if the collection was not added to the sensitive detector manager
    // when the scorer was registered with G4MultiFunctionalDetector
    assert (id >= 0);
    
    fMap[id] = new G4THitsMap<G4double>(detectorName, scorer->GetName());
  }
}

ElectronRun::~ElectronRun()
{
  // Important to clean up the map
  std::map<G4int, G4THitsMap<G4double>* >::iterator iter = fMap.begin();

  while (iter != fMap.end()) {
    delete iter->second;			
    iter++;
  }
}

void ElectronRun::RecordEvent(const G4Event* anEvent)
{
  numberOfEvent++;

  // Get the hits collection
  G4HCofThisEvent* eventHitCollection = anEvent->GetHCofThisEvent();
 
  if (!eventHitCollection) return;
 
  // Update our private fMap
  std::map< G4int, G4THitsMap<G4double>* >::iterator iter = fMap.begin();
 
  while (iter != fMap.end()) {
    G4int id = iter->first;
   
    // Get the hit collection corresponding to "id"
    G4THitsMap<G4double>* eventHitsMap
      = dynamic_cast< G4THitsMap<G4double>* >(eventHitCollection->GetHC(id));
 
    // Expect this to exist
    assert (0 != eventHitsMap);
 
    // Accumulate event data into our G4THitsMap<G4double> map
    *(iter->second) += *eventHitsMap;
 
    iter++;
  }
}

void ElectronRun::DumpData(G4String &outputFileSpec) const
{
  // Titles
  std::vector<G4String> title;
  title.push_back("Radius");

  // Output map - energy binning on x axis, theta on y
  std::map< G4int, std::vector<G4double> > output;

  G4int nThetaBins = 233;

  // Energy bins depends on the number of scorers
  G4int nEnergyBins = fMap.size();

  G4int i(0), j(0);

  // Initialise current to 0 in all bins
  for (i=0; i<nThetaBins; i++) {
    for (j=0; j<nEnergyBins; j++) {
      output[i].push_back(0);
    }
  }
  
  i=0;

  // Fill the output map
  std::map< G4int, G4THitsMap<G4double>* >::const_iterator iter = fMap.begin();

  while (iter != fMap.end()) {
    G4THitsMap<G4double>* hitMap = iter->second;
    
    title.push_back(hitMap->GetName());

    std::map<G4int,G4double*>* myMap = hitMap->GetMap();

    for (j=0; j<nThetaBins; j++) {
      G4double* current = (*myMap)[j];
      if (0 != current) output[j][i] = (*current);
    }
    
    i++;
    iter++;
  }
  
  Print(title, output, outputFileSpec);
}

void ElectronRun::Print(const std::vector<G4String>& title,
			const std::map< G4int, std::vector<G4double> >&myMap,
			G4String &outputFileSpec) const
{
  // Print to G4cout and an output file
  std::ofstream outFile(outputFileSpec);

  // Print title vector
  std::vector<G4String>::const_iterator titleIter = title.begin();
  
  while (titleIter != title.end()) {
    G4cout << std::setw(8)<<*titleIter<<" ";
    titleIter++;
  }
  
  G4cout<<G4endl;

  // Print current data
  std::map< G4int, std::vector<G4double> >::const_iterator iter = myMap.begin();
  
  while (iter != myMap.end()) {
    G4cout << std::setw(8)<<std::setprecision(3)<< iter->first<<" ";
    
    std::vector<G4double>::const_iterator energyBinIter = iter->second.begin();
    
	G4int ringNum = iter->first;
	G4double area = (ringNum+1)*(ringNum+1) - ringNum*ringNum; 
    G4int counter = 0;
 
    while (energyBinIter != iter->second.end()) {
      G4double value = *energyBinIter;
      if (counter < 2) value = value/area;
	  G4cout << std::setw(10)<<std::setprecision(5)<< value*mm*mm<<" ";
      outFile << value*mm*mm;
      if (counter < 3) outFile <<",";
      counter++;
      energyBinIter++;
    }
    outFile<<G4endl;
    G4cout<<G4endl;
    iter++;
  }
}
