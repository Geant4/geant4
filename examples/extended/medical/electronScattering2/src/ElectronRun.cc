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
/// \file medical/electronScattering2/src/ElectronRun.cc
/// \brief Implementation of the ElectronRun class

#include "ElectronRun.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4SDManager.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4SystemOfUnits.hh"
#include <assert.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ElectronRun::ElectronRun()
: G4Run(), fMap()
{
    fMap[0] = new G4THitsMap<G4double>("MyDetector", "cell flux");
    fMap[1] = new G4THitsMap<G4double>("MyDetector", "e cell flux");
    fMap[2] = new G4THitsMap<G4double>("MyDetector", "population");
    fMap[3] = new G4THitsMap<G4double>("MyDetector", "e population");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ElectronRun::~ElectronRun()
{
    // Important to clean up the map
    std::map<G4int, G4THitsMap<G4double>* >::iterator iter = fMap.begin();
    
    while (iter != fMap.end()) {
        delete iter->second;
        iter++;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ElectronRun::RecordEvent(const G4Event* anEvent)
{
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
    
    G4Run::RecordEvent(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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
        
        // The built-in scorers do not automatically account for the area of the
        // cylinder replica rings. We must account for this now by multiplying our result
        // by the ratio of the area of the full cylinder end over the area of the actual
        // scoring ring.
        // In this ratio, PI divides out, as does the width of the scoring rings.
        // Left with just the number of rings squared divided by the ring index plus
        // 1 squared minus ring index squared.
        G4int ringNum = iter->first;
        G4double areaCorrection = 233.*233. /
        ( (ringNum+1)*(ringNum+1) - ringNum*ringNum );
        G4int counter = 0;
        
        while (energyBinIter != iter->second.end()) {
            G4double value = *energyBinIter;
            if (counter < 2) value = value*areaCorrection;
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ElectronRun::Merge(const G4Run* aRun)
{
    // This method is called at the end of the run for each worker thread.
    // It accumulates the worker's results into global results.
    const ElectronRun* localRun = static_cast<const ElectronRun*>(aRun);
    const std::map< G4int, G4THitsMap<G4double>* >& localMap = localRun->fMap;
    std::map< G4int, G4THitsMap<G4double>* >::const_iterator iter = localMap.begin();
    for ( ; iter != localMap.end() ; ++iter)
        (*(fMap[iter->first])) += (*(iter->second));
    
    // This call lets Geant4 maintain overall summary information.
    G4Run::Merge(aRun);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
