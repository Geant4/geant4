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
// $Id: Run.cc,v 1.1 2007-05-24 21:57:03 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// Jane Tinslay, May 2007. Creation, based on 
//                         BeamTestRun by T. Aso
//

#include "Run.hh"

#include "ConfigData.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4SDManager.hh"
#include "G4VPrimitiveScorer.hh"


Run::Run()
{
  G4MultiFunctionalDetector* detector = 
    dynamic_cast<G4MultiFunctionalDetector*>(G4SDManager::GetSDMpointer()->FindSensitiveDetector(ConfigData::GetDetectorName()));
  
  assert (0 != detector);
  
  // Loop over the registered primitive scorers.
  for (G4int icol=0; icol<detector->GetNumberOfPrimitives(); ++icol) {
    
    // Get Primitive Scorer object.
    G4VPrimitiveScorer* scorer = detector->GetPrimitive(icol);
    
    // collection name and collectionID for HitsCollection,
    // where type of HitsCollection is G4THitsMap in case of primitive scorer.
    // The collection name is given by <MFD name>/<Primitive Scorer name>.
    G4String collectionName = scorer->GetName();
    G4String fullCollectionName = ConfigData::GetDetectorName() + "/" + collectionName;
    G4int    collectionID = G4SDManager::GetSDMpointer()->GetCollectionID(fullCollectionName);
    
    assert (collectionID >= 0);
    
    // Store obtained HitsCollection information into data members.
    // And, creates new G4THitsMap for accumulating quantities during RUN.
    theCollName.push_back(fullCollectionName);
    theCollID.push_back(collectionID);
    theRunMap.push_back(new G4THitsMap<G4double>(ConfigData::GetDetectorName(), collectionName));
  }
}

Run::~Run()
{
  // Cleanup
  G4int Nmap = theRunMap.size();
  for ( G4int i = 0; i < Nmap; i++){
    if(theRunMap[i] ) theRunMap[i]->clear();
  }
  theCollName.clear();
  theCollID.clear();
  theRunMap.clear();
}

void Run::RecordEvent(const G4Event* aEvent)
{
  numberOfEvent++; 

  G4HCofThisEvent* HCE = aEvent->GetHCofThisEvent();
  if (!HCE) return;


  G4int Ncol = theCollID.size();
  for (G4int i=0; i<Ncol; ++i) {
    assert (theCollID[i] >= 0 );

    G4THitsMap<G4double>* EvtMap = static_cast<G4THitsMap<G4double>*>(HCE->GetHC(theCollID[i]));
    *theRunMap[i] += *EvtMap;
  }
}
