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
// $Id: G4FissionStore.cc,v 1.17 2010-12-15 07:41:05 gunter Exp $
//
// 20100728  Move ::addConfig() implementation to .cc file

#include "G4FissionStore.hh"
#include "G4FissionConfiguration.hh"
#include <cmath>

G4FissionStore::G4FissionStore() : verboseLevel(0) {
  if (verboseLevel > 1) 
    G4cout << " >>> G4FissionStore::G4FissionStore" << G4endl;
}

void G4FissionStore::addConfig(G4double a, G4double z, G4double ez, 
			       G4double ek, G4double ev) {
  G4FissionConfiguration config(a, z, ez, ek, ev);
  configurations.push_back(config);
  if (verboseLevel > 2) config.print();
}

G4FissionConfiguration G4FissionStore::generateConfiguration(G4double amax, 
							     G4double rand) const {
  if (verboseLevel > 1)
    G4cout << " >>> G4FissionStore::generateConfiguration" << G4endl;
  
  const G4double small = -30.0;

  G4double totProb = 0.0;
  std::vector<G4double> probs(size());

  if (verboseLevel > 3)
    G4cout << " amax " << amax << " ic " << size() << G4endl;

  for (G4int i = 0; i < size(); i++) {
    G4double ez = configurations[i].ezet;
    G4double pr = ez - amax;

    if (pr < small) pr = small;
    pr = std::exp(pr); 
    //  configurations[i].print();
    //  G4cout << " probability " << pr << G4endl; 
    totProb += pr;
    probs[i] = totProb;  
  };

  G4double st = totProb * rand;
  G4int igen = 0;

  while (probs[igen] <= st && igen < size()) igen++;

  if (verboseLevel > 3) G4cout << " igen " << igen << G4endl;

  return configurations[igen];
}		    
