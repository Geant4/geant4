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
// $Id: G4FissionStore.cc 66241 2012-12-13 18:34:42Z gunter $
//
// 20100728  Move ::addConfig() implementation to .cc file
// 20110801  Make configuration probs a data member array, reduce memory churn
// 20110922  Replace config::print() with stream output.
// 20150608  M. Kelsey -- Label all while loops as terminating.
// 20150619  M. Kelsey -- Replace std::exp with G4Exp

#include "G4FissionStore.hh"
#include "G4FissionConfiguration.hh"
#include "G4Exp.hh"


G4FissionStore::G4FissionStore() : verboseLevel(0) {
  if (verboseLevel > 1) 
    G4cout << " >>> G4FissionStore::G4FissionStore" << G4endl;
}

void G4FissionStore::addConfig(G4double a, G4double z, G4double ez, 
			       G4double ek, G4double ev) {
  G4FissionConfiguration config(a, z, ez, ek, ev);
  configurations.push_back(config);
  if (verboseLevel > 2) G4cout << config << G4endl;
}

G4FissionConfiguration G4FissionStore::generateConfiguration(G4double amax, 
							     G4double rand) const {
  if (verboseLevel > 1)
    G4cout << " >>> G4FissionStore::generateConfiguration" << G4endl;
  
  const G4double small = -30.0;

  G4double totProb = 0.0;
  configProbs.resize(size(),0.);

  if (verboseLevel > 3)
    G4cout << " amax " << amax << " ic " << size() << G4endl;

  for (size_t i = 0; i < size(); i++) {
    G4double ez = configurations[i].ezet;
    G4double pr = ez - amax;

    if (pr < small) pr = small;
    pr = G4Exp(pr);
    if (verboseLevel > 2) {
      G4cout << configurations[i] << "\n probability " << pr << G4endl; 
    }
    totProb += pr;
    configProbs[i] = totProb;  
  };

  G4double st = totProb * rand;

  size_t igen = 0;
  /* Loop checking 08.06.2015 MHK */
  while (configProbs[igen] <= st && igen < size()) igen++;

  if (verboseLevel > 3) G4cout << " igen " << igen << G4endl;

  return configurations[igen];
}		    
