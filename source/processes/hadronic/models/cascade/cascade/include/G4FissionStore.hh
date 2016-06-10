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
// $Id: G4FissionStore.hh 66241 2012-12-13 18:34:42Z gunter $
//
// 20100728  Move ::addConfig() to .cc file, add setVerboseLevel(), clear()
// 20110801  Make configuration probs a data member array, reduce memory churn

#ifndef G4FISSION_STORE_HH
#define G4FISSION_STORE_HH

#include "G4FissionConfiguration.hh"
#include <vector>

class G4FissionStore {
public:
  G4FissionStore();

  void setVerboseLevel(G4int verbose=1) { verboseLevel = verbose; }

  void addConfig(G4double a, G4double z, G4double ez, G4double ek, G4double ev);

  void clear() { configurations.clear(); }

  size_t size() const { return configurations.size(); }

  G4FissionConfiguration generateConfiguration(G4double amax, 
					       G4double rand) const;

private:
  G4int verboseLevel;
  std::vector<G4FissionConfiguration> configurations;
  mutable std::vector<G4double> configProbs;
};

#endif // G4FISSION_STORE_HH 



