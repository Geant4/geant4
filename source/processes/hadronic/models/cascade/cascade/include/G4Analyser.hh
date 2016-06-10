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
// $Id: G4Analyser.hh 66241 2012-12-13 18:34:42Z gunter $
//
// 20101010  M. Kelsey -- Migrate to integer A and Z

#ifndef G4ANALYSER_HH
#define G4ANALYSER_HH

#define WITH_NUCLEI

#include "G4CollisionOutput.hh"
//#include "G4InuclParticle.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclNuclei.hh"
#include "G4NuclWatcher.hh"
//#include "G4ExitonConfiguration.hh"

#include <vector>

class G4Analyser {

public:

  G4Analyser();
  void setInelCsec(G4double csec, G4bool withn);
  void setWatchers(const std::vector<G4NuclWatcher>& watchers);
  void try_watchers(G4int a, G4int z, G4bool if_nucl);
  void analyse(const G4CollisionOutput& output);
  void printResults();
  void printResultsSimple();
  void handleWatcherStatistics();
  void printResultsNtuple();

private: 

  G4int verboseLevel;
  G4double eventNumber;
  G4double averageMultiplicity;
  G4double averageProtonNumber;
  G4double averageNeutronNumber;
  G4double averagePionNumber;
  G4double averageNucleonKinEnergy;
  G4double averageProtonKinEnergy;
  G4double averageNeutronKinEnergy;
  G4double averagePionKinEnergy;
  G4double averageExitationEnergy;
  G4double averageOutgoingNuclei;
  G4double fissy_prob;
  G4double averagePionPl;
  G4double averagePionMin;
  G4double averagePion0;
  G4double averageA;
  G4double averageZ;
  std::vector<G4NuclWatcher> ana_watchers;
  G4double inel_csec;
  G4bool withNuclei;
};        

#endif // G4ANALYSER_HH
