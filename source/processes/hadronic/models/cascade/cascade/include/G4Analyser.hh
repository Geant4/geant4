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
#define G4ANALYSER_HH

#include "G4CollisionOutput.hh"
//#include "G4InuclParticle.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclNuclei.hh"
#include "G4NuclWatcher.hh"
//#include "G4ExitonConfiguration.hh"

#include "g4std/vector"

class G4Analyser {

public:

  G4Analyser();
  void setInelCsec(G4double csec, 
		   G4bool withn);
  void setWatchers(const G4std::vector<G4NuclWatcher>& watchers);
  void try_watchers(G4double a, 
		    G4double z, 
		    G4bool if_nucl);
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
  G4double averageNucleiFragments;
  G4double fissy_prob;
  G4double averagePionPl;
  G4double averagePionMin;
  G4double averagePion0;
  G4double averageA;
  G4double averageZ;
  G4std::vector<G4NuclWatcher> ana_watchers;
  G4double inel_csec;
  G4bool withNuclei;
};        

#endif // G4ANALYSER_HH
