#include "G4VPhysicalVolume.hh"
#include "G4RunManager.hh"
#include "g4std/set"
#include "g4std/iomanip"

#include "B02DetectorConstruction.hh"
#include "B02PhysicsList.hh"
#include "B02PrimaryGeneratorAction.hh"

#include "B02ScoringDetectorConstruction.hh"

// Files specific for scoring 
#include "B02Scorer.hh"
#include "G4Sigma.hh"
#include "G4ParallelScoreManager.hh"

// helper function for print out
G4std::string FillString(const G4std::string &name, char c, G4int n, G4bool back = true);


int main(int argc, char **argv) {
  

  ostream *myout = &G4cout;
  G4int numberOfEvent = 1000;

  G4String random_status_out_file, random_status_in_file;
  G4long myseed = 345354;

  HepRandom::setTheSeed(myseed);

  G4RunManager *runManager = new G4RunManager;
  
  // create the   "tracking"    detector      -----------------
  runManager->SetUserInitialization(new B02DetectorConstruction);
  //  ---------------------------------------------------
  runManager->SetUserInitialization(new B02PhysicsList);
  runManager->SetUserAction(new B02PrimaryGeneratorAction);
  runManager->Initialize();

  // create the detector for scoring
  B02ScoringDetectorConstruction scoringdetector;
  
  // create scorer and manager to score neutrons in the "scoring" detector
  B02Scorer pScorer;
  G4ParallelScoreManager pmgr(*(scoringdetector.Construct()), 
			      "neutron", pScorer);
  pmgr.Initialize();

  runManager->BeamOn(numberOfEvent);

  // ======= after running ============================

  // print all the numbers calculated from the scorer
  *myout << "output pScorer, scoring detector, neutron" << G4endl;
  *myout << pScorer << G4endl;
  *myout << "----------------------------------------------"  << G4endl;

  // print some exclusive numbers

  // head line
  G4int FieldName = 25;
  G4int FieldValue = 12;
  G4std::string vname = FillString("Volume name", ' ', FieldName+1);
  *myout << vname << '|';
  vname = FillString(" AV E/Track ", ' ', FieldValue+1, false);
  *myout << vname << '|';
  vname = FillString(" sigma", ' ', FieldValue+1);
  *myout << vname << '|';
  vname = FillString("Coll_Ent.Tr", ' ', FieldValue+1, false);
  *myout << vname << '|';
  *myout << G4endl;



  const G4PMapPtkTallys &m = pScorer.GetMapPtkTallys();
  for (G4PMapPtkTallys::const_iterator mit = m.begin();
       mit != m.end(); mit++) {
    G4PTouchableKey ptk = mit->first; // get a key identifying a volume
    G4PMapNameTally mtallies = mit->second; // get tallies of the volume
    G4String name(ptk.fVPhysiclaVolume->GetName()); // print volume name
    G4int nEnteringTracks = 0;
    G4double colli_EnteringTrack = 0;
    G4double meanTrackEnergy = 0, sigmaTrackEnergy = 0;
    for (G4PMapNameTally::iterator mt = mtallies.begin();
	 mt != mtallies.end(); mt++) {
      G4String tmp(mt->first);
      if (tmp == "HistorysEntering") {
	nEnteringTracks = mt->second.GetXsum();
      }
      if (tmp == "EnergyEnteringHistory") {
	meanTrackEnergy =  mt->second.GetMean();
	sigmaTrackEnergy = mt->second.GetSigma();
      }
      if (tmp == "Collisions") {
	if (!nEnteringTracks) {
	  G4cout << "exampleB01: Error nEnteringTracks=0" <<G4endl;
	}
	else {
	  colli_EnteringTrack =  mt->second.GetXsum() / nEnteringTracks;
	}
      }
    }


    // print values

    G4std::string fname = FillString(name, '.', FieldName);
    *myout << fname << " |";
    *myout << G4std::setw(FieldValue) << meanTrackEnergy << " |"; 
    *myout << G4std::setw(FieldValue) << sigmaTrackEnergy << " |";
    *myout << G4std::setw(FieldValue) << colli_EnteringTrack << " |";
    *myout << G4endl;
  }


  return 0;
}

G4std::string FillString(const G4std::string &name, char c, G4int n, bool back){
  G4std::string fname;
  G4int k = n - name.size();
  if (k > 0) {
    if (back) {
      fname = name;
      fname += G4std::string(k,c);
    }
    else {
      fname = G4std::string(k,c);
      fname += name;
    }
  }
  else {
    fname = name;
  }
  return fname;
}
