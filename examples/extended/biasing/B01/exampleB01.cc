#include "G4VPhysicalVolume.hh"
#include "G4RunManager.hh"
#include <unistd.h>
#define _GNU_SOURCE
#include <getopt.h>
#include "g4std/set"

#include <iomanip>

#include "B01DetectorConstruction.hh"
#include "B01PhysicsList.hh"
#include "B01PrimaryGeneratorAction.hh"

// Files specific for scoring 
#include "G4PScorer.hh"
#include "G4Sigma.hh"
#include "G4MassScoreManager.hh"

// helper function for print out
std::string FillString(const std::string &name, char c, int n, bool back = true);

int main(int argc, char **argv) {
  

  ostream *myout = &G4cout;
  int numberOfEvent = 1000;

  G4String random_status_out_file, random_status_in_file;
  G4long myseed = 345354;

  HepRandom::setTheSeed(myseed);

  G4RunManager *runManager = new G4RunManager;
  
  // create the detector      ---------------------------
  runManager->
    SetUserInitialization(new B01DetectorConstruction);
  //  ---------------------------------------------------
  runManager->SetUserInitialization(new B01PhysicsList);
  runManager->SetUserAction(new B01PrimaryGeneratorAction);
  runManager->Initialize();

  // create scorer and manager to score neutrons in the detector
  G4PScorer mScorer;
  G4MassScoreManager msm(mScorer, "neutron"); // to be don after 
  msm.Initialize();                           // runManager->Initialize()

  runManager->BeamOn(numberOfEvent);

  // ======= after running ============================

  // print all the numbers calculated from the scorer
  *myout << "output mScorer, mass world, neutron" << G4endl;
  *myout << mScorer << G4endl;
  *myout << "----------------------------------------------"  << G4endl;

  // print some exclusive numbers

  // head line
  int FieldName = 25;
  int FieldValue = 12;
  std::string vname = FillString("Volume name", ' ', FieldName+1);
  *myout << vname << '|';
  vname = FillString(" AV weight ", ' ', FieldValue+1, false);
  *myout << vname << '|';
  vname = FillString(" sigma", ' ', FieldValue+1);
  *myout << vname << '|';
  vname = FillString(" AV E/Track ", ' ', FieldValue+1, false);
  *myout << vname << '|';
  vname = FillString(" sigma", ' ', FieldValue+1);
  *myout << vname << '|';
  *myout << G4endl;



  const G4PMapPtkTallys &m = mScorer.GetMapPtkTallys();
  for (G4PMapPtkTallys::const_iterator mit = m.begin();
       mit != m.end(); mit++) {
    G4PTouchableKey ptk = mit->first; // get a key identifying a volume
    G4PMapNameTally mtallies = mit->second; // get tallies of the volume
    G4String name(ptk.fVPhysiclaVolume->GetName()); // print volume name
    G4double avweight = 0, avsig = 0;
    G4double meanTrackEnergy = 0, sigmaTrackEnergy = 0;
    for (G4PMapNameTally::iterator mt = mtallies.begin();
	 mt != mtallies.end(); mt++) {
      G4String tmp(mt->first);
      if (tmp == "WeighteOfHistorysEntering") {
	avweight =  mt->second.GetMean();
	avsig = mt->second.GetSigma();
      }
      if (tmp == "WeightedEnergyEnteringHistory") {
	meanTrackEnergy =  mt->second.GetMean();
	sigmaTrackEnergy = mt->second.GetSigma();
      }
    }


    // print values

    std::string fname = FillString(name, '.', FieldName);
    *myout << fname << " |";
    *myout << std::setw(FieldValue) << avweight << " |"; 
    *myout << std::setw(FieldValue) << avsig << " |";
    *myout << std::setw(FieldValue) << meanTrackEnergy << " |";
    *myout << std::setw(FieldValue) << sigmaTrackEnergy << " |"; 
    *myout << G4endl;
  }

  return 0;
}

std::string FillString(const std::string &name, char c, int n, bool back){
  std::string fname;
  int k = n - name.size();
  if (k > 0) {
    if (back) {
      fname = name;
      fname += std::string(k,c);
    }
    else {
      fname = std::string(k,c);
      fname += name;
    }
  }
  else {
    fname = name;
  }
  return fname;
}
