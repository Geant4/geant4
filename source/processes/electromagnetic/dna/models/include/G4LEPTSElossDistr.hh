#ifndef G4LEPTSElossDistr_h
#define G4LEPTSElossDistr_h 1

#include "globals.hh"
#include "Randomize.hh"
#include <string.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <map> 
using namespace std;
class G4LEPTSDistribution;
typedef std::map<G4double, G4LEPTSDistribution*> mddist;
typedef std::map<G4double, mddist > mdmddist;


class G4LEPTSElossDistr {

public:

  G4LEPTSElossDistr( string);    // Constructor
  void ReadFile();          // Read file
  G4double Sample( G4double, G4double);

  G4bool IsFileFound() const {
    return bFileFound;
  }

private:
  mdmddist theDistributions; // Energy , angle , distribution
  G4int theNDistributions;
  string fileName;
  int    NoBins;
#define NMAX 20000
  G4bool bFileFound;
};


#endif
