#ifndef G4LEPTSDistribution_h
#define G4LEPTSDistribution_h 1

#include "globals.hh"
#include "Randomize.hh"
#include <string.h>
#include <stdio.h>
#include <string>
#include <iostream> 
using namespace std;


class G4LEPTSDistribution {

public:

  G4LEPTSDistribution();    // Constructor
  void ReadFile( G4String fileName );          // Read file
  G4bool ReadFile( FILE* fp, G4int nData );
  G4double Sample( G4double, G4double);

  G4bool IsFileFound() const {
    return bFileFound;
  }

private:
  int    NoBins;
#define NMAX 20000
  G4double E[NMAX], f[NMAX], F[NMAX], eF[NMAX];
  G4bool bFileFound;
};


#endif
