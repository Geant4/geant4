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
//
#include "G4ios.hh"
#include <time.h>
#include <vector>

// Test speed of pow(x, 2) compared to x * x

int main() {
  G4double y[3]= {1.0, 1.3, 1.2}; 
  y[1]=1.2;
  G4int LOOPS = 20000000; // Set test parameters
  G4double x = 1.2;
  clock_t startTime;
  clock_t endTime;

  startTime = clock();

  for(G4int i = 1; i < LOOPS; i++){
    G4double ans = std::pow(y[2], 2);
  };
  endTime = clock();
  G4double firstTime = (G4double)(endTime - startTime) /
    (CLOCKS_PER_SEC * 1000000.0);
  cout << "pow(x, 2) time: " << firstTime  << endl;

  startTime = clock();
  for(G4int j = 1; j < LOOPS; j++){
    G4double ans = y[2] * y[2];
  };

  endTime = clock();
  G4double secondTime = (G4double)(endTime - startTime) / 
    (CLOCKS_PER_SEC * 1000000.0);
  cout << "x * x time: " << secondTime << endl;
  cout << "pow / * speed ratio = " << firstTime / secondTime << endl;
}
