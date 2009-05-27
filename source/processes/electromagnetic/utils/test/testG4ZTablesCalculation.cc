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

#include <fstream>
#include <iostream>
#include "G4Timer.hh"
#include "G4Power.hh"
#include "G4String.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char** argv)
{
 // Control on input

  if(argc < 2) {
    G4cout << "Input parameters are not specified! Exit" << G4endl;
    exit(1);
  }

  G4Timer* timer1 = new G4Timer();

  G4String fname = argv[1];
  fname += "_time_comparison_double.out";

  std::ofstream asciiFile;
  asciiFile.open(fname.c_str(), std::ios::out);
  if(asciiFile.is_open()) {
    asciiFile << "Time comparison of " << argv[1] << " calculation" <<G4endl;
    asciiFile << "(loop  Z = 1-256 over 1000000 bins)" <<G4endl;
  } else {
    G4cout << "ERROR file <" << fname << "> is not opened" << G4endl;
    exit(1);
  }

  G4int zmax = 256;

  G4double d;
  G4double zm = 256.;
  G4double nx = (zm - 1.)/1000000.;
  G4Power* power = G4Power::Power();

  timer1->Start();
  for (G4double x = 1.; x < zm; ) {
    //d = power->Z13(x);
    d = power->LogZ(x);
    //    asciiFile << "d[" << x << "] = " << d << G4endl;
    x += nx; 
  } 

  timer1->Stop();
  G4String s1 = "spline: ";
  asciiFile << s1 << " " << *timer1 << G4endl;

  delete timer1;

  G4Timer* timer2 = new G4Timer();

  G4double onethird = 1.0/3.0;

  timer2->Start();
  for (G4double x = 1.; x < zm; ) {
    //d = std::pow((x),onethird);
    d = std::log(x);
    //    asciiFile << "d[" << x << "] = " << d << G4endl;
    x += nx; 
  } 

  timer2->Stop();
  G4String s2 = "formula: ";
  asciiFile << s2 << " " << *timer2 << G4endl;

  delete timer2;
  asciiFile.close();
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
