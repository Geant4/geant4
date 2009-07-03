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

#include "globals.hh"
#include "G4ios.hh"
#include "G4Timer.hh"
#include "G4Pow.hh"

#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char** argv)
{
  if(argc < 2) {
    G4cout << "Name is not specified! Exit" << G4endl;
    exit(1);
  }

  G4double onethird = 1.0/3.0;
  G4int zmax = 255;
  const G4int n = 2551;
  G4int nn = 10000;

  G4cout.setf( std::ios::scientific, std::ios::floatfield );

  G4cout << "Loop  Z = 1-256 over " << n << " points in " << nn << " seria " << G4endl;
  G4cout << "Time comparison of G4Pow calculation of " << argv[1] << G4endl;

  G4int i, j;
  G4double x[n];
  G4double y[n];
  G4double z[n];

  G4double del = G4double(zmax)/G4double(n-1);
  G4double t   = 0.5;
  for (i=0; i<n; i++) {
    x[i] = t;
    t += del;
  }

  G4Pow* power = G4Pow::GetInstance();
  G4Timer* timer1 = new G4Timer();

  timer1->Start();
  for(j=0; j<nn; j++) {
    for (i=0; i<n; i++) {
      t = x[i];
      //y[i] = power->A13(t);
      // y[i] = power->A23(t);
      y[i] = power->logA(t);
    } 
  }

  timer1->Stop();
  G4cout << "======= G4pow   ";
  G4cout << *timer1 << std::endl;

  timer1->Start();
  for(j=0; j<nn; j++) {
    for (i=0; i<n; i++) {
      t = x[i];
      //z[i] = std::pow(t,onethird);
      z[i] = std::log(t);
      //    G4cout << "d[" << x << "] = " << d << G4endl;
    } 
  }

  timer1->Stop();
  G4cout << std::endl;
  G4cout << "======= std:   ";
  G4cout << *timer1 << std::endl;
  
  nn *= 10;

  timer1->Start();
  for(j=0; j<nn; j++) {
    for (i=1; i<256; i++) {
      //y[i] = power->Z13(i);
      // y[i] = power->Z23(i);
      y[i] = power->logZ(i);
    } 
  }

  timer1->Stop();
  G4cout << std::endl;
  G4cout << "======= G4pow  Z ";
  G4cout << *timer1 << std::endl;

  timer1->Start();
  for(j=0; j<nn; j++) {
    for (i=1; i<256; i++) {
      //z[i] = std::pow(G4double(i),onethird);
      z[i] = std::log(G4double(i));
      //    G4cout << "d[" << x << "] = " << d << G4endl;
    } 
  }

  timer1->Stop();
  G4cout << std::endl;
  G4cout << "======= std:   ";
  G4cout << *timer1 << std::endl;

  delete timer1;
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
