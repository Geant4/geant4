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
