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
// G4LEPTSDistribution
//
// Author: Pedro Arce (CIEMAT), 2014
// --------------------------------------------------------------------
#ifndef G4LEPTSDistribution_hh
#define G4LEPTSDistribution_hh 1

#include "globals.hh"
#include "Randomize.hh"

class G4LEPTSDistribution
{

 public:

  G4LEPTSDistribution() = default;
  void ReadFile( const G4String& fileName );
  G4bool ReadFile( FILE* fp, G4int nData );
  G4double Sample( G4double, G4double );

  G4bool IsFileFound() const { return bFileFound; }

 private:

  G4int NoBins;
#define NMAX 20000
  G4double E[NMAX], f[NMAX], F[NMAX], eF[NMAX];
  G4bool bFileFound;
};

#endif
