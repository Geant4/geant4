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
// G4LEPTSElossDistr
//
// Author: Pedro Arce (CIEMAT), 2014
// --------------------------------------------------------------------
#ifndef G4LEPTSElossDistr_hh
#define G4LEPTSElossDistr_hh 1

#include "globals.hh"
#include "Randomize.hh"

#include <map> 

class G4LEPTSDistribution;

class G4LEPTSElossDistr
{
 public:

  G4LEPTSElossDistr(const G4String&);

  void ReadFile();
  G4double Sample( G4double, G4double );
  G4bool IsFileFound() const { return bFileFound; }

 private:

  std::map<G4double, std::map<G4double, G4LEPTSDistribution *> > theDistributions;
    // Energy , angle , distribution

  G4int theNDistributions;
  G4String fileName;
  G4int NoBins;
  G4bool bFileFound;
};

#endif
