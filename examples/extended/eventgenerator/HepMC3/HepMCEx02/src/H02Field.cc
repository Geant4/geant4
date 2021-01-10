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
/// \file eventgenerator/HepMC/HepMCEx02/src/H02Field.cc
/// \brief Implementation of the H02Field class
//
//   H02Field.hh
//
#include "G4SystemOfUnits.hh"
#include "H02Field.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void H02Field::GetFieldValue(const G4double Point[3], G4double* Bfield) const
{
  const G4double Bz= 3.0*tesla;
  const G4double rmax_sq = sqr(1.*m);
  const G4double zmax = 1.*m;

  Bfield[0]= 0.;
  Bfield[1] = 0.;
  if(std::abs(Point[2])<zmax && (sqr(Point[0])+sqr(Point[1]))<rmax_sq) {
    Bfield[2]= Bz;
  } else {
    Bfield[2]= 0.;
  }
}
