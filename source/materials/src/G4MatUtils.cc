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

#include "G4MatUtils.hh"
#include "G4PhysicsFreeVector.hh"
#include <fstream>
#include <sstream>

namespace
{
  constexpr G4int nmax = 40;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ExtendedPhysicsVector*
G4MatUtils::BuildExtendedVector(const G4String& dirpath,
				const G4String& filename,
				const G4int nxsections,
				const G4int length,
				const G4double unitE,
				const G4double unitS)
{
  std::ostringstream ss;
  ss << dirpath << "/" << filename;
  std::ifstream infile(ss.str(), std::ios::in);

  if (nxsections < 1 || length <= 0 || nxsections > nmax) {
    G4ExceptionDescription ed;
    ed << " Wrong data size: nsections=" << nxsections << "  length="
       << length << "  nmax=" << nmax;
    G4Exception("G4MatUtils::BuildExtendedVector(..)","mat004",
                FatalException, ed, "Check input values");
    return nullptr;
  }

  // file is not opened
  if (!infile.is_open()) {
    G4ExceptionDescription ed;
    ed << " Fail to open file: <" << ss.str() << ">";
    G4Exception("G4MatUtils::BuildExtendedVector(..)","mat004",
                FatalException, ed, "Check file path");
    return nullptr;
  }
  // file is opened
  auto vtot = new G4PhysicsFreeVector(false);
  G4ExtendedPhysicsVector* v = new G4ExtendedPhysicsVector(vtot, nxsections);
  v->SetDataLength(length);

  G4double y[nmax];
  G4double e, yy, sum;
  auto nn = (std::size_t)length;
  
  for (std::size_t i=0; i<nn; ++i) {
    sum = 0.0;
    infile >> e;
    e *= unitE;
    for (G4int j=0; j<nxsections; ++j) {
      infile >> yy;
      yy *= unitS;
      sum += yy;
      y[j] = yy;
    }
    vtot->PutValues(i, e, sum);
    v->PutPartialXSData(i, y);
  }
  return v;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
