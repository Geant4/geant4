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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software 
// and the DNA geometry given in the Geom_DNA example 
// shall cite the following Geant4-DNA collaboration publications:
// [1] NIM B 298 (2013) 47-54
// [2] Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
//
/// \file ChromosomeParameterisation.cc
/// \brief Implementation of the ChromosomeParameterisation class

#include <ChromosomeParameterisation.hh>

// STD
#include <fstream>

// G4
#include <CLHEP/Units/SystemOfUnits.h>
#include <G4VPhysicalVolume.hh>

using namespace std;
using CLHEP::nanometer;
using CLHEP::degree;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ChromosomeParameterisation::ChromosomeParameterisation(const char* filename):
    G4VPVParameterisation()
{
  ifstream f(filename, ios::in);
  if (!f)
    return;

  fPositions.reserve(100);
  fRotations.reserve(100);

  while (!f.eof())
  {
    double x, y, z, rot;
    f >> x >> y >> z >> rot;
    fPositions.push_back(new G4ThreeVector(x * nanometer,
                                           y * nanometer,
                                           z * nanometer));
    fRotations.push_back(new G4RotationMatrix(0, 0, rot * degree));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ChromosomeParameterisation::~ChromosomeParameterisation()
{
  unsigned int i;
  for (i = 0; i < fPositions.size(); i++)
    delete fPositions[i];
  for (i = 0; i < fRotations.size(); i++)
    delete fRotations[i];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChromosomeParameterisation::ComputeDimensions(
    G4Tubs& /*rosette*/,
    const G4int /*copyNo*/,
    const G4VPhysicalVolume* /*physVol*/) const
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChromosomeParameterisation::ComputeTransformation(
    const G4int copyNo, G4VPhysicalVolume* physVol) const
{
  // see passive vs active method to specify
  // this transformation
  physVol->SetTranslation(*fPositions[copyNo]);
  physVol->SetRotation(fRotations[copyNo]);
}
