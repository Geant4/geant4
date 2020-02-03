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
/// \file runAndEvent/RE01/src/RE01CalorimeterParametrisation.cc
/// \brief Implementation of the RE01CalorimeterParametrisation class
//
//

#include "RE01CalorimeterParametrisation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE01CalorimeterParametrisation::RE01CalorimeterParametrisation()
  :G4VPVParameterisation()
{
#include "RE01DetectorParameterDef.icc"
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE01CalorimeterParametrisation::~RE01CalorimeterParametrisation()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RE01CalorimeterParametrisation::ComputeTransformation
(const G4int,G4VPhysicalVolume *physVol) const
{
  G4ThreeVector origin;
  physVol->SetTranslation(origin);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RE01CalorimeterParametrisation::ComputeDimensions
(G4Tubs & calorimeterLayer, const G4int copyNo, const G4VPhysicalVolume*) const
{
  G4double innerRad = fCaloTubs_rmin
              + copyNo*(fAbsorber_thick+fScinti_thick);
  calorimeterLayer.SetInnerRadius(innerRad);
  calorimeterLayer.SetOuterRadius(innerRad+fAbsorber_thick);
  calorimeterLayer.SetZHalfLength(fCaloTubs_dz);
  calorimeterLayer.SetStartPhiAngle(fCaloTubs_sphi);
  calorimeterLayer.SetDeltaPhiAngle(fCaloTubs_dphi);
}
