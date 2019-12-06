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
//
/// \file ExN04CalorimeterParametrisation.cc
/// \brief Implementation of the ExN04CalorimeterParametrisation class
//

#include "ExN04CalorimeterParametrisation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4SystemOfUnits.hh"

ExN04CalorimeterParametrisation::ExN04CalorimeterParametrisation()
{
#include "ExN04DetectorParameterDef.icc"
}

ExN04CalorimeterParametrisation::~ExN04CalorimeterParametrisation()
{;}

void ExN04CalorimeterParametrisation::ComputeTransformation
(const G4int,G4VPhysicalVolume *physVol) const
{
  G4ThreeVector origin;
  physVol->SetTranslation(origin);
}

void ExN04CalorimeterParametrisation::ComputeDimensions
(G4Tubs & calorimeterLayer, const G4int copyNo, const G4VPhysicalVolume*) const
{
  G4double innerRad = caloTubs_rmin
              + copyNo*(absorber_thick+scinti_thick);
  calorimeterLayer.SetInnerRadius(innerRad);
  calorimeterLayer.SetOuterRadius(innerRad+absorber_thick);
  calorimeterLayer.SetZHalfLength(caloTubs_dz);
  calorimeterLayer.SetStartPhiAngle(caloTubs_sphi);
  calorimeterLayer.SetDeltaPhiAngle(caloTubs_dphi);
}
