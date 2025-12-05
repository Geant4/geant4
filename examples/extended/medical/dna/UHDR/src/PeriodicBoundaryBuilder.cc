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
/// \file PeriodicBoundaryBuilder.cc
/// \brief Implementation of the PeriodicBoundaryBuilder class

#include "PeriodicBoundaryBuilder.hh"

#include "LogicalVolumePeriodic.hh"

#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4VisAttributes.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* PeriodicBoundaryBuilder::Construct(G4LogicalVolume* logical_world)
{
  auto world = dynamic_cast<G4Box*>(logical_world->GetSolid());

  if (world == nullptr) {
    G4ExceptionDescription ed;
    ed << " PeriodicBoundaryBuilder::Construct: "
       << " Unsupported period boundary for this solid : " << logical_world->GetName() << G4endl;
    G4Exception("G4PeriodicBoundaryProcess::G4PeriodicBoundaryProcess", "Builder01", FatalException,
                ed, "Unsupported period boundary for this solid");
  }
  G4double buffer = 0.01 * nanometer;

  G4double periodic_world_hx = world->GetXHalfLength();
  G4double periodic_world_hy = world->GetYHalfLength();
  G4double periodic_world_hz = world->GetZHalfLength();

  world->SetXHalfLength(world->GetXHalfLength() + buffer / 2);
  world->SetYHalfLength(world->GetYHalfLength() + buffer / 2);
  world->SetZHalfLength(world->GetZHalfLength() + buffer / 2);

  auto* periodic_world = new G4Box("PBC", periodic_world_hx, periodic_world_hy, periodic_world_hz);

  fLogicalPeriodic = new LogicalVolumePeriodic(periodic_world, logical_world->GetMaterial(), "PBC");

  fLogicalPeriodic->SetVisAttributes(G4Color::Magenta());

  new G4PVPlacement(nullptr, G4ThreeVector(), fLogicalPeriodic, "PBC", logical_world, false, 0,
                    true);

  return fLogicalPeriodic;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
