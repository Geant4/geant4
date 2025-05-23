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
/// \file VG01DetectorConstruction.cc
/// \brief Implementation of the VG01DetectorConstruction class

//  Authors:   J. Apostolakis & S. Wenzel (CERN)   2018-2021
//
//  Derived from FullCMS code by Mihaly Novak  (2017-18)

#include "VG01DetectorConstruction.hh"

#include "G4VecGeomConverter.h"
#include "VG01DetectorMessenger.hh"

#include "G4FieldManager.hh"
#include "G4PropagatorInField.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4TransportationManager.hh"
#include "G4UniformMagField.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"
#include "globals.hh"

// NOTE: VecGeom versioning was introduced in 1.2.4, the minimum required to
// build this file
#include <VecGeom/base/Version.h>
#include <VecGeom/management/GeoManager.h>
#include <VecGeom/management/LoggerTypes.h>

#if VECGEOM_VERSION >= 0x010208
// NOTE: Logger was introduced in 1.2.8; earlier version log messages will go
// to stdout via printf/cout
#  include <VecGeom/management/Logger.h>
#  define VECGEOM_HAS_LOGGER 1
#endif

#ifdef VECGEOM_HAS_LOGGER
namespace
{
using vecgeom::LogLevel;
using vecgeom::Provenance;

// This logger adapts VecGeom messages, warnings, and errors to Geant4 io/exceptions.
void G4VecGeomLogger(Provenance prov, LogLevel lev, std::string str)
{
  std::ostringstream where;
  if (!prov.file.empty()) {
    // Only print out the last path component of the originating file
    auto pos = prov.file.rfind('/') + 1;
    where << prov.file.substr(pos);
  }
  if (prov.line != 0) {
    where << ':' << prov.line;
  }

  if (lev < LogLevel::warning) {
    G4cout << "VecGeom: " << where.str() << ": " << vecgeom::to_cstring(lev) << ": " << str << std::endl;
    return;
  }

  G4Exception(where.str().c_str(), "VecGeom",
              lev == LogLevel::warning ? JustWarning : FatalException, str.c_str());
}
}  // namespace
#endif

G4double VG01DetectorConstruction::fglobFieldValue = 0.0;

VG01DetectorConstruction::VG01DetectorConstruction()
  : G4VUserDetectorConstruction(),
    fWorld(nullptr),
    fFieldMgr(nullptr),
    fUniformMagField(nullptr),
    fDetectorMessenger(nullptr)
{
  fGDMLFileName = "TestNTST.gdml";
  fDetectorMessenger = new VG01DetectorMessenger(this);
#ifdef VECGEOM_HAS_LOGGER
  // Have the global VecGeom logger redirect to G4cout rather than the default std:clog
  auto orig_level = vecgeom::logger().level();
  vecgeom::Logger new_logger(&G4VecGeomLogger);
  new_logger.level(orig_level);
  vecgeom::logger() = std::move(new_logger);
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

VG01DetectorConstruction::~VG01DetectorConstruction()
{
  delete fDetectorMessenger;
  if (fUniformMagField) {
    delete fUniformMagField;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* VG01DetectorConstruction::Construct()
{
  //  parser.SetOverlapCheck(true);
  G4TransportationManager* trMgr = G4TransportationManager::GetTransportationManager();
  assert(trMgr);

  fParser.Read(fGDMLFileName, false);
  fWorld = (G4VPhysicalVolume*)fParser.GetWorldVolume();

  // enabling 'check' mode in G4 Navigator
  G4Navigator* nav = trMgr->GetNavigatorForTracking();
  assert(nav);
  nav->CheckMode(true);
  // nav->SetVerboseLevel(1);
  std::cout << "Enabled Check mode in G4Navigator";

  // write back
  // fParser.Write("out.gdml", fWorld);

  fFieldMgr = trMgr->GetFieldManager();
  fWorld->GetLogicalVolume()->SetVisAttributes(G4VisAttributes::GetInvisible());
  if (fWorld == nullptr) {
    G4ExceptionDescription ed;
    ed << "World volume not set properly check your setup selection criteria"
       << "or GDML input!" << G4endl;
    G4Exception("VG01DetectorConstruction::Construct()", "G4VecGeomNavExtExample_0001",
                FatalException, ed);
  }
  CreateMagFieldAndIntegrator();

  if (fUseVecGeom) {
    // This is converting the geometry to VecGeom and implicitely also setting
    // the navigator; We should pull this out really
    G4VecGeomConverter::Instance().SetVerbose(1);
    G4VecGeomConverter::Instance().ConvertG4Geometry(fWorld);
    G4cout << vecgeom::GeoManager::Instance().getMaxDepth() << "\n";
  }
  return fWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void VG01DetectorConstruction::CreateMagFieldAndIntegrator()
{
  delete fUniformMagField;

  if (std::abs(fglobFieldValue) > 0.0) {
    // Apply a global uniform magnetic field along the Z axis.
    // Note:  only when the magnetic field (pointer) is NOT zero,
    //        does the Geant4 transportion in field get activated.
    fUniformMagField = new G4UniformMagField(G4ThreeVector(0.0, 0.0, fglobFieldValue));
    fFieldMgr->SetDetectorField(fUniformMagField);
    fFieldMgr->CreateChordFinder(fUniformMagField);
    G4cout << G4endl << " *** SETTING MAGNETIC FIELD in Z direction : fieldValue = "
           << fglobFieldValue / tesla << " tesla *** " << G4endl << G4endl;
  }
  else {
    fFieldMgr->SetDetectorField(nullptr);
    G4cout << G4endl << " *** NO MAGNETIC FIELD SET  *** " << G4endl << G4endl;
  }
}
