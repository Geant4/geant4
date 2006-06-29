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
#ifndef Tst34DetectorConstruction_h
#define Tst34DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class GFlashHomoShowerParameterisation;
class GFlashHitMaker;
class GFlashShowerModel;
class GFlashParticleBounds;

class Tst34DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    Tst34DetectorConstruction();
    ~Tst34DetectorConstruction();

    G4VPhysicalVolume* Construct();

  private:

    G4LogicalVolume* m_experimentalHall_log;
    G4LogicalVolume* m_calo_log;

    G4VPhysicalVolume* m_experimentalHall_phys;
    G4VPhysicalVolume* m_calo_phys;

    G4Box *m_experimentalHall_box;

    G4double m_experimentalHall_x;
    G4double m_experimentalHall_y;
    G4double m_experimentalHall_z;

    G4double m_calo_xside;
    G4double m_calo_yside;
    G4double m_calo_zside;

  // Gflash members

    GFlashHomoShowerParameterisation *m_theParametrisation;
    GFlashHitMaker *m_theHMaker;
    GFlashParticleBounds *m_theParticleBounds;
    GFlashShowerModel* m_theFastShowerModel;
    G4Region * aRegion;
};

#endif
