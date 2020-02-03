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
/// \file runAndEvent/RE04/include/RE04DetectorConstruction.hh
/// \brief Definition of the RE04DetectorConstruction class
//
//
#ifndef RE04DetectorConstruction_h
#define RE04DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4Material;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Region;

//
/// User detector construction class
///
/// - void DefineMaterials()
///     defines materials
/// - void SetupGeometry()
///     creates the world volume and detector geometries
//
class RE04DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    RE04DetectorConstruction();
    virtual ~RE04DetectorConstruction();

  public:
    virtual G4VPhysicalVolume* Construct();
  
  private:
    void DefineMaterials();
    void SetupGeometry();

  private:
    G4Material* fAir;
    G4Material* fWater;
    G4Material* fPb;
    G4VPhysicalVolume* fWorldPhys;
    G4bool fConstructed;

};

#endif

