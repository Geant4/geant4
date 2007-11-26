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
#ifndef LXeDetectorConstruction_H
#define LXeDetectorConstruction_H 1

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Box;
class G4Tubs;
class LXePMTSD;
class LXeScintSD;
class G4Sphere;

#include "G4Material.hh"
#include "LXeDetectorMessenger.hh"
#include "G4VisAttributes.hh"
#include "G4RotationMatrix.hh"
#include "G4VUserDetectorConstruction.hh"

class LXeDetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    LXeDetectorConstruction();
    ~LXeDetectorConstruction();

    G4VPhysicalVolume* Construct();

  //Functions to modify the geometry
  void SetDimensions(G4ThreeVector dims);
  void SetHousingThickness(G4double d_mtl);
  void SetNX(G4int nx);
  void SetNY(G4int ny);
  void SetNZ(G4int nz);
  void SetPMTRadius(G4double outerRadius_pmt);
  void SetDefaults();

  //Get values
  G4double GetScintX(){return scint_x;}
  G4double GetScintY(){return scint_y;}
  G4double GetScintZ(){return scint_z;}
  G4double GetHousingThickness(){return d_mtl;}
  G4int GetNX(){return nx;}
  G4int GetNY(){return ny;}
  G4int GetNZ(){return nz;}
  G4double GetPMTRadius(){return outerRadius_pmt;}
  G4double GetSlabZ(){return slab_z;}
    
  //rebuild the geometry based on changes. must be called
  void UpdateGeometry();
  G4bool GetUpdated(){return updated;}

  void SetSphereOn(G4bool b){sphereOn=b;updated=true;}
  static G4bool GetSphereOn(){return sphereOn;}

  void SetHousingReflectivity(G4double r){refl=r;updated=true;}
  G4double GetHousingReflectivity(){return refl;}

  void SetWLSSlabOn(G4bool b){WLSslab=b;updated=true;}
  G4bool GetWLSSlabOn(){return WLSslab;}

  void SetMainVolumeOn(G4bool b){mainVolume=b;updated=true;}
  G4bool GetMainVolumeOn(){return mainVolume;}

  void SetNFibers(G4int n){nfibers=n;updated=true;}
  G4int GetNFibers(){return nfibers;}

  void SetMainScintYield(G4double y);
  void SetWLSScintYield(G4double y);

private:

  void DefineMaterials();
  G4VPhysicalVolume* ConstructDetector();

  LXeDetectorMessenger* detectorMessenger;

  G4bool updated;
   
  G4Box* experimentalHall_box;
  G4LogicalVolume* experimentalHall_log;
  G4VPhysicalVolume* experimentalHall_phys;

  //Materials & Elements
  G4Material* LXe;
  G4Material* Al;
  G4Element* N;
  G4Element* O;
  G4Material* Air;
  G4Material* Vacuum;
  G4Element* C;
  G4Element* H;
  G4Material* Glass;
  G4Material* Pstyrene;
  G4Material* PMMA;
  G4Material* Pethylene;
  G4Material* fPethylene; 



  //Geometry
  G4double scint_x;
  G4double scint_y;
  G4double scint_z;
  G4double d_mtl;
  G4int nx;
  G4int ny;
  G4int nz;
  G4double outerRadius_pmt;
  G4int nfibers;
  static G4bool sphereOn;
  G4double refl;
  G4bool WLSslab;
  G4bool mainVolume;
  G4double slab_z;

  G4MaterialPropertiesTable* LXe_mt;
  G4MaterialPropertiesTable* MPTPStyrene;

};

#endif




