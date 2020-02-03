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
/// \file optical/LXe/include/LXeDetectorConstruction.hh
/// \brief Definition of the LXeDetectorConstruction class
//
//
#ifndef LXeDetectorConstruction_H
#define LXeDetectorConstruction_H 1

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Box;
class G4Tubs;
class LXeMainVolume;
class G4Sphere;

#include "G4Material.hh"
#include "LXeDetectorMessenger.hh"
#include "G4VisAttributes.hh"
#include "G4RotationMatrix.hh"

#include "LXeScintSD.hh"
#include "LXePMTSD.hh"

#include "G4VUserDetectorConstruction.hh"
#include "G4Cache.hh"

class LXeDetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    LXeDetectorConstruction();
    virtual ~LXeDetectorConstruction();

    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

    //Functions to modify the geometry
    void SetDimensions(G4ThreeVector );
    void SetHousingThickness(G4double );
    void SetNX(G4int );
    void SetNY(G4int );
    void SetNZ(G4int );
    void SetPMTRadius(G4double );
    void SetDefaults();
    void SetSaveThreshold(G4int );

    //Get values
    G4int GetNX() const {return fNx;};
    G4int GetNY() const {return fNy;};
    G4int GetNZ() const {return fNz;};
    G4int GetSaveThreshold() const {return fSaveThreshold;};
    G4double GetScintX() const {return fScint_x;}
    G4double GetScintY() const {return fScint_y;}
    G4double GetScintZ() const {return fScint_z;}
    G4double GetHousingThickness() const {return fD_mtl;}
    G4double GetPMTRadius() const {return fOuterRadius_pmt;}
    G4double GetSlabZ() const {return fSlab_z;}
 
    void SetSphereOn(G4bool );
    static G4bool GetSphereOn(){return fSphereOn;}

    void SetHousingReflectivity(G4double );
    G4double GetHousingReflectivity() const {return fRefl;}

    void SetWLSSlabOn(G4bool b);
    G4bool GetWLSSlabOn() const {return fWLSslab;}

    void SetMainVolumeOn(G4bool b);
    G4bool GetMainVolumeOn() const {return fMainVolumeOn;}

    void SetNFibers(G4int n);
    G4int GetNFibers() const {return fNfibers;}

    void SetMainScintYield(G4double );
    void SetWLSScintYield(G4double );

  private:

    void DefineMaterials();

    LXeDetectorMessenger* fDetectorMessenger;

    G4Box* fExperimentalHall_box;
    G4LogicalVolume* fExperimentalHall_log;
    G4VPhysicalVolume* fExperimentalHall_phys;

    //Materials & Elements
    G4Material* fLXe;
    G4Material* fAl;
    G4Element* fN;
    G4Element* fO;
    G4Material* fAir;
    G4Material* fVacuum;
    G4Element* fC;
    G4Element* fH;
    G4Material* fGlass;
    G4Material* fPstyrene;
    G4Material* fPMMA;
    G4Material* fPethylene1;
    G4Material* fPethylene2;

    //Geometry
    G4double fScint_x;
    G4double fScint_y;
    G4double fScint_z;
    G4double fD_mtl;
    G4int fNx;
    G4int fNy;
    G4int fNz;
    G4int fSaveThreshold;
    G4double fOuterRadius_pmt;
    G4int fNfibers;
    static G4bool fSphereOn;
    G4double fRefl;
    G4bool fWLSslab;
    G4bool fMainVolumeOn;
    G4double fSlab_z;

    LXeMainVolume* fMainVolume;

    G4MaterialPropertiesTable* fLXe_mt;
    G4MaterialPropertiesTable* fMPTPStyrene;

    //Sensitive Detectors
    G4Cache<LXeScintSD*> fScint_SD;
    G4Cache<LXePMTSD*> fPmt_SD;

};

#endif
