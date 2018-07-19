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
/// \file optical/OpNovice2/include/DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4OpticalSurface.hh"

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4RunManager.hh"


class DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction();
    virtual ~DetectorConstruction();

    G4VPhysicalVolume* GetTank() {return fTank;}
    G4double GetTankXSize() {return fTank_x;}

    G4OpticalSurface* GetSurface(void) {return fSurface;}

    void SetSurfaceFinish(const G4OpticalSurfaceFinish finish) {
      fSurface->SetFinish(finish);
      G4RunManager::GetRunManager()->GeometryHasBeenModified();
    }
    G4OpticalSurfaceFinish GetSurfaceFinish(void) 
      {return fSurface->GetFinish();}

    void SetSurfaceType(const G4SurfaceType type) {
      fSurface->SetType(type);
      G4RunManager::GetRunManager()->GeometryHasBeenModified();
    }
    
    void SetSurfaceModel(const G4OpticalSurfaceModel model) {
      fSurface->SetModel(model);
      G4RunManager::GetRunManager()->GeometryHasBeenModified();
    }
    G4OpticalSurfaceModel GetSurfaceModel(void) 
      {return fSurface->GetModel();}

    void SetSurfaceSigmaAlpha(G4double v);

    void AddBoxMPV(const char* c, G4MaterialPropertyVector* mpv);
    void AddBoxMPCV(const char* c, G4double v);
    G4MaterialPropertiesTable* GetBoxMaterialPropertiesTable() 
      {return fBoxMPT;}

    void AddWorldMPV(const char* c, G4MaterialPropertyVector* mpv);
    void AddWorldMPCV(const char* c, G4double v);
    G4MaterialPropertiesTable* GetWorldMaterialPropertiesTable() 
      {return fWorldMPT;}

    void AddSurfaceMPV(const char* c, G4MaterialPropertyVector* mpv);
    G4MaterialPropertiesTable* GetSurfaceMaterialPropertiesTable() 
      {return fSurfaceMPT;}


    virtual G4VPhysicalVolume* Construct();

  private:
    G4double fExpHall_x;
    G4double fExpHall_y;
    G4double fExpHall_z;

    G4VPhysicalVolume* fTank;

    G4double fTank_x;
    G4double fTank_y;
    G4double fTank_z;

    G4OpticalSurface* fSurface;

    DetectorMessenger* fDetectorMessenger;

    G4MaterialPropertiesTable* fBoxMPT;
    G4MaterialPropertiesTable* fWorldMPT;
    G4MaterialPropertiesTable* fSurfaceMPT;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*DetectorConstruction_h*/
