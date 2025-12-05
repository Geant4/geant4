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
// gpaterno, October 2025
//
/// \file DetectorConstruction.hh
/// \brief Description of the DetectorConstruction class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4ios.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include <vector>

#include "G4Region.hh"
#include "G4PVPlacement.hh"

#include "DetectorConstructionMessenger.hh"
#include "G4ChannelingFastSimModel.hh"

#define NSpheresMax 10000

class G4VPhysicalVolume;
class G4LogicalVolume;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// Detector construction class to define materials and geometry.

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction();
    ~DetectorConstruction() override = default;

    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;
    
    //method to get the scoring volumes
    std::vector<G4LogicalVolume*> GetScoringVolume() const {
        return fScoringVolume;}        
    
    //method to set if it is a hybrid source or not
    void SetHybridSource(G4bool val) {fHybridSource = val;}
   
    //methods to set the Crystal (Radiator) features
    void SetCrystalMaterial(G4String val) {fCrystalMaterialStr = val;}
    void SetCrystalSize(G4ThreeVector val) {fCrystalSize = val;}
    void SetCrystalBendingAngle(G4double val) {fBendingAngle = val;}    
    void SetCrystalLattice(G4String val) {fLattice = val;}
    void SetCrystalAngleX(G4double val) {fAngleX = val;}
    void SetCrystalAngleY(G4double val) {fAngleY = val;}
    G4double GetCrystalZ() const {return fCrystalZ;}
    void SetRadiationModel(G4bool val) {fActivateRadiationModel = val;}
    void SetOCeffects(G4bool val) {fActivateOCeffects = val;}
    G4bool GetOCeffects() const {return fActivateOCeffects;}
    G4LogicalVolume* GetCrystalVolume() const {return fCrystalLogic;}
    void SetPotentialPath(const G4String path){fPotentialPath = path;}
    
    //method to set/get the Converter (Target) features
    void SetRadiatorConverterSepDistance(G4double val) {
        fRadiatorConverterSepDistance = val;}
    G4double GetRadiatorConverterSepDistance() const {
        return fRadiatorConverterSepDistance;}
    void SetConverterSize(G4ThreeVector val) {fConverterSize = val;}
    void SetConverterMaterial(G4String val) {fConverterMaterialStr = val;}
    void SetGranularConverter(G4bool val) {fGranularConverter = val;}
    void SetSphereRadius(G4double val) {fSphereRadius = val;}
    G4int GetNSpheres() const {return fNSpheres;}
    G4LogicalVolume* GetConverterVolume() const {return fConverterLogic;}
       
    //methods to set the Magnetic field features
    void SetMagneticField(G4bool val) {fSetMagneticField = val;}
    void SetFieldValue(G4double val) {fFieldValue = val;}
    void SetFieldRegionLength(G4double val) {fFieldRegionLength = val;}
       
    //methods to set the Collimator features
    void SetCollimator(G4bool val) {fSetCollimator = val;}
    void SetCollimatorHole(G4String val) {fCollimatorHole = val;}
    void SetCollimatorAperture(G4double val) {fCollimatorAperture = val;}
    void SetCollimatorThickness(G4double val) {fCollimatorThickness = val;}
    void SetCollimatorSide(G4double val) {fCollimatorSide = val;}
    void SetRadiatorCollimatorSepDistance(G4double val) {
        fRadiatorCollimatorSepDistance = val;}
    G4double GetRadiatorCollimatorSepDistance() const {
        return fRadiatorCollimatorSepDistance;}
       
    //methods to set/Get the Virtual Detector features
    void SetVirtualDetectorSize(G4ThreeVector val) {fVirtualDetectorSize = val;}
    std::vector<G4ThreeVector> GetVirtualDetectorPositionVector() const {
        return fVirtualDetectorPositionVector;}
       
    //methods to set and get ScoreCrystalExit (27/09/2024)
    void SetScoringCrystalExit(G4bool bval) {fScoringCrystalExit = bval;} 
    G4bool GetScoringCrystalExit() const {return fScoringCrystalExit;}
        
protected:
  std::vector<G4LogicalVolume*> fScoringVolume; //for spheres only

private:
    DetectorConstructionMessenger* fMessenger;  
    
    G4bool fHybridSource = true;
        
    G4Region* fCrystalRegion{nullptr};
    G4LogicalVolume* fCrystalLogic{nullptr};
    G4String fCrystalMaterialStr = "W";
    G4Material* fCrystalMaterial{nullptr};
    G4ThreeVector fCrystalSize = G4ThreeVector(7.*mm, 7.*mm, 2.*mm);
    G4double fBendingAngle = 0.e-6; //rad
    G4String fLattice = "<111>";  
    G4double fAngleX = 0.e-6; //rad
    G4double fAngleY = 0.e-6; //rad
    G4double fCrystalZ = 0.;
    G4bool fActivateRadiationModel = true;
    G4bool fActivateOCeffects = true;
    G4String fPotentialPath = "";
    
    G4double fRadiatorConverterSepDistance = 60.*cm;    
    G4ThreeVector fConverterSize = G4ThreeVector(199.75*mm, 199.75*mm, 11.6*mm);
    G4double fConverterZ = 0.;
    G4LogicalVolume* fConverterLogic{nullptr};
    G4bool fGranularConverter = false;
    G4String fConverterMaterialStr = "W";
    G4Material* fConverterMaterial{nullptr};
    G4double fSphereRadius = 1.1*mm;
    G4LogicalVolume* fSphereLogic[NSpheresMax];
    G4int fNSpheres = 0;
    G4bool fConverter = true;
       
    G4bool fSetMagneticField = false;
    G4double fFieldValue = 100.*tesla;
    G4double fFieldRegionLength = 90.*cm;
    G4LogicalVolume* fMFlogic{nullptr}; 
    
    G4bool fSetCollimator = false;
    G4double fCollimatorAperture = 2.*mm;
    G4String fCollimatorHole = "squared";
    G4double fCollimatorThickness = 50.*cm;
    G4double fCollimatorSide = 2.5*m;
    G4double fRadiatorCollimatorSepDistance = 5.*cm;
    G4LogicalVolume* fCollimatorLogic{nullptr}; 
      
    G4ThreeVector fVirtualDetectorSize = G4ThreeVector(40.*cm, 40.*cm, 0.01*mm);
    std::vector<G4ThreeVector> fVirtualDetectorPositionVector;
    G4LogicalVolume* fVirtualDetectorLogic0{nullptr};
    G4LogicalVolume* fVirtualDetectorLogic1{nullptr};
    G4LogicalVolume* fVirtualDetectorLogic2{nullptr};
            
    G4bool fScoringCrystalExit = false;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
