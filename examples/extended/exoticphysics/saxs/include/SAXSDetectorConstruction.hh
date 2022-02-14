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
/// \file SAXSDetectorConstruction.hh
/// \brief Implementation of the SAXSDetectorConstruction class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef SAXSDetectorConstruction_h
#define SAXSDetectorConstruction_h 1
#endif

#include "SAXSDetectorConstructionMessenger.hh"

#include "G4VUserDetectorConstruction.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"

#include "G4RunManager.hh"
#include "globals.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4NistManager.hh"
#include "G4ExtendedMaterial.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// Detector construction.

class SAXSDetectorConstruction : public G4VUserDetectorConstruction
{
public:   
  SAXSDetectorConstruction();
  ~SAXSDetectorConstruction();
  
  void DefineMaterials();     
  void SetGeometricalVariables();        
  
  virtual G4VPhysicalVolume* Construct();

  G4LogicalVolume* GetSensitiveVolume() const {return fSensitiveVolume;}        
  G4LogicalVolume* GetPhantom() const {return fPhantomLogic;}
        
protected:
    G4LogicalVolume* fSensitiveVolume;        

private:     
    virtual void ConstructSDandField();
    
    SAXSDetectorConstructionMessenger* fMessenger;   
    
    G4bool fIWantSlits;
      
        //-----------------definition of Volumes-------------------
        //World
    G4Box* fWorldSolid;
    G4LogicalVolume* fWorldLogic;
    G4VPhysicalVolume* fWorldPhysical; 
        
    //Phantom
    G4LogicalVolume* fPhantomLogic;
    G4VPhysicalVolume* fPhantomPhysical;  
    G4Material* fPhantomMaterial;
    G4int fPhantomMaterialIndex;        
                
    //Slits
    G4LogicalVolume* fSlit1Logic;
    G4LogicalVolume* fSlit2Logic;
    G4LogicalVolume* fSlit3Logic;
    G4LogicalVolume* fSlit4Logic;
    G4VPhysicalVolume* fSlit1Physical;
    G4VPhysicalVolume* fSlit2Physical;
    G4VPhysicalVolume* fSlit3Physical;
    G4VPhysicalVolume* fSlit4Physical;
        
        //Detector
        G4LogicalVolume* fDetectorLogic;
    G4VPhysicalVolume* fDetectorPhysical;
    
    //Shielding
        G4LogicalVolume* fShieldingLogic;
    G4VPhysicalVolume* fShieldingPhysical;
    G4LogicalVolume* fShieldingBackLogic;
    G4VPhysicalVolume* fShieldingBackPhysical;
    
    //--------------------definition of Materials--------------
    //materials for MIFF study
    G4Material* fFat; 
    G4Material* fWater; 
    G4Material* fBoneMatrix; 
    G4Material* fMineral;
    G4Material* fMedMat;
    G4Material* fPMMA;
    G4Material* fAdipose; 
    G4Material* fGlandular; 
    G4Material* fBreast5050; 
    G4Material* fcarcinoma; 
    G4Material* fLexan; 
    G4Material* fKapton; 
    G4Material* fNylon; 
    G4Material* fPolyethylene; 
    G4Material* fPolystyrene; 
    G4Material* fGrayMatter; 
    G4Material* fWhiteMatter;
        G4Material* fbeefBlood;
        G4Material* fFormaline;
        G4Material* fAcetone;
        G4Material* fHperoxide;
        G4Material* fCIRS3070;
        G4Material* fCIRS5050;
        G4Material* fCIRS7030;
        G4Material* fRMI454;        
        G4Material* fBone;
        G4Material* ffatLowX;
        G4Material* fbonematrixLowX;
        G4Material* fdryBoneLowX;
        G4Material* fliver;
        G4Material* fkidney;

        //custom material
        G4ExtendedMaterial* fCustomMat;
        G4double fCustomMatDensity;
        G4double fCustomMatHmassfract;
        G4double fCustomMatCmassfract;
        G4double fCustomMatNmassfract;
        G4double fCustomMatOmassfract;
        G4double fCustomMatNamassfract;
        G4double fCustomMatPmassfract;
        G4double fCustomMatSmassfract;
        G4double fCustomMatClmassfract;
        G4double fCustomMatKmassfract;
        G4double fCustomMatCamassfract;
        G4String fCustomMatFF;
        
        //definitions of variables for MedMat composition
        G4double fComp0,fComp1,fComp2,fComp3;
        
        //other materials
        G4Material* fAir;          
        G4Material* fTungsten;
        G4Material* fLead;
        G4Material* fGe;
    
    //-------------definition of Geometrical Variables---------
    //World
    G4double fWorldSize;    
          
    //Phantom
    G4double fPhantomDiameter;
    G4double fPhantomHeight;
    G4double fPhantomZ;     
        
    //setup angle (rad)
    G4double fthetaSetup;
    
    //Slits
    G4double fSlitSize;
    G4double fSlit1Thickness;        
    G4double fSlit2Thickness;        
    G4double fSlit3Thickness;        
    G4double fSlit4Thickness;        
    G4double fSlit1SampleDistance;                //center-center
    G4double fSlit2SampleDistance;
    G4double fSlit3SampleDistance;
    G4double fSlit4SampleDistance;
    G4double fSlit1xAperture;
    G4double fSlit2xAperture;
    G4double fSlit3xAperture;
    G4double fSlit4xAperture;
    G4double fSlit1yAperture;
    G4double fSlit2yAperture;
    G4double fSlit3yAperture;
    G4double fSlit4yAperture;

    //Detector
        G4double fDetectorThickness;
        G4double fDetectorSize;
        G4double fDetectorSampleDistance;        //center-center
        
        //Shielding
        G4double fShieldingThickness; 

        //-------------set methods for the messenger---------------
public:   
        void SetCustomMatFF(const G4String& ffname) {fCustomMatFF = ffname;}

        void SetCustomMatDensity(G4double csd) {fCustomMatDensity = csd;}
        void SetCustomMatHmassfract(G4double csHmf) {fCustomMatHmassfract = csHmf;}
        void SetCustomMatCmassfract(G4double csCmf) {fCustomMatCmassfract = csCmf;}
        void SetCustomMatNmassfract(G4double csNmf) {fCustomMatNmassfract = csNmf;}
        void SetCustomMatOmassfract(G4double csOmf) {fCustomMatOmassfract = csOmf;}
        void SetCustomMatNamassfract(G4double csNamf) {fCustomMatNamassfract = csNamf;}
        void SetCustomMatPmassfract(G4double csPmf) {fCustomMatPmassfract = csPmf;}
        void SetCustomMatSmassfract(G4double csSmf) {fCustomMatSmassfract = csSmf;}
        void SetCustomMatClmassfract(G4double csClmf) {fCustomMatClmassfract = csClmf;}
        void SetCustomMatKmassfract(G4double csKmf) {fCustomMatKmassfract = csKmf;}
        void SetCustomMatCamassfract(G4double csCamf) {fCustomMatCamassfract = csCamf;}
                   
        void SetPhantomMaterial(G4int mat) {fPhantomMaterialIndex = mat;}
        
        void SetPhantomDiameter(G4double diam) {fPhantomDiameter = diam;}
    void SetPhantomHeight(G4double ht) {fPhantomHeight = ht;}
    void SetPhantomZ(G4double PhZ) {fPhantomZ = PhZ;}
           
    void SetComp0(G4double c0) {fComp0 = c0;}
    void SetComp1(G4double c1) {fComp1 = c1;}
    void SetComp2(G4double c2) {fComp2 = c2;}
    void SetComp3(G4double c3) {fComp3 = c3;}
    
    void SetThetaSetup(G4double theta) {fthetaSetup = theta;}
    
    void SetSlits(G4bool bslits) {fIWantSlits = bslits;} 
    void SetSlit1Thickness(G4double sl1th) {fSlit1Thickness = sl1th;}
    void SetSlit2Thickness(G4double sl2th) {fSlit2Thickness = sl2th;}
    void SetSlit3Thickness(G4double sl3th) {fSlit3Thickness = sl3th;}
    void SetSlit4Thickness(G4double sl4th) {fSlit4Thickness = sl4th;}
    void SetSlit1SampleDistance(G4double slSampleDist1) 
            {fSlit1SampleDistance = slSampleDist1;}
    void SetSlit2SampleDistance(G4double slSampleDist2) 
            {fSlit2SampleDistance = slSampleDist2;}
    void SetSlit3SampleDistance(G4double slSampleDist3) 
            {fSlit3SampleDistance = slSampleDist3;}
    void SetSlit4SampleDistance(G4double slSampleDist4) 
            {fSlit4SampleDistance = slSampleDist4;}
    void SetSlit1xAperture(G4double aperture1x) {fSlit1xAperture = aperture1x;}
    void SetSlit2xAperture(G4double aperture2x) {fSlit2xAperture = aperture2x;}
    void SetSlit3xAperture(G4double aperture3x) {fSlit3xAperture = aperture3x;}
    void SetSlit4xAperture(G4double aperture4x) {fSlit4xAperture = aperture4x;}
    void SetSlit1yAperture(G4double aperture1y) {fSlit1yAperture = aperture1y;}
    void SetSlit2yAperture(G4double aperture2y) {fSlit2yAperture = aperture2y;}
    void SetSlit3yAperture(G4double aperture3y) {fSlit3yAperture = aperture3y;}
    void SetSlit4yAperture(G4double aperture4y) {fSlit4yAperture = aperture4y;}
    
    void SetDetectorSize(G4double detSize) {fDetectorSize = detSize;}
    void SetDetectorThickness(G4double detTh) {fDetectorThickness = detTh;}
    void SetDetectorSampleDistance(G4double detDist) 
            {fDetectorSampleDistance = detDist;}
    
};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

