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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "SensitiveDetector.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Transform3D.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4AnalysisManager.hh"
#include "G4RegionStore.hh"
#include "G4VisAttributes.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

#include "G4SubtractionSolid.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4PSDoseDeposit.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4SDParticleFilter.hh"

#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include <G4ChordFinder.hh>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
{
    //instantiate the messenger
    fMessenger = new DetectorConstructionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    G4cout << G4endl << "### DetectorConstruction::Construct() ###" << G4endl 
           << G4endl;
    
    
    //sanity check
    if (fHybridSource) {
        G4cout << "This is a hybrid positron source!" << G4endl << G4endl;
    } else {
        fConverter = false;
        fRadiatorConverterSepDistance = 0.;
        G4cout << "This is a single-crystal positron source!" << G4endl 
               << G4endl;    
    }
    

    //check overlap option
    G4bool checkOverlaps = true;


    //Materials
    G4NistManager* nist = G4NistManager::Instance();
    G4Material* Silicon = nist->FindOrBuildMaterial("G4_Si");  
    G4Material* PWO = nist->FindOrBuildMaterial("G4_PbWO4");
    G4Material* Diamond = nist->FindOrBuildMaterial("G4_C");    
    G4Material* Tungsten = nist->FindOrBuildMaterial("G4_W");
    G4Material* Iridium = nist->FindOrBuildMaterial("G4_W");
    G4Material* Germanium = nist->FindOrBuildMaterial("G4_Ge");
    G4Element* elBi = nist->FindOrBuildElement("Bi");    
    G4Element* elGe = nist->FindOrBuildElement("Ge");
    G4Element* elO = nist->FindOrBuildElement("O");
    G4Material* BGO = new G4Material("G4_BGO", 7.13*g/cm3, 3);
    BGO->AddElement(elBi, 0.671054);
    BGO->AddElement(elGe, 0.17482);
    BGO->AddElement(elO, 0.154126);
  
    //Vacuum
    G4double z = 7.;
    G4double a = 14.007 * CLHEP::g/CLHEP::mole;
    G4double density = CLHEP::universe_mean_density;
    G4double pressure = 1.E-6 * 1.E-3 * CLHEP::bar; //10-6 mbar
    G4double temperature = 300. * CLHEP::kelvin; //300 K
    G4Material* Vacuum = new G4Material("Vacuum", 
                                        z, 
                                        a, 
                                        density, 
                                        kStateGas, 
                                        temperature, 
                                        pressure);
  

    // ------------------- World -------------------------
    G4Box* solidWorld = new G4Box("World", 3.*m, 3.*m, 20.*m);
    
    G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, 
                                                      Vacuum, 
                                                      "World");
                                                      
    logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());  
    
    G4VPhysicalVolume* physWorld = new G4PVPlacement
                                   (0,              // no rotation
                                   G4ThreeVector(), // centre position
                                   logicWorld,      // its logical volume
                                   "World",         // its name
                                   0,               // its mother volume
                                   false,           // no boolean operation
                                   0,               // copy number
                                   checkOverlaps);  // overlaps checking


    // --------------- Crystal (Radiator) -------------------------    
    //set visualization attributes
    G4VisAttributes* CrystalVisAttribute = 
        new G4VisAttributes(G4Colour(0., 0., 1., 1.));  
    CrystalVisAttribute->SetForceSolid(true);
    
    //select crystal material
    if (fCrystalMaterialStr == "PWO") {
        fCrystalMaterial = PWO;
    } else if (fCrystalMaterialStr == "BGO") {
        fCrystalMaterial = BGO;
    } else if (fCrystalMaterialStr == "C") {
        fCrystalMaterial = Diamond;
    } else if (fCrystalMaterialStr == "W") {
        fCrystalMaterial = Tungsten;
    } else if (fCrystalMaterialStr == "Ir") {
        fCrystalMaterial = Iridium;
    } else if (fCrystalMaterialStr == "Ge") {
        fCrystalMaterial = Germanium;
    } else {
        fCrystalMaterial = Silicon;
    }
        
    //Crystal solid and logic volumes
    G4Box* crystalSolid = new G4Box("Crystal",
                                    fCrystalSize.x()*0.5,
                                    fCrystalSize.y()*0.5,
                                    fCrystalSize.z()*0.5);
    
    fCrystalLogic = new G4LogicalVolume(crystalSolid,
                                        fCrystalMaterial,
                                        "Crystal"); 
                                  
    //Set actually or not the radiatior crystal
    fCrystalLogic->SetVisAttributes(CrystalVisAttribute);
    
    //Crystal position
    G4double tollfCrystalZ = 0.*mm;
    if (fAngleX > 0 || fAngleY > 0) {
        tollfCrystalZ = 0.15*mm; 
        //this allows us to tolerate misalignements up to 30 mrad on 1 cm thick crystals
    }
    fCrystalZ = -fCrystalSize.z()*0.5 - tollfCrystalZ;
    G4ThreeVector posCrystal = G4ThreeVector(0.*mm, 0.*mm, fCrystalZ);       
       
    //Crystal rotation angle (also the angle of crystal planes vs the beam)
    G4RotationMatrix* crystalRotationMatrix = new G4RotationMatrix;
    crystalRotationMatrix->rotateY(-fAngleX);
    crystalRotationMatrix->rotateX(-fAngleY);

    //Crystal placement
    new G4PVPlacement(crystalRotationMatrix, 
                      posCrystal,
                      fCrystalLogic,
                      "Crystal",
                      logicWorld,
                      false,
                      0,
                      checkOverlaps);
    
    //Crystal region (necessary for the FastSim model)
    fCrystalRegion = new G4Region("Crystal");
    fCrystalRegion->AddRootLogicalVolume(fCrystalLogic);
        
    //Print Crystal info
    G4cout << "Radiator Crystal set!" << G4endl;
    G4cout << "Crystal material: " << fCrystalMaterial->GetName() << G4endl;
    G4cout << "Crystal size: " << fCrystalSize.x()/mm 
                               << "x" << fCrystalSize.y()/mm
           << "x" << fCrystalSize.z()/mm << " mm3" << G4endl;
    G4cout << "RadiatorZ: " << fCrystalZ/mm << " mm" << G4endl;
    G4cout << G4endl;    
    

    // -------- volume for the magnetic field ------------
    if (fSetMagneticField && 
        fRadiatorConverterSepDistance > fFieldRegionLength + 0.1*mm) {
        G4double magFieldRegionZ = fRadiatorConverterSepDistance*0.5;
  
        G4Tubs* tub1 = new G4Tubs("MFvolume",
                                  0.*cm,
                                  fFieldRegionLength*0.5,
                                  fFieldRegionLength*0.5,
                                  0.*deg,
                                  360.*deg);
  
        G4RotationMatrix *xRot = new G4RotationMatrix;
        xRot->rotateX(90.*deg);
        
        fMFlogic = new G4LogicalVolume(tub1, 
                                       Vacuum, 
                                       "MFvolume");
                                                       
        G4VisAttributes* MFregionVisAttribute = 
            new G4VisAttributes(G4Colour(1., 0., 0., 0.35));  
        MFregionVisAttribute->SetForceSolid(true);
        fMFlogic->SetVisAttributes(MFregionVisAttribute);
        new G4PVPlacement(xRot,
                          G4ThreeVector(0., 0., magFieldRegionZ),
                          fMFlogic,
                          "MFvolume",
                          logicWorld,
                          false,
                          0,
                          checkOverlaps);
                          
        G4cout << "Magnetic Field set!" << G4endl;
        G4cout << "Field value: " << fFieldValue/tesla << " T" << G4endl;
        G4cout << "Field Region Diameter: " << fFieldRegionLength/mm 
               << " mm" << G4endl;
        G4cout << "Field Region Length: " << fFieldRegionLength/mm 
               << " mm" << G4endl;
        G4cout << G4endl; 
        
        fSetCollimator = false;
    }
    

    // ------------------ Collimator ---------------------
    else if (fSetCollimator && 
             fRadiatorConverterSepDistance > fCollimatorThickness + 0.1*mm) {
        G4double CollimatorZ = fRadiatorCollimatorSepDistance +
                               fCollimatorThickness*0.5;
          
        G4Box* outerCollimSolid = new G4Box("outerCollimSolid",
                                            fCollimatorSide*0.5,
                                            fCollimatorSide*0.5,
                                            fCollimatorThickness*0.5);
        
        G4VSolid* innerCollimSolid;
        if (fCollimatorHole == "circular") {
            innerCollimSolid = new G4Tubs("innerCollimSolid",
                                          0.*cm,
                                          fCollimatorAperture*0.5,
                                          fCollimatorThickness*0.51,
                                          0.*deg,
                                          360.*deg);
        } else {
            innerCollimSolid = new G4Box("innerCollimSolid",
                                           fCollimatorAperture*0.5,
                                         fCollimatorAperture*0.5,
                                         fCollimatorThickness*0.6);        
        }
                                                   
        G4SubtractionSolid* CollimSolid = new G4SubtractionSolid("Collimator", 
                                                                 outerCollimSolid, 
                                                                 innerCollimSolid);
                                                                 
        fCollimatorLogic = new G4LogicalVolume(CollimSolid, 
                                               Tungsten, 
                                               "Collimator");        
                                                                 
        G4VisAttributes* CollimatorVisAttribute = 
            new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));  
        CollimatorVisAttribute->SetForceSolid(true);
        fCollimatorLogic->SetVisAttributes(CollimatorVisAttribute);
        new G4PVPlacement(0,
                          G4ThreeVector(0., 0., CollimatorZ),
                          fCollimatorLogic,
                          "Collimator",
                          logicWorld,
                          false,
                          0,
                          checkOverlaps);
                          
        G4cout << "Collimator set!" << G4endl;
        G4cout << "Collimator aperture: " << fCollimatorAperture/mm 
               << " mm" << G4endl;
        G4cout << "Collimator hole shape: " << fCollimatorHole << G4endl;
        G4cout << "Collimator thickness: " << fCollimatorThickness/mm 
               << " mm" << G4endl;
        G4cout << "CollimatorZ: " << CollimatorZ/mm << " mm" << G4endl;
        G4cout << G4endl;            
    } 
    

    // -------- Converter (Target) randomly oriented -----
    //visualization attributes
    G4VisAttributes* ConverterVisAttribute = 
        new G4VisAttributes(G4Colour(0., 0.2, 0.8, 0.8));  
    
    //select Converter material
    G4Material* sphereMaterial;
    if (fConverterMaterialStr == "PWO") {
        fConverterMaterial = PWO;
        sphereMaterial = PWO;
    } else if (fConverterMaterialStr == "BGO") {
        fConverterMaterial = BGO;
        sphereMaterial = BGO;        
    } else if (fConverterMaterialStr == "Ir") {
        fConverterMaterial = Iridium;
        sphereMaterial = Iridium;
    } else {
        fConverterMaterial = Tungsten; //defualt is W
        sphereMaterial = Tungsten;
    }    
    
    if (fGranularConverter) {
        fConverterMaterial = Vacuum; 
        //Converter will be a vacuum box filled with spheres (of fConverterMaterial)
        ConverterVisAttribute->SetForceSolid(false);
    } else {
        ConverterVisAttribute->SetForceSolid(true);
    }
    
    //set Converter size
    G4double ConverterWidth = fConverterSize.x();
    G4double ConverterHeight = fConverterSize.y();
    G4double ConverterThickness = fConverterSize.z();
    
    //Converter position
    G4double toll = fVirtualDetectorSize.z();    
    if (fRadiatorConverterSepDistance > 0) {
        toll = 2.*fVirtualDetectorSize.z();
    }    
    fConverterZ = fRadiatorConverterSepDistance + ConverterThickness*0.5 +
                  fVirtualDetectorSize.z() + toll;
    G4ThreeVector posConverter = G4ThreeVector(0.*mm, 0.*mm, fConverterZ);
                                                 
    //Converter volume
    G4Box* ConverterSolid = new G4Box("Converter",
                                      ConverterWidth*0.5, 
                                      ConverterHeight*0.5, 
                                      ConverterThickness*0.5);
                                   
    fConverterLogic = new G4LogicalVolume(ConverterSolid, 
                                          fConverterMaterial,
                                          "Converter");
                                                                                            
    fConverterLogic->SetVisAttributes(ConverterVisAttribute);
    
    if (fConverter) {
        new G4PVPlacement(0, 
                          posConverter,
                          fConverterLogic,
                          "Converter",
                          logicWorld,
                          false,
                          0,
                          checkOverlaps);                   
        
        //print Converter info
        G4cout << "RadiatorConverterSepDistance: " << fRadiatorConverterSepDistance 
               << " mm " << G4endl;
        G4cout << "ConverterZ: " << fConverterZ/mm << " mm" << G4endl;  
        G4cout << "Converter material: " << fConverterMaterial->GetName() << G4endl;
        G4cout << "Converter size: " << ConverterWidth/mm 
               << "x" << ConverterHeight/mm
               << "x" << ConverterThickness/mm << " mm3" << G4endl;    
        G4cout << G4endl;        
            
        //Is it a glanular converter?
        if (fGranularConverter) {
            G4cout << "Granular Converter set!" << G4endl;
            G4cout << "sphere material: " << sphereMaterial->GetName() << G4endl;
            G4cout << "sphere radius: " << fSphereRadius/mm << " mm" << G4endl;
            
            G4int numLayers = G4int(ConverterThickness/(std::sqrt(2.)*fSphereRadius));
            G4int NintX = G4int(ConverterWidth/(2.*fSphereRadius)) - 0;
            G4int NintY = G4int(ConverterHeight/(2.*fSphereRadius)) - 0;    
            G4double spacing = 2*fSphereRadius;
            
            //numLayers = 10;
            //NintX = NintY = 10;
            
            G4cout << "numLayers: " << numLayers 
                   << ", NintX: " << NintX 
                   << ", NintY: " << NintY << G4endl;
            
            G4Sphere* sphereSolid = new G4Sphere("Sphere", 
                                                 0, 
                                                 fSphereRadius, 
                                                 0, 
                                                 360*deg, 
                                                 0, 
                                                 180*deg);
                                                      
            G4VisAttributes* sphereVisAtt = 
                new G4VisAttributes(G4Colour(0., 0., 1., 0.7)); 
            sphereVisAtt->SetForceSolid(true);
            
            G4int k = 0;
            for (int layer = 0; layer < numLayers; layer++) {                
                //alternating number of spheres in each layer so as
                //they are even for odd layers and odd for even layers.
                G4int numSpheresX = (layer % 2 == 0) ? 
                    (NintX % 2 == 0) ? NintX-1 : NintX : 
                    (NintX % 2 == 0) ? NintX : NintX-1; 
                G4int numSpheresY = (layer % 2 == 0) ?
                    (NintY % 2 == 0) ? NintY-1 : NintY : 
                    (NintY % 2 == 0) ? NintY : NintY-1;  
                    
                //calculate the total width and height of the sphere arrangement in each layer
                G4double totalWidth = numSpheresX*(2.*fSphereRadius);
                G4double totalHeight = numSpheresY*(2.*fSphereRadius);        

                //calculate the starting position for the first sphere in each layer
                G4double startX = -totalWidth/2. + fSphereRadius;
                G4double startY = -totalHeight/2. + fSphereRadius;
                
                //sphere positioning
                for (int i = 0; i < numSpheresX; ++i) {
                    for (int j = 0; j < numSpheresY; ++j) {            
                        std::string sphereName = "Sphere_" + std::to_string(layer) + 
                                                       "_" + std::to_string(i) + 
                                                       "_" + std::to_string(j);
                                    
                        fSphereLogic[k] = new G4LogicalVolume(sphereSolid, 
                                                              sphereMaterial,
                                                              sphereName);
                        fSphereLogic[k]->SetVisAttributes(sphereVisAtt);
                        fScoringVolume.push_back(fSphereLogic[k]);

                        G4double xPosition = startX + i*spacing;
                        G4double yPosition = startY + j*spacing;
                        G4double zPosition = -ConverterThickness/2. + fSphereRadius 
                                             + layer*std::sqrt(2)*fSphereRadius;
                        G4ThreeVector position(xPosition, yPosition, zPosition);
                        
                        new G4PVPlacement(0, 
                                          position, 
                                          fSphereLogic[k], 
                                          sphereName, 
                                          fConverterLogic, 
                                          false, 
                                          k, 
                                          checkOverlaps);
                        
                        k++;
                    }
                }
            }
                    
            fNSpheres = k;
            G4cout << "Nspheres positioned: " << fNSpheres << G4endl << G4endl;
            
        }
    
    }
    
    
    // --------------- virtual Detectors -----------------
    //position
    G4double VirtualDetector0Z = fConverterZ - ConverterThickness*0.5 
                                             - fVirtualDetectorSize.z()*0.5;      
    G4ThreeVector posVirtualDetector0 = G4ThreeVector(0, 0, VirtualDetector0Z);    
    G4ThreeVector frontVirtualDetector0 = G4ThreeVector(0, 0, VirtualDetector0Z - 
                                                              fVirtualDetectorSize.z()*0.5);
    G4cout << "VirtualDetector0Z: " << VirtualDetector0Z/mm << " mm" << G4endl;
    fVirtualDetectorPositionVector.push_back(frontVirtualDetector0);
    
    G4double VirtualDetector1Z = fConverterZ + ConverterThickness*0.5 
                                             + fVirtualDetectorSize.z()*0.5;
    G4ThreeVector posVirtualDetector1 = G4ThreeVector(0, 0, VirtualDetector1Z);
    G4ThreeVector frontVirtualDetector1 = G4ThreeVector(0, 0, VirtualDetector1Z - 
                                                              fVirtualDetectorSize.z()*0.5);
    if (fHybridSource) {
        G4cout << "VirtualDetector1Z: " << VirtualDetector1Z/mm << " mm" << G4endl;
        fVirtualDetectorPositionVector.push_back(frontVirtualDetector1);      
    }
                                                  
    G4double VirtualDetector2Z = fVirtualDetectorSize.z()*0.5;
    G4ThreeVector posVirtualDetector2 = G4ThreeVector(0, 0, VirtualDetector2Z);
    G4ThreeVector frontVirtualDetector2 = G4ThreeVector(0, 0, VirtualDetector2Z - 
                                                              fVirtualDetectorSize.z()*0.5);
    if (fRadiatorConverterSepDistance > 0) {
        G4cout << "VirtualDetector2Z: " << VirtualDetector2Z/mm << " mm" << G4endl;
        fVirtualDetectorPositionVector.push_back(frontVirtualDetector2);    
    }

    //virtual Detector volume
    G4Box* VirtualDetectorSolid = new G4Box("VirtualDetector",
                                            fVirtualDetectorSize.x()*0.5, 
                                            fVirtualDetectorSize.y()*0.5, 
                                            fVirtualDetectorSize.z()*0.5);
    
    fVirtualDetectorLogic0 = new G4LogicalVolume(VirtualDetectorSolid, 
                                                 Vacuum,
                                                 "VirtualDetector0");
                                                                
    fVirtualDetectorLogic1 = new G4LogicalVolume(VirtualDetectorSolid, 
                                                 Vacuum,
                                                 "VirtualDetector1");
    
    fVirtualDetectorLogic2 = new G4LogicalVolume(VirtualDetectorSolid, 
                                                 Vacuum,
                                                 "VirtualDetector2");
                                                                
    G4VisAttributes* VirtualDetectorVisAttribute = 
        new G4VisAttributes(G4Colour(1., 1., 1.));  
    VirtualDetectorVisAttribute->SetForceSolid(false);
    fVirtualDetectorLogic0->SetVisAttributes(VirtualDetectorVisAttribute);
    fVirtualDetectorLogic1->SetVisAttributes(VirtualDetectorVisAttribute);
    fVirtualDetectorLogic2->SetVisAttributes(VirtualDetectorVisAttribute);
       
    new G4PVPlacement(0, 
                      posVirtualDetector0, 
                      fVirtualDetectorLogic0, 
                      "VirtualDetector0", 
                      logicWorld, 
                      false, 
                      0, 
                      checkOverlaps);
                      
    if (fHybridSource) {                 
        new G4PVPlacement(0, 
                          posVirtualDetector1, 
                          fVirtualDetectorLogic1, 
                          "VirtualDetector1", 
                          logicWorld, 
                          false, 
                          1, 
                          checkOverlaps);      
    }
  
    if (fRadiatorConverterSepDistance > 0) {             
        new G4PVPlacement(0, 
                          posVirtualDetector2, 
                          fVirtualDetectorLogic2, 
                          "VirtualDetector2", 
                          logicWorld, 
                          false, 
                          2, 
                          checkOverlaps);
    }
       

    //always return the physical World
    return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
    //activate OC effects through FastSim model
    if (fActivateOCeffects) {    
        G4RegionStore* regionStore = G4RegionStore::GetInstance();
        G4Region* RegionCh = regionStore->GetRegion("Crystal");

        G4ChannelingFastSimModel* ChannelingModel = 
            new G4ChannelingFastSimModel("ChannelingModel", RegionCh);

        ChannelingModel->Input(fCrystalLogic->GetMaterial(), fLattice, fPotentialPath);
        ChannelingModel->GetCrystalData()->SetBendingAngle(fBendingAngle, fCrystalLogic);

        G4double fParticleLEth = 1.*GeV; //deafult 200.*MeV (5.*GeV -> much faster)
        G4double fLindhardAngles = 10; //default 100
        ChannelingModel->SetLowKineticEnergyLimit(fParticleLEth, "e-"); 
        ChannelingModel->SetLowKineticEnergyLimit(fParticleLEth, "e+");
        ChannelingModel->SetLindhardAngleNumberHighLimit(fLindhardAngles, "e-");
        ChannelingModel->SetLindhardAngleNumberHighLimit(fLindhardAngles, "e+"); 
        
        G4cout << G4endl;
        G4cout << "Oriented Crystal effects set through FastSim model" << G4endl;
        G4cout << "Crystal bending angle: " << fBendingAngle << " rad" << G4endl;
        G4cout << "Crystal Lattice: " << fLattice << G4endl;
        G4cout << "Crystal AngleX: " << fAngleX << " rad" << G4endl;   
        G4cout << "Crystal AngleY: " << fAngleY << " rad" << G4endl;
        G4cout << "fParticleLEth: " << fParticleLEth/MeV << " MeV" << G4endl;
        G4cout << "fLindhardAngles: " << fLindhardAngles << G4endl;
        G4cout << "ActivateRadiationModel: " << fActivateRadiationModel << G4endl;
        
        if (fActivateRadiationModel) {
            ChannelingModel->RadiationModelActivate();
            G4int fSamplingPhotonsNumber = 150; //default 150
            G4int fNSmallTrajectorySteps = 10000; //default 10000
            G4double fRadiactionAngleFactor = 4.; //deafult 4
            G4double fSinglePhotonRadProbLimit = 0.25; //default 0.25
            G4double fLEthreshold = 1.*MeV;
            ChannelingModel->GetRadiationModel()->
                SetSamplingPhotonsNumber(fSamplingPhotonsNumber);
            ChannelingModel->GetRadiationModel()->
                SetNSmallTrajectorySteps(fNSmallTrajectorySteps);
            ChannelingModel->GetRadiationModel()->
                SetRadiationAngleFactor(fRadiactionAngleFactor);
            ChannelingModel->GetRadiationModel()
                ->SetSinglePhotonRadiationProbabilityLimit(fSinglePhotonRadProbLimit);
            ChannelingModel->GetRadiationModel()
                ->SetSpectrumEnergyRange(fLEthreshold, 20.*GeV, 100);
            
            G4cout << "SamplingPhotonsNumber: " 
                   << fSamplingPhotonsNumber << G4endl;
            G4cout << "NSmallTrajectorySteps: " 
                   << fNSmallTrajectorySteps << G4endl;                    
            G4cout << "fRadiactionAngleFactor: " 
                   << fRadiactionAngleFactor << G4endl;
            G4cout << "fSinglePhotonRadProbLimit: " 
                   << fSinglePhotonRadProbLimit << G4endl;    
            G4cout << "Low Eenergy threshold to emit photons and record their energy: " 
                   << fLEthreshold/MeV << " MeV" << G4endl << G4endl;               
        } else {
            G4cout << G4endl;
        }   
    }
    

    //built-in Edep Scorer in the Radiator Crystal
    G4MultiFunctionalDetector* multisd = new G4MultiFunctionalDetector("multisd");
    G4VPrimitiveScorer* edepscorer = new G4PSEnergyDeposit("edep");
    multisd->RegisterPrimitive(edepscorer);
    SetSensitiveDetector(fCrystalLogic->GetName(), multisd);
    G4SDManager::GetSDMpointer()->AddNewDetector(multisd);

    
    //Sensitive Volumes (Virtual Detectors)
    G4VSensitiveDetector* vDetector = new SensitiveDetector("det");
    G4SDManager::GetSDMpointer()->AddNewDetector(vDetector);   
    fVirtualDetectorLogic0->SetSensitiveDetector(vDetector);
    if (fHybridSource) { 
        fVirtualDetectorLogic1->SetSensitiveDetector(vDetector);
    }
    if (fRadiatorConverterSepDistance > 0) {
        fVirtualDetectorLogic2->SetSensitiveDetector(vDetector);
    }        

    
    //Magnetic field
    if (fSetMagneticField & fHybridSource) {
        G4UniformMagField* myField = 
            new G4UniformMagField(G4ThreeVector(0., fFieldValue, 0.));
        G4FieldManager* localfieldMgr = new G4FieldManager(myField);
        localfieldMgr->CreateChordFinder(myField);
        fMFlogic->SetFieldManager(localfieldMgr, true);
    }
    
    
    G4cout << "### End of DetectorConstruction ###" << G4endl << G4endl << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

