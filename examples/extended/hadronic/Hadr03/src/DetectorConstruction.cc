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
/// \file hadronic/Hadr03/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//

//
// $Id: DetectorConstruction.cc,v 1.8 2007-11-12 15:48:58 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UniformMagField.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:fPBox(0), fLBox(0), fMaterial(0), fMagField(0)
{
  fBoxSize = 10*m;
  DefineMaterials();
  SetMaterial("Molybdenum98");  
  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ delete fDetectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
 // define a Material from isotopes
 //
 G4int ncomponents;
 G4double abundance, massfraction;

 // Boron10
 // 
 G4Isotope* B10 = new G4Isotope("B10", 5, 10);
 //
 G4Element* b10  = new G4Element("Boron10","B10",ncomponents=1);
 b10->AddIsotope(B10, abundance= 100.*perCent);
 //
 G4Material* boron10 
         = new G4Material("Boron10", 2.46*g/cm3, ncomponents=1);
 boron10->AddElement(b10, massfraction=100.*perCent);

 // Boron11
 // 
 G4Isotope* B11 = new G4Isotope("B11", 5, 11);
 //
 G4Element* b11  = new G4Element("Boron11","B11",ncomponents=1);
 b11->AddIsotope(B11, abundance= 100.*perCent);
 //
 G4Material* boron11 
         = new G4Material("Boron11", 2.46*g/cm3, ncomponents=1);
 boron11->AddElement(b11, massfraction=100.*perCent);
 
 // Carbon12
 // 
 ////G4Isotope* C12 = new G4Isotope("C12", 6, 12);
 //
 ////G4Element* c12  = new G4Element("Carbon12","C12",ncomponents=1);
 ////c12->AddIsotope(C12, abundance= 100.*perCent);
 //
 ////G4Material* carbon12 
 ////        = new G4Material("Carbon12", 2.27*g/cm3, ncomponents=1);
 ////carbon12->AddElement(c12, massfraction=100.*perCent);
  
 // Oxygen16
 // 
 G4Isotope* O16 = new G4Isotope("O16", 8, 16);
 //
 G4Element* o16  = new G4Element("Oxygen16","O16",ncomponents=1);
 o16->AddIsotope(O16, abundance= 100.*perCent);
 //
 G4Material* oxygen16
         = new G4Material("Oxygen16", 1.429*g/cm3, ncomponents=1);
 oxygen16->AddElement(o16, massfraction=100.*perCent);
 
 // Calcium40
 // 
 G4Isotope* Ca40 = new G4Isotope("Ca40", 20, 40);
 //
 G4Element* ca40  = new G4Element("Calcium40","Ca40",ncomponents=1);
 ca40->AddIsotope(Ca40, abundance= 100.*perCent);
 //
 G4Material* calcium40
         = new G4Material("Calcium40", 1.55*g/cm3, ncomponents=1);
 calcium40->AddElement(ca40, massfraction=100.*perCent);
 
 // Calcium48
 // 
 G4Isotope* Ca48 = new G4Isotope("Ca48", 20, 48);
 //
 G4Element* ca48  = new G4Element("Calcium48","Ca48",ncomponents=1);
 ca48->AddIsotope(Ca48, abundance= 100.*perCent);
 //
 G4Material* calcium48
         = new G4Material("Calcium48", 1.55*g/cm3, ncomponents=1);
 calcium48->AddElement(ca48, massfraction=100.*perCent);
   
 // Zirconium90
 // 
 G4Isotope* Zr90 = new G4Isotope("Zr90", 40, 90);
 //
 G4Element* zr90  = new G4Element("Zirconium90","Zr90",ncomponents=1);
 zr90->AddIsotope(Zr90, abundance= 100.*perCent);
 //
 G4Material* zirconium90 
         = new G4Material("Zirconium90", 6.51*g/cm3, ncomponents=1);
 zirconium90->AddElement(zr90, massfraction=100.*perCent); 
  
 // Molybdenum98
 // 
 G4Isotope* Mo98 = new G4Isotope("Mo98", 42, 98);
 //
 G4Element* mo98  = new G4Element("Molybdenum98","Mo98",ncomponents=1);
 mo98->AddIsotope(Mo98, abundance= 100.*perCent);
 //
 G4Material* molybden98 
         = new G4Material("Molybdenum98", 10.28*g/cm3, ncomponents=1);
 molybden98->AddElement(mo98, massfraction=100.*perCent);
 
 // Molybdenum100
 //  
 G4Isotope* Mo100 = new G4Isotope("Mo100", 42, 100);
 //
 G4Element* mo100  = new G4Element("Molybdenum100","Mo100",ncomponents=1);
 mo100->AddIsotope(Mo100, abundance= 100.*perCent);
 //
 G4Material* molybden100 
         = new G4Material("Molybdenum100", 10.28*g/cm3, ncomponents=1);
 molybden100->AddElement(mo100, massfraction=100.*perCent);
 
 // Lead208
 //  
 G4Isotope* Pb208 = new G4Isotope("Pb208", 82, 208);
 //
 G4Element* pb208  = new G4Element("Lead208","Pb208",ncomponents=1);
 pb208->AddIsotope(Pb208, abundance= 100.*perCent);
 //
 G4Material* lead208 
         = new G4Material("Lead208", 11.34*g/cm3, ncomponents=1);
 lead208->AddElement(pb208, massfraction=100.*perCent); 
 
 // Uranium 235
 //  
 G4Isotope* U235 = new G4Isotope("U235", 92, 235);
 //
 G4Element* u235  = new G4Element("Uranium235","U235",ncomponents=1);
 u235->AddIsotope(U235, abundance= 100.*perCent);
 //
 G4Material* uranium235 
         = new G4Material("Uranium235", 19.05*g/cm3, ncomponents=1);
 uranium235->AddElement(u235, massfraction=100.*perCent);
    
 // Uranium 238
 //  
 G4Isotope* U238 = new G4Isotope("U238", 92, 238);
 //
 G4Element* u238  = new G4Element("Uranium238","U238",ncomponents=1);
 u238->AddIsotope(U238, abundance= 100.*perCent);
 //
 G4Material* uranium238 
         = new G4Material("Uranium238", 19.05*g/cm3, ncomponents=1);
 uranium238->AddElement(u238, massfraction=100.*perCent);
    
 // or use G4-NIST materials data base
 //
 G4NistManager* man = G4NistManager::Instance();
 man->FindOrBuildMaterial("G4_B");

 ///G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  G4Box*
  sBox = new G4Box("Container",                         //its name
                   fBoxSize/2,fBoxSize/2,fBoxSize/2);   //its dimensions

  fLBox = new G4LogicalVolume(sBox,                     //its shape
                             fMaterial,                 //its material
                             fMaterial->GetName());     //its name

  fPBox = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(),            //at (0,0,0)
                            fLBox,                      //its logical volume                           
                            fMaterial->GetName(),       //its name
                            0,                          //its mother  volume
                            false,                      //no boolean operation
                            0);                         //copy number
                           
  PrintParameters();
  
  //always return the root volume
  //
  return fPBox;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters()
{
  G4cout << "\n The Box is " << G4BestUnit(fBoxSize,"Length")
         << " of " << fMaterial->GetName() 
         << "\n \n" << fMaterial << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaterial(G4String materialChoice)
{
  // search the material by its name
  ////G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
  G4Material* pttoMaterial = 
     G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);
  
  if (pttoMaterial) { fMaterial = pttoMaterial;
    } else {
    G4cout << "\n--> warning from DetectorConstruction::SetMaterial : "
           << materialChoice << " not found" << G4endl;
  }              
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSize(G4double value)
{
  fBoxSize = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

void DetectorConstruction::SetMagField(G4double fieldValue)
{
  //apply a global uniform magnetic field along Z axis
  G4FieldManager* fieldMgr
   = G4TransportationManager::GetTransportationManager()->GetFieldManager();

  if (fMagField) delete fMagField;        //delete the existing magn field

  if (fieldValue!=0.)                        // create a new one if non nul
    {
      fMagField = new G4UniformMagField(G4ThreeVector(0.,0.,fieldValue));
      fieldMgr->SetDetectorField(fMagField);
      fieldMgr->CreateChordFinder(fMagField);
    }
   else
    {
      fMagField = 0;
      fieldMgr->SetDetectorField(fMagField);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4RunManager.hh"

void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructVolumes());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
