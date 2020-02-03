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
/// \file persistency/gdml/G03/src/G03DetectorConstruction.cc
/// \brief Implementation of the G03DetectorConstruction class
//
//
//
// Class G03DetectorConstruction implementation
//
// ----------------------------------------------------------------------------

#include "G03DetectorConstruction.hh"

// Geant4 includes
//
#include "globals.hh"
#include "G4Material.hh"
#include "G4VPhysicalVolume.hh"

// Messenger
//
#include "G03DetectorMessenger.hh"

// Color extension include for reading
//
#include "G03ColorReader.hh"

// Color extension include for writing
//
#include "G03ColorWriter.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G03DetectorConstruction::G03DetectorConstruction()
  : G4VUserDetectorConstruction(),
    fAir(0), fAluminum(0), fPb(0), fXenon(0),
    fReader(0), fWriter(0), fParser(0),
    fDetectorMessenger(0)
{  
  fReadFile = "color_extension.gdml";
  fWriteFile = "color_extension_test.gdml";
  fWritingChoice = 1;

  fDetectorMessenger = new G03DetectorMessenger( this );

  fReader = new G03ColorReader;
  fWriter = new G03ColorWriter;
  fParser = new G4GDMLParser(fReader, fWriter);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G03DetectorConstruction::~G03DetectorConstruction()
{
  delete fDetectorMessenger;
  delete fReader;
  delete fWriter;
  delete fParser;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* G03DetectorConstruction::Construct()
{ 
  // Reading of Geometry from GDML

  G4VPhysicalVolume* fWorldPhysVol;

  fParser->Read(fReadFile,false);
    //
    // 2nd Boolean argument "Validate" set to false.
    // Disabling Schema validation for reading extended GDML file.
     
  // Prints the material information
  //
  G4cout << *(G4Material::GetMaterialTable() ) << G4endl;

  // Giving World Physical Volume from GDML Parser
  //
  fWorldPhysVol = fParser->GetWorldVolume();     

  if(fWritingChoice!=0)
  {
    fParser->Write(fWriteFile, fWorldPhysVol, true,
                  "./SimpleExtensionSchema/SimpleExtension.xsd");
  }
      
  return fWorldPhysVol;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G03DetectorConstruction::ListOfMaterials()
{
  G4double a;  // atomic mass
  G4double z;  // atomic number
  G4double density,temperature,pressure;
  G4double fractionmass;
  G4String name, symbol;
  G4int ncomponents;

  // Elements needed for the materials

  a = 14.01*g/mole;
  G4Element* elN = new G4Element(name="Nitrogen", symbol="N", z=7., a);

  a = 16.00*g/mole;
  G4Element* elO = new G4Element(name="Oxygen", symbol="O", z=8., a);
          
  a = 26.98*g/mole;
  G4Element* elAl = new G4Element(name="Aluminum", symbol="Al", z=13., a);
        
  // Print the Element information
  //
  G4cout << *(G4Element::GetElementTable()) << G4endl;

  // Air
  //
  density = 1.29*mg/cm3;
  fAir = new G4Material(name="Air", density, ncomponents=2);
  fAir->AddElement(elN, fractionmass=0.7);
  fAir->AddElement(elO, fractionmass=0.3);

  // Aluminum
  //
  density = 2.70*g/cm3;
  fAluminum = new G4Material(name="Aluminum", density, ncomponents=1);
  fAluminum->AddElement(elAl, fractionmass=1.0);

  // Lead
  //
  fPb = new G4Material("Lead", z=82., a= 207.19*g/mole, density= 11.35*g/cm3);

  // Xenon gas
  //
  fXenon = new G4Material("XenonGas", z=54., a=131.29*g/mole,
                         density= 5.458*mg/cm3, kStateGas,
                         temperature= 293.15*kelvin, pressure= 1*atmosphere);

  // Prints the material information
  //
  G4cout << *(G4Material::GetMaterialTable() ) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G03DetectorConstruction::SetReadFile( const G4String& fname )
{
  fReadFile=fname;
  fWritingChoice=0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G03DetectorConstruction::SetWriteFile( const G4String& fname )
{
  fWriteFile=fname;
  fWritingChoice=1;
}
