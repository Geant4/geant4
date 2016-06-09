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
// $Id: DetectorConstruction.cc,v 1.3 2009-04-15 13:26:26 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Class DetectorConstruction implementation
//
// ----------------------------------------------------------------------------

#include "DetectorConstruction.hh"

// Geant4 includes
//
#include "globals.hh"
#include "G4Material.hh"
#include "G4VPhysicalVolume.hh"

// Messenger
//
#include "DetectorMessenger.hh"

// Color extension include for reading
//
#include "ColorReader.hh"

// Color extension include for writing
//
#include "ColorWriter.hh"

// ----------------------------------------------------------------------------
//
// Constructor
//
DetectorConstruction::DetectorConstruction()
  : Air(0), Aluminum(0), Pb(0), Xenon(0)
{  
  fReadFile = "color_extension.gdml";
  fWriteFile = "color_extension_test.gdml";
  writingChoice = 1;

  detectorMessenger = new DetectorMessenger( this );

  reader = new ColorReader;
  writer = new ColorWriter;
  parser = new G4GDMLParser(reader, writer);
}

// ----------------------------------------------------------------------------
//
// Destructor
//
DetectorConstruction::~DetectorConstruction()
{
  delete detectorMessenger;
  delete reader;
  delete writer;
  delete parser;
}

// ----------------------------------------------------------------------------
//
// Constructs geometries and materials
//
G4VPhysicalVolume* DetectorConstruction::Construct()
{ 
  // Reading of Geometry from GDML

  G4VPhysicalVolume* fWorldPhysVol;

  parser->Read(fReadFile,false);
    //
    // 2nd Boolean argument "Validate" set to false.
    // Disabling Schema validation for reading extended GDML file.
     
  // Prints the material information
  //
  G4cout << *(G4Material::GetMaterialTable() ) << G4endl;

  // Giving World Physical Volume from GDML Parser
  //
  fWorldPhysVol = parser->GetWorldVolume();     

  if(writingChoice!=0)
  {
    parser->Write(fWriteFile, fWorldPhysVol, true,
                  "./SimpleExtensionSchema/SimpleExtension.xsd");
  }
      
  return fWorldPhysVol;
}

// ----------------------------------------------------------------------------
//
// Utility to build and list necessary materials
//
void DetectorConstruction::ListOfMaterials()
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
  Air = new G4Material(name="Air", density, ncomponents=2);
  Air->AddElement(elN, fractionmass=0.7);
  Air->AddElement(elO, fractionmass=0.3);

  // Aluminum
  //
  density = 2.70*g/cm3;
  Aluminum = new G4Material(name="Aluminum", density, ncomponents=1);
  Aluminum->AddElement(elAl, fractionmass=1.0);

  // Lead
  //
  Pb = new G4Material("Lead", z=82., a= 207.19*g/mole, density= 11.35*g/cm3);

  // Xenon gas
  //
  Xenon = new G4Material("XenonGas", z=54., a=131.29*g/mole,
                         density= 5.458*mg/cm3, kStateGas,
                         temperature= 293.15*kelvin, pressure= 1*atmosphere);

  // Prints the material information
  //
  G4cout << *(G4Material::GetMaterialTable() ) << G4endl;
}


// ----------------------------------------------------------------------------
//
// SetReadFile
//
void DetectorConstruction::SetReadFile( const G4String& fname )
{
  fReadFile=fname;
  writingChoice=0;
}

// ----------------------------------------------------------------------------
//
// SetWriteFile
//
void DetectorConstruction::SetWriteFile( const G4String& fname )
{
  fWriteFile=fname;
  writingChoice=1;
}
