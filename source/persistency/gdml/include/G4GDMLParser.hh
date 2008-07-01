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
// $Id: G4GDMLParser.hh,v 1.44 2008-07-01 08:18:02 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4GDMLParser
//
// Class description:
//
// GDML main parser.

// History:
// - Created.                                  Zoltan Torzsok, November 2007
// -------------------------------------------------------------------------
 
#ifndef _G4GDMLPARSER_INCLUDED_
#define _G4GDMLPARSER_INCLUDED_

// #include "G4STRead.hh"

#include "G4GDMLReadStructure.hh"
#include "G4GDMLWriteStructure.hh"

#define G4GDML_DEFAULT_SCHEMALOCATION G4String("http://cern.ch/service-spi/app/releases/GDML/GDML_2_10_0/src/GDMLSchema/gdml.xsd")

class G4GDMLParser {
   G4GDMLReadStructure reader;
   G4GDMLWriteStructure writer;
public:
   G4GDMLParser() { xercesc::XMLPlatformUtils::Initialize(); }
   ~G4GDMLParser() { xercesc::XMLPlatformUtils::Terminate(); }

   void Read(const G4String& name,bool Validate=true) { 
      
      const bool IsModule = false;
      reader.Read(name,Validate,IsModule); 
   }

   void Write(const G4String& name,const G4LogicalVolume* const logvol,const G4String& SchemaLocation=G4GDML_DEFAULT_SCHEMALOCATION) { 
   
      const G4int depth = 0;
      writer.Write(name,logvol,SchemaLocation,depth);
   }
/*
   G4VPhysicalVolume* ReadST(const G4String& name,G4Material* medium,G4Material* solid) {
   
      G4STRead reader;
      return reader.Read(name,medium,solid);
   }
*/
   G4double GetConstant(const G4String& name) { return reader.getConstant(name); }
   G4double GetVariable(const G4String& name) { return reader.getVariable(name); }
   G4double GetQuantity(const G4String& name) { return reader.getQuantity(name); }
   G4ThreeVector GetPosition(const G4String& name) { return reader.getPosition(name); }
   G4ThreeVector GetRotation(const G4String& name) { return reader.getRotation(name); }
   G4ThreeVector GetScale(const G4String& name) { return reader.getScale(name); }
   G4GDMLMatrix GetMatrix(const G4String& name) { return reader.getMatrix(name); }
   G4LogicalVolume* GetVolume(const G4String& name) { return reader.getVolume(name); }
   G4VPhysicalVolume* GetWorldVolume(const G4String& setupName="Default") { return reader.GetWorldVolume(setupName); }
   G4GDMLAuxListType GetVolumeAuxiliaryInformation(const G4LogicalVolume* const logvol) { return reader.getVolumeAuxiliaryInformation(logvol); }
   void AddModule(const G4VPhysicalVolume* const physvol) { writer.AddModule(physvol); }
   void AddModule(const G4long depth) { writer.AddModule(depth); }
   void SetAddPointerToName(bool set) { writer.SetAddPointerToName(set); }
};

#endif
