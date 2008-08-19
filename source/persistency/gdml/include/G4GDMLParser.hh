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
// $Id: G4GDMLParser.hh,v 1.52 2008-08-19 15:03:17 gcosmo Exp $
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

#include "G4GDMLReadStructure.hh"
#include "G4GDMLWriteStructure.hh"
#include "G4STRead.hh"

#include "G4TransportationManager.hh"  // Used for the writing the whole
                                       // geometrical structure
#define G4GDML_DEFAULT_SCHEMALOCATION G4String("http://service-spi.web.cern.ch/service-spi/app/releases/GDML/GDML_2_10_0/src/GDMLSchema/gdml.xsd")

class G4GDMLParser
{
  public:  // with description

   G4GDMLParser() { xercesc::XMLPlatformUtils::Initialize(); }
   ~G4GDMLParser() { xercesc::XMLPlatformUtils::Terminate(); }

   inline void Read(const G4String& filename, G4bool Validate=true);

   inline void Write(const G4String& filename,
                     const G4VPhysicalVolume* const pvol = 0,
                           G4bool storeReferences = true,
                     const G4String& SchemaLocation = G4GDML_DEFAULT_SCHEMALOCATION);

   inline G4LogicalVolume* ParseST(const G4String& name,
                                   G4Material* medium, G4Material* solid);

   // Methods for Reader
   //
   inline G4double GetConstant(const G4String& name);
   inline G4double GetVariable(const G4String& name);
   inline G4double GetQuantity(const G4String& name);
   inline G4ThreeVector GetPosition(const G4String& name);
   inline G4ThreeVector GetRotation(const G4String& name);
   inline G4ThreeVector GetScale(const G4String& name);
   inline G4GDMLMatrix GetMatrix(const G4String& name);
   inline G4LogicalVolume* GetVolume(const G4String& name);
   inline G4VPhysicalVolume* GetWorldVolume(const G4String& setupName="Default");
   inline G4GDMLAuxListType GetVolumeAuxiliaryInformation(const G4LogicalVolume* const logvol);

   // Methods for Writer
   //
   inline void AddModule(const G4VPhysicalVolume* const physvol);
   inline void AddModule(const G4int depth);
   inline void SetAddPointerToName(G4bool set);

  private:

   G4GDMLReadStructure reader;
   G4GDMLWriteStructure writer;

};

#include "G4GDMLParser.icc"

#endif
