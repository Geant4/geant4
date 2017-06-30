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
// $Id: G4GDMLParser.hh 104688 2017-06-12 12:16:28Z gcosmo $
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
#include "G4GDMLMessenger.hh"
#include "G4GDMLEvaluator.hh"

#include "G4TransportationManager.hh"
#include "G4Navigator.hh"


#define G4GDML_DEFAULT_SCHEMALOCATION G4String("http://service-spi.web.cern.ch/service-spi/app/releases/GDML/schema/gdml.xsd")

class G4GDMLParser
{
  public:  // with description

   G4GDMLParser();
   G4GDMLParser(G4GDMLReadStructure*);
   G4GDMLParser(G4GDMLReadStructure*, G4GDMLWriteStructure*);
  ~G4GDMLParser();
     //
     // Parser constructors & destructor

   inline void Read(const G4String& filename, G4bool Validate=true);
     //
     // Imports geometry with world-volume, specified by the GDML filename
     // in input. Validation against schema is activated by default.

   inline void ReadModule(const G4String& filename, G4bool Validate=true);
     //
     // Imports a single GDML module, specified by the GDML filename
     // in input. Validation against schema is activated by default.

   inline void Write(const G4String& filename,
                     const G4VPhysicalVolume* pvol = 0,
                           G4bool storeReferences = true,
                     const G4String& SchemaLocation = G4GDML_DEFAULT_SCHEMALOCATION);
     //
     // Exports on a GDML file, specified by 'filename' a geometry tree
     // starting from 'pvol' as top volume. Uniqueness of stored entities
     // is guaranteed by storing pointer-references by default.
     // Alternative path for the schema location can be specified; by default
     // the URL to the GDML web site is used.

   inline void Write(const G4String& filename,
                     const G4LogicalVolume* lvol,
                           G4bool storeReferences = true,
                     const G4String& SchemaLocation = G4GDML_DEFAULT_SCHEMALOCATION);
     //
     // Exports on a GDML file, specified by 'filename' a geometry tree
     // starting from 'pvol' as top volume. Uniqueness of stored entities
     // is guaranteed by storing pointer-references by default.
     // Alternative path for the schema location can be specified; by default
     // the URL to the GDML web site is used. Same as method above except
     // that the logical volume must be provided here.

   inline G4LogicalVolume* ParseST(const G4String& name,
                                         G4Material* medium,
                                         G4Material* solid);
     //
     // Imports a tessellated geometry stored as STEP-Tools files
     // 'name.geom' and 'name.tree'. It returns a pointer of a generated
     // mother volume with 'medium' material associated, including the
     // imported tessellated geometry with 'solid' material associated.

   // Methods for Reader
   //
   inline G4bool IsValid(const G4String& name) const;
   inline G4double GetConstant(const G4String& name) const;
   inline G4double GetVariable(const G4String& name) const;
   inline G4double GetQuantity(const G4String& name) const;
   inline G4ThreeVector GetPosition(const G4String& name) const;
   inline G4ThreeVector GetRotation(const G4String& name) const;
   inline G4ThreeVector GetScale(const G4String& name) const;
   inline G4GDMLMatrix GetMatrix(const G4String& name) const;
   inline G4LogicalVolume* GetVolume(const G4String& name) const;
   inline G4VPhysicalVolume* GetWorldVolume(const G4String& setupName="Default") const;
   inline G4GDMLAuxListType GetVolumeAuxiliaryInformation(G4LogicalVolume* lvol) const;
   inline const G4GDMLAuxMapType* GetAuxMap() const;
   inline const G4GDMLAuxListType* GetAuxList() const;
   inline void AddAuxiliary(G4GDMLAuxStructType myaux);
   inline void StripNamePointers() const;
   inline void SetStripFlag(G4bool);
   inline void SetOverlapCheck(G4bool);
   inline void SetRegionExport(G4bool);
   inline void SetEnergyCutsExport(G4bool);
   inline void SetSDExport(G4bool);

   inline G4int GetMaxExportLevel() const;  // Manage max number of levels
   inline void SetMaxExportLevel(G4int);    // to export

   inline void Clear();                  // Clears the evaluator

   // Methods for Writer
   //
   inline void AddModule(const G4VPhysicalVolume* const physvol);
   inline void AddModule(const G4int depth);
   inline void SetAddPointerToName(G4bool set);
   inline void AddVolumeAuxiliary(G4GDMLAuxStructType myaux, const G4LogicalVolume* const lvol);

  private:

   void ImportRegions();
   void ExportRegions(G4bool storeReferences = true);

  private:

   G4GDMLEvaluator eval;
   G4GDMLReadStructure* reader;
   G4GDMLWriteStructure* writer;
   G4GDMLAuxListType *rlist, *ullist;
   G4GDMLMessenger* messenger;
   G4bool urcode, uwcode, strip, rexp;

};

#include "G4GDMLParser.icc"

#endif
