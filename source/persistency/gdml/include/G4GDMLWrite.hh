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
// $Id: G4GDMLWrite.hh,v 1.37 2008/08/19 15:03:17 gcosmo Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
//
// class G4GDMLWrite
//
// Class description:
//
// GDML writer.

// History:
// - Created.                                  Zoltan Torzsok, November 2007
// -------------------------------------------------------------------------

#ifndef _G4GDMLWRITE_INCLUDED_
#define _G4GDMLWRITE_INCLUDED_

#include <sys/stat.h>
#include <iostream>

#include <xercesc/dom/DOM.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/framework/LocalFileFormatTarget.hpp>

#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4PVDivision.hh"

class G4GDMLWrite
{
   typedef std::map<const G4LogicalVolume*,G4Transform3D> VolumeMapType;
   typedef std::map<const G4VPhysicalVolume*,G4String> PhysVolumeMapType;
   typedef std::map<G4int,G4int> DepthMapType;

 public:  // without description

   G4Transform3D Write(const G4String& filename,
                       const G4LogicalVolume* const topLog,
                       const G4String& schemaPath,
                       const G4int depth, G4bool storeReferences=true);
   void AddModule(const G4VPhysicalVolume* const topVol);
   void AddModule(const G4int depth);
   static void SetAddPointerToName(G4bool);

 protected:

   G4String SchemaLocation;

   VolumeMapType& VolumeMap();

   G4String GenerateName(const G4String&,const void* const);
   xercesc::DOMAttr* NewAttribute(const G4String&, const G4String&);
   xercesc::DOMAttr* NewAttribute(const G4String&, const G4double&);
   xercesc::DOMElement* NewElement(const G4String&);
   virtual void DefineWrite(xercesc::DOMElement*)=0;
   virtual void MaterialsWrite(xercesc::DOMElement*)=0;
   virtual void SolidsWrite(xercesc::DOMElement*)=0;
   virtual void StructureWrite(xercesc::DOMElement*)=0;
   virtual G4Transform3D TraverseVolumeTree(const G4LogicalVolume* const,
                                            const G4int)=0;
   virtual void SetupWrite(xercesc::DOMElement*,
                           const G4LogicalVolume* const)=0;
   G4String Modularize(const G4VPhysicalVolume* const topvol,
                       const G4int depth);

 private:

   G4bool FileExists(const G4String&) const;
   PhysVolumeMapType& PvolumeMap();
   DepthMapType& DepthMap();

 private:

   static G4bool addPointerToName;
   xercesc::DOMDocument* doc;
   XMLCh tempStr[100];

};

#endif
