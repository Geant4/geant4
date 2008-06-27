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
// Class description:
//
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

class G4GDMLWrite {
private:
   typedef std::vector<const G4VPhysicalVolume*> VolumeListType;
   typedef std::vector<G4long> DepthListType;
   static bool addPointerToName;
   xercesc::DOMDocument* doc;
   XMLCh tempStr[100];

   bool FileExists(const G4String&) const;
   VolumeListType& volumeList();
   DepthListType& depthList();
protected:
   G4String SchemaLocation;

   G4String GenerateName(const G4String&,const void* const);
   xercesc::DOMAttr* newAttribute(const G4String&,const G4String&);
   xercesc::DOMAttr* newAttribute(const G4String&,const G4double&);
   xercesc::DOMElement* newElement(const G4String&);
   virtual void defineWrite(xercesc::DOMElement*)=0;
   virtual void materialsWrite(xercesc::DOMElement*)=0;
   virtual void solidsWrite(xercesc::DOMElement*)=0;
   virtual void structureWrite(xercesc::DOMElement*)=0;
   virtual G4Transform3D TraverseVolumeTree(const G4LogicalVolume* const,G4long)=0;
   virtual void setupWrite(xercesc::DOMElement*,const G4LogicalVolume* const)=0;
   bool Modularize(const G4VPhysicalVolume* const,const G4long);
public:
   G4Transform3D Write(const G4String&,const G4LogicalVolume* const,const G4String&,const G4long);
   void AddModule(const G4VPhysicalVolume* const);
   void AddModule(const G4long);
   static void SetAddPointerToName(bool);
};

#endif
