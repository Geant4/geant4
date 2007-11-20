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
// $Id: G4GDMLSetup.cc,v 1.8 2007-11-20 09:37:11 gcosmo Exp $
// GEANT4 tag $ Name:$
//
// class G4GDMLSetup Implementation
//
// Original author: Zoltan Torzsok, November 2007
//
// --------------------------------------------------------------------

#include "G4GDMLSetup.hh"

G4VPhysicalVolume *G4GDMLSetup::Get(const G4String& ref) {

   if (setupMap.find(ref) == setupMap.end()) {
   
      G4cout << "GDML: Error! Referenced setup '" << ref << "' was not found!" << G4endl;
      return 0;
   }

   G4String worldref = setupMap[ref];

   G4LogicalVolume *volume = G4LogicalVolumeStore::GetInstance()->GetVolume(worldref,false);

   if (!volume) {
   
      G4cout << "G4GDML ERROR! volume '" << worldref << "' referenced in setup '" << ref << "' was not found!" << G4endl;
      return 0;
   }

   volume->SetVisAttributes(G4VisAttributes::Invisible);

   return new G4PVPlacement(0,G4ThreeVector(),volume,ref,0,0,0);
}

G4bool G4GDMLSetup::Read(const xercesc::DOMElement* const element,const G4String& newModule) {

   module = newModule;

   G4String name;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name") { name = attribute_value; }
   }

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag != "world") continue;

      XMLCh *ref_attr = xercesc::XMLString::transcode("ref");
      G4String ref = xercesc::XMLString::transcode(child->getAttribute(ref_attr));
      xercesc::XMLString::release(&ref_attr);

      setupMap[name] = module+ref;

      return true;
   }

   return false;
}
