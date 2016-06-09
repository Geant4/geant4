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
// $Id: G4GDMLBase.cc,v 1.3 2007/11/30 11:58:46 ztorzsok Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
//
// class G4GDMLBase
//
// Class description:
//
// History:
// - Created.                                  Zoltan Torzsok, November 2007
// -------------------------------------------------------------------------

#include "G4GDMLBase.hh"

G4GDMLBase::G4GDMLBase() {

   try {

      xercesc::XMLPlatformUtils::Initialize();
   }
   catch(xercesc::XMLException& e) {

      char* message = xercesc::XMLString::transcode(e.getMessage());
      G4cerr << "XML toolkit initialization error: " << message << G4endl;
      xercesc::XMLString::release(&message);
   }

   parser = new xercesc::XercesDOMParser;

   parser->setValidationScheme(xercesc::XercesDOMParser::Val_Always);
   parser->setDoNamespaces(true);
   parser->setDoSchema(true);
   parser->setValidationSchemaFullChecking(true);
}

G4GDMLBase::~G4GDMLBase() {

   if (parser) delete parser;

   try {

      xercesc::XMLPlatformUtils::Terminate();
   }
   catch(xercesc::XMLException& e) {
    
      char* message = xercesc::XMLString::transcode(e.getMessage());
      G4cerr << "XML toolkit termination error: " << message << G4endl;
      xercesc::XMLString::release(&message);
   }
}

G4String G4GDMLBase::GenerateName(const G4String& in) {

   std::string out(prename);
   
   std::string::size_type open = in.find("[",0);

   out.append(in,0,open);
   
   while (open != std::string::npos) {
   
      std::string::size_type close = in.find("]",open);

      if (close == std::string::npos) G4Exception("Bracket mismatch in loop!");
   
      std::string expr = in.substr(open+1,close-open-1);

      std::stringstream stream;
      
      stream << "[" << eval.EvaluateInteger(expr) << "]";
   
      out.append(stream.str());

      open = in.find("[",close);
   }

   return out;
}

void G4GDMLBase::Parse(const G4String& fileName) {

   prename = fileName + "_";

   try {

      parser->parse(fileName.c_str());
   }
   catch (const xercesc::XMLException &e) {
   
      char* message = xercesc::XMLString::transcode(e.getMessage());
      G4cout << "XML: " << message << G4endl;
      xercesc::XMLString::release(&message);
   }
   catch (const xercesc::DOMException &e) {
   
      char* message = xercesc::XMLString::transcode(e.getMessage());
      G4cout << "DOM: " << message << G4endl;
      xercesc::XMLString::release(&message);
   }

   xercesc::DOMDocument* doc = parser->getDocument();

   if (!doc) G4Exception("GDML: Unable to open document: "+fileName);

   xercesc::DOMElement* element = doc->getDocumentElement();

   if (!element) G4Exception("GDML: Empty document!");

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="define"   ) defineRead(child); else
      if (tag=="materials") materialsRead(child); else
      if (tag=="solids"   ) solidsRead(child); else
      if (tag=="setup"    ) setupRead(child); else
      if (tag=="structure") structureRead(child);
   }
}

G4PVPlacement* G4GDMLBase::getTopVolume(const G4String& setupName) {

   G4LogicalVolume* volume = getVolume(getSetup(setupName));

   volume->SetVisAttributes(G4VisAttributes::Invisible);

   return new G4PVPlacement(0,G4ThreeVector(),volume,"",0,0,0);
}
