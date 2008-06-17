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
// History:
// - Created.                                  Zoltan Torzsok, November 2007
// -------------------------------------------------------------------------

#include "G4GDMLRead.hh"

G4String G4GDMLRead::GenerateName(const G4String& nameIn) {

   G4String nameOut(nameIn);

   if (InLoop>0) nameOut = eval.SolveBrackets(nameOut);

   return G4String(ModuleName + nameOut);
}

void G4GDMLRead::GeneratePhysvolName(const G4String& nameIn,G4VPhysicalVolume* physvol) {

   G4String nameOut = nameIn;

   if (nameOut.empty()) {
      
      std::stringstream stream;
      stream << physvol->GetLogicalVolume()->GetName() << "_in_" << physvol->GetMotherLogical()->GetName() << physvol;
      nameOut = stream.str();
   } else {
   
     nameOut = ModuleName + nameOut;
   }
   
   physvol->SetName(eval.SolveBrackets(nameOut));
}

void G4GDMLRead::loopRead(const xercesc::DOMElement* const element,void(G4GDMLRead::*func)(const xercesc::DOMElement* const)) {

   G4String var;
   G4String from;
   G4String to;
   G4String step;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="for") var = attribute_value; else
      if (attribute_name=="from") from = attribute_value; else
      if (attribute_name=="to") to = attribute_value; else
      if (attribute_name=="step") step = attribute_value;
   }

   if (var.empty()) G4Exception("G4GDML: ERROR! No variable is determined for loop!");

   if (!eval.isVariable(var)) G4Exception("G4GDML: ERROR! Variable is not defined in loop!");

   G4int _var = eval.EvaluateInteger(var);
   G4int _from = eval.EvaluateInteger(from);
   G4int _to = eval.EvaluateInteger(to);
   G4int _step = eval.EvaluateInteger(step);
   
   if (!from.empty()) _var = _from;

   if (_from == _to) G4Exception("G4GDML: ERROR! Empty loop!");
   if (_from < _to && _step <= 0) G4Exception("G4GDML: ERROR! Infinite loop!");
   if (_from > _to && _step >= 0) G4Exception("G4GDML: ERROR! Infinite loop!");
   
   InLoop++;
   
   while (_var <= _to) {
   
      eval.setVariable(var,_var);
      (this->*func)(element);

      _var += _step;
   }

   InLoop--;
}

void G4GDMLRead::Read(const G4String& fileName,bool SetValidate,bool IsModule) {

   if (IsModule) G4cout << "G4GDML: Reading module '" << fileName << "'..." << G4endl;
   else G4cout << "G4GDML: Reading '" << fileName << "'..." << G4endl;

   InLoop = 0;
   ModuleName.clear();
   Validate = SetValidate;

   if (IsModule) { 
   
      ModuleName = fileName;
      ModuleName.remove(ModuleName.length()-5,5); // remove ".gdml"
      ModuleName += "_";
   }

   xercesc::ErrorHandler* handler = new MyErrorHandler(!Validate);
   xercesc::XercesDOMParser* parser = new xercesc::XercesDOMParser;

   parser->setValidationScheme(xercesc::XercesDOMParser::Val_Always);
   parser->setValidationSchemaFullChecking(true);
   parser->setCreateEntityReferenceNodes(false);   // Entities will be automatically resolved by Xerces
   parser->setDoNamespaces(true);
   parser->setDoSchema(true);
   parser->setErrorHandler(handler);

   try {

      parser->parse(fileName.c_str());
   }
   catch (const xercesc::XMLException &e) {
   
      char* message = xercesc::XMLString::transcode(e.getMessage());
      G4cout << "G4GDML: " << message << G4endl;
      xercesc::XMLString::release(&message);
   }
   catch (const xercesc::DOMException &e) {
   
      char* message = xercesc::XMLString::transcode(e.getMessage());
      G4cout << "G4GDML: " << message << G4endl;
      xercesc::XMLString::release(&message);
   }

   xercesc::DOMDocument* doc = parser->getDocument();

   if (!doc) G4Exception("G4GDML: ERROR! Unable to open document: "+fileName);

   xercesc::DOMElement* element = doc->getDocumentElement();

   if (!element) G4Exception("G4GDML: ERROR! Empty document!");

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="define") defineRead(child); else
      if (tag=="materials") materialsRead(child); else
      if (tag=="solids") solidsRead(child); else
      if (tag=="setup") setupRead(child); else
      if (tag=="structure") structureRead(child); else
      G4Exception("G4GDML: ERROR! Unknown tag in gdml: "+tag);
   }

   if (parser) delete parser;
   if (handler) delete handler;

   if (IsModule) G4cout << "G4GDML: Reading module '" << fileName << "' done!" << G4endl;
   else G4cout << "G4GDML: Reading '" << fileName << "' done!" << G4endl;
}
