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
// $Id: G4GDMLRead.cc 110108 2018-05-15 11:46:54Z gcosmo $
//
// class G4GDMLRead Implementation
//
// History:
// - Created.                                  Zoltan Torzsok, November 2007
// -------------------------------------------------------------------------

#include "globals.hh"

#include "G4GDMLRead.hh"

#include "G4UnitsTable.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"

G4GDMLRead::G4GDMLRead()
  : validate(true), check(false), dostrip(true), inLoop(0), loopCount(0)
{
   G4UnitDefinition::BuildUnitsTable();
}

G4GDMLRead::~G4GDMLRead()
{
}

G4String G4GDMLRead::Transcode(const XMLCh* const toTranscode)
{
   char* char_str = xercesc::XMLString::transcode(toTranscode);
   G4String my_str(char_str);
   xercesc::XMLString::release(&char_str);
   return my_str;
}

void G4GDMLRead::OverlapCheck(G4bool flag)
{
   check = flag;
}

G4String G4GDMLRead::GenerateName(const G4String& nameIn, G4bool strip)
{
   G4String nameOut(nameIn);

   if (inLoop>0)
   {
     nameOut = eval.SolveBrackets(nameOut);
   }
   if (strip) { StripName(nameOut); }
   
   return nameOut;
}

void G4GDMLRead::GeneratePhysvolName(const G4String& nameIn,
                                     G4VPhysicalVolume* physvol)
{
   G4String nameOut(nameIn);

   if (nameIn.empty())
   {
     std::stringstream stream;
     stream << physvol->GetLogicalVolume()->GetName() << "_PV";
     nameOut = stream.str();
   }
   nameOut = eval.SolveBrackets(nameOut);

   physvol->SetName(nameOut);
}

G4String G4GDMLRead::Strip(const G4String& name) const
{
  G4String sname(name);
  return sname.remove(sname.find("0x"));
}

void G4GDMLRead::StripName(G4String& name) const
{
  name.remove(name.find("0x"));
}

void G4GDMLRead::StripNames() const
{
  // Strips off names of volumes, solids elements and materials from possible
  // reference pointers or IDs attached to their original identifiers.

  G4PhysicalVolumeStore* pvols = G4PhysicalVolumeStore::GetInstance();
  G4LogicalVolumeStore* lvols = G4LogicalVolumeStore::GetInstance();
  G4SolidStore* solids = G4SolidStore::GetInstance();
  const G4ElementTable* elements = G4Element::GetElementTable();
  const G4MaterialTable* materials = G4Material::GetMaterialTable();

  G4cout << "Stripping off GDML names of materials, solids and volumes ..."
         << G4endl;

  G4String sname;
  size_t i;

  // Solids...
  //
  for (i=0; i<solids->size(); ++i)
  {
    G4VSolid* psol = (*solids)[i];
    sname = psol->GetName();
    StripName(sname);
    psol->SetName(sname);
  }

  // Logical volumes...
  //
  for (i=0; i<lvols->size(); ++i)
  {
    G4LogicalVolume* lvol = (*lvols)[i];
    sname = lvol->GetName();
    StripName(sname);
    lvol->SetName(sname);
  }

  // Physical volumes...
  //
  for (i=0; i<pvols->size(); ++i)
  {
    G4VPhysicalVolume* pvol = (*pvols)[i];
    sname = pvol->GetName();
    StripName(sname);
    pvol->SetName(sname);
  }

  // Materials...
  //
  for (i=0; i<materials->size(); ++i)
  {
    G4Material* pmat = (*materials)[i];
    sname = pmat->GetName();
    StripName(sname);
    pmat->SetName(sname);
  }

  // Elements...
  //
  for (i=0; i<elements->size(); ++i)
  {
    G4Element* pelm = (*elements)[i];
    sname = pelm->GetName();
    StripName(sname);
    pelm->SetName(sname);
  }
}

void G4GDMLRead::LoopRead(const xercesc::DOMElement* const element,
     void(G4GDMLRead::*func)(const xercesc::DOMElement* const))
{
   G4String var;
   G4String from;
   G4String to;
   G4String step;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;
        attribute_index<attributeCount;attribute_index++)
   {
      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE)
      { continue; }

      const xercesc::DOMAttr* const attribute
            = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   
      if (!attribute)
      {
        G4Exception("G4GDMLRead::LoopRead()", "InvalidRead",
                    FatalException, "No attribute found!");
        return;
      }
      const G4String attribute_name = Transcode(attribute->getName());
      const G4String attribute_value = Transcode(attribute->getValue());

      if (attribute_name=="for")  { var = attribute_value; }  else
      if (attribute_name=="from") { from = attribute_value; } else
      if (attribute_name=="to")   { to = attribute_value; }   else
      if (attribute_name=="step") { step = attribute_value; }
   }

   if (var.empty())
   {
     G4Exception("G4GDMLRead::loopRead()", "InvalidRead",
                 FatalException, "No variable is determined for loop!");
   }

   if (!eval.IsVariable(var))
   {
     G4Exception("G4GDMLRead::loopRead()", "InvalidRead",
                 FatalException, "Variable is not defined in loop!");
   }

   G4int _var = eval.EvaluateInteger(var);
   G4int _from = eval.EvaluateInteger(from);
   G4int _to = eval.EvaluateInteger(to);
   G4int _step = eval.EvaluateInteger(step);
   
   if (!from.empty()) { _var = _from; }

   if (_from == _to)
   {
     G4Exception("G4GDMLRead::loopRead()", "InvalidRead",
                 FatalException, "Empty loop!");
   }
   if ((_from < _to) && (_step <= 0))
   {
     G4Exception("G4GDMLRead::loopRead()", "InvalidRead",
                 FatalException, "Infinite loop!");
   }
   if ((_from > _to) && (_step >= 0))
   {
     G4Exception("G4GDMLRead::loopRead()", "InvalidRead",
                 FatalException, "Infinite loop!");
   }

   inLoop++;

   while (_var <= _to)
   {
      eval.SetVariable(var,_var);
      (this->*func)(element);
      _var += _step;
      loopCount++;
   }

   inLoop--;
   if (!inLoop) { loopCount = 0; }
}

G4GDMLAuxStructType G4GDMLRead::
AuxiliaryRead(const xercesc::DOMElement* const auxiliaryElement)
{
   G4GDMLAuxStructType auxstruct = {"","","",0};
   G4GDMLAuxListType* auxList=0;
  
   const xercesc::DOMNamedNodeMap* const attributes
         = auxiliaryElement->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;
        attribute_index<attributeCount; attribute_index++)
   {
      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE)
      { continue; }

      const xercesc::DOMAttr* const attribute
            = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   
      if (!attribute)
      {
        G4Exception("G4GDMLRead::AuxiliaryRead()",
                    "InvalidRead", FatalException, "No attribute found!");
        return auxstruct;
      }
      const G4String attName = Transcode(attribute->getName());
      const G4String attValue = Transcode(attribute->getValue());

      if (attName=="auxtype") { auxstruct.type = attValue; } else
      if (attName=="auxvalue") { auxstruct.value = attValue; } else
      if (attName=="auxunit") { auxstruct.unit = attValue; }
   }

   for (xercesc::DOMNode* iter = auxiliaryElement->getFirstChild();
        iter != 0; iter = iter->getNextSibling())
   {
      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE)  { continue; }
       
      const xercesc::DOMElement* const child
        = dynamic_cast<xercesc::DOMElement*>(iter);
      if (!child)
      {
         G4Exception("G4GDMLRead::AuxiliaryRead()",
                     "InvalidRead", FatalException, "No child found!");
         break;
      }
      const G4String tag = Transcode(child->getTagName());
       
      if (tag=="auxiliary")
      {
         if(!auxList) { auxList = new G4GDMLAuxListType; }
         auxList->push_back(AuxiliaryRead(child));
      }
   }
   
   if (auxList) { auxstruct.auxList = auxList; }

   return auxstruct;
}

void G4GDMLRead::UserinfoRead(const xercesc::DOMElement* const userinfoElement)
{
#ifdef G4VERBOSE
   G4cout << "G4GDML: Reading userinfo..." << G4endl;
#endif
   for (xercesc::DOMNode* iter = userinfoElement->getFirstChild();
        iter != 0; iter = iter->getNextSibling())
   {
      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE)  { continue; }

      const xercesc::DOMElement* const child
            = dynamic_cast<xercesc::DOMElement*>(iter);
      if (!child)
      {
        G4Exception("G4GDMLRead::UserinfoRead()",
                    "InvalidRead", FatalException, "No child found!");
        return;
      }
      const G4String tag = Transcode(child->getTagName());

      if (tag=="auxiliary")
      {
        auxGlobalList.push_back(AuxiliaryRead(child));
      }
      else
      {
        G4String error_msg = "Unknown tag in structure: " + tag;
        G4Exception("G4GDMLRead::UserinfoRead()",
                    "ReadError", FatalException, error_msg);
      }
   }
}

void G4GDMLRead::ExtensionRead(const xercesc::DOMElement* const)
{
   G4String error_msg = "No handle to user-code for parsing extensions!";
   G4Exception("G4GDMLRead::ExtensionRead()",
               "NotImplemented", JustWarning, error_msg);
}

void G4GDMLRead::Read(const G4String& fileName,
                            G4bool validation,
                            G4bool isModule,
                            G4bool strip)
{
   dostrip = strip;
#ifdef G4VERBOSE
   if (isModule)
   {
      G4cout << "G4GDML: Reading module '" << fileName << "'..." << G4endl;
   }
   else
   {
      G4cout << "G4GDML: Reading '" << fileName << "'..." << G4endl;
   }
#endif
   inLoop = 0;
   validate = validation;

   xercesc::ErrorHandler* handler = new G4GDMLErrorHandler(!validate);
   xercesc::XercesDOMParser* parser = new xercesc::XercesDOMParser;

   if (validate)
   {
     parser->setValidationScheme(xercesc::XercesDOMParser::Val_Always);
   }
   parser->setValidationSchemaFullChecking(validate);
   parser->setCreateEntityReferenceNodes(false); 
     // Entities will be automatically resolved by Xerces

   parser->setDoNamespaces(true);
   parser->setDoSchema(validate);
   parser->setErrorHandler(handler);

   try { parser->parse(fileName.c_str()); }
   catch (const xercesc::XMLException &e)
     { G4cout << "G4GDML: " << Transcode(e.getMessage()) << G4endl; }
   catch (const xercesc::DOMException &e)
     { G4cout << "G4GDML: " << Transcode(e.getMessage()) << G4endl; }

   xercesc::DOMDocument* doc = parser->getDocument();

   if (!doc)
   {
     G4String error_msg = "Unable to open document: " + fileName;
     G4Exception("G4GDMLRead::Read()", "InvalidRead",
                 FatalException, error_msg);
     return;
   }
   xercesc::DOMElement* element = doc->getDocumentElement();

   if (!element)
   {
     std::ostringstream message;
     message << "ERROR - Empty document or unable to validate schema!" << G4endl
             << "        Check Internet connection is ON in case of schema"
             << G4endl
             << "        validation enabled and location defined as URL in"
             << G4endl
             << "        the GDML file - " << fileName << " - being imported!"
             << G4endl
             << "        Otherwise, verify GDML schema server is reachable!";
     G4Exception("G4GDMLRead::Read()", "InvalidRead", FatalException, message);
     return;
   }

   for (xercesc::DOMNode* iter = element->getFirstChild();
        iter != 0; iter = iter->getNextSibling())
   {
      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE)  { continue; }

      const xercesc::DOMElement* const child
            = dynamic_cast<xercesc::DOMElement*>(iter);
      if (!child)
      {
        G4Exception("G4GDMLRead::Read()", "InvalidRead",
                    FatalException, "No child found!");
        return;
      }
      const G4String tag = Transcode(child->getTagName());

      if (tag=="define")    { DefineRead(child);    } else
      if (tag=="materials") { MaterialsRead(child); } else
      if (tag=="solids")    { SolidsRead(child);    } else
      if (tag=="setup")     { SetupRead(child);     } else
      if (tag=="structure") { StructureRead(child); } else
      if (tag=="userinfo")  { UserinfoRead(child);  } else
      if (tag=="extension") { ExtensionRead(child); }
      else
      {
        G4String error_msg = "Unknown tag in gdml: " + tag;
        G4Exception("G4GDMLRead::Read()", "InvalidRead",
                    FatalException, error_msg);
      }
   }

   delete parser;
   delete handler;

   if (isModule)
   {
#ifdef G4VERBOSE
      G4cout << "G4GDML: Reading module '" << fileName << "' done!" << G4endl;
#endif
   }
   else
   {
      G4cout << "G4GDML: Reading '" << fileName << "' done!" << G4endl;
      if (strip)  { StripNames(); }
   }
}

const G4GDMLAuxListType* G4GDMLRead::GetAuxList() const
{
   return &auxGlobalList;
}
