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
/// \file persistency/gdml/G03/src/G03ColorReader.cc
/// \brief Implementation of the G03ColorReader class
//
//
// --------------------------------------------------------------------

#include "G03ColorReader.hh"

#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G03ColorReader::G03ColorReader()
  : G4GDMLReadStructure()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G03ColorReader::~G03ColorReader()
{
  std::map<G4String, G4VisAttributes*>::iterator pos;
  for (pos=fAttribs.begin(); pos!=fAttribs.end(); pos++)
  {  delete pos->second;  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G03ColorReader::ExtensionRead(const xercesc::DOMElement* const extElement)
{
   G4cout << "G4GDML: Reading GDML extension..." << G4endl;

   for (xercesc::DOMNode* iter = extElement->getFirstChild();
        iter != 0; iter = iter->getNextSibling())
   {
      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE)  { continue; }

      const xercesc::DOMElement* const child
            = dynamic_cast<xercesc::DOMElement*>(iter);
      const G4String tag = Transcode(child->getTagName());

      if (tag=="color") { ColorRead(child); }
      else
      {
        G4String error_msg = "Unknown tag in structure: " + tag;
        G4Exception("G03ColorReader::ExtensionRead()",
                    "ReadError", FatalException, error_msg);
      }
   }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G03ColorReader::VolumeRead(const xercesc::DOMElement* const volumeElement)
{
   G4VSolid* solidPtr = 0;
   G4Material* materialPtr = 0;
   G4VisAttributes* attrPtr = 0;
   G4GDMLAuxListType auxList;
   
   XMLCh *name_attr = xercesc::XMLString::transcode("name");
   const G4String name = Transcode(volumeElement->getAttribute(name_attr));
   xercesc::XMLString::release(&name_attr);

   for (xercesc::DOMNode* iter = volumeElement->getFirstChild();
        iter != 0; iter = iter->getNextSibling())
   {
      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE)  { continue; }

      const xercesc::DOMElement* const child
            = dynamic_cast<xercesc::DOMElement*>(iter);
      const G4String tag = Transcode(child->getTagName());

      if (tag=="auxiliary")
        { auxList.push_back(AuxiliaryRead(child)); } else
      if (tag=="materialref")
        { materialPtr = GetMaterial(GenerateName(RefRead(child),true)); } else
      if (tag=="solidref")
        { solidPtr = GetSolid(GenerateName(RefRead(child))); } else
      if (tag == "colorref")
        { attrPtr = GetVisAttribute(GenerateName(RefRead(child))); }
   }

   pMotherLogical = new G4LogicalVolume(solidPtr,materialPtr,
                                        GenerateName(name),0,0,0);
   pMotherLogical->SetVisAttributes(attrPtr);

   if (!auxList.empty()) { auxMap[pMotherLogical] = auxList; }

   Volume_contentRead(volumeElement);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G03ColorReader::ColorRead(const xercesc::DOMElement* const colorElement)
{
   G4String name;
   G4VisAttributes* color = 0;
   G4double r=0., g=0., b=0., a=0.;

   const xercesc::DOMNamedNodeMap* const attributes
         = colorElement->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;
        attribute_index<attributeCount; attribute_index++)
   {
      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE)
        { continue; }

      const xercesc::DOMAttr* const attribute
            = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   
      const G4String attName = Transcode(attribute->getName());
      const G4String attValue = Transcode(attribute->getValue());

      if (attName=="name")
        { name = GenerateName(attValue); } else
      if (attName=="R")
        { r = eval.Evaluate(attValue); } else
      if (attName=="G")
        { g = eval.Evaluate(attValue); } else
      if (attName=="B")
        { b = eval.Evaluate(attValue); } else
      if (attName=="A")
        { a = eval.Evaluate(attValue); }
   }

   G4cout << "Color attribute (R,G,B,A) is: "
          << r << ", " << g << ", " << b << ", " << a << " !" << G4endl;
   color = new G4VisAttributes(G4Color(r,g,b,a));
   fAttribs.insert(std::make_pair(name,color));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VisAttributes* G03ColorReader::GetVisAttribute(const G4String& ref)
{
  G4VisAttributes* col = 0;
  std::map<G4String, G4VisAttributes*>::iterator pos = fAttribs.find(ref);

  if (pos != fAttribs.end())
  {
    col = pos->second;
  }
  else
  {
    G4String err_mess = "Attribute: " + ref + " NOT found !";
    G4Exception("G03ColorReader::GetVisAttribute()",
                "ReadError", FatalException, err_mess);
  }
  return col;
}
