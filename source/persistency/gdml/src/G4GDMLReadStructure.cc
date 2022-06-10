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
// G4GDMLReadStructure implementation
//
// Author: Zoltan Torzsok, November 2007
// --------------------------------------------------------------------

#include "G4GDMLReadStructure.hh"

#include "G4UnitsTable.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4AssemblyVolume.hh"
#include "G4ReflectionFactory.hh"
#include "G4PVDivisionFactory.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4VisAttributes.hh"

// --------------------------------------------------------------------
G4GDMLReadStructure::G4GDMLReadStructure()
  : G4GDMLReadParamvol()
{
}

// --------------------------------------------------------------------
G4GDMLReadStructure::~G4GDMLReadStructure()
{
}

// --------------------------------------------------------------------
void G4GDMLReadStructure::BorderSurfaceRead(
  const xercesc::DOMElement* const bordersurfaceElement)
{
  G4String name;
  G4VPhysicalVolume* pv1  = nullptr;
  G4VPhysicalVolume* pv2  = nullptr;
  G4SurfaceProperty* prop = nullptr;
  G4int index             = 0;

  const xercesc::DOMNamedNodeMap* const attributes =
    bordersurfaceElement->getAttributes();
  XMLSize_t attributeCount = attributes->getLength();

  for(XMLSize_t attribute_index = 0; attribute_index < attributeCount;
      ++attribute_index)
  {
    xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

    if(attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE)
    {
      continue;
    }

    const xercesc::DOMAttr* const attribute =
      dynamic_cast<xercesc::DOMAttr*>(attribute_node);
    if(attribute == nullptr)
    {
      G4Exception("G4GDMLReadStructure::BorderSurfaceRead()", "InvalidRead",
                  FatalException, "No attribute found!");
      return;
    }
    const G4String attName  = Transcode(attribute->getName());
    const G4String attValue = Transcode(attribute->getValue());

    if(attName == "name")
    {
      name = GenerateName(attValue);
    }
    else if(attName == "surfaceproperty")
    {
      prop = GetSurfaceProperty(GenerateName(attValue));
    }
  }

  for(xercesc::DOMNode* iter = bordersurfaceElement->getFirstChild();
                        iter != nullptr; iter = iter->getNextSibling())
  {
    if(iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE)
    {
      continue;
    }

    const xercesc::DOMElement* const child =
      dynamic_cast<xercesc::DOMElement*>(iter);
    if(child == nullptr)
    {
      G4Exception("G4GDMLReadStructure::BorderSurfaceRead()", "InvalidRead",
                  FatalException, "No child found!");
      return;
    }
    const G4String tag = Transcode(child->getTagName());

    if(tag != "physvolref")
    {
      continue;
    }

    if(index == 0)
    {
      pv1 = GetPhysvol(GenerateName(RefRead(child)));
      ++index;
    }
    else if(index == 1)
    {
      pv2 = GetPhysvol(GenerateName(RefRead(child)));
      ++index;
    }
    else
      break;
  }

  new G4LogicalBorderSurface(Strip(name), pv1, pv2, prop);
}

// --------------------------------------------------------------------
void G4GDMLReadStructure::DivisionvolRead(
  const xercesc::DOMElement* const divisionvolElement)
{
  G4String name;
  G4double unit           = 1.0;
  G4double width          = 0.0;
  G4double offset         = 0.0;
  G4int number            = 0;
  EAxis axis              = kUndefined;
  G4LogicalVolume* logvol = nullptr;

  const xercesc::DOMNamedNodeMap* const attributes =
    divisionvolElement->getAttributes();
  XMLSize_t attributeCount = attributes->getLength();
  G4String unitname;

  for(XMLSize_t attribute_index = 0; attribute_index < attributeCount;
      ++attribute_index)
  {
    xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

    if(attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE)
    {
      continue;
    }

    const xercesc::DOMAttr* const attribute =
      dynamic_cast<xercesc::DOMAttr*>(attribute_node);
    if(attribute == nullptr)
    {
      G4Exception("G4GDMLReadStructure::DivisionvolRead()", "InvalidRead",
                  FatalException, "No attribute found!");
      return;
    }
    const G4String attName  = Transcode(attribute->getName());
    const G4String attValue = Transcode(attribute->getValue());

    if(attName == "name")
    {
      name = attValue;
    }
    else if(attName == "unit")
    {
      unit     = G4UnitDefinition::GetValueOf(attValue);
      unitname = G4UnitDefinition::GetCategory(attValue);
    }
    else if(attName == "width")
    {
      width = eval.Evaluate(attValue);
    }
    else if(attName == "offset")
    {
      offset = eval.Evaluate(attValue);
    }
    else if(attName == "number")
    {
      number = eval.EvaluateInteger(attValue);
    }
    else if(attName == "axis")
    {
      if(attValue == "kXAxis")
      {
        axis = kXAxis;
      }
      else if(attValue == "kYAxis")
      {
        axis = kYAxis;
      }
      else if(attValue == "kZAxis")
      {
        axis = kZAxis;
      }
      else if(attValue == "kRho")
      {
        axis = kRho;
      }
      else if(attValue == "kPhi")
      {
        axis = kPhi;
      }
    }
  }

  if(((axis == kXAxis || axis == kYAxis || axis == kZAxis) &&
      unitname != "Length") ||
     ((axis == kRho || axis == kPhi) && unitname != "Angle"))
  {
    G4Exception("G4GDMLReadStructure::DivisionvolRead()", "InvalidRead",
                FatalException, "Invalid unit!");
  }

  width *= unit;
  offset *= unit;

  for(xercesc::DOMNode* iter = divisionvolElement->getFirstChild();
                        iter != nullptr; iter = iter->getNextSibling())
  {
    if(iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE)
    {
      continue;
    }

    const xercesc::DOMElement* const child =
      dynamic_cast<xercesc::DOMElement*>(iter);
    if(child == nullptr)
    {
      G4Exception("G4GDMLReadStructure::DivisionvolRead()", "InvalidRead",
                  FatalException, "No child found!");
      return;
    }
    const G4String tag = Transcode(child->getTagName());

    if(tag == "volumeref")
    {
      logvol = GetVolume(GenerateName(RefRead(child)));
    }
  }

  if(logvol == nullptr)
  {
    return;
  }

  G4PVDivisionFactory::GetInstance();
  G4PhysicalVolumesPair pair;

  G4String pv_name = logvol->GetName() + "_div";
  if((number != 0) && (width == 0.0))
  {
    pair = G4ReflectionFactory::Instance()->Divide(
      pv_name, logvol, pMotherLogical, axis, number, offset);
  }
  else if((number == 0) && (width != 0.0))
  {
    pair = G4ReflectionFactory::Instance()->Divide(
      pv_name, logvol, pMotherLogical, axis, width, offset);
  }
  else
  {
    pair = G4ReflectionFactory::Instance()->Divide(
      pv_name, logvol, pMotherLogical, axis, number, width, offset);
  }

  if(pair.first != nullptr)
  {
    GeneratePhysvolName(name, pair.first);
  }
  if(pair.second != nullptr)
  {
    GeneratePhysvolName(name, pair.second);
  }
}

// --------------------------------------------------------------------
G4LogicalVolume* G4GDMLReadStructure::FileRead(
  const xercesc::DOMElement* const fileElement)
{
  G4String name;
  G4String volname;

  const xercesc::DOMNamedNodeMap* const attributes =
    fileElement->getAttributes();
  XMLSize_t attributeCount = attributes->getLength();

  for(XMLSize_t attribute_index = 0; attribute_index < attributeCount;
      ++attribute_index)
  {
    xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

    if(attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE)
    {
      continue;
    }

    const xercesc::DOMAttr* const attribute =
      dynamic_cast<xercesc::DOMAttr*>(attribute_node);
    if(attribute == nullptr)
    {
      G4Exception("G4GDMLReadStructure::FileRead()", "InvalidRead",
                  FatalException, "No attribute found!");
      return nullptr;
    }
    const G4String attName  = Transcode(attribute->getName());
    const G4String attValue = Transcode(attribute->getValue());

    if(attName == "name")
    {
      name = attValue;
    }
    else if(attName == "volname")
    {
      volname = attValue;
    }
  }

  const G4bool isModule = true;
  G4GDMLReadStructure structure;
  structure.Read(name, validate, isModule);

  // Register existing auxiliar information defined in child module
  //
  const G4GDMLAuxMapType* aux = structure.GetAuxMap();
  if(!aux->empty())
  {
    for(auto pos = aux->cbegin(); pos != aux->cend(); ++pos)
    {
      auxMap.insert(std::make_pair(pos->first, pos->second));
    }
  }

  // Return volume structure from child module
  //
  if(volname.empty())
  {
    return structure.GetVolume(structure.GetSetup("Default"));
  }
  else
  {
    return structure.GetVolume(structure.GenerateName(volname));
  }
}

// --------------------------------------------------------------------
void G4GDMLReadStructure::PhysvolRead(
  const xercesc::DOMElement* const physvolElement, G4AssemblyVolume* pAssembly)
{
  G4String name;
  G4LogicalVolume* logvol    = nullptr;
  G4AssemblyVolume* assembly = nullptr;
  G4ThreeVector position(0.0, 0.0, 0.0);
  G4ThreeVector rotation(0.0, 0.0, 0.0);
  G4ThreeVector scale(1.0, 1.0, 1.0);
  G4int copynumber = 0;

  const xercesc::DOMNamedNodeMap* const attributes =
    physvolElement->getAttributes();
  XMLSize_t attributeCount = attributes->getLength();

  for(XMLSize_t attribute_index = 0; attribute_index < attributeCount;
      ++attribute_index)
  {
    xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

    if(attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE)
    {
      continue;
    }

    const xercesc::DOMAttr* const attribute =
      dynamic_cast<xercesc::DOMAttr*>(attribute_node);
    if(attribute == nullptr)
    {
      G4Exception("G4GDMLReadStructure::PhysvolRead()", "InvalidRead",
                  FatalException, "No attribute found!");
      return;
    }
    const G4String attName  = Transcode(attribute->getName());
    const G4String attValue = Transcode(attribute->getValue());

    if(attName == "name")
    {
      name = attValue;
    }
    if(attName == "copynumber")
    {
      copynumber = eval.EvaluateInteger(attValue);
    }
  }

  for(xercesc::DOMNode* iter = physvolElement->getFirstChild(); iter != nullptr;
      iter                   = iter->getNextSibling())
  {
    if(iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE)
    {
      continue;
    }

    const xercesc::DOMElement* const child =
      dynamic_cast<xercesc::DOMElement*>(iter);
    if(child == nullptr)
    {
      G4Exception("G4GDMLReadStructure::PhysvolRead()", "InvalidRead",
                  FatalException, "No child found!");
      return;
    }
    const G4String tag = Transcode(child->getTagName());

    if(tag == "volumeref")
    {
      const G4String& child_name = GenerateName(RefRead(child));
      assembly                   = GetAssembly(child_name);
      if(assembly == nullptr)
      {
        logvol = GetVolume(child_name);
      }
    }
    else if(tag == "file")
    {
      logvol = FileRead(child);
    }
    else if(tag == "position")
    {
      VectorRead(child, position);
    }
    else if(tag == "rotation")
    {
      VectorRead(child, rotation);
    }
    else if(tag == "scale")
    {
      VectorRead(child, scale);
    }
    else if(tag == "positionref")
    {
      position = GetPosition(GenerateName(RefRead(child)));
    }
    else if(tag == "rotationref")
    {
      rotation = GetRotation(GenerateName(RefRead(child)));
    }
    else if(tag == "scaleref")
    {
      scale = GetScale(GenerateName(RefRead(child)));
    }
    else
    {
      G4String error_msg = "Unknown tag in physvol: " + tag;
      G4Exception("G4GDMLReadStructure::PhysvolRead()", "ReadError",
                  FatalException, error_msg);
      return;
    }
  }

  G4Transform3D transform(GetRotationMatrix(rotation).inverse(), position);
  transform = transform * G4Scale3D(scale.x(), scale.y(), scale.z());

  if(pAssembly != nullptr)  // Fill assembly structure
  {
    if(assembly != nullptr)  // Case of recursive assemblies
    {
      pAssembly->AddPlacedAssembly(assembly, transform);
    }
    if(logvol == nullptr)
    {
      return;
    }
    pAssembly->AddPlacedVolume(logvol, transform);
  }
  else  // Generate physical volume tree or do assembly imprint
  {
    if(assembly != nullptr)
    {
      assembly->MakeImprint(pMotherLogical, transform, 0, check);
    }
    else
    {
      if(logvol == nullptr)
      {
        return;
      }
      G4String pv_name           = logvol->GetName() + "_PV";
      G4PhysicalVolumesPair pair = G4ReflectionFactory::Instance()->Place(
        transform, pv_name, logvol, pMotherLogical, false, copynumber, check);

      if(pair.first != nullptr)
      {
        GeneratePhysvolName(name, pair.first);
      }
      if(pair.second != nullptr)
      {
        GeneratePhysvolName(name, pair.second);
      }
    }
  }
}

// --------------------------------------------------------------------
void G4GDMLReadStructure::ReplicavolRead(
  const xercesc::DOMElement* const replicavolElement, G4int number)
{
  G4LogicalVolume* logvol = nullptr;
  for(xercesc::DOMNode* iter = replicavolElement->getFirstChild();
                        iter != nullptr; iter = iter->getNextSibling())
  {
    if(iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE)
    {
      continue;
    }

    const xercesc::DOMElement* const child =
      dynamic_cast<xercesc::DOMElement*>(iter);
    if(child == nullptr)
    {
      G4Exception("G4GDMLReadStructure::ReplicavolRead()", "InvalidRead",
                  FatalException, "No child found!");
      return;
    }
    const G4String tag = Transcode(child->getTagName());

    if(tag == "volumeref")
    {
      logvol = GetVolume(GenerateName(RefRead(child)));
    }
    else if(tag == "replicate_along_axis")
    {
      if(logvol)
      {
        ReplicaRead(child, logvol, number);
      }
    }
    else
    {
      G4String error_msg = "Unknown tag in ReplicavolRead: " + tag;
      G4Exception("G4GDMLReadStructure::ReplicavolRead()", "ReadError",
                  FatalException, error_msg);
    }
  }
}

// --------------------------------------------------------------------
void G4GDMLReadStructure::ReplicaRead(
  const xercesc::DOMElement* const replicaElement, G4LogicalVolume* logvol,
  G4int number)
{
  G4double width  = 0.0;
  G4double offset = 0.0;
  G4ThreeVector position(0.0, 0.0, 0.0);
  G4ThreeVector rotation(0.0, 0.0, 0.0);
  EAxis axis = kUndefined;
  G4String name;

  for(xercesc::DOMNode* iter = replicaElement->getFirstChild(); iter != nullptr;
      iter                   = iter->getNextSibling())
  {
    if(iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE)
    {
      continue;
    }

    const xercesc::DOMElement* const child =
      dynamic_cast<xercesc::DOMElement*>(iter);
    if(child == nullptr)
    {
      G4Exception("G4GDMLReadStructure::ReplicaRead()", "InvalidRead",
                  FatalException, "No child found!");
      return;
    }
    const G4String tag = Transcode(child->getTagName());

    if(tag == "position")
    {
      VectorRead(child, position);
    }
    else if(tag == "rotation")
    {
      VectorRead(child, rotation);
    }
    else if(tag == "positionref")
    {
      position = GetPosition(GenerateName(RefRead(child)));
    }
    else if(tag == "rotationref")
    {
      rotation = GetRotation(GenerateName(RefRead(child)));
    }
    else if(tag == "direction")
    {
      axis = AxisRead(child);
    }
    else if(tag == "width")
    {
      width = QuantityRead(child);
    }
    else if(tag == "offset")
    {
      offset = QuantityRead(child);
    }
    else
    {
      G4String error_msg = "Unknown tag in ReplicaRead: " + tag;
      G4Exception("G4GDMLReadStructure::ReplicaRead()", "ReadError",
                  FatalException, error_msg);
    }
  }

  G4String pv_name           = logvol->GetName() + "_PV";
  G4PhysicalVolumesPair pair = G4ReflectionFactory::Instance()->Replicate(
    pv_name, logvol, pMotherLogical, axis, number, width, offset);

  if(pair.first != nullptr)
  {
    GeneratePhysvolName(name, pair.first);
  }
  if(pair.second != nullptr)
  {
    GeneratePhysvolName(name, pair.second);
  }
}

// --------------------------------------------------------------------
EAxis G4GDMLReadStructure::AxisRead(
  const xercesc::DOMElement* const axisElement)
{
  EAxis axis = kUndefined;

  const xercesc::DOMNamedNodeMap* const attributes =
    axisElement->getAttributes();
  XMLSize_t attributeCount = attributes->getLength();

  for(XMLSize_t attribute_index = 0; attribute_index < attributeCount;
      ++attribute_index)
  {
    xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

    if(attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE)
    {
      continue;
    }

    const xercesc::DOMAttr* const attribute =
      dynamic_cast<xercesc::DOMAttr*>(attribute_node);
    if(attribute == nullptr)
    {
      G4Exception("G4GDMLReadStructure::AxisRead()", "InvalidRead",
                  FatalException, "No attribute found!");
      return axis;
    }
    const G4String attName  = Transcode(attribute->getName());
    const G4String attValue = Transcode(attribute->getValue());
    if(attName == "x")
    {
      if(eval.Evaluate(attValue) == 1.)
      {
        axis = kXAxis;
      }
    }
    else if(attName == "y")
    {
      if(eval.Evaluate(attValue) == 1.)
      {
        axis = kYAxis;
      }
    }
    else if(attName == "z")
    {
      if(eval.Evaluate(attValue) == 1.)
      {
        axis = kZAxis;
      }
    }
    else if(attName == "rho")
    {
      if(eval.Evaluate(attValue) == 1.)
      {
        axis = kRho;
      }
    }
    else if(attName == "phi")
    {
      if(eval.Evaluate(attValue) == 1.)
      {
        axis = kPhi;
      }
    }
  }

  return axis;
}

// --------------------------------------------------------------------
G4double G4GDMLReadStructure::QuantityRead(
  const xercesc::DOMElement* const readElement)
{
  G4double value = 0.0;
  G4double unit  = 0.0;
  const xercesc::DOMNamedNodeMap* const attributes =
    readElement->getAttributes();
  XMLSize_t attributeCount = attributes->getLength();

  for(XMLSize_t attribute_index = 0; attribute_index < attributeCount;
      ++attribute_index)
  {
    xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

    if(attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE)
    {
      continue;
    }
    const xercesc::DOMAttr* const attribute =
      dynamic_cast<xercesc::DOMAttr*>(attribute_node);
    if(attribute == nullptr)
    {
      G4Exception("G4GDMLReadStructure::QuantityRead()", "InvalidRead",
                  FatalException, "No attribute found!");
      return value;
    }
    const G4String attName  = Transcode(attribute->getName());
    const G4String attValue = Transcode(attribute->getValue());

    if(attName == "unit")
    {
      unit = G4UnitDefinition::GetValueOf(attValue);
      if(G4UnitDefinition::GetCategory(attValue) != "Length" &&
         G4UnitDefinition::GetCategory(attValue) != "Angle")
      {
        G4Exception("G4GDMLReadStructure::QuantityRead()", "InvalidRead",
                    FatalException,
                    "Invalid unit for length or angle (width, offset)!");
      }
    }
    else if(attName == "value")
    {
      value = eval.Evaluate(attValue);
    }
  }

  return value * unit;
}

// --------------------------------------------------------------------
void G4GDMLReadStructure::VolumeRead(
  const xercesc::DOMElement* const volumeElement)
{
  G4VSolid* solidPtr      = nullptr;
  G4Material* materialPtr = nullptr;
  G4GDMLAuxListType auxList;

  XMLCh* name_attr    = xercesc::XMLString::transcode("name");
  const G4String name = Transcode(volumeElement->getAttribute(name_attr));
  xercesc::XMLString::release(&name_attr);

  for(xercesc::DOMNode* iter = volumeElement->getFirstChild(); iter != nullptr;
      iter                   = iter->getNextSibling())
  {
    if(iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE)
    {
      continue;
    }

    const xercesc::DOMElement* const child =
      dynamic_cast<xercesc::DOMElement*>(iter);
    if(child == nullptr)
    {
      G4Exception("G4GDMLReadStructure::VolumeRead()", "InvalidRead",
                  FatalException, "No child found!");
      return;
    }
    const G4String tag = Transcode(child->getTagName());

    if(tag == "auxiliary")
    {
      auxList.push_back(AuxiliaryRead(child));
    }
    else if(tag == "materialref")
    {
      materialPtr = GetMaterial(GenerateName(RefRead(child), true));
    }
    else if(tag == "solidref")
    {
      solidPtr = GetSolid(GenerateName(RefRead(child)));
    }
  }

  pMotherLogical =
    new G4LogicalVolume(solidPtr, materialPtr, GenerateName(name), 0, 0, 0);

  if(!auxList.empty())
  {
    auxMap[pMotherLogical] = auxList;
  }

  Volume_contentRead(volumeElement);
}

// --------------------------------------------------------------------
void G4GDMLReadStructure::AssemblyRead(
  const xercesc::DOMElement* const assemblyElement)
{
  XMLCh* name_attr    = xercesc::XMLString::transcode("name");
  const G4String name = Transcode(assemblyElement->getAttribute(name_attr));
  xercesc::XMLString::release(&name_attr);

  G4AssemblyVolume* pAssembly = new G4AssemblyVolume();
  auto aName = GenerateName(name);
  if(reverseSearch)
  {
    assemblyMap.insert_or_assign(aName, pAssembly);
  }
  else
  { 
    assemblyMap.insert(std::make_pair(aName, pAssembly));
  }

  for(xercesc::DOMNode* iter = assemblyElement->getFirstChild();
                        iter != nullptr; iter = iter->getNextSibling())
  {
    if(iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE)
    {
      continue;
    }
    const xercesc::DOMElement* const child =
      dynamic_cast<xercesc::DOMElement*>(iter);
    if(child == nullptr)
    {
      G4Exception("G4GDMLReadStructure::AssemblyRead()", "InvalidRead",
                  FatalException, "No child found!");
      return;
    }
    const G4String tag = Transcode(child->getTagName());

    if(tag == "physvol")
    {
      PhysvolRead(child, pAssembly);
    }
    else
    {
      G4cout << "Unsupported GDML tag '" << tag
             << "' for Geant4 assembly structure !" << G4endl;
    }
  }
}

// --------------------------------------------------------------------
void G4GDMLReadStructure::SkinSurfaceRead(
  const xercesc::DOMElement* const skinsurfaceElement)
{
  G4String name;
  G4LogicalVolume* logvol = nullptr;
  G4SurfaceProperty* prop = nullptr;

  const xercesc::DOMNamedNodeMap* const attributes =
    skinsurfaceElement->getAttributes();
  XMLSize_t attributeCount = attributes->getLength();

  for(XMLSize_t attribute_index = 0; attribute_index < attributeCount;
      ++attribute_index)
  {
    xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

    if(attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE)
    {
      continue;
    }

    const xercesc::DOMAttr* const attribute =
      dynamic_cast<xercesc::DOMAttr*>(attribute_node);
    if(attribute == nullptr)
    {
      G4Exception("G4GDMLReadStructure::SkinsurfaceRead()", "InvalidRead",
                  FatalException, "No attribute found!");
      return;
    }
    const G4String attName  = Transcode(attribute->getName());
    const G4String attValue = Transcode(attribute->getValue());

    if(attName == "name")
    {
      name = GenerateName(attValue);
    }
    else if(attName == "surfaceproperty")
    {
      prop = GetSurfaceProperty(GenerateName(attValue));
    }
  }

  for(xercesc::DOMNode* iter = skinsurfaceElement->getFirstChild();
                        iter != nullptr; iter = iter->getNextSibling())
  {
    if(iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE)
    {
      continue;
    }

    const xercesc::DOMElement* const child =
      dynamic_cast<xercesc::DOMElement*>(iter);
    if(child == nullptr)
    {
      G4Exception("G4GDMLReadStructure::SkinsurfaceRead()", "InvalidRead",
                  FatalException, "No child found!");
      return;
    }
    const G4String tag = Transcode(child->getTagName());

    if(tag == "volumeref")
    {
      logvol = GetVolume(GenerateName(RefRead(child)));
    }
    else
    {
      G4String error_msg = "Unknown tag in skinsurface: " + tag;
      G4Exception("G4GDMLReadStructure::SkinsurfaceRead()", "ReadError",
                  FatalException, error_msg);
    }
  }

  new G4LogicalSkinSurface(Strip(name), logvol, prop);
}

// --------------------------------------------------------------------
void G4GDMLReadStructure::Volume_contentRead(
  const xercesc::DOMElement* const volumeElement)
{
  for(xercesc::DOMNode* iter = volumeElement->getFirstChild(); iter != nullptr;
      iter                   = iter->getNextSibling())
  {
    if(iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE)
    {
      continue;
    }

    const xercesc::DOMElement* const child =
      dynamic_cast<xercesc::DOMElement*>(iter);
    if(child == nullptr)
    {
      G4Exception("G4GDMLReadStructure::Volume_contentRead()", "InvalidRead",
                  FatalException, "No child found!");
      return;
    }
    const G4String tag = Transcode(child->getTagName());

    if((tag == "auxiliary") || (tag == "materialref") || (tag == "solidref"))
    {
      // These are already processed in VolumeRead()
    }
    else if(tag == "paramvol")
    {
      ParamvolRead(child, pMotherLogical);
    }
    else if(tag == "physvol")
    {
      PhysvolRead(child);
    }
    else if(tag == "replicavol")
    {
      G4int number = 1;
      const xercesc::DOMNamedNodeMap* const attributes = child->getAttributes();
      XMLSize_t attributeCount = attributes->getLength();
      for(XMLSize_t attribute_index = 0; attribute_index < attributeCount;
          ++attribute_index)
      {
        xercesc::DOMNode* attribute_node = attributes->item(attribute_index);
        if(attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE)
        {
          continue;
        }
        const xercesc::DOMAttr* const attribute =
          dynamic_cast<xercesc::DOMAttr*>(attribute_node);
        if(attribute == nullptr)
        {
          G4Exception("G4GDMLReadStructure::Volume_contentRead()",
                      "InvalidRead", FatalException, "No attribute found!");
          return;
        }
        const G4String attName  = Transcode(attribute->getName());
        const G4String attValue = Transcode(attribute->getValue());
        if(attName == "number")
        {
          number = eval.EvaluateInteger(attValue);
        }
      }
      ReplicavolRead(child, number);
    }
    else if(tag == "divisionvol")
    {
      DivisionvolRead(child);
    }
    else if(tag == "loop")
    {
      LoopRead(child, &G4GDMLRead::Volume_contentRead);
    }
    else
    {
      G4cout << "Treating unknown GDML tag in volume '" << tag
             << "' as GDML extension..." << G4endl;
    }
  }
}

// --------------------------------------------------------------------
void G4GDMLReadStructure::StructureRead(
  const xercesc::DOMElement* const structureElement)
{
#ifdef G4VERBOSE
  G4cout << "G4GDML: Reading structure..." << G4endl;
#endif
  for(xercesc::DOMNode* iter = structureElement->getFirstChild();
                        iter != nullptr; iter = iter->getNextSibling())
  {
    if(iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE)
    {
      continue;
    }

    const xercesc::DOMElement* const child =
      dynamic_cast<xercesc::DOMElement*>(iter);
    if(child == nullptr)
    {
      G4Exception("G4GDMLReadStructure::StructureRead()", "InvalidRead",
                  FatalException, "No child found!");
      return;
    }
    const G4String tag = Transcode(child->getTagName());

    if(tag == "bordersurface")
    {
      BorderSurfaceRead(child);
    }
    else if(tag == "skinsurface")
    {
      SkinSurfaceRead(child);
    }
    else if(tag == "volume")
    {
      VolumeRead(child);
    }
    else if(tag == "assembly")
    {
      AssemblyRead(child);
    }
    else if(tag == "loop")
    {
      LoopRead(child, &G4GDMLRead::StructureRead);
    }
    else
    {
      G4String error_msg = "Unknown tag in structure: " + tag;
      G4Exception("G4GDMLReadStructure::StructureRead()", "ReadError",
                  FatalException, error_msg);
    }
  }
}

// --------------------------------------------------------------------
G4VPhysicalVolume* G4GDMLReadStructure::GetPhysvol(const G4String& ref) const
{
  G4VPhysicalVolume* physvolPtr
    = G4PhysicalVolumeStore::GetInstance()->GetVolume(ref,false,reverseSearch);

  if(physvolPtr == nullptr)
  {
    G4String error_msg = "Referenced physvol '" + ref + "' was not found!";
    G4Exception("G4GDMLReadStructure::GetPhysvol()", "ReadError",
                FatalException, error_msg);
  }

  return physvolPtr;
}

// --------------------------------------------------------------------
G4LogicalVolume* G4GDMLReadStructure::GetVolume(const G4String& ref) const
{
  G4LogicalVolume* volumePtr
    = G4LogicalVolumeStore::GetInstance()->GetVolume(ref,false,reverseSearch);

  if(volumePtr == nullptr)
  {
    G4String error_msg = "Referenced volume '" + ref + "' was not found!";
    G4Exception("G4GDMLReadStructure::GetVolume()", "ReadError", FatalException,
                error_msg);
  }

  return volumePtr;
}

// --------------------------------------------------------------------
G4AssemblyVolume* G4GDMLReadStructure::GetAssembly(const G4String& ref) const
{
  auto pos = assemblyMap.find(ref);
  if(pos != assemblyMap.cend())
  {
    return pos->second;
  }
  return nullptr;
}

// --------------------------------------------------------------------
G4GDMLAuxListType G4GDMLReadStructure::GetVolumeAuxiliaryInformation(
  G4LogicalVolume* logvol) const
{
  auto pos = auxMap.find(logvol);
  if(pos != auxMap.cend())
  {
    return pos->second;
  }
  else
  {
    return G4GDMLAuxListType();
  }
}

// --------------------------------------------------------------------
G4VPhysicalVolume* G4GDMLReadStructure::GetWorldVolume(
  const G4String& setupName)
{
  G4String sname = GetSetup(setupName);
  if(sname == "")
  {
    return nullptr;
  }

  G4LogicalVolume* volume = GetVolume(GenerateName(sname, dostrip));
  volume->SetVisAttributes(G4VisAttributes::GetInvisible());

  G4VPhysicalVolume* pvWorld = nullptr;

  if(setuptoPV[setupName])
  {
    pvWorld = setuptoPV[setupName];
  }
  else
  {
    pvWorld = new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), volume,
                                volume->GetName() + "_PV", 0, 0, 0);
    setuptoPV[setupName] = pvWorld;
  }
  return pvWorld;
}

// --------------------------------------------------------------------
void G4GDMLReadStructure::Clear()
{
  eval.Clear();
  setuptoPV.clear();
  auxMap.clear();
}
