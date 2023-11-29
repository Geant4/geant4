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
// G4GDMLReadParamvol implementation
//
// Author: Zoltan Torzsok, November 2007
// --------------------------------------------------------------------

#include "G4GDMLReadParamvol.hh"
#include "G4GDMLReadSolids.hh"
#include "G4LogicalVolume.hh"
#include "G4PVParameterised.hh"
#include "G4PVPlacement.hh"
#include "G4VPhysicalVolume.hh"
#include "G4UnitsTable.hh"

// --------------------------------------------------------------------
G4GDMLReadParamvol::G4GDMLReadParamvol()
  : G4GDMLReadSetup()
{
}

G4GDMLReadParamvol::~G4GDMLReadParamvol()
{
}

// --------------------------------------------------------------------
void G4GDMLReadParamvol::Box_dimensionsRead(
  const xercesc::DOMElement* const element,
  G4GDMLParameterisation::PARAMETER& parameter)
{
  G4double lunit = 1.0;

  const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
  XMLSize_t attributeCount                         = attributes->getLength();

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
      G4Exception("G4GDMLReadParamvol::Box_dimensionsRead()", "InvalidRead",
                  FatalException, "No attribute found!");
      return;
    }
    const G4String attName  = Transcode(attribute->getName());
    const G4String attValue = Transcode(attribute->getValue());

    if(attName == "lunit")
    {
      lunit = G4UnitDefinition::GetValueOf(attValue);
      if(G4UnitDefinition::GetCategory(attValue) != "Length")
      {
        G4Exception("G4GDMLReadParamvol::Box_dimensionsRead()", "InvalidRead",
                    FatalException, "Invalid unit for length!");
      }
    }
    else if(attName == "x")
    {
      parameter.dimension[0] = eval.Evaluate(attValue);
    }
    else if(attName == "y")
    {
      parameter.dimension[1] = eval.Evaluate(attValue);
    }
    else if(attName == "z")
    {
      parameter.dimension[2] = eval.Evaluate(attValue);
    }
  }

  parameter.dimension[0] *= 0.5 * lunit;
  parameter.dimension[1] *= 0.5 * lunit;
  parameter.dimension[2] *= 0.5 * lunit;
}

// --------------------------------------------------------------------
void G4GDMLReadParamvol::Trd_dimensionsRead(
  const xercesc::DOMElement* const element,
  G4GDMLParameterisation::PARAMETER& parameter)
{
  G4double lunit = 1.0;

  const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
  XMLSize_t attributeCount                         = attributes->getLength();

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
      G4Exception("G4GDMLReadParamvol::Trd_dimensionsRead()", "InvalidRead",
                  FatalException, "No attribute found!");
      return;
    }
    const G4String attName  = Transcode(attribute->getName());
    const G4String attValue = Transcode(attribute->getValue());

    if(attName == "lunit")
    {
      lunit = G4UnitDefinition::GetValueOf(attValue);
      if(G4UnitDefinition::GetCategory(attValue) != "Length")
      {
        G4Exception("G4GDMLReadParamvol::Trd_dimensionsRead()", "InvalidRead",
                    FatalException, "Invalid unit for length!");
      }
    }
    else if(attName == "x1")
    {
      parameter.dimension[0] = eval.Evaluate(attValue);
    }
    else if(attName == "x2")
    {
      parameter.dimension[1] = eval.Evaluate(attValue);
    }
    else if(attName == "y1")
    {
      parameter.dimension[2] = eval.Evaluate(attValue);
    }
    else if(attName == "y2")
    {
      parameter.dimension[3] = eval.Evaluate(attValue);
    }
    else if(attName == "z")
    {
      parameter.dimension[4] = eval.Evaluate(attValue);
    }
  }

  parameter.dimension[0] *= 0.5 * lunit;
  parameter.dimension[1] *= 0.5 * lunit;
  parameter.dimension[2] *= 0.5 * lunit;
  parameter.dimension[3] *= 0.5 * lunit;
  parameter.dimension[4] *= 0.5 * lunit;
}

// --------------------------------------------------------------------
void G4GDMLReadParamvol::Trap_dimensionsRead(
  const xercesc::DOMElement* const element,
  G4GDMLParameterisation::PARAMETER& parameter)
{
  G4double lunit = 1.0;
  G4double aunit = 1.0;

  const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
  XMLSize_t attributeCount                         = attributes->getLength();

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
      G4Exception("G4GDMLReadParamvol::Trap_dimensionsRead()", "InvalidRead",
                  FatalException, "No attribute found!");
      return;
    }
    const G4String attName  = Transcode(attribute->getName());
    const G4String attValue = Transcode(attribute->getValue());

    if(attName == "lunit")
    {
      lunit = G4UnitDefinition::GetValueOf(attValue);
      if(G4UnitDefinition::GetCategory(attValue) != "Length")
      {
        G4Exception("G4GDMLReadParamvol::Trap_dimensionsRead()", "InvalidRead",
                    FatalException, "Invalid unit for length!");
      }
    }
    else if(attName == "aunit")
    {
      aunit = G4UnitDefinition::GetValueOf(attValue);
      if(G4UnitDefinition::GetCategory(attValue) != "Angle")
      {
        G4Exception("G4GDMLReadParamvol::Trap_dimensionsRead()", "InvalidRead",
                    FatalException, "Invalid unit for angle!");
      }
    }
    else if(attName == "z")
    {
      parameter.dimension[0] = eval.Evaluate(attValue);
    }
    else if(attName == "theta")
    {
      parameter.dimension[1] = eval.Evaluate(attValue);
    }
    else if(attName == "phi")
    {
      parameter.dimension[2] = eval.Evaluate(attValue);
    }
    else if(attName == "y1")
    {
      parameter.dimension[3] = eval.Evaluate(attValue);
    }
    else if(attName == "x1")
    {
      parameter.dimension[4] = eval.Evaluate(attValue);
    }
    else if(attName == "x2")
    {
      parameter.dimension[5] = eval.Evaluate(attValue);
    }
    else if(attName == "alpha1")
    {
      parameter.dimension[6] = eval.Evaluate(attValue);
    }
    else if(attName == "y2")
    {
      parameter.dimension[7] = eval.Evaluate(attValue);
    }
    else if(attName == "x3")
    {
      parameter.dimension[8] = eval.Evaluate(attValue);
    }
    else if(attName == "x4")
    {
      parameter.dimension[9] = eval.Evaluate(attValue);
    }
    else if(attName == "alpha2")
    {
      parameter.dimension[10] = eval.Evaluate(attValue);
    }
  }

  parameter.dimension[0] *= 0.5 * lunit;
  parameter.dimension[1] *= aunit;
  parameter.dimension[2] *= aunit;
  parameter.dimension[3] *= 0.5 * lunit;
  parameter.dimension[4] *= 0.5 * lunit;
  parameter.dimension[5] *= 0.5 * lunit;
  parameter.dimension[6] *= aunit;
  parameter.dimension[7] *= 0.5 * lunit;
  parameter.dimension[8] *= 0.5 * lunit;
  parameter.dimension[9] *= 0.5 * lunit;
  parameter.dimension[10] *= aunit;
}

// --------------------------------------------------------------------
void G4GDMLReadParamvol::Tube_dimensionsRead(
  const xercesc::DOMElement* const element,
  G4GDMLParameterisation::PARAMETER& parameter)
{
  G4double lunit = 1.0;
  G4double aunit = 1.0;

  const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
  XMLSize_t attributeCount                         = attributes->getLength();

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
      G4Exception("G4GDMLReadParamvol::Tube_dimensionsRead()", "InvalidRead",
                  FatalException, "No attribute found!");
      return;
    }
    const G4String attName  = Transcode(attribute->getName());
    const G4String attValue = Transcode(attribute->getValue());

    if(attName == "lunit")
    {
      lunit = G4UnitDefinition::GetValueOf(attValue);
      if(G4UnitDefinition::GetCategory(attValue) != "Length")
      {
        G4Exception("G4GDMLReadParamvol::Tube_dimensionsRead()", "InvalidRead",
                    FatalException, "Invalid unit for length!");
      }
    }
    else if(attName == "aunit")
    {
      aunit = G4UnitDefinition::GetValueOf(attValue);
      if(G4UnitDefinition::GetCategory(attValue) != "Angle")
      {
        G4Exception("G4GDMLReadParamvol::Tube_dimensionsRead()", "InvalidRead",
                    FatalException, "Invalid unit for angle!");
      }
    }
    else if(attName == "InR")
    {
      parameter.dimension[0] = eval.Evaluate(attValue);
    }
    else if(attName == "OutR")
    {
      parameter.dimension[1] = eval.Evaluate(attValue);
    }
    else if(attName == "hz")
    {
      parameter.dimension[2] = eval.Evaluate(attValue);
    }
    else if(attName == "StartPhi")
    {
      parameter.dimension[3] = eval.Evaluate(attValue);
    }
    else if(attName == "DeltaPhi")
    {
      parameter.dimension[4] = eval.Evaluate(attValue);
    }
  }

  parameter.dimension[0] *= lunit;
  parameter.dimension[1] *= lunit;
  parameter.dimension[2] *= 0.5 * lunit;
  parameter.dimension[3] *= aunit;
  parameter.dimension[4] *= aunit;
}

// --------------------------------------------------------------------
void G4GDMLReadParamvol::Cone_dimensionsRead(
  const xercesc::DOMElement* const element,
  G4GDMLParameterisation::PARAMETER& parameter)
{
  G4double lunit = 1.0;
  G4double aunit = 1.0;

  const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
  XMLSize_t attributeCount                         = attributes->getLength();

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
      G4Exception("G4GDMLReadParamvol::Cone_dimensionsRead()", "InvalidRead",
                  FatalException, "No attribute found!");
      return;
    }
    const G4String attName  = Transcode(attribute->getName());
    const G4String attValue = Transcode(attribute->getValue());

    if(attName == "lunit")
    {
      lunit = G4UnitDefinition::GetValueOf(attValue);
      if(G4UnitDefinition::GetCategory(attValue) != "Length")
      {
        G4Exception("G4GDMLReadParamvol::Cone_dimensionsRead()", "InvalidRead",
                    FatalException, "Invalid unit for length!");
      }
    }
    else if(attName == "aunit")
    {
      aunit = G4UnitDefinition::GetValueOf(attValue);
      if(G4UnitDefinition::GetCategory(attValue) != "Angle")
      {
        G4Exception("G4GDMLReadParamvol::Cone_dimensions()", "InvalidRead",
                    FatalException, "Invalid unit for angle!");
      }
    }
    else if(attName == "rmin1")
    {
      parameter.dimension[0] = eval.Evaluate(attValue);
    }
    else if(attName == "rmax1")
    {
      parameter.dimension[1] = eval.Evaluate(attValue);
    }
    else if(attName == "rmin2")
    {
      parameter.dimension[2] = eval.Evaluate(attValue);
    }
    else if(attName == "rmax2")
    {
      parameter.dimension[3] = eval.Evaluate(attValue);
    }
    else if(attName == "z")
    {
      parameter.dimension[4] = eval.Evaluate(attValue);
    }
    else if(attName == "startphi")
    {
      parameter.dimension[5] = eval.Evaluate(attValue);
    }
    else if(attName == "deltaphi")
    {
      parameter.dimension[6] = eval.Evaluate(attValue);
    }
  }

  parameter.dimension[0] *= lunit;
  parameter.dimension[1] *= lunit;
  parameter.dimension[2] *= lunit;
  parameter.dimension[3] *= lunit;
  parameter.dimension[4] *= 0.5 * lunit;
  parameter.dimension[5] *= aunit;
  parameter.dimension[6] *= aunit;
}

// --------------------------------------------------------------------
void G4GDMLReadParamvol::Sphere_dimensionsRead(
  const xercesc::DOMElement* const element,
  G4GDMLParameterisation::PARAMETER& parameter)
{
  G4double lunit = 1.0;
  G4double aunit = 1.0;

  const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
  XMLSize_t attributeCount                         = attributes->getLength();

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
      G4Exception("G4GDMLReadParamvol::Sphere_dimensionsRead()", "InvalidRead",
                  FatalException, "No attribute found!");
      return;
    }
    const G4String attName  = Transcode(attribute->getName());
    const G4String attValue = Transcode(attribute->getValue());

    if(attName == "lunit")
    {
      lunit = G4UnitDefinition::GetValueOf(attValue);
      if(G4UnitDefinition::GetCategory(attValue) != "Length")
      {
        G4Exception("G4GDMLReadParamvol::Sphere_dimensionsRead()",
                    "InvalidRead", FatalException, "Invalid unit for length!");
      }
    }
    else if(attName == "aunit")
    {
      aunit = G4UnitDefinition::GetValueOf(attValue);
      if(G4UnitDefinition::GetCategory(attValue) != "Angle")
      {
        G4Exception("G4GDMLReadParamvol::Sphere_dimensionsRead()",
                    "InvalidRead", FatalException, "Invalid unit for angle!");
      }
    }
    else if(attName == "rmin")
    {
      parameter.dimension[0] = eval.Evaluate(attValue);
    }
    else if(attName == "rmax")
    {
      parameter.dimension[1] = eval.Evaluate(attValue);
    }
    else if(attName == "startphi")
    {
      parameter.dimension[2] = eval.Evaluate(attValue);
    }
    else if(attName == "deltaphi")
    {
      parameter.dimension[3] = eval.Evaluate(attValue);
    }
    else if(attName == "starttheta")
    {
      parameter.dimension[4] = eval.Evaluate(attValue);
    }
    else if(attName == "deltatheta")
    {
      parameter.dimension[5] = eval.Evaluate(attValue);
    }
  }

  parameter.dimension[0] *= lunit;
  parameter.dimension[1] *= lunit;
  parameter.dimension[2] *= aunit;
  parameter.dimension[3] *= aunit;
  parameter.dimension[4] *= aunit;
  parameter.dimension[5] *= aunit;
}

// --------------------------------------------------------------------
void G4GDMLReadParamvol::Orb_dimensionsRead(
  const xercesc::DOMElement* const element,
  G4GDMLParameterisation::PARAMETER& parameter)
{
  G4double lunit = 1.0;

  const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
  XMLSize_t attributeCount                         = attributes->getLength();

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
      G4Exception("G4GDMLReadParamvol::Orb_dimensionsRead()", "InvalidRead",
                  FatalException, "No attribute found!");
      return;
    }
    const G4String attName  = Transcode(attribute->getName());
    const G4String attValue = Transcode(attribute->getValue());

    if(attName == "lunit")
    {
      lunit = G4UnitDefinition::GetValueOf(attValue);
    }
    if(G4UnitDefinition::GetCategory(attValue) != "Length")
    {
      G4Exception("G4GDMLReadParamvol::Orb_dimensionsRead()", "InvalidRead",
                  FatalException, "Invalid unit for length!");
    }
    else if(attName == "r")
    {
      parameter.dimension[0] = eval.Evaluate(attValue);
    }
  }

  parameter.dimension[0] *= lunit;
}

// --------------------------------------------------------------------
void G4GDMLReadParamvol::Torus_dimensionsRead(
  const xercesc::DOMElement* const element,
  G4GDMLParameterisation::PARAMETER& parameter)
{
  G4double lunit = 1.0;
  G4double aunit = 1.0;

  const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
  XMLSize_t attributeCount                         = attributes->getLength();

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
      G4Exception("G4GDMLReadParamvol::Torus_dimensionsRead()", "InvalidRead",
                  FatalException, "No attribute found!");
      return;
    }
    const G4String attName  = Transcode(attribute->getName());
    const G4String attValue = Transcode(attribute->getValue());

    if(attName == "lunit")
    {
      lunit = G4UnitDefinition::GetValueOf(attValue);
      if(G4UnitDefinition::GetCategory(attValue) != "Length")
      {
        G4Exception("G4GDMLReadParamvol::Torus_dimensionsRead()", "InvalidRead",
                    FatalException, "Invalid unit for length!");
      }
    }
    else if(attName == "aunit")
    {
      aunit = G4UnitDefinition::GetValueOf(attValue);
      if(G4UnitDefinition::GetCategory(attValue) != "Angle")
      {
        G4Exception("G4GDMLReadParamvol::Torus_dimensionsRead()", "InvalidRead",
                    FatalException, "Invalid unit for angle!");
      }
    }
    else if(attName == "rmin")
    {
      parameter.dimension[0] = eval.Evaluate(attValue);
    }
    else if(attName == "rmax")
    {
      parameter.dimension[1] = eval.Evaluate(attValue);
    }
    else if(attName == "rtor")
    {
      parameter.dimension[2] = eval.Evaluate(attValue);
    }
    else if(attName == "startphi")
    {
      parameter.dimension[3] = eval.Evaluate(attValue);
    }
    else if(attName == "deltaphi")
    {
      parameter.dimension[4] = eval.Evaluate(attValue);
    }
  }

  parameter.dimension[0] *= lunit;
  parameter.dimension[1] *= lunit;
  parameter.dimension[2] *= lunit;
  parameter.dimension[3] *= aunit;
  parameter.dimension[4] *= aunit;
}

// --------------------------------------------------------------------
void G4GDMLReadParamvol::Ellipsoid_dimensionsRead(
  const xercesc::DOMElement* const element,
  G4GDMLParameterisation::PARAMETER& parameter)
{
  G4double lunit                                   = 1.0;
  const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
  XMLSize_t attributeCount                         = attributes->getLength();

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
      G4Exception("G4GDMLReadParamvol::Ellipsoid_dimensionsRead()",
                  "InvalidRead", FatalException, "No attribute found!");
      return;
    }
    const G4String attName  = Transcode(attribute->getName());
    const G4String attValue = Transcode(attribute->getValue());

    if(attName == "lunit")
    {
      lunit = G4UnitDefinition::GetValueOf(attValue);
      if(G4UnitDefinition::GetCategory(attValue) != "Length")
      {
        G4Exception("G4GDMLReadParamvol::Ellipsoid_dimensionsRead()",
                    "InvalidRead", FatalException, "Invalid unit for length!");
      }
    }
    else if(attName == "ax")
    {
      parameter.dimension[0] = eval.Evaluate(attValue);
    }
    else if(attName == "by")
    {
      parameter.dimension[1] = eval.Evaluate(attValue);
    }
    else if(attName == "cz")
    {
      parameter.dimension[2] = eval.Evaluate(attValue);
    }
    else if(attName == "zcut1")
    {
      parameter.dimension[3] = eval.Evaluate(attValue);
    }
    else if(attName == "zcut2")
    {
      parameter.dimension[4] = eval.Evaluate(attValue);
    }
  }

  parameter.dimension[0] *= lunit;
  parameter.dimension[1] *= lunit;
  parameter.dimension[2] *= lunit;
  parameter.dimension[3] *= lunit;
  parameter.dimension[4] *= lunit;
}

// --------------------------------------------------------------------
void G4GDMLReadParamvol::Para_dimensionsRead(
  const xercesc::DOMElement* const element,
  G4GDMLParameterisation::PARAMETER& parameter)
{
  G4double lunit = 1.0;
  G4double aunit = 1.0;

  const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
  XMLSize_t attributeCount                         = attributes->getLength();

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
      G4Exception("G4GDMLReadParamvol::Para_dimensionsRead()", "InvalidRead",
                  FatalException, "No attribute found!");
      return;
    }
    const G4String attName  = Transcode(attribute->getName());
    const G4String attValue = Transcode(attribute->getValue());

    if(attName == "lunit")
    {
      lunit = G4UnitDefinition::GetValueOf(attValue);
      if(G4UnitDefinition::GetCategory(attValue) != "Length")
      {
        G4Exception("G4GDMLReadParamvol::Para_dimensionsRead()", "InvalidRead",
                    FatalException, "Invalid unit for length!");
      }
    }
    else if(attName == "aunit")
    {
      aunit = G4UnitDefinition::GetValueOf(attValue);
      if(G4UnitDefinition::GetCategory(attValue) != "Angle")
      {
        G4Exception("G4GDMLReadParamvol::Para_dimensionsRead()", "InvalidRead",
                    FatalException, "Invalid unit for angle!");
      }
    }
    else if(attName == "x")
    {
      parameter.dimension[0] = eval.Evaluate(attValue);
    }
    else if(attName == "y")
    {
      parameter.dimension[1] = eval.Evaluate(attValue);
    }
    else if(attName == "z")
    {
      parameter.dimension[2] = eval.Evaluate(attValue);
    }
    else if(attName == "alpha")
    {
      parameter.dimension[3] = eval.Evaluate(attValue);
    }
    else if(attName == "theta")
    {
      parameter.dimension[4] = eval.Evaluate(attValue);
    }
    else if(attName == "phi")
    {
      parameter.dimension[5] = eval.Evaluate(attValue);
    }
  }

  parameter.dimension[0] = 0.5 * lunit;
  parameter.dimension[1] = 0.5 * lunit;
  parameter.dimension[2] = 0.5 * lunit;
  parameter.dimension[3] = aunit;
  parameter.dimension[4] = aunit;
  parameter.dimension[5] = aunit;
}

// --------------------------------------------------------------------
void G4GDMLReadParamvol::Hype_dimensionsRead(
  const xercesc::DOMElement* const element,
  G4GDMLParameterisation::PARAMETER& parameter)
{
  G4double lunit = 1.0;
  G4double aunit = 1.0;

  const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
  XMLSize_t attributeCount                         = attributes->getLength();

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
      G4Exception("G4GDMLReadParamvol::Hype_dimensionsRead()", "InvalidRead",
                  FatalException, "No attribute found!");
      return;
    }
    const G4String attName  = Transcode(attribute->getName());
    const G4String attValue = Transcode(attribute->getValue());

    if(attName == "lunit")
    {
      lunit = G4UnitDefinition::GetValueOf(attValue);
      if(G4UnitDefinition::GetCategory(attValue) != "Length")
      {
        G4Exception("G4GDMLReadParamvol::Hype_dimensionsRead()", "InvalidRead",
                    FatalException, "Invalid unit for length!");
      }
    }
    else if(attName == "aunit")
    {
      aunit = G4UnitDefinition::GetValueOf(attValue);
      if(G4UnitDefinition::GetCategory(attValue) != "Angle")
      {
        G4Exception("G4GDMLReadParamvol::Hype_dimensionsRead()", "InvalidRead",
                    FatalException, "Invalid unit for angle!");
      }
    }
    else if(attName == "rmin")
    {
      parameter.dimension[0] = eval.Evaluate(attValue);
    }
    else if(attName == "rmax")
    {
      parameter.dimension[1] = eval.Evaluate(attValue);
    }
    else if(attName == "inst")
    {
      parameter.dimension[2] = eval.Evaluate(attValue);
    }
    else if(attName == "outst")
    {
      parameter.dimension[3] = eval.Evaluate(attValue);
    }
    else if(attName == "z")
    {
      parameter.dimension[4] = eval.Evaluate(attValue);
    }
  }

  parameter.dimension[0] = lunit;
  parameter.dimension[1] = lunit;
  parameter.dimension[2] = aunit;
  parameter.dimension[3] = aunit;
  parameter.dimension[4] = 0.5 * lunit;
}

// --------------------------------------------------------------------
void G4GDMLReadParamvol::Polycone_dimensionsRead(
  const xercesc::DOMElement* const element,
  G4GDMLParameterisation::PARAMETER& parameter)
{
  G4double lunit = 1.0;
  G4double aunit = 1.0;

  std::vector<zplaneType> zplaneList;

  const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
  XMLSize_t attributeCount                         = attributes->getLength();

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
      G4Exception("G4GDMLReadParamvol::Polycone_dimensionsRead()",
                  "InvalidRead", FatalException, "No attribute found!");
      return;
    }
    const G4String attName  = Transcode(attribute->getName());
    const G4String attValue = Transcode(attribute->getValue());

    if(attName == "lunit")
    {
      lunit = G4UnitDefinition::GetValueOf(attValue);
      if(G4UnitDefinition::GetCategory(attValue) != "Length")
      {
        G4Exception("G4GDMLReadParamvol::Polycone_dimensionsRead()",
                    "InvalidRead", FatalException, "Invalid unit for length!");
      }
    }
    else if(attName == "aunit")
    {
      aunit = G4UnitDefinition::GetValueOf(attValue);
      if(G4UnitDefinition::GetCategory(attValue) != "Angle")
      {
        G4Exception("G4GDMLReadParamvol::Polycone_dimensionsRead()",
                    "InvalidRead", FatalException, "Invalid unit for angle!");
      }
    }
    else if(attName == "startPhi")
    {
      parameter.dimension[0] = eval.Evaluate(attValue);
    }
    else if(attName == "openPhi")
    {
      parameter.dimension[1] = eval.Evaluate(attValue);
    }
    else if(attName == "numRZ")
    {
      parameter.dimension[2] = eval.Evaluate(attValue);
    }
  }

  parameter.dimension[0] *= aunit;
  parameter.dimension[1] *= aunit;

  for(xercesc::DOMNode* iter = element->getFirstChild(); iter != nullptr;
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
      G4Exception("G4GDMLReadParamVol::Polycone_dimensionsRead()",
                  "InvalidRead", FatalException, "No child found!");
      return;
    }
    const G4String tag = Transcode(child->getTagName());

    if(tag == "zplane")
    {
      zplaneList.push_back(ZplaneRead(child));
    }
  }

  std::size_t numZPlanes = zplaneList.size();
  for(std::size_t i = 0; i < numZPlanes; ++i)
  {
    parameter.dimension[3 + i * 3] = zplaneList[i].rmin * lunit;
    parameter.dimension[4 + i * 3] = zplaneList[i].rmax * lunit;
    parameter.dimension[5 + i * 3] = zplaneList[i].z * lunit;
  }
}

// --------------------------------------------------------------------
void G4GDMLReadParamvol::Polyhedra_dimensionsRead(
  const xercesc::DOMElement* const element,
  G4GDMLParameterisation::PARAMETER& parameter)
{
  G4double lunit = 1.0;
  G4double aunit = 1.0;

  std::vector<zplaneType> zplaneList;

  const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
  XMLSize_t attributeCount                         = attributes->getLength();

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
      G4Exception("G4GDMLReadParamvol::Polycone_dimensionsRead()",
                  "InvalidRead", FatalException, "No attribute found!");
      return;
    }
    const G4String attName  = Transcode(attribute->getName());
    const G4String attValue = Transcode(attribute->getValue());

    if(attName == "lunit")
    {
      lunit = G4UnitDefinition::GetValueOf(attValue);
      if(G4UnitDefinition::GetCategory(attValue) != "Length")
      {
        G4Exception("G4GDMLReadParamvol::Polyhedra_dimensionsRead()",
                    "InvalidRead", FatalException, "Invalid unit for length!");
      }
    }
    else if(attName == "aunit")
    {
      aunit = G4UnitDefinition::GetValueOf(attValue);
      if(G4UnitDefinition::GetCategory(attValue) != "Angle")
      {
        G4Exception("G4GDMLReadParamvol::Polyhedra_dimensionsRead()",
                    "InvalidRead", FatalException, "Invalid unit for angle!");
      }
    }
    else if(attName == "startPhi")
    {
      parameter.dimension[0] = eval.Evaluate(attValue);
    }
    else if(attName == "openPhi")
    {
      parameter.dimension[1] = eval.Evaluate(attValue);
    }
    else if(attName == "numRZ")
    {
      parameter.dimension[2] = eval.Evaluate(attValue);
    }
    else if(attName == "numSide")
    {
      parameter.dimension[3] = eval.Evaluate(attValue);
    }
  }

  parameter.dimension[0] *= aunit;
  parameter.dimension[1] *= aunit;

  for(auto iter = element->getFirstChild(); iter != nullptr;
           iter = iter->getNextSibling())
  {
    if(iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE)
    {
      continue;
    }

    const xercesc::DOMElement* const child =
      dynamic_cast<xercesc::DOMElement*>(iter);
    if(child == nullptr)
    {
      G4Exception("G4GDMLReadParamvo::PolyhedraRead()", "InvalidRead",
                  FatalException, "No child found!");
      return;
    }
    const G4String tag = Transcode(child->getTagName());

    if(tag == "zplane")
    {
      zplaneList.push_back(ZplaneRead(child));
    }
  }

  std::size_t numZPlanes = zplaneList.size();
  for(std::size_t i = 0; i < numZPlanes; ++i)
  {
    parameter.dimension[4 + i * 3] = zplaneList[i].rmin * lunit;
    parameter.dimension[5 + i * 3] = zplaneList[i].rmax * lunit;
    parameter.dimension[6 + i * 3] = zplaneList[i].z * lunit;
  }
}

// --------------------------------------------------------------------
void G4GDMLReadParamvol::ParametersRead(
  const xercesc::DOMElement* const element)
{
  G4ThreeVector rotation(0.0, 0.0, 0.0);
  G4ThreeVector position(0.0, 0.0, 0.0);

  G4GDMLParameterisation::PARAMETER parameter;

  for(xercesc::DOMNode* iter = element->getFirstChild(); iter != nullptr;
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
      G4Exception("G4GDMLReadParamvol::ParametersRead()", "InvalidRead",
                  FatalException, "No child found!");
      return;
    }
    const G4String tag = Transcode(child->getTagName());
    if(tag == "rotation")
    {
      VectorRead(child, rotation);
    }
    else if(tag == "position")
    {
      VectorRead(child, position);
    }
    else if(tag == "positionref")
    {
      position = GetPosition(GenerateName(RefRead(child)));
    }
    else if(tag == "rotationref")
    {
      rotation = GetRotation(GenerateName(RefRead(child)));
    }
    else if(tag == "box_dimensions")
    {
      Box_dimensionsRead(child, parameter);
    }
    else if(tag == "trd_dimensions")
    {
      Trd_dimensionsRead(child, parameter);
    }
    else if(tag == "trap_dimensions")
    {
      Trap_dimensionsRead(child, parameter);
    }
    else if(tag == "tube_dimensions")
    {
      Tube_dimensionsRead(child, parameter);
    }
    else if(tag == "cone_dimensions")
    {
      Cone_dimensionsRead(child, parameter);
    }
    else if(tag == "sphere_dimensions")
    {
      Sphere_dimensionsRead(child, parameter);
    }
    else if(tag == "orb_dimensions")
    {
      Orb_dimensionsRead(child, parameter);
    }
    else if(tag == "torus_dimensions")
    {
      Torus_dimensionsRead(child, parameter);
    }
    else if(tag == "ellipsoid_dimensions")
    {
      Ellipsoid_dimensionsRead(child, parameter);
    }
    else if(tag == "para_dimensions")
    {
      Para_dimensionsRead(child, parameter);
    }
    else if(tag == "polycone_dimensions")
    {
      Polycone_dimensionsRead(child, parameter);
    }
    else if(tag == "polyhedra_dimensions")
    {
      Polyhedra_dimensionsRead(child, parameter);
    }
    else if(tag == "hype_dimensions")
    {
      Hype_dimensionsRead(child, parameter);
    }
    else
    {
      G4String error_msg = "Unknown tag in parameters: " + tag;
      G4Exception("G4GDMLReadParamvol::ParametersRead()", "ReadError",
                  FatalException, error_msg);
    }
  }

  parameter.pRot = new G4RotationMatrix();

  parameter.pRot->rotateX(rotation.x());
  parameter.pRot->rotateY(rotation.y());
  parameter.pRot->rotateZ(rotation.z());

  parameter.position = position;

  parameterisation->AddParameter(parameter);
}

// --------------------------------------------------------------------
void G4GDMLReadParamvol::ParameterisedRead(
  const xercesc::DOMElement* const element)
{
  for(xercesc::DOMNode* iter = element->getFirstChild(); iter != nullptr;
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
      G4Exception("G4GDMLReadParamvol::ParameterisedRead()", "InvalidRead",
                  FatalException, "No child found!");
      return;
    }
    const G4String tag = Transcode(child->getTagName());

    if(tag == "parameters")
    {
      const xercesc::DOMNamedNodeMap* const attributes =
        element->getAttributes();
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
          G4Exception("G4GDMLReadParamvol::ParameterisedRead()", "InvalidRead",
                      FatalException, "No attribute found!");
          return;
        }
        const G4String attName  = Transcode(attribute->getName());
        const G4String attValue = Transcode(attribute->getValue());

        if(attName == "number")
        {
          eval.Evaluate(attValue);
        }
      }
      ParametersRead(child);
    }
    else
    {
      if(tag == "loop")
      {
        LoopRead(child, &G4GDMLRead::Paramvol_contentRead);
      }
    }
  }
}

// --------------------------------------------------------------------
void G4GDMLReadParamvol::Paramvol_contentRead(
  const xercesc::DOMElement* const element)
{
  for(xercesc::DOMNode* iter = element->getFirstChild(); iter != nullptr;
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
      G4Exception("G4GDMLReadParamvol::Paramvol_contentRead()", "InvalidRead",
                  FatalException, "No child found!");
      return;
    }
    const G4String tag = Transcode(child->getTagName());
    if(tag == "parameterised_position_size")
    {
      ParameterisedRead(child);
    }
    else if(tag == "loop")
    {
      LoopRead(child, &G4GDMLRead::Paramvol_contentRead);
    }
  }
}

// --------------------------------------------------------------------
void G4GDMLReadParamvol::ParamvolRead(const xercesc::DOMElement* const element,
                                      G4LogicalVolume* mother)
{
  G4String volumeref;

  parameterisation = new G4GDMLParameterisation();
  for(xercesc::DOMNode* iter = element->getFirstChild(); iter != nullptr;
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
      G4Exception("G4GDMLReadParamvol::ParamvolRead()", "InvalidRead",
                  FatalException, "No child found!");
      return;
    }
    const G4String tag = Transcode(child->getTagName());

    if(tag == "volumeref")
    {
      volumeref = RefRead(child);
    }
  }

  Paramvol_contentRead(element);

  G4LogicalVolume* logvol = GetVolume(GenerateName(volumeref));

  if(parameterisation->GetSize() == 0)
  {
    G4Exception("G4GDMLReadParamvol::ParamvolRead()", "ReadError",
                FatalException,
                "No parameters are defined in parameterised volume!");
  }
  G4String pv_name = logvol->GetName() + "_param";
  new G4PVParameterised(pv_name, logvol, mother, kUndefined,
                        parameterisation->GetSize(), parameterisation, check);
}
