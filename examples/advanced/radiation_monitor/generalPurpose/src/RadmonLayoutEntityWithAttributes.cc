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
// File name:     RadmonLayoutEntityWithAttributes.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonLayoutEntityWithAttributes.cc,v 1.5 2006/06/29 16:14:45 gunter Exp $
// Tag:           $Name: geant4-08-01 $
//

// Include files
#include "RadmonLayoutEntityWithAttributes.hh"
#include "RadmonDumpStyle.hh"
#include "RadmonMessenger.hh"
#include "G4UIcommand.hh"
#include <iomanip>



G4int                                           RadmonLayoutEntityWithAttributes :: GetNAttributes(void) const
{
 return attributesVector.size();
}



const G4String &                                RadmonLayoutEntityWithAttributes :: GetAttributeName(G4int index) const
{
 return attributesVector[index].first;
}


 
G4String                                        RadmonLayoutEntityWithAttributes :: GetAttribute(const G4String & attributeName, const G4String & defaultValue) const
{
 AttributesVector::const_iterator i(attributesVector.begin());
 AttributesVector::const_iterator end(attributesVector.end());
 
 while (i!=end)
 {
  if (i->first==attributeName)
   return i->second;
 
  i++;
 }
 
 return defaultValue;
}



G4bool                                          RadmonLayoutEntityWithAttributes :: ExistsAttribute(const G4String & attributeName) const
{
 AttributesVector::const_iterator i(attributesVector.begin());
 AttributesVector::const_iterator end(attributesVector.end());
 
 while (i!=end)
 {
  if (i->first==attributeName)
   return true;
 
  i++;
 }
 
 return false;
}



void                                            RadmonLayoutEntityWithAttributes :: SetAttribute(const G4String & attributeName, const G4String & value)
{
 AttributesVector::iterator i(attributesVector.begin());
 AttributesVector::iterator end(attributesVector.end());
 
 while (i!=end)
 {
  if (i->first==attributeName)
  {
   i->second=value;
   return;
  }
 
  i++;
 }

 size_t n(attributesVector.size());
 attributesVector.resize(n+1);
 attributesVector[n].first=attributeName;
 attributesVector[n].second=value;
}



void                                            RadmonLayoutEntityWithAttributes :: ClearAttribute(const G4String & attributeName)
{
 AttributesVector::iterator i(attributesVector.begin());
 AttributesVector::iterator end(attributesVector.end());
 
 while (i!=end)
 {
  if (i->first==attributeName)
  {
   attributesVector.erase(i);
   return;
  }
 
  i++;
 }
}



void                                            RadmonLayoutEntityWithAttributes :: ClearAllAttributes(void)
{
 attributesVector.clear();
}





G4double                                        RadmonLayoutEntityWithAttributes :: GetAttributeAsDouble(const G4String & attributeName, double defaultValue) const
{
 G4String str;
 
 str=GetAttribute(attributeName, "#");
 if (str=="#")
  return defaultValue;
  
 G4String args[1];
 if (!RadmonMessenger::ProcessArguments(str, 1, args))
  return defaultValue;
  
 return G4UIcommand::ConvertToDouble(args[0]);
}



G4double                                        RadmonLayoutEntityWithAttributes :: GetAttributeAsMeasure(const G4String & attributeName, const char * category, double defaultValue) const
{
 G4String str;
 
 str=GetAttribute(attributeName, "#");
 if (str=="#")
  return defaultValue;
 
 G4String args[2];
 if (!RadmonMessenger::ProcessArguments(str, 2, args))
  return defaultValue;

 G4double unit(RadmonMessenger::GetUnit(args[1], category));
 if (unit<=0.)
  return defaultValue;
  
 return G4UIcommand::ConvertToDouble(args[0])*unit;
}



G4int                                           RadmonLayoutEntityWithAttributes :: GetAttributeAsInteger(const G4String & attributeName, G4int defaultValue) const
{
 G4String str;
 
 str=GetAttribute(attributeName, "#");
 if (str=="#")
  return defaultValue;
  
 G4String args[1];
 if (!RadmonMessenger::ProcessArguments(str, 1, args))
  return defaultValue;
  
 return G4UIcommand::ConvertToInt(args[0]);
}





G4ThreeVector                                   RadmonLayoutEntityWithAttributes :: GetAttributeAsThreeVector(const G4String & attributeName, const G4ThreeVector & defaultValue) const
{
 G4String str;
 
 str=GetAttribute(attributeName, "#");
 if (str=="#")
  return defaultValue;
  
 G4String args[3];
 if (!RadmonMessenger::ProcessArguments(str, 3, args))
  return defaultValue;

 return G4ThreeVector(G4UIcommand::ConvertToDouble(args[0]), G4UIcommand::ConvertToDouble(args[1]), G4UIcommand::ConvertToDouble(args[2])); 
}





G4ThreeVector                                   RadmonLayoutEntityWithAttributes :: GetAttributeAsThreeVectorWithMeasure(const G4String & attributeName, const char * category, const G4ThreeVector & defaultValue) const
{
 G4String str;
 
 str=GetAttribute(attributeName, "#");
 if (str=="#")
  return defaultValue;
  
 G4String args[4];
 if (!RadmonMessenger::ProcessArguments(str, 4, args))
  return defaultValue;

 G4double unit(RadmonMessenger::GetUnit(args[3], category));
 if (unit<=0.)
  return defaultValue;
  
 return G4ThreeVector(G4UIcommand::ConvertToDouble(args[0])*unit, G4UIcommand::ConvertToDouble(args[1])*unit, G4UIcommand::ConvertToDouble(args[2])*unit);
}





G4ThreeVector                                   RadmonLayoutEntityWithAttributes :: GetAttributeAsDirection(const G4String & attributeName, const G4ThreeVector & defaultValue) const
{
 G4String str;
 
 str=GetAttribute(attributeName, "#");
 if (str=="#")
  return defaultValue;
  
 G4String args[3];
 if (!RadmonMessenger::ProcessArguments(str, 3, args))
  return defaultValue;

 G4double theta(RadmonMessenger::GetUnit(args[2], "Angle"));
 if (theta<=0.)
  return defaultValue;

 G4double phi(theta*G4UIcommand::ConvertToDouble(args[1]));
 theta*=G4UIcommand::ConvertToDouble(args[0]);
  
 G4ThreeVector axis;
 axis.setRThetaPhi(1., theta/rad, phi/rad);
 
 return axis;
}



G4RotationMatrix                                RadmonLayoutEntityWithAttributes :: GetAttributeAsRotationMatrix(const G4String & attributeName, const G4RotationMatrix & defaultValue) const
{
 G4String str;
 
 str=GetAttribute(attributeName, "#");
 if (str=="#")
  return defaultValue;
  
 G4String args[4];
 if (!RadmonMessenger::ProcessArguments(str, 4, args))
  return defaultValue;

 G4double theta(RadmonMessenger::GetUnit(args[3], "Angle"));
 if (theta<=0.)
  return defaultValue;

 G4double phi(theta*G4UIcommand::ConvertToDouble(args[1]));
 G4double delta(theta*G4UIcommand::ConvertToDouble(args[2]));
 theta*=G4UIcommand::ConvertToDouble(args[0]);
  
 G4ThreeVector axis;
 axis.setRThetaPhi(1., theta/rad, phi/rad);
 
 return G4RotationMatrix(axis, delta/rad);
}





void                                            RadmonLayoutEntityWithAttributes :: DumpAttributesLayout(std::ostream & out, const G4String & indent) const
{
 AttributesVector::const_iterator i(attributesVector.begin());
 AttributesVector::const_iterator end(attributesVector.end());
 
 if (i==end)
 {
  out << indent << "No attributes defined.\n";
  return;
 }

 G4int width(RADMONDUMP_INDENT_WIDTH-1-indent.length());
 if (width<0)
  width=0;
  
 G4int width2;
 
 while (i!=end)
 {
  width2=width-i->first.length();
  if (width2<0)
   width2=0;
   
  out << indent << '\"' << i->first << std::setw(width2); 
  out.setf(std::ostream::left, std::ostream::adjustfield); 
  out << "\""; 
  out.setf(std::ostream::right, std::ostream::adjustfield); 
  out << " = \"" << i->second << "\"\n";
  
  i++;
 }
}





void                                            RadmonLayoutEntityWithAttributes :: CopyFrom(const RadmonLayoutEntityWithAttributes & copy)
{
 attributesVector=copy.attributesVector;
}

