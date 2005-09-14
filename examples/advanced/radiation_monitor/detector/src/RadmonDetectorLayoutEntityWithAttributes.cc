//
// File name:     RadmonDetectorLayoutEntityWithAttributes.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorLayoutEntityWithAttributes.cc,v 1.3 2005-09-14 12:28:31 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonDetectorLayoutEntityWithAttributes.hh"
#include "RadmonDetectorDumpStyle.hh"
#include <iomanip>



G4String                                        RadmonDetectorLayoutEntityWithAttributes :: GetAttribute(const G4String & attributeName, const G4String & defaultValue) const
{
 AttributesMap::const_iterator i(attributesMap.find(attributeName));
 
 if (i==attributesMap.end())
  return defaultValue;
 
 return i->second;
}



G4bool                                          RadmonDetectorLayoutEntityWithAttributes :: ExistsAttribute(const G4String & attributeName) const
{
 AttributesMap::const_iterator i(attributesMap.find(attributeName));
 
 return (i!=attributesMap.end());
}



void                                            RadmonDetectorLayoutEntityWithAttributes :: SetAttribute(const G4String & attributeName, const G4String & value)
{
 attributesMap[attributeName]=value;
}



void                                            RadmonDetectorLayoutEntityWithAttributes :: ClearAttribute(const G4String & attributeName)
{
 AttributesMap::iterator i(attributesMap.find(attributeName));
 
 if (i==attributesMap.end())
  return;
 
 attributesMap.erase(i);
}



void                                            RadmonDetectorLayoutEntityWithAttributes :: ClearAllAttributes(void)
{
 attributesMap.clear();
}





void                                            RadmonDetectorLayoutEntityWithAttributes :: DumpAttributesLayout(std::ostream & out, const G4String & indent) const
{
 AttributesMap::const_iterator i(attributesMap.begin());
 AttributesMap::const_iterator end(attributesMap.end());
 
 if (i==end)
 {
  out << indent << "No attributes defined.\n";
  return;
 }

 G4int width(RADMONDETECTORDUMP_INDENT_WIDTH-1-indent.length());
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





void                                            RadmonDetectorLayoutEntityWithAttributes :: CopyFrom(const RadmonDetectorLayoutEntityWithAttributes & copy)
{
 attributesMap=copy.attributesMap;
}

