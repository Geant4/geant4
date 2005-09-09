//
// File name:     RadmonDetectorLayoutEntityWithAttributes.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorLayoutEntityWithAttributes.cc,v 1.1 2005-09-09 08:26:24 capra Exp $
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





void                                            RadmonDetectorLayoutEntityWithAttributes :: DumpAttributesLayout(std::ostream & out) const
{
 AttributesMap::const_iterator i(attributesMap.begin());
 AttributesMap::const_iterator end(attributesMap.end());
 
 while (i!=end)
 {
  out << '\"' << std::setw(RADMONDETECTORDUMPWIDTH-2) << i->first << "\" = \"" << i->second << '\"' << std::endl;
  i++;
 }
}





void                                            RadmonDetectorLayoutEntityWithAttributes :: CopyFrom(const RadmonDetectorLayoutEntityWithAttributes & copy)
{
 attributesMap=copy.attributesMap;
}

