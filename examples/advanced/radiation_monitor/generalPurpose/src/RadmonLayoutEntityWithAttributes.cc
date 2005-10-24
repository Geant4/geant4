//
// File name:     RadmonLayoutEntityWithAttributes.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonLayoutEntityWithAttributes.cc,v 1.1 2005-10-24 14:51:36 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonLayoutEntityWithAttributes.hh"
#include "RadmonDumpStyle.hh"
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

