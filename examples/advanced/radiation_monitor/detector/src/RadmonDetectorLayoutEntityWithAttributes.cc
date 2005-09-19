//
// File name:     RadmonDetectorLayoutEntityWithAttributes.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorLayoutEntityWithAttributes.cc,v 1.4 2005-09-19 19:42:13 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonDetectorLayoutEntityWithAttributes.hh"
#include "RadmonDetectorDumpStyle.hh"
#include <iomanip>



G4int                                           RadmonDetectorLayoutEntityWithAttributes :: GetNAttributes(void) const
{
 return attributesVector.size();
}



const G4String &                                RadmonDetectorLayoutEntityWithAttributes :: GetAttributeName(G4int index) const
{
 return attributesVector[index].first;
}


 
G4String                                        RadmonDetectorLayoutEntityWithAttributes :: GetAttribute(const G4String & attributeName, const G4String & defaultValue) const
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



G4bool                                          RadmonDetectorLayoutEntityWithAttributes :: ExistsAttribute(const G4String & attributeName) const
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



void                                            RadmonDetectorLayoutEntityWithAttributes :: SetAttribute(const G4String & attributeName, const G4String & value)
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



void                                            RadmonDetectorLayoutEntityWithAttributes :: ClearAttribute(const G4String & attributeName)
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



void                                            RadmonDetectorLayoutEntityWithAttributes :: ClearAllAttributes(void)
{
 attributesVector.clear();
}





void                                            RadmonDetectorLayoutEntityWithAttributes :: DumpAttributesLayout(std::ostream & out, const G4String & indent) const
{
 AttributesVector::const_iterator i(attributesVector.begin());
 AttributesVector::const_iterator end(attributesVector.end());
 
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
 attributesVector=copy.attributesVector;
}

