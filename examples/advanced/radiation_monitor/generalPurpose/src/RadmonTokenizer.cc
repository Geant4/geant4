//
// File name:     RadmonTokenizer.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonTokenizer.cc,v 1.1 2005-09-14 12:30:15 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonTokenizer.hh"

G4String                                        RadmonTokenizer :: operator()(const char * separators, const char * quotations)
{
 size_t i;
 char c;

 while (position<data.size())
 {
  i=0;
  c=data[position];
  
  while (separators[i]!=0)
  {
   if (c==separators[i])
    break;

   i++;
  }
  
  if (separators[i]==0)
   break;
   
  position++;
 }

 if (position>=data.size())
  return G4String();
  
 i=0;
 while (quotations[i]!=0)
 {
  if (c==quotations[i])
   break;
  
  i++;
 }

 char pattern[2];

 if (quotations[i]!=0)
 {
  pattern[0]=quotations[i];
  pattern[1]=0;
  
  separators=pattern;
  position++;
 }
 
 str_size start(position);
 
 while (position<data.size())
 {
  i=0;
  c=data[position];
  
  while (separators[i]!=0)
  {
   if (c==separators[i])
    break;

   i++;
  }
  
  if (separators[i]!=0)
   break;
  
  position++;
 }
 
 str_size stop(position);
 position++;
 
 return data(start, stop-start);
}
