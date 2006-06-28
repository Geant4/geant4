//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// File name:     RadmonTokenizer.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonTokenizer.cc,v 1.2 2006-06-28 13:52:38 gunter Exp $
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
