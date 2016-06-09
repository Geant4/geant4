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
// File name:     RadmonTokenizer.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonTokenizer.cc,v 1.3 2006/06/29 16:14:49 gunter Exp $
// Tag:           $Name: geant4-08-01 $
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
