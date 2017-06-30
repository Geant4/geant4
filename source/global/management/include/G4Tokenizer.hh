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
// $Id: G4Tokenizer.hh 102049 2016-12-19 09:04:20Z gcosmo $
//
// 
//---------------------------------------------------------------
//  GEANT 4 class header file
//
//  G4Tokenizer
//
//  Class description:
//
//  String tokenizer.
//  It derives from the implementation of the Rogue Wave
//  RWTokenizer. It intrinsically uses STL string.

//---------------------------------------------------------------

#ifndef __G4Tokenizer
#define __G4Tokenizer

#include "G4String.hh"

class G4Tokenizer 
{
public:
  G4Tokenizer(const G4String& stn):string2tokenize(stn),actual(0){}

  G4String operator()(const char* str=" \t\n",size_t l=0)
    {
      size_t i,j,tmp;
      G4bool hasws=false;
      if(l==0) l=strlen(str);
      //Skip leading delimeters
      while(actual<string2tokenize.size())
	{
	  
	  for(i=0;i<l;i++)
	    if(string2tokenize[actual]==str[i]) hasws=true;
	  if(hasws)
	    {
	      actual++;
	      hasws=false;
	    }
	  else
	    break;
	}
	  
      for(j=actual;j<string2tokenize.size();j++)
	{
	  for(i=0;i<l;i++)
	    if(string2tokenize[j]==str[i]) break;
	  if(i<l) break;
	}
      if(j!=string2tokenize.size())
	{
	  tmp=actual;
	  actual=j+1;
	  return string2tokenize(tmp,j-tmp);
	}
      else
	{
	  tmp=actual;
	  actual=j;
	  return string2tokenize(tmp,j-tmp);
	}
    } 

private:

  G4String string2tokenize;
  size_t actual;

};

#endif

