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
// $Id: G4Tokenizer.hh,v 1.1 2001-10-11 14:04:05 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
  G4Tokenizer(const G4String& s):string2tokenize(s),actual(0){}

  G4SubString operator()(const char* str=" \t\n",size_t l=0)
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

