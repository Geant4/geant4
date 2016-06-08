// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ctoken.h,v 1.5 1999/11/29 10:17:41 gcosmo Exp $
// GEANT4 tag $Name: geant4-01-01 $
//
// 
//---------------------------------------------------------------
//  GEANT 4 class header file
//
//  G4Tokenizer
//
//  Class description:
//
//  STL wrapper class for String tokenizer.
//  It implements Rogue Wave RWTokenizer signature but
//  intrinsically using STL string.

//---------------------------------------------------------------

#ifndef __ctoken
#define __ctoken

#include "g4rw/defs.h"
#include "g4rw/cstring.h"

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

