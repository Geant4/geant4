#include <stdlib.h>
#include <string.h>
#include "Conversions.hh"
#include "Boolean.h"
#include "heapwr.hh"
#include "String.hh"

String& Conversion(char* s,String& t) 
{
  return t = s;
}

int& Conversion(char* string,int& n)
{
   return ( n = atoi(string) );
}

long& Conversion(char* string,long& n)
{
   return ( n = atol(string) );
}

double& Conversion(char* string,double& n)
{
   return ( n = atof(string) );
}

Boolean& Conversion(char* string,Boolean& n)
{
  if (string)
    if (strcmp(string,"true") )
      n = Boolean::True;
    else
      if (strcmp(string,"false") )
	n = Boolean::False;
      else
	throw ConversionImpossible("Boolean",string);
  else
    n = Boolean::True;
  return n;
}

char*& Conversion(char* string,char*& s)
{
  s = strdup(string);
  return s;
}

ConversionImpossible::ConversionImpossible(char* t,char* s)
{
   strcpy(Type = NEW char[strlen(t)+1],t);
   strcpy(String = NEW char[strlen(s)+1],s);
}

ConversionImpossible::~ConversionImpossible()
{
   delete [] Type;
   delete [] String;
}

void ConversionImpossible::writeMessage(G4std::ostream& o) const
{
   o << "Cannot convert '" << String << "' to type " << Type << G4endl;
}

