#include "Conversions.hh"
#include <string.h>

template<class t>
ArgumentEntry<t>::ArgumentEntry(char* string,const t& def,int n) 
  : ArgumentEntryBase(string,n),Fallible<t>(def)
{
  
}

template<class t>
ArgumentEntry<t>::ArgumentEntry(char* string,const t* def,int n) 
  : ArgumentEntryBase(string,n),Fallible<t>(*def)
{
}

template<class t>
inline 
void ArgumentEntry<t>::convertChar(char* s,int n)
{
  if ( n<getNumber() || !getNumber() ) {
    t value;
    Conversion(s,value);
    validate(value);
  }
  else
      throw IndexOutOfRange(n,getNumber()-1);
}
/*
inline 
void ArgumentEntry<char*>::convertChar(char* s,int n)
{
  if ( n<getNumber() || !getNumber() ) {
    char* value = strdup(s);
    validate(value);
  }
  else
      throw IndexOutOfRange(n,getNumber()-1);
}

*/


