#ifndef __ctoken
#define __ctoken

#include <string>
#include "rw/defs.h"
#include "rw/cstring.h"

class RWCTokenizer 
{
public:
  RWCTokenizer(const RWCString& s):string2tokenize(s),actual(0){}

  RWCSubString operator()(const char* str=" \t\n",size_t l=0)
    {
      size_t i,j,tmp;
      RWBoolean hasws=false;
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

  RWCString string2tokenize;
  size_t actual;

};

#endif

