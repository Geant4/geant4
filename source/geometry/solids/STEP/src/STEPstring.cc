// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: STEPstring.cc,v 1.1 1999-01-07 16:08:14 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

/*
* NIST STEP Core Class Library
* clstepcore/STEPstring.cc
* May 1995
* KC Morris

* Development of this software was funded by the United States Government,
* and is not subject to copyright.
*/

/*  */

#include <STEPstring.h>
#ifdef WIN32
#  include <Strstrea.h>
#else
#  include <strstream.h>
#endif

void 
SdaiString::STEPwrite (ostream& out) const
{
    const char *str = 0;
// strings that exist but do not contain any chars should be written as '',
// not $ --Josh L, 4/28/95
//    if (is_null ())
    if (is_undefined ())
	out << "$";
    else
    {
	out << "\'";
	str = chars ();
	while (*str)
	{
	    if(*str == STRING_DELIM)
		out << STRING_DELIM;
	    out << *str;
	    str++;
	}
	out << "\'";
    }
}

void 
SdaiString::STEPwrite (SCLstring &s) const
{
    const char *str = 0;
// null strings should be represented by '', not $ --Josh L, 4/28/95
//    if (is_null ())
    if (is_undefined ())
    {
//	s.set_null(); // this would free up space? nope
	s = "$";
    }
    else
    {
	s = "\'";
	str = chars ();
	while (*str)
	{
	    if(*str == STRING_DELIM)
		s.Append(STRING_DELIM);
	    s.Append(*str);
	    str++;
	}
	s.Append(STRING_DELIM);
    }
}

Severity
SdaiString::StrToVal (const char * s)
{
  operator= (s);
  if (! strcmp (chars (),  s))  return SEVERITY_NULL ; 
  else return SEVERITY_INPUT_ERROR; 
}

//  STEPread reads a string in exchange file format
//  starting with a single quote
Severity 
SdaiString::STEPread (istream& in, ErrorDescriptor *err)
{
    int foundEndQuote = 0; // need so this string is not ok: 'hi''
    set_null ();  // clear the old string
    char c;
    in >> ws; // skip white space
    in >> c;

	// remember the current format state to restore the previous settings
    long int flags = in.flags();
    in.unsetf(ios::skipws);

    if (c == STRING_DELIM)
    {
	while( (c != '\0') && in.good() && in.get(c) )
	{
	    if (c == STRING_DELIM)   {
		in.get (c);
		if( ! in.good() )
		{	// it is the final quote and no extra char was read
		    foundEndQuote = 1;
		    c = '\0';
		}
		else if( !(c == STRING_DELIM) )
		{	// it is the final quote and extra char was read
		    in.putback (c); // put back non-quote extra char
		    foundEndQuote = 1;
		    c = '\0';
		}
		// else { ; } // do nothing it is an embedded quote
	    } 
	    Append (c);
	}
	Append ('\0');

	if(foundEndQuote)
	    return SEVERITY_NULL;
	else
	{    // non-recoverable error
	    err->AppendToDetailMsg("Missing closing quote on string value.\n");
	    err->AppendToUserMsg("Missing closing quote on string value.\n");
	    err->GreaterSeverity(SEVERITY_INPUT_ERROR);
	    return SEVERITY_INPUT_ERROR;
	}
    }
    //  otherwise there was not a quote
    in.putback (c);
    in.flags(flags); // set the format state back to previous settings

    set_undefined();

    return err -> GreaterSeverity (SEVERITY_INCOMPLETE);  
}

Severity 
SdaiString::STEPread (const char *s, ErrorDescriptor *err)
{
    istrstream in((char *)s);
    return STEPread (in, err);
}
