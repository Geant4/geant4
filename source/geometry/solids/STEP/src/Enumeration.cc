

//



//
// $Id: Enumeration.cc,v 1.2 1999-05-21 20:20:47 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include <Enumeration.h>

/*
* NIST STEP Core Class Library
* clstepcore/Enumeration.cc
* May 1995
* K. C. Morris
* David Sauder

* Development of this software was funded by the United States Government,
* and is not subject to copyright.
*/

/*   */

//static char rcsid[] ="";

#ifdef WIN32
#  include <Strstrea.h>
#else
#  include <strstream.h>
#endif

// Josh L, 3/31/95
// These constants have to be initialized outside the SDAI struct.  They
// are initialized here instead of in the header file in order to avoid
// mulitple inclusions when building SCL applications.
const Logical SDAI::UNKNOWN = sdaiUNKNOWN;
const Boolean SDAI::SdaiTRUE = sdaiTRUE;
const Boolean SDAI::SdaiFALSE = sdaiFALSE;

///////////////////////////////////////////////////////////////////////////////
// class Logical
///////////////////////////////////////////////////////////////////////////////

Logical::operator LOGICAL () const  {
  switch (v) {
  case sdaiUNKNOWN:  return sdaiUNKNOWN;
  case sdaiFALSE: return sdaiFALSE;
  case sdaiTRUE: return sdaiTRUE;
  default: return sdaiUNKNOWN;
}}

const char * 
Logical::element_at (int n)  const {
  switch (n)  {
  case sdaiUNKNOWN:  return "U";
  case sdaiFALSE: return "F";
  case sdaiTRUE: return "T";
  default: return "";
  }
}

///////////////////////////////////////////////////////////////////////////////
// class Boolean  29-Sep-1994
///////////////////////////////////////////////////////////////////////////////

Boolean::Boolean (LOGICAL val)  {
  if (val == sdaiUNKNOWN) return;
  set_value (val);
}

Boolean::operator LOGICAL () const  {
  switch (v) {
  case sdaiFALSE: return sdaiFALSE;
  case sdaiTRUE: return sdaiTRUE;
  default: return sdaiUNKNOWN;
}}

const char * 
Boolean::element_at (int n)  const {
  switch (n)  {
  case sdaiFALSE: return "F";
  case sdaiTRUE: return "T";
  default: return "";
  }
}

Boolean::operator int () const  {
  if (v == sdaiFALSE)  return 0;
  else return 1;
}

int Boolean::operator ==( LOGICAL log ) const
{
  if( v == log )
      return 1;
  return 0;
}

///////////////////////////////////////////////////////////////////////////////

/******************************************************************
 ** Procedure:  DebugDisplay
 ** Parameters:  ostream& out
 ** Returns:  
 ** Description:  prints out some information on the enumerated 
 **               item for debugging purposes
 ** Side Effects:  
 ** Status:  ok 2/1/91
 ******************************************************************/
void
STEPenumeration::DebugDisplay (ostream& out) const  {
    SCLstring tmp;
    out << "Value  :  " <<  v  << "  " 
	<< asStr (tmp) 
	    << ": for type " << Name () << "\n";
    out << "valid values are  ";
    int i =0;
    while (element_at (i))  
      {
	  out << element_at (i++) << "  ";
      }
    STEPwrite (out);
    out << "\n";
	
}	

// Read an Enumeration value 
// ENUMERATION = "." UPPER { UPPER | DIGIT } "."
// *note* UPPER is defined as alpha or underscore.
// returns: Severity of the error.
// error message and error Severity is written to ErrorDescriptor *err.
// int AssignVal is: 
// true => value is assigned to the STEPenumeration; 
// true or false => value is read and appropriate error info is set and 
// 	returned.
// int needDelims is: 
// false => absence of the period delimiters is not an error; 
// true => delimiters must be valid; 
// true or false => non-matching delimiters are flagged as an error

Severity 
STEPenumeration::ReadEnum(istream& in, ErrorDescriptor *err, int AssignVal,
			  int needDelims)
{
    if(AssignVal)
	set_null();

    SCLstring str;
    char c;
    char messageBuf[512];
    messageBuf[0] = '\0';

    int validDelimiters = 1;

    in >> ws; // skip white space

    if( in.good() )
    {
	in.get(c);
	if( c == '.' || isalpha(c))
	{
	    if( c == '.' )
	    {
		in.get(c); // push past the delimiter
		// since found a valid delimiter it is now invalid until the 
		//   matching ending delim is found
		validDelimiters = 0;
	    }

	    // look for UPPER
	    if( in.good() && ( isalpha(c) || c == '_' ) )
	    {
		str.Append(c);
		in.get(c);
	    }

	    // look for UPPER or DIGIT
	    while( in.good() && ( isalnum(c) || c == '_' ) )
	    {
		str.Append(c);
		in.get(c);
	    }
	    // if character is not the delimiter unread it
	    if(in.good() && (c != '.') )
		in.putback(c);

	    // a value was read
	    if(str.Length() > 0)
	    {
		int i =0;
		const char *strval = str.chars();
		SCLstring tmp;
		while( (i < no_elements ())  &&  
		       (strcmp( (char *)StrToUpper( strval, tmp ), 
					     element_at (i) ) != 0) )
		    ++i;
		if ( no_elements () == i)
		{	//  exhausted all the possible values 
		    err->GreaterSeverity(SEVERITY_WARNING);
		    err->AppendToDetailMsg("Invalid Enumeration value.\n");
		    err->AppendToUserMsg("Invalid Enumeration value.\n");
		}
		else
		{
		    if(AssignVal)
			v = i;
		}

		// now also check the delimiter situation
		if(c == '.') // if found ending delimiter
		{
		    // if expecting delim (i.e. validDelimiter == 0)
		    if(!validDelimiters) 
		    {
			validDelimiters = 1; // everything is fine
		    }
		    else // found ending delimiter but no initial delimiter
		    {
			validDelimiters = 0;
		    }
		}
		// didn't find any delimiters at all and need them.
		else if(needDelims) 
		{
		    validDelimiters = 0;
		}

		if (!validDelimiters)
		{	
		    err->GreaterSeverity(SEVERITY_WARNING);
		    if(needDelims)
			sprintf(messageBuf, 
			  "Enumerated value has invalid period delimiters.\n");
		    else
			sprintf(messageBuf, 
			   "Mismatched period delimiters for enumeration.\n");
		    err->AppendToDetailMsg(messageBuf);
		    err->AppendToUserMsg(messageBuf);
		}
		return err->severity();
	    }
	    // found valid or invalid delimiters with no associated value 
	    else if( (c == '.') || !validDelimiters)
	    {
		err->GreaterSeverity(SEVERITY_WARNING);
		err->AppendToDetailMsg(
	   "Enumerated has valid or invalid period delimiters with no value.\n"
				      );
		err->AppendToUserMsg(
	   "Enumerated has valid or invalid period delimiters with no value.\n"
				      );
		return err->severity();
	    }
	    else // no delims and no value
		err->GreaterSeverity(SEVERITY_INCOMPLETE);

	}
	else if( (c == ',') || (c == ')') )
	{
	    in.putback(c);
	    err->GreaterSeverity(SEVERITY_INCOMPLETE);
	}
	else
	{
	    in.putback(c);
	    err->GreaterSeverity(SEVERITY_WARNING);
	    sprintf(messageBuf, "Invalid enumeration value.\n");
	    err->AppendToDetailMsg(messageBuf);
	    err->AppendToUserMsg(messageBuf);
	}
    }
    else // Hit eof (assuming there was no error state for istream passed in)
    {
	err->GreaterSeverity(SEVERITY_INCOMPLETE);
    }
    return err->severity();
}

/*
Severity 
STEPenumeration::StrToVal (const char * s)
{
    put (s);
    return SEVERITY_NULL;
}
*/

Severity 
STEPenumeration::StrToVal (const char * s, ErrorDescriptor *err, int optional)
{
    istrstream in ((char *)s); // sz defaults to length of s

    Severity sev = ReadEnum(in, err, 1, 0);
    if( (err->severity() == SEVERITY_INCOMPLETE) && optional)
	err->severity(SEVERITY_NULL);

    return err->severity();
}

// reads an enumerated value in STEP file format 
Severity
STEPenumeration::STEPread (const char *s, ErrorDescriptor *err, int optional)
{
    istrstream in((char *)s);
    return STEPread (in, err, optional);
}

// reads an enumerated value in STEP file format 
Severity
STEPenumeration::STEPread (istream& in, ErrorDescriptor *err, int optional)
{
    Severity sev = ReadEnum(in, err, 1, 1);
    if( (err->severity() == SEVERITY_INCOMPLETE) && optional)
	err->severity(SEVERITY_NULL);

    return err->severity();
}


const char * 
STEPenumeration::asStr (SCLstring &s) const  {
    if (v != ENUM_NULL) 
    {
//	s = elements[v];
	return s = element_at (v);
//	return s.chars();
    }
    else return "";
}

void 
STEPenumeration::STEPwrite (ostream& out)  const  {
    if( is_null() )
	out << '$';
    else
    {
	SCLstring tmp;
	out << "." <<  asStr (tmp) << ".";
    }
}

const char * 
STEPenumeration::STEPwrite (SCLstring &s) const
{
    if( is_null() )
    {
	s.set_null();
    }
    else
    {
	SCLstring tmp;
	s = ".";
	s.Append(asStr(tmp));
	s.Append('.');
    }
    return s.chars();
}

//STEPenumeration::STEPenumeration (const char * const e)
//:  elements (e)
//{  
//}

/******************************************************************
 ** Procedure:  set_elements
 ** Parameters:  
 ** Returns:  
 ** Description:  
 ** Side Effects:  
 ** Status:  
 ******************************************************************/
#ifdef OBSOLETE
void
STEPenumeration::set_elements (const char * const e [])  {
    elements = e;
}
#endif
Severity 
STEPenumeration::EnumValidLevel(istream &in, ErrorDescriptor *err,
				int optional, char *tokenList, 
				int needDelims, int clearError)
{
    if(clearError)
	err->ClearErrorMsg();

    in >> ws; // skip white space
    char c = ' '; 
    c = in.peek();
    if(c == '$' || in.eof())
    {
	if(!optional)
	    err->GreaterSeverity(SEVERITY_INCOMPLETE);
	if(in)
	    in >> c;
	CheckRemainingInput(in, err, "enumeration", tokenList);
	return err->severity();
    }
    else
    {
	ErrorDescriptor error;

	ReadEnum(in, &error, 0, needDelims);
	CheckRemainingInput(in, &error, "enumeration", tokenList);

	Severity sev = error.severity();
	if(sev < SEVERITY_INCOMPLETE)
	{
	    err->AppendToDetailMsg(error.DetailMsg());
	    err->AppendToUserMsg(error.UserMsg());
	    err->GreaterSeverity(error.severity());
	}
	else if(sev == SEVERITY_INCOMPLETE && !optional)
	    err->GreaterSeverity(SEVERITY_INCOMPLETE);
    }
    return err->severity();
}

Severity 
STEPenumeration::EnumValidLevel(const char *value, ErrorDescriptor *err,
				int optional, char *tokenList, 
				int needDelims, int clearError)
{
    istrstream in((char *)value);
    return EnumValidLevel (in, err, optional, tokenList, needDelims,
			   clearError);
/*

    char messageBuf[BUFSIZ];
    messageBuf[0] = '\0';

    if(attrValue)
    {
	int len = strlen (attrValue);
	char *valstart = new char [len + 1];
	char *val = valstart;

	int numFound = sscanf(attrValue," %s", val);
	if(numFound != EOF)
	{
	    int i = 0;
	    if(val [0] == '.')  // strip the delims
	    {

		val++;
		char * pos = strchr(val, '.');
		if (pos) 
		    *pos = '\0';
		else
		{
		    err->AppendToDetailMsg(
		    "Missing ending period delimiter for enumerated value.\n");
		    err->AppendToUserMsg(
		    "Missing ending period delimiter for enumerated value.\n");
		    err->GreaterSeverity(SEVERITY_WARNING);
		}
	    }

	    SCLstring tmp;
	    while((i < no_elements() ) && 
	    (strcmp( (char *)StrToUpper(val, tmp), element_at (i) ) != 0))
		++i;
	    if(no_elements() == i)	// exhausted all the possible values 
	    {
		err->GreaterSeverity(SEVERITY_WARNING);
		sprintf(messageBuf, 
			"attribute %s: Invalid enumeration value: '%s'",
			Name(), val);
		err->AppendToUserMsg(messageBuf);
		err->AppendToDetailMsg(messageBuf);
//		DebugDisplay ();
		return SEVERITY_WARNING;
	    }
	    err->GreaterSeverity(SEVERITY_NULL);
	    return SEVERITY_NULL;
	}
	delete [] valstart;
    }
    if(optional) 
    {
	err->GreaterSeverity(SEVERITY_NULL);
	return SEVERITY_NULL;
    }
    else
    {
	err->GreaterSeverity(SEVERITY_INCOMPLETE);
	return SEVERITY_INCOMPLETE;
    }
*/
}

/******************************************************************
 ** Procedure:  set_value
 ** Parameters:  char * n  OR  in i  -- value to be set
 ** Returns:  value set 
 ** Description:  sets the value of an enumerated attribute
 **     case is not important in the character based version
 **     if value is not acceptable, a warning is printed and 
 **     processing continues
 ** Side Effects:  
 ** Status:  ok 2.91
 ******************************************************************/
int
STEPenumeration::set_value (const char * n)  {  
    //  assigns the appropriate value based on n
    if  ( !n || (!strcmp (n, "")) )  return set_value (ENUM_NULL);
	
    int i =0;
    SCLstring tmp;
    while ((i < no_elements ())  &&  
	   (strcmp ( (char *)StrToUpper( n, tmp ),  element_at (i)) != 0 ) )
	++i;
    if ( no_elements () == i)  {  //  exhausted all the possible values 
	return set_value (ENUM_NULL);
    }
    v = i;	
    return v;
    
}

//  set_value is the same function as put
int
STEPenumeration::set_value (const int i)  {  
    if (i == ENUM_NULL)  {
	v = ENUM_NULL;
	return ENUM_NULL;
    }
    const char *tmp = element_at( i );
    if ( tmp[0] != '\0' )  return (v =i);
    // otherwise 
    G4cerr << "(OLD Warning:) invalid enumeration value " << i
	<< " for " <<  Name () << "\n";
    DebugDisplay ();
    return  ENUM_NULL ;
    
}

STEPenumeration& 
STEPenumeration::operator= (const int i)
{
    put (i);
    return *this;
}

STEPenumeration&
STEPenumeration::operator= (const STEPenumeration& Senum)
{
    put (Senum.asInt());
    return *this;
}

ostream &operator<< ( ostream& out, const STEPenumeration& a )
{
    SCLstring tmp;
    out << a.asStr( tmp );
    return out;

}


#ifdef pojoldStrToValNstepRead

Severity 
STEPenumeration::StrToVal (const char * s, ErrorDescriptor *err, int optional)
{
    const char *sPtr = s;
    while(isspace(*sPtr)) sPtr++;
    if(*sPtr == '\0')
    {
	if(optional) 
	{
	    err->GreaterSeverity(SEVERITY_NULL);
	    return SEVERITY_NULL;
	}
	else
	{
	    err->GreaterSeverity(SEVERITY_INCOMPLETE);
	    return SEVERITY_INCOMPLETE;
	}
    }
    else if(*sPtr == '.') // look for initial period delimiter
    {
	return STEPread(sPtr, err);
    }
    else
    {
		// look for ending period delimiter (an error)
	char *periodLoc = strchr(sPtr, '.');
	if (periodLoc)
	{	// found an ending period w/out initial period
	    char *tmp = new char[strlen(sPtr) + 1];
	    strcpy(tmp, sPtr);
	    tmp[periodLoc - sPtr] = '\0'; // write over ending period
	    err->GreaterSeverity(SEVERITY_WARNING);
	    err->AppendToDetailMsg(
		"Ending period delimiter without initial period delimiter.\n");
	    err->AppendToUserMsg(
		"Ending period delimiter without initial period delimiter.\n");
	    delete [] tmp;
	    if( ValidLevel(sPtr, err, optional) )
	    { // remaining value is valid so assign it
		put(tmp);
		return SEVERITY_WARNING;
	    }
	    else
	    {
		err->AppendToDetailMsg("Invalid Enumerated value.\n");
		err->AppendToUserMsg("Invalid Enumerated value.\n");
		return SEVERITY_WARNING;
	    }
	}
	// no surrounding delimiters
	else if( ValidLevel(sPtr, err, optional) )
	{ // value is valid so assign it
	    put (sPtr);  
	    return SEVERITY_NULL;
	}
	else
	{
	    err->AppendToDetailMsg("Invalid Enumerated value.\n");
	    err->AppendToUserMsg("Invalid Enumerated value.\n");
	    return SEVERITY_WARNING;
	}
    }
}


Severity
STEPenumeration::STEPread (istream& in, ErrorDescriptor *err, int optional)
{
    char enumValue [BUFSIZ];
    char c;
    char errStr[BUFSIZ];
    errStr[0] = '\0';

    err->severity(SEVERITY_NULL); // assume ok until error happens
    in >> c;
    switch (c)  
      {
	case '.':
	  in.getline (enumValue, BUFSIZ, '.');// reads but does not Store the .
/*
  // gcc 2.3.3 - It does and should read the . It doesn't Store it DAS 4/27/93
	  char * pos = index(enumValue, '.');
	  if (pos) *pos = '\0';
	  //  NON-STANDARD (GNUism)  getline should not Retrieve .
	  //  function gcount is unavailable
*/
	  if(in.fail())
	  {
	      err->GreaterSeverity(SEVERITY_WARNING);
	      err->AppendToUserMsg(
		    "Missing ending period delimiter for enumerated value.\n");
	      err->AppendToDetailMsg(
		    "Missing ending period delimiter for enumerated value.\n");
	  }
	  if(ValidLevel(enumValue, err, optional) == SEVERITY_NULL)
	      set_value (enumValue);
	  else
	  {
	      err->AppendToDetailMsg("Invalid enumerated value.\n");
	      err->GreaterSeverity(SEVERITY_WARNING);
	      set_value(ENUM_NULL);
	  }	      
	  break;
	  
	case ',':	// for next attribute or next aggregate value?
	case ')':	// for end of aggregate value?
	default:
	  in.putback (c);
	  set_value (ENUM_NULL);
	  if(optional) err->GreaterSeverity(SEVERITY_NULL);
	  else	       err->GreaterSeverity(SEVERITY_INCOMPLETE);
	  break;
	  
/*
	default:
	  set_value (ENUM_NULL);
		// read didn't know what to do
	  err->GreaterSeverity(SEVERITY_INPUT_ERROR); 
	  sprintf(errStr,
	       "STEPenumeration::STEPread(): warning : poorly delimited %s %s",
	        Name(), "enumerated value was ignored.");
	  err->AppendToDetailMsg(errStr);
*/
      }  
    return err->severity();
}

#endif
