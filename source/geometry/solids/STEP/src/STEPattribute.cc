// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: STEPattribute.cc,v 1.1 1999-01-07 16:08:12 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
/*
* NIST STEP Core Class Library
* clstepcore/STEPattribute.cc
* May 1995
* K. C. Morris
* David Sauder

* Development of this software was funded by the United States Government,
* and is not subject to copyright.
*/

/*  */

//static char rcsid[] ="";


#include <read_func.h>
#include <STEPattribute.h>
#include <sdai.h>
#include <instmgr.h>
#include <STEPundefined.h>
#include <STEPaggregate.h>
#include <STEPentity.h>
#include <STEPselect.h>
#include <SdaiBinary.h>

// REAL_NUM_PRECISION is defined in STEPattribute.h, and is also used
// in aggregate real handling (STEPaggregate.cc)  -- IMS 6 Jun 95
const int Real_Num_Precision = REAL_NUM_PRECISION;

/******************************************************************
**	  Functions for manipulating attribute

**  KNOWN BUGS:  
**    -- error reporting does not include line number information
**    -- null attributes are only handled through the attribute pointer class
**       direct access of a null attribute will not indicate whether the 
**       attribute is null in all cases (although a null value is available
**       for entities and enumerations.) 
**/


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// STEPattribute Functions
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/******************************************************************
// the value of the attribute is assigned from the supplied string
 ******************************************************************/

Severity 
STEPattribute::StrToVal (const char* s, InstMgr * instances, int addFileId)
{
    _error.ClearErrorMsg();	// also sets Severity to SEVERITY_NULL

    //  set the value to be null (reinitialize the attribute value)
    set_null();

    int nullable = (aDesc->Optional().asInt() == sdaiTRUE);

    // an empty str gets assigned NULL
    if (!s )
    {
	if( nullable || IsDerived() )	// if it is derived it doesn't 
	    return SEVERITY_NULL;	// matter if it is null DAS
	else {
	  _error.severity(SEVERITY_INCOMPLETE);
	  return SEVERITY_INCOMPLETE;
	}
    }
    if (s[0] == '\0')
    {
	if(NonRefType() == STRING_TYPE)
	{ // this is interpreted as a string with no value i.e. "".
	    *(ptr.S) = s; // using string class - don't need to declare space
	    return SEVERITY_NULL;
	}
	if( nullable || IsDerived() )	// if it is derived it doesn't 
	    return SEVERITY_NULL;	// matter if it is null DAS
	else  {
	  _error.severity(SEVERITY_INCOMPLETE);
	  return SEVERITY_INCOMPLETE;
	}
    }

    //  an overridden attribute always has a \'*\' value
    if( IsDerived() )	// check to see if value contains: optional space,
    {			// followed by *, followed by optional space
	const char * tmpSptr = s;
	while( isspace( *tmpSptr ) ) tmpSptr++;
	if(*tmpSptr == '*')
	{
	    tmpSptr++;
	    char tmpC;
	    int charsFound = sscanf(tmpSptr, "%c", &tmpC);
	    if(charsFound == EOF) // no non-white chars followed the *
		return SEVERITY_NULL;
	}
	_error.AppendToDetailMsg (
			 "Derived attribute must have \'*\' for its value.\n");
	return _error.severity (SEVERITY_INPUT_ERROR);
    }

    istrstream in ((char *)s); // sz defaults to length of s
    int valAssigned = 0;

    // read in value for attribute
    switch (NonRefType())
    {
      case INTEGER_TYPE:
      {
	  valAssigned = ReadInteger(*(ptr.i), s, &_error, 0);
	  break;
      }
      case REAL_TYPE:
      {
	  valAssigned = ReadReal(*(ptr.r), s, &_error, 0);
	  break;
      }
      case NUMBER_TYPE:
      {
	  valAssigned = ReadNumber(*(ptr.r), s, &_error, 0);
	  break;
      }

      case ENTITY_TYPE:
      {
	  STEPentity *se = ReadEntityRef(s, &_error, 0, instances, addFileId);
	  if( se != S_ENTITY_NULL )
	  {
	      if(EntityValidLevel(se, aDesc->NonRefTypeDescriptor(), 
				  &_error) 
		 == SEVERITY_NULL)
		  *(ptr.c) = se;
	      else
		  *(ptr.c) = S_ENTITY_NULL;
	  }
	  else
	      *(ptr.c) = S_ENTITY_NULL;
	  break;
      }

      case BINARY_TYPE:
      {
	  ptr.b->StrToVal(s, &_error); // call class SdaiBinary::StrToVal()
	  break;
      }
      case STRING_TYPE:
      {
	  *(ptr.S) = s; // using string class - don't need to declare space
	  break;
      }

      case BOOLEAN_TYPE:
      case LOGICAL_TYPE:
      case ENUM_TYPE:
	{
	    ptr.e->StrToVal (s, &_error, nullable);
	    break;
	}

      case AGGREGATE_TYPE:
      case ARRAY_TYPE:		// DAS
      case BAG_TYPE:		// DAS
      case SET_TYPE:		// DAS
      case LIST_TYPE:		// DAS
        ptr.a -> StrToVal (s, &_error, 
			   aDesc -> AggrElemTypeDescriptor (), 
			   instances, addFileId);
	break;

      case SELECT_TYPE:
	if (_error.severity( ptr.sh->STEPread(in, &_error, instances) ) 
		!= SEVERITY_NULL)
	_error.AppendToDetailMsg (ptr.sh ->Error ());
        break;

//      case UNKNOWN_TYPE: // should never have this case
      case GENERIC_TYPE:
      default:
	// other cases are the same for StrToVal and file
	return STEPread(in, instances, addFileId);
      }
    return _error.severity();
}


/******************************************************************
 ** Procedure:  STEPread
 ** Returns:  Severity, which indicates success, or failure
 **         value >= SEVERITY_WARNING means program can continue parsing input,
 **         value <= SEVERITY_INPUT_ERROR  is fatal read error
 ** Description:  the value of the attribute is read from istream
 **               the format expected is STEP exchange file
 **          modified to accept '$' as value for OPTIONAL ATTRIBUTE,
 **             (since this function is used for reading files in both
 **              formats, the function accepts either '$' (new format),
 **              or nothing (old format)) 31-Jul-1992
 ******************************************************************/
// does not read the delimiter separating individual attributes (i.e. ',') or 
// the delim separating the last attribute from the end of the entity (')').

Severity
STEPattribute::STEPread (istream& in, InstMgr * instances, int addFileId)
{

    char errStr[BUFSIZ];
    errStr[0] = '\0';
    char c ='\0';

    _error.ClearErrorMsg();	// also sets Severity to SEVERITY_NULL

    //  set the value to be null (reinitialize the attribute value)
    set_null();

    in >> ws; // skip whitespace
    in >> c;
    in.putback (c);    //  leave input stream alone

    //  check for NULL or derived attribute value, return if either
    switch (c)
    {
      case '$':
      case ',':
      case ')':
	if (c == '$') 
	{
	    in.get(c);
	    CheckRemainingInput(in, &_error, aDesc->TypeName(), ",)");
	}
	if(Nullable())  {
	    _error.severity(SEVERITY_NULL);
	}
	else {
	    _error.severity(SEVERITY_INCOMPLETE);
	    sprintf(errStr, " missing and required\n");
//		    " Warning: attribute '%s : %s' is missing and required.\n",
//		    Name(), TypeName());
	    _error.AppendToDetailMsg(errStr);
	}
	return _error.severity();

      case '*':
	in.get(c);  // take * off the istream
	if (IsDerived ())  
	    _error.severity(SEVERITY_NULL);
	else {
	    _error.severity(SEVERITY_INCOMPLETE);
	    sprintf(errStr, " attribute not derived\n");
	    _error.AppendToDetailMsg(errStr);
	}
	CheckRemainingInput(in, &_error, aDesc->TypeName(), ",)");
	return _error.severity();
    }

    PrimitiveType attrBaseType = NonRefType();
    switch( attrBaseType )
    {
      case INTEGER_TYPE:
  	{
	    int valAssigned = ReadInteger(*(ptr.i), in, &_error, ",)");
	    return _error.severity();
	}
      case REAL_TYPE:
  	{
	    int valAssigned = ReadReal(*(ptr.r), in, &_error, ",)");
	    return _error.severity();
	}
      case NUMBER_TYPE:
  	{
	    int valAssigned = ReadNumber(*(ptr.r), in, &_error, ",)");
	    return _error.severity();
	}
      case STRING_TYPE: 
  	{
	    ptr.S->STEPread (in, &_error);
	    CheckRemainingInput(in, &_error, "string", ",)");
	    return _error.severity();
	}
      case BINARY_TYPE:
	{	// call class SdaiBinary::STEPread()
	    ptr.b->STEPread(in, &_error);
	    CheckRemainingInput(in, &_error, "binary", ",)");
	    return _error.severity();
	}
      case BOOLEAN_TYPE:
  	{
//	    int nullable = (aDesc->Optional().asInt() == sdaiTRUE);
	    ptr.e->STEPread (in, &_error,  Nullable());
	    CheckRemainingInput(in, &_error, "boolean", ",)");
	    return _error.severity();
	}
      case LOGICAL_TYPE:
  	{
	    ptr.e->STEPread (in, &_error,  Nullable());
	    CheckRemainingInput(in, &_error, "logical", ",)");
	    return _error.severity();
	}
      case ENUM_TYPE:
  	{
	    ptr.e->STEPread (in, &_error,  Nullable());
	    CheckRemainingInput(in, &_error, "enumeration", ",)");
	    return _error.severity();
	}
      case AGGREGATE_TYPE:
      case ARRAY_TYPE:		// DAS
      case BAG_TYPE:		// DAS
      case SET_TYPE:		// DAS
      case LIST_TYPE:		// DAS
  	{
	    ptr.a->STEPread(in, &_error, 
			    aDesc->AggrElemTypeDescriptor(), 
			    instances, addFileId);

	    // cannot recover so give up and let STEPentity recover
	    if(_error.severity() < SEVERITY_WARNING)
		return _error.severity();

	    // check for garbage following the aggregate
	    CheckRemainingInput(in, &_error, "aggregate", ",)");
	    return _error.severity();
	}
      case ENTITY_TYPE:
  	{
	    STEPentity *se = ReadEntityRef(in, &_error, ",)", instances, 
					   addFileId);
	    if( se != S_ENTITY_NULL )
	    {
		if(EntityValidLevel(se, 
				    aDesc->NonRefTypeDescriptor(), 
				    &_error) == SEVERITY_NULL)
		    *(ptr.c) = se;
		else
		{
		    *(ptr.c) = S_ENTITY_NULL;
		}
	    }
	    else
		*(ptr.c) = S_ENTITY_NULL;
//	    CheckRemainingInput(in, &_error, "enumeration", ",)");
	    return _error.severity();

	}
      case SELECT_TYPE:
	    if (_error.severity( ptr.sh->STEPread(in, &_error, instances, addFileId) ) 
		  != SEVERITY_NULL)
		_error.AppendToDetailMsg (ptr.sh ->Error ());
	    CheckRemainingInput(in, &_error, "select", ",)");
	    return _error.severity();

      case GENERIC_TYPE:
  	{
	   G4cerr << "Internal error:  " << __FILE__ <<  __LINE__
	       << "\n" << _POC_ "\n";
	   _error.GreaterSeverity (SEVERITY_BUG);
	    return _error.severity();
	}

      case UNKNOWN_TYPE:
      case REFERENCE_TYPE:
      default:
  	{
	    // bug
	   G4cerr << "Internal error:  " << __FILE__ <<  __LINE__
	       << "\n" << _POC_ "\n";
	   _error.GreaterSeverity (SEVERITY_BUG);
	    return _error.severity();
	}
    }
}

/******************************************************************
 ** Procedure:  asStr
 ** Parameters:  
 ** Returns:  
 ** Description:  the value of the attribute is returned as a string.
 ** Side Effects:  
 ** Status:  complete 3/91
 ******************************************************************/
const char *
STEPattribute::asStr (SCLstring& str) const
{
//  char attrVal[BUFSIZ];
//  attrVal[0] = '\0';

    str.set_null();
    if (IsDerived ())  { str = "*"; return str.chars(); }
    if (is_null ())  { str = ""; return str.chars(); }

    switch (NonRefType()) {
      case INTEGER_TYPE:
      str.Append (*(ptr.i));
//	sprintf ( attrVal,"%ld",*(ptr.i));
	break;

      case NUMBER_TYPE:
      case REAL_TYPE:
      str.Append (*(ptr.r));
//	sprintf (attrVal, "%.*g", (int) Real_Num_Precision, *(ptr.r));
	break;	  

      case ENTITY_TYPE:
	// Print instance id only if not empty pointer
	// and has value assigned
	if ((*(ptr.c) == S_ENTITY_NULL) || (*(ptr.c) == 0))
	    break;
	else
	{
	    (*(ptr.c))->STEPwrite_reference (str);
//	    strncpy(attrVal, str.chars(), str.Length()+1);
	}
	break;	  
	
      case BINARY_TYPE:
	if( !( (ptr.b)->is_null() ) )
	    (ptr.b) -> STEPwrite (str);
//	    strncpy(attrVal, str.chars(), str.Length()+1);
	break;

      case STRING_TYPE:
	if( !( (ptr.S)->is_null() ) )
	    return (ptr.S) -> asStr (str);
	break;
	
      case AGGREGATE_TYPE:
      case ARRAY_TYPE:		// DAS
      case BAG_TYPE:		// DAS
      case SET_TYPE:		// DAS
      case LIST_TYPE:		// DAS
      return  ptr.a->asStr(str) ;
//	sprintf ( attrVal, "%s", ptr.a->asStr() );

      case ENUM_TYPE:
      case BOOLEAN_TYPE:
      case LOGICAL_TYPE:
	    return ptr.e -> asStr(str);

      case SELECT_TYPE:
	  ptr.sh -> STEPwrite (str);
	  return str;
//	  strncpy (attrVal, tmp.chars(), BUFSIZ);

      case REFERENCE_TYPE:
      case GENERIC_TYPE:
	   G4cerr << "Internal error:  " << __FILE__ <<  __LINE__
	       << "\n" << _POC_ "\n";
	   return 0;

      case UNKNOWN_TYPE:
      default:
	return (ptr.u -> asStr (str));
    }
//    return (attrVal);
    return str;
}

// The value of the attribute is printed to the output stream specified by out.
// The output is in physical file format.

void
STEPattribute::STEPwrite (ostream& out) 
{ 
    if (IsDerived ())  {
      out << "*";
      return;
    }

    if (is_null ()) {
	out << "$";
	return;
    }

    switch (NonRefType())
    {
      case INTEGER_TYPE:
	out << *(ptr.i);
	break;

      case NUMBER_TYPE:
      case REAL_TYPE:
      {
	real tmp = *(ptr.r);

        // IMS (6 Jun 95):
        // Previously, this code was explicitly rounding reals with
        // Magnitude smaller than 1e-9, by doing the following:

	// if(tmp > -.00000001 && tmp < .00000001) out << 0.0;
        // else...

        // No one's quite sure why this was happening, because the %g
        // style output used below should shift such things to e-notation
        // automagically.
        
#ifndef WIN32
        out << form("%.*g", (int) Real_Num_Precision,tmp);
#else
        out << (int) Real_Num_Precision << " " << tmp;
#endif
	break;	  
      }

      case ENTITY_TYPE:
	// Print instance id only if not empty pointer
	if ((*(ptr.c) == 0) ||  
	    // no value was assigned  <-- this would be a BUG
	    (*(ptr.c) == S_ENTITY_NULL) )
	{
	    out << "$";
	    G4cerr << "Internal error:  " << __FILE__ <<  __LINE__
		 << "\n" << _POC_ "\n";

	    char errStr[BUFSIZ];
	    errStr[0] = '\0';
	    _error.GreaterSeverity(SEVERITY_BUG);
	    sprintf(errStr,
		   " Warning: attribute '%s : %s' is null and shouldn't be.\n",
		    Name(), TypeName());
	    _error.AppendToUserMsg(errStr);
	    _error.AppendToDetailMsg(errStr);
	}
	else  (*(ptr.c)) -> STEPwrite_reference (out);
	break;	  

      case STRING_TYPE:
	// if null pointer or pointer to a string of length zero
	if(ptr.S)
	    (ptr.S) -> STEPwrite (out);
	else
	{
	    out << "$";
	    G4cerr << "Internal error:  " << __FILE__ <<  __LINE__
		 << "\n" << _POC_ "\n";

	    char errStr[BUFSIZ];
	    errStr[0] = '\0';
	    _error.GreaterSeverity(SEVERITY_BUG);
	    sprintf(errStr,
		   " Warning: attribute '%s : %s' should be pointing at %s",
		    Name(), TypeName(), "an SdaiString.\n");
	    _error.AppendToUserMsg(errStr);
	    _error.AppendToDetailMsg(errStr);
	}
	break;

      case BINARY_TYPE:	
	// if null pointer or pointer to a string of length zero
	if(ptr.b)
	    (ptr.b) -> STEPwrite (out);
	else
	{
	    out << "$";
	    G4cerr << "Internal error:  " << __FILE__ << ": " <<  __LINE__
		 << "\n" << _POC_ "\n";

	    char errStr[BUFSIZ];
	    errStr[0] = '\0';
	    _error.GreaterSeverity(SEVERITY_BUG);
	    sprintf(errStr,
		   " Warning: attribute '%s : %s' should be pointing at %s",
		    Name(), TypeName(), "an SdaiBinary.\n");
	    _error.AppendToUserMsg(errStr);
	    _error.AppendToDetailMsg(errStr);
	}
	break;

      case AGGREGATE_TYPE:
      case ARRAY_TYPE:		// DAS
      case BAG_TYPE:		// DAS
      case SET_TYPE:		// DAS
      case LIST_TYPE:		// DAS
	ptr.a -> STEPwrite (out);
	break;

      case ENUM_TYPE:
      case BOOLEAN_TYPE:
      case LOGICAL_TYPE:
	if(ptr.e)
	    ptr.e -> STEPwrite (out);
	else
	{
	    out << "$";
	    G4cerr << "Internal error:  " << __FILE__ << ": " <<  __LINE__
		 << "\n" << _POC_ "\n";

	    char errStr[BUFSIZ];
	    errStr[0] = '\0';
	    _error.GreaterSeverity(SEVERITY_BUG);
	    sprintf(errStr,
		   " Warning: attribute '%s : %s' should be pointing at %s",
		    Name(), TypeName(), "a STEPenumeration class.\n");
	    _error.AppendToUserMsg(errStr);
	    _error.AppendToDetailMsg(errStr);
	}
	break;

      case SELECT_TYPE:
	if(ptr.sh)
	    ptr.sh -> STEPwrite (out);
	else
	{
	    out << "$";
	    G4cerr << "Internal error:  " << __FILE__ <<  __LINE__
		 << "\n" << _POC_ "\n";

	    char errStr[BUFSIZ];
	    errStr[0] = '\0';
	    _error.GreaterSeverity(SEVERITY_BUG);
	    sprintf(errStr,
		   " Warning: attribute '%s : %s' should be pointing at %s",
		    Name(), TypeName(), "a STEPselect class.\n");
	    _error.AppendToUserMsg(errStr);
	    _error.AppendToDetailMsg(errStr);
	}
	break;

      case REFERENCE_TYPE:
      case GENERIC_TYPE:
	G4cerr << "Internal error:  " << __FILE__ << ": " <<  __LINE__
	     << "\n" << _POC_ "\n";
	_error.GreaterSeverity (SEVERITY_BUG);
	return;

      case UNKNOWN_TYPE:
      default:
	ptr.u -> STEPwrite (out);
	break;
	  
    }
}


BOOLEAN 
STEPattribute::ShallowCopy(STEPattribute *sa)
{
    switch(sa->NonRefType())
    {
      case INTEGER_TYPE:
	    *ptr.i = *(sa->ptr.i);
	    break;
      case BINARY_TYPE:
	    *(ptr.b) = *(sa->ptr.b);
	    break;
      case STRING_TYPE:
	    *(ptr.S) = *(sa->ptr.S);
	    break;
      case REAL_TYPE:
      case NUMBER_TYPE:
	    *ptr.r = *(sa->ptr.r);
	    break;
      case ENTITY_TYPE:
	    *ptr.c = *(sa->ptr.c);
	    break;
      case AGGREGATE_TYPE:
      case ARRAY_TYPE:		// DAS
      case BAG_TYPE:		// DAS
      case SET_TYPE:		// DAS
      case LIST_TYPE:		// DAS
	    ptr.a -> ShallowCopy (*(sa -> ptr.a));
	    break;
      case SELECT_TYPE:
	    *ptr.sh = *(sa->ptr.sh);
	    break;
      case ENUM_TYPE:
      case BOOLEAN_TYPE:
      case LOGICAL_TYPE:
	    ptr.e->put(sa->ptr.e->asInt());
	    break;

	  default:
	    *ptr.u = *(sa->ptr.u);
	    break;
    }
    return 1;
}

// for a string attribute this means make it not exist i.e. STEPstring will 
// exist in member variable ptr but STEPstring will be told to report as not
// containing a value (even a value of no chars).

Severity 
STEPattribute::set_null()
{  
    switch (NonRefType()) {
      case INTEGER_TYPE:
	*(ptr.i) = S_INT_NULL;
	break;

      case NUMBER_TYPE:
	*(ptr.r) = S_NUMBER_NULL;
	break;	  

      case REAL_TYPE:
	*(ptr.r) = S_REAL_NULL;
	break;	  

      case ENTITY_TYPE:
	*(ptr.c) = S_ENTITY_NULL;
	break;
	
      case STRING_TYPE:
	ptr.S -> set_undefined ();
	break;

      case BINARY_TYPE:
	ptr.b -> set_null ();
	break;

      case AGGREGATE_TYPE:
      case ARRAY_TYPE:		// DAS
      case BAG_TYPE:		// DAS
      case SET_TYPE:		// DAS
      case LIST_TYPE:		// DAS
      {
	ptr.a -> Empty ();
	break;
      }

      case ENUM_TYPE:
      case BOOLEAN_TYPE:
      case LOGICAL_TYPE:
	ptr.e -> set_null();
	break;

      case SELECT_TYPE:
	ptr.sh -> set_null();
	break;

      case REFERENCE_TYPE:  
      case GENERIC_TYPE:	
	G4cerr << "Internal error:  " << __FILE__ <<  __LINE__
	     << "\n" << _POC_ "\n";
	return SEVERITY_BUG;

      case UNKNOWN_TYPE:
      default:
      {
	ptr.u -> set_null();
	char errStr[BUFSIZ];
	errStr[0] = '\0';
	sprintf(errStr, " Warning: attribute '%s : %s : %d' - %s.\n",
		Name(), TypeName(), Type(),
		"Don't know how to make attribute NULL");
	_error.AppendToDetailMsg(errStr);
	_error.GreaterSeverity(SEVERITY_WARNING);
	return SEVERITY_WARNING;
      }
    }
    if(Nullable())
	return SEVERITY_NULL;
    else
	return SEVERITY_INCOMPLETE;
}

// For a string value this reports whether the string exists (as reported by 
// STEPstring) not whether STEPstring contains a null string.

BOOLEAN
STEPattribute::is_null ()  const
{
    switch ( NonRefType() )  
      {
	case INTEGER_TYPE:
	  return (*(ptr.i) == S_INT_NULL);

	case NUMBER_TYPE:
	  return  (*(ptr.r) == S_NUMBER_NULL);

	case REAL_TYPE:
	  return  (*(ptr.r) == S_REAL_NULL);

	case ENTITY_TYPE:
	  return  (*(ptr.c) == S_ENTITY_NULL);

	case STRING_TYPE:
	  return ptr.S -> is_undefined ();

	case BINARY_TYPE:
	  return ptr.b -> is_null ();

	case AGGREGATE_TYPE:
	case ARRAY_TYPE:		// DAS
	case BAG_TYPE:		// DAS
	case SET_TYPE:		// DAS
	case LIST_TYPE:		// DAS
	{
	  return  ( ptr.a -> is_null() );
	}

	case ENUM_TYPE:
	case BOOLEAN_TYPE:
	case LOGICAL_TYPE:
	  return ( ptr.e -> is_null() );

	case SELECT_TYPE:
	  return( ptr.sh->is_null() );

	case REFERENCE_TYPE:
	case GENERIC_TYPE:
	  G4cerr << "Internal error:  " << __FILE__ << ": " <<  __LINE__
	       << "\n" << _POC_ "\n";
	  return SEVERITY_BUG;

        case UNKNOWN_TYPE:
	default:
	  return (ptr.u -> is_null());
      }
}    

/******************************************************************
 ** Procedure:  operator == 
 ** Parameters:  STEPattribute & a1 and a2
 ** Returns:  int -- if 0 => not equal
 ** Description:  evaluates the equality of two attributes
 ** Side Effects:  none
 ** Status:  stub -- needs alot of work
 ******************************************************************/

//	equality for STEPattribute
int operator == (STEPattribute &a1, STEPattribute &a2)
{
  G4cerr << "Internal error:  " << __FILE__ << ": " <<  __LINE__
       << "\n" << _POC_ "\n";
    if (a1.aDesc->NonRefType() == a2.aDesc->NonRefType()) 
	;
    return 0;
}	



///////////////////////////////////////////////////////////////////////////////
// returns the severity level that the parameter attrValue would pass if it 
// was the value for this attribute.
// *note* for string values - (attrValue = 0) => string value does not exist,
//        attrValue exists it is valid.
///////////////////////////////////////////////////////////////////////////////

Severity 
STEPattribute::ValidLevel (const char *attrValue, ErrorDescriptor *error, 
			  InstMgr *im, int clearError)
{
    if(clearError)
	ClearErrorMsg();

    int optional = Nullable();

    if( !attrValue )
    {
	if(optional)
	    return error->severity();
	else
	{
	    error->GreaterSeverity(SEVERITY_INCOMPLETE);
	    return SEVERITY_INCOMPLETE;
	}
    }
    if( attrValue[0] == '\0' )
    {
	if(NonRefType() == STRING_TYPE)
	{ // this is interpreted as a string with no value i.e. "".
	  // Thus if it exists it has to be valid.
	    return SEVERITY_NULL;
	}
	if(optional)
	    return error->severity();
	else
	{
	    error->GreaterSeverity(SEVERITY_INCOMPLETE);
	    return SEVERITY_INCOMPLETE;
	}
    }

    //  an overridden attribute always has a \'*\' value
    if (IsDerived ())  {
	if (!strcmp (attrValue, "*")) return SEVERITY_NULL;
	else {
	  _error.AppendToDetailMsg ("attr is derived - value not permitted\n");
	  return _error.severity (SEVERITY_INPUT_ERROR);
	}
    }

    switch (NonRefType() ) {
      case INTEGER_TYPE:
	    return IntValidLevel(attrValue, error, clearError, optional, 0);

      case STRING_TYPE:
	    // if a value exists (checked above) then that is the string value
	    return SEVERITY_NULL;
	
      case REAL_TYPE:
	    return RealValidLevel(attrValue, error, clearError, optional, 0);

      case NUMBER_TYPE:
	    return NumberValidLevel(attrValue, error, clearError, optional, 0);

      case ENTITY_TYPE:
	    return EntityValidLevel(attrValue, 
				    aDesc->NonRefTypeDescriptor(), 
				    error, im, 0);
      case BINARY_TYPE:
	    return ptr.b->BinaryValidLevel(attrValue, &_error, optional, 0);

      case AGGREGATE_TYPE:
      case ARRAY_TYPE:		// DAS
      case BAG_TYPE:		// DAS
      case SET_TYPE:		// DAS
      case LIST_TYPE:		// DAS
      {
	return ptr.a->AggrValidLevel(attrValue, error, 
			      aDesc->AggrElemTypeDescriptor(), im, 
			      optional, 0, 0, 0);
      }
      case ENUM_TYPE:
      case BOOLEAN_TYPE:
      case LOGICAL_TYPE:
	    return ptr.e->EnumValidLevel(attrValue, error, optional, 0, 0, 1);
      case SELECT_TYPE:
	    return ptr.sh->SelectValidLevel(attrValue, error, im, 0);

      default:
	    G4cerr << "Internal error:  " << __FILE__ <<  __LINE__
		 << "\n" << _POC_ "\n";
	    return error->severity();
    }
}
  
/******************************************************************
 ** Procedure:  operator <<
 ** Parameters:  ostream & out -- output stream
 **              STEPattribute & a -- attribute to output
 ** Returns:  ostream &
 ** Description:  overloads the output operator to Print an attribute 
 ******************************************************************/

ostream &operator<< ( ostream& out, STEPattribute& a )
{
    a.STEPwrite (out);
    return out;

}


///////////////////////////// AddErrorInfo() //////////////////////////////////
// This adds prepends attribute information to the detailed error msg.  This 
// is intended to add information to the error msgs written by Enumerations, 
// Aggregates, and STEPstrings which don't know they are a STEPattribute value.
///////////////////////////////////////////////////////////////////////////////

void 
STEPattribute::AddErrorInfo()
{
    char errStr[BUFSIZ];
    errStr[0] = '\0';

    if(SEVERITY_INPUT_ERROR < _error.severity() &&
       _error.severity() < SEVERITY_NULL)
    {
	sprintf(errStr, " Warning: ATTRIBUTE '%s : %s : %d' - ",
		Name(), TypeName(), Type());
	_error.PrependToDetailMsg(errStr);
    }
    else if(_error.severity() == SEVERITY_INPUT_ERROR)
    {
	sprintf(errStr, " Error: ATTRIBUTE '%s : %s : %d' - ",
		Name(), TypeName(), Type());
	_error.PrependToDetailMsg(errStr);
    }
    else if(_error.severity() <= SEVERITY_BUG)
    {
	sprintf(errStr, " BUG: ATTRIBUTE '%s : %s : %d' - ",
		Name(), TypeName(), Type());
	_error.PrependToDetailMsg(errStr);
    }
}

///////////////////////////////////////////////////////////////////////////////

// this function reads until it hits eof or one of the StopChars.
// if it hits one of StopChars it puts it back.
// RETURNS: the last char it read.
char
STEPattribute::SkipBadAttr(istream& in, char *StopChars)
{
#ifndef WIN32
    int sk = in.skip(0);  //turn skipping whitespace off
#endif

	// read bad data until end of this attribute or entity.
    char *foundCh = 0;
    char c = '\0';
    char errStr[BUFSIZ];
    errStr[0] = '\0';

    _error.GreaterSeverity(SEVERITY_WARNING);
    in >> c;
    while( !in.eof() && !(foundCh = strchr(StopChars, c)) )
	  in >> c;
    if(in.eof())
    {
	_error.GreaterSeverity(SEVERITY_INPUT_ERROR);
	sprintf(errStr, " Error: attribute '%s : %s : %d' - %s.\n",
		Name(), TypeName(), Type(),
		"Unexpected EOF when skipping bad attr value");
	_error.AppendToDetailMsg(errStr);
    }
    else
    {
	sprintf(errStr, " Error: attribute '%s : %s : %d' - %s.\n",
		Name(), TypeName(), Type(), "Invalid value");
	_error.AppendToDetailMsg(errStr);
    }
    in.putback (c);
#ifndef WIN32
    in.skip(sk);     //set skip whitespace to its original state
#endif
    return c;
}


/////////////////// READ 

#ifdef OBSOLETE

// appends char 'c' to char * 's' starting at index 'index'.  'index' is 
// incremented by one and returned.  If assigning 'c' will cause 's' to be 
// larger than size 'sSize' then 's' is deleted and reallocated and 's' and 
// 'sSize' are changed and returned appropriately.

void
AppendChar(char c, int& index, char *&s, int& sSize)
{
    if(index >= sSize - 1)
    {
	char *tmp = s;
	s = new char[sSize * 2];
	strcpy(s, tmp);
	delete [] tmp;
	sSize = sSize * 2;
    }
    s[index] = c;
    index++;
    s[index] = '\0';
}


/////////////////// STEPread_reference()

STEPentity *
STEPread_reference (const char * s, ErrorDescriptor *err, InstMgr * instances, 
		     int addFileId)
{
    char errStr[BUFSIZ];
    errStr[0] = '\0';

    int fileId = -1;
    int numfound;
    if( numfound = sscanf((char *)s, " #%d ", &fileId) )
	;
    else ( numfound = sscanf((char *)s, " @%d ", &fileId) )
	;
    if( (numfound != EOF) && (numfound > 0) )
    {
	fileId = fileId + addFileId;

	if (!instances)
	{
	    sprintf(errStr, "STEPread_reference(): %s - entity #%d %s.\n",
		    "BUG - cannot read reference without the InstMgr", fileId, 
		    "is unknown to attribute");
	    err->AppendToDetailMsg(errStr);
	    err->GreaterSeverity(SEVERITY_BUG);
	    return S_ENTITY_NULL;
	}
	
	//  lookup which object has id as its instance fileId
	STEPentity* inst;
	/* If there is a ManagerNode it should have a STEPentity */
	MgrNode* mn = 0;
	mn = instances->FindFileId(fileId);
	if (mn)
	{
	    inst =  mn->GetSTEPentity() ;
	    if (inst) { return (inst); }
	    else
	    {
		sprintf(errStr, 
			"%s - entity #%d %s.\n",
		        "BUG - MgrNode::GetSTEPentity returned NULL pointer",
			fileId, "is unknown to attribute");
		err->AppendToDetailMsg(errStr);
		err->GreaterSeverity(SEVERITY_BUG);
		return S_ENTITY_NULL;
	    }
	}
	else
	{
	    sprintf(errStr,"Reference to non-existent ENTITY #%d.\n", fileId);
	    err->AppendToDetailMsg(errStr);
	    err->GreaterSeverity(SEVERITY_WARNING);
	    return S_ENTITY_NULL;
	}
    }
    else
    {
	numfound = sscanf((char *)s, " %*s ");
	if(numfound == EOF)
	{
	    err->GreaterSeverity(SEVERITY_INCOMPLETE);
	}
	else 
	{
	    sprintf(errStr,"  Invalid entity identifier # %s.\n", s);
	    err->AppendToDetailMsg(errStr);
	    err->GreaterSeverity(SEVERITY_WARNING);
	}
	return S_ENTITY_NULL;
    }
}

#endif
