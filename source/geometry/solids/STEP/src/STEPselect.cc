

//



//
// $Id: STEPselect.cc,v 1.2 1999-05-21 20:20:52 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

/*
* NIST STEP Core Class Library
* clstepcore/STEPselect.cc
* May 1995
* Dave Helfrick
* KC Morris

* Development of this software was funded by the United States Government,
* and is not subject to copyright.
*/

/*  */

#include <STEPselect.h>
#include <stdio.h> // to get the BUFSIZ #define
#include <ExpDict.h>
#ifdef WIN32
#  include <Strstrea.h>
#else
#  include <strstream.h>
#endif
#include <STEPentity.h>

/**********
	(member) functions for the select class SdaiSelect
**********/
/*
SdaiSelect::SdaiSelect()
{
}
SdaiSelect::SdaiSelect( const SdaiSelect& s )
{
}
*/

Severity SdaiSelect::severity() const
{		/**  fn Severity  **/
   return _error.severity();
}		/**  end function  **/

Severity  SdaiSelect::severity( Severity s )
{		/**  fn Severity  **/
   return _error.severity( s );
}		/**  end function  **/

const char * SdaiSelect::Error()
{		/**  fn Error  **/
   return _error.DetailMsg();
}		/**  end function  **/

void SdaiSelect::Error( char * e )
{		/**  fn Error  **/
   _error.DetailMsg( e );
}		/**  end function  **/

void  SdaiSelect::ClearError()
{		/**  end function  **/
   _error.ClearErrorMsg();
}		/**  end function  **/

const TypeDescriptor *
SdaiSelect::CanBe (const char * n) const  
{
    const TypeDescLinkNode * tdn = 
		(const TypeDescLinkNode *) _type -> GetElements ().GetHead ();
    const TypeDescriptor *td = tdn -> TypeDesc ();

    SCLstring tmp1;
    SCLstring tmp2;
  
    while (tdn)  {
	td = tdn -> TypeDesc ();
	if (strcmp( (char *)StrToUpper( n, tmp1 ), 
		    (char *)StrToUpper( td->Name(), tmp2 ) ))
	    ;
	else 
	    return td;  // they are the same 
	tdn = (TypeDescLinkNode *) (tdn -> NextNode ());
    }
    return 0;
}


const TypeDescriptor *
SdaiSelect::CanBe (BASE_TYPE bt) const  
{
  const TypeDescLinkNode * tdn = (const TypeDescLinkNode *) _type -> GetElements ().GetHead ();
  const TypeDescriptor *td = tdn -> TypeDesc ();

  while (tdn)  {
    td = tdn -> TypeDesc ();    
    if (td -> NonRefType () == bt)
      return td;  // they are the same 
    tdn = (TypeDescLinkNode *) (tdn -> NextNode ());
  }
  return 0;
}

const TypeDescriptor *
SdaiSelect::CanBe (const TypeDescriptor * td) const  
{
  return _type -> CanBe (td);
}

SdaiString SdaiSelect::UnderlyingTypeName () const  
{
  return underlying_type-> Name ();
}

const TypeDescriptor *  SdaiSelect::CurrentUnderlyingType() const 
{
  return underlying_type;
}

const TypeDescriptor * 
SdaiSelect::SetUnderlyingType (const TypeDescriptor * td)
{
  //  don\'t do anything if the descriptor is bad
  if (!td || !(_type -> CanBe (td))) return 0;

  base_type = td -> NonRefType ();
  return underlying_type = td;
}

int SdaiSelect::exists() const 
{
  return (int) underlying_type;
}

void SdaiSelect::nullify() 
{
  underlying_type = 0;
}

Severity 
SdaiSelect::SelectValidLevel(const char *attrValue, ErrorDescriptor *err, 
			     InstMgr *im, int clearError)
{
  SdaiSelect * tmp = NewSelect();
  Severity s = SEVERITY_NULL;

//  COMMENTED OUT TO ALLOW FOR TYPE RESOLUTION WHEN TYPE CHANGES
//  if (CurrentUnderlyingType ())
//    s = tmp -> StrToVal 
//      (attrValue, CurrentUnderlyingType () -> Name (), err, im);
//  else
  istrstream paska((char *)attrValue);
  s = tmp -> STEPread (paska, err, im);
  delete tmp;
  return s;
}

Severity 
SdaiSelect::StrToVal(const char *Val, const char *selectType, 
		     ErrorDescriptor *err, InstMgr * instances)
{
    if (SetUnderlyingType (CanBe (selectType)))  

	//  the underlying type is set to a valid type
	// call read on underlying type in subclass

	switch (base_type)  {
	  case ENTITY_TYPE:
	  {  
	    STEPentity * tmp = 
	      ReadEntityRef (Val, err, ",)", instances, 0);
	    if( tmp && ( tmp != ENTITY_NULL ) && AssignEntity(tmp) ) 
	      return SEVERITY_NULL;  
	    else
	    {
		err->AppendToDetailMsg( 
		"Reference to entity that is not a valid type for SELECT.\n" );
		nullify ();	   
		err->GreaterSeverity( SEVERITY_WARNING );
		return SEVERITY_WARNING;
	    }	
	  }

	    // call StrToVal on the contents
	  case STRING_TYPE:
	  case AGGREGATE_TYPE:
	  case ARRAY_TYPE:		// DAS
	  case BAG_TYPE:		// DAS
	  case SET_TYPE:		// DAS
	  case LIST_TYPE:		// DAS
	  case ENUM_TYPE:
	  case SELECT_TYPE:
	  case BOOLEAN_TYPE:
	  case LOGICAL_TYPE:
	  {
	      err->GreaterSeverity( StrToVal_content (Val, instances) );
	      if(_error.severity() != SEVERITY_NULL)
		  err->AppendFromErrorArg(&_error);
	      return err->severity();
	  }

	    // do the same as STEPread_content
	  case BINARY_TYPE:
	  case NUMBER_TYPE:
	  case REAL_TYPE:
	  case INTEGER_TYPE:
	  default:
	  {
	    istrstream paska((char *)Val);
	      err->GreaterSeverity( STEPread_content( paska ) );
	      if(_error.severity() != SEVERITY_NULL)
		  err->AppendFromErrorArg(&_error);
	      return err->severity();
	  }
	}
    return SEVERITY_INPUT_ERROR;
}

Severity
SdaiSelect::STEPread  (istream& in, ErrorDescriptor *err, InstMgr * instances, 
		       int addFileId)
{
    char c ='\0';
    SCLstring tmp;

    // find out what case we have
    //  NOTE case C falls out of recursive calls in cases A and B

    in >> ws;
    in >> c;
//    c = in.peek();
    if( isalpha(c) ) //  case B 
    {
	//  token is a type name - get the type
      int eot =0;  // end of token flag
	while( (c != '(') && in.good() )
	{
	    if (!eot && ! (eot = isspace (c)) ) 
	      // as long as eot hasn\'t been reached keep appending
		tmp.Append(c);
	    in >> c; 
	}

	//  check for valid type and  set the underlying type 

	if( SetUnderlyingType( CanBe( tmp.chars() ) ) )
	{  //  assign the value to the underlying type
	    in >> ws; // skip white space
	    STEPread_content (in, instances, addFileId);
//  STEPread_content uses the ErrorDesc data member from the STEPselect class
	    err->AppendToDetailMsg (Error ());
	    err->GreaterSeverity( severity () );
	    in >> ws >> c;  
	    if( c != ')' )
	    {
		err->AppendToDetailMsg( 
				"Bad data or missing closing ')' for SELECT type.\n" );
		err->GreaterSeverity( SEVERITY_WARNING );
		in.putback (c);
		return SEVERITY_WARNING;
	    }
	    return err->severity();
	}
	else // ERROR  -- the type wasn't one of the choices
	{
	    if( !in.good() )
	    {
		err->GreaterSeverity( SEVERITY_INPUT_ERROR );
		return SEVERITY_INPUT_ERROR;
	    }
	    else
	    {
		err->AppendToDetailMsg (
			"The type name for the SELECT type is not valid.\n");
		err->GreaterSeverity (SEVERITY_WARNING);
		return SEVERITY_WARNING;
	    }
	}
/*	if (!in.good())
	{
	    err->GreaterSeverity( SEVERITY_INPUT_ERROR );
	    return SEVERITY_INPUT_ERROR;
	}
*/
    }
    else // case A
    {    //  the type can be determined from the value 
      if (_type && ! _type -> UniqueElements ())  {
	err->AppendToDetailMsg("Type for value of SELECT is ambiguous.\n" );
	err->GreaterSeverity( SEVERITY_WARNING );
      }
	switch (c)
	{
	  case '$':  
	    nullify ();
	    err->GreaterSeverity (SEVERITY_INCOMPLETE);
	    return SEVERITY_INCOMPLETE;

	  case ',':
	  case '\0':
	    // ERROR  IN INPUT
	    in.putback (c);
	    err->AppendToDetailMsg( "No value found for SELECT type.\n" );
	    err->GreaterSeverity( SEVERITY_WARNING );
	    return SEVERITY_WARNING;

	  case '.': // assign enum
	    base_type = ENUM_TYPE;
	    break;
	    // set the underlying type
	    // call STEPread
	    // return

	  case '\'': // assign string
	    base_type = STRING_TYPE;
	    break;

	  case '"': // assign string
	    base_type = BINARY_TYPE;
	    break;

	  case '#':
	    base_type = ENTITY_TYPE;
	    break;
	    // call STEPread_reference
	    // set the underlying type

	    // assign entity
	    // read the reference
	    // match type to underlying type
	    // assign the value
	    // set the underlying type

	  case '(':
	  {
	    char n;
	    in >> n;
	    in.putback(n);
	    if( isalpha(n) )  
		base_type = SELECT_TYPE;
	    else 
		base_type = AGGREGATE_TYPE;
	    in.putback( n );
	    break;
	  }

	  case '0':
	  case '1':
	  case '2':
	  case '3':
	  case '4':
	  case '5':
	  case '6':
	  case '7':
	  case '8':
	  case '9':
	  case '-':
	    if( CanBe( REAL_TYPE ) ) 
		base_type = REAL_TYPE;
	    else 	
		base_type = INTEGER_TYPE;
	    break;

	  default:
	    // ambiguous - ERROR:  underlying type should have been set
	    //	     STEPread_error ();
	    err->AppendToDetailMsg( 
		    "type for SELECT could not be determined from value.\n" );
	    nullify ();
	    in.putback (c);
	    err->GreaterSeverity( SEVERITY_WARNING );
	    return SEVERITY_WARNING;
	}

	in.putback (c);

	// now the type descriptor should be derivable from the base_type
	if( base_type == ENTITY_TYPE ) 
	{  // you don\'t know if this is an ENTITY or a SELECT
	    // have to do this here - not in STEPread_content
	    STEPentity * tmp = 
		ReadEntityRef (in, err, ",)", instances, addFileId);
	    if( tmp && ( tmp != ENTITY_NULL ) && AssignEntity(tmp) ) 
		return SEVERITY_NULL;  
	    else
	    {
		err->AppendToDetailMsg( 
		"Reference to entity that is not a valid type for SELECT.\n" );
		nullify ();	   
		err->GreaterSeverity( SEVERITY_WARNING );
		return SEVERITY_WARNING;
	    }	
	}
	else if ( SetUnderlyingType( CanBe(base_type) ) )
	    STEPread_content( in, instances, addFileId );

	else // ERROR  -- the type wasn\'t one of the choices
	{
	    err->AppendToDetailMsg( 
			    "The type of the SELECT type is not valid.\n");
	    err->GreaterSeverity( SEVERITY_WARNING );
	    return SEVERITY_WARNING;
	}
    }
//    if (!in.good())	severity (SEVERITY_INPUT_ERROR);
    return err->severity ();
}


void
SdaiSelect::STEPwrite(ostream& out)  const
{
    if (!exists ()) {  
	out << "$";
	return;
    }
    if (_type -> UniqueElements ())
	STEPwrite_content( out ); 
    else
    {
	SCLstring tmp;
	out << StrToUpper( CurrentUnderlyingType()->Name(), tmp ) << "(";
	STEPwrite_content (out); 
	out << ")";
    }
}

const char *
SdaiSelect::STEPwrite(SCLstring& s)  const
{
  strstream buf;
  STEPwrite (buf);
  buf << ends;  // have to add the terminating \0 char
  char * tmp;
  tmp = buf.str ();
  s = tmp;
  delete tmp;
  return s;
}


//SdaiSelect& 
//SdaiSelect::operator= (SdaiSelect& x)  
//{
//    return *this;
//}

#ifdef OBSOLETE
char* 
SdaiSelect::asStr()  const
{
  strstream buf;
  STEPwrite (buf);
/*  STEPwrite_content (buf);*/
  buf << ends;  // have to add the terminating \0 char
  return buf.str ();  // need to delete this space
}
#endif

int
SdaiSelect::set_null ()  
{
  nullify ();
    return 1;
}

int
SdaiSelect::is_null ()  
{
  return (!exists ());
}
