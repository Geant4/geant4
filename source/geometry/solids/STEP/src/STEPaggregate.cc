

//



//
// $Id: STEPaggregate.cc,v 1.2 1999-05-21 20:20:48 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

/*
* NIST STEP Core Class Library
* clstepcore/STEPaggregate.cc
* May 1995
* K. C. Morris
* David Sauder

* Development of this software was funded by the United States Government,
* and is not subject to copyright.
*/

/*  */

#include <stdio.h> 

#include <read_func.h>
#include <STEPaggregate.h>
#include <STEPattribute.h>
#include <STEPentity.h>
#include <instmgr.h>

const int Real_Num_Precision = REAL_NUM_PRECISION; // from STEPattribute.h

#define STRING_DELIM '\''

static char rcsid[] = "";

/******************************************************************
**	  Functions for manipulating aggregate attributes

**  KNOWN BUGs:  
**     -- treatment of aggregates of reals or ints is inconsistent with
**        other aggregates (there's no classes for these)
**     -- no two- dimensional aggregates are implemented
**/

STEPaggregate NilSTEPaggregate;


///////////////////////////////////////////////////////////////////////////////
// STEPaggregate
///////////////////////////////////////////////////////////////////////////////

STEPaggregate::STEPaggregate  ()  
{
    _null = 1;
}

STEPaggregate::~STEPaggregate  ()  
{
}

STEPaggregate& 
STEPaggregate::ShallowCopy (const STEPaggregate& a) 
{
  G4cerr << "Internal error:  " << __FILE__ << ": " <<  __LINE__
       << "\n" << _POC_ "\n";
  G4cerr << "function:  STEPaggregate::ShallowCopy \n" << "\n";
  return *this;
}

    // do not require exchange file format
Severity 
STEPaggregate::AggrValidLevel(const char *value, ErrorDescriptor *err, 
			      const TypeDescriptor *elem_type, InstMgr *insts, 
			      int optional, char *tokenList, int addFileId, 
			      int clearError)
{
    if(clearError)
	err->ClearErrorMsg();

    istrstream in ((char *)value); // sz defaults to length of s

    ReadValue(in, err, elem_type, insts, addFileId, 0, 0);
    CheckRemainingInput(in, err, elem_type->AttrTypeName(), tokenList);
    if( optional && (err->severity() == SEVERITY_INCOMPLETE) )
	err->severity(SEVERITY_NULL);
    return err->severity();
}

    // require exchange file format
Severity 
STEPaggregate::AggrValidLevel(istream &in, ErrorDescriptor *err, 
			      const TypeDescriptor *elem_type, InstMgr *insts, 
			      int optional, char *tokenList, int addFileId, 
			      int clearError)
{
    if(clearError)
	err->ClearErrorMsg();

    ReadValue(in, err, elem_type, insts, addFileId, 0, 1);
    CheckRemainingInput(in, err, elem_type->AttrTypeName(), tokenList);
    if( optional && (err->severity() == SEVERITY_INCOMPLETE) )
	err->severity(SEVERITY_NULL);
    return err->severity();
}

// if exchangeFileFormat == 1 then paren delims are required.

Severity 
STEPaggregate::ReadValue(istream &in, ErrorDescriptor *err, 
			 const TypeDescriptor *elem_type, 
			 InstMgr *insts, int addFileId, 
			 int assignVal, int exchangeFileFormat)
{
    ErrorDescriptor errdesc;
    char errmsg[BUFSIZ];
    int value_cnt = 0;

    if(assignVal)
	Empty ();  // read new values and discard existing ones

    char c;
    int validDelims = 1;

    in >> ws; // skip white space

    c = in.peek(); // does not advance input

    if(in.eof() || c == '$')
    {
	_null = 1;
	err->GreaterSeverity(SEVERITY_INCOMPLETE);
	return SEVERITY_INCOMPLETE;
    }

    if(c == '(')
    {
	in.get(c);
	validDelims = 0; // signal expectation for end paren delim
    }
    else if(exchangeFileFormat)
    {	// error did not find opening delim
	    // cannot recover so give up and let STEPattribute recover
	err->GreaterSeverity(SEVERITY_INPUT_ERROR);
	return SEVERITY_INPUT_ERROR;
    }
    else if(!in.good())
    {// this should actually have been caught by skipping white space above
	err->GreaterSeverity(SEVERITY_INCOMPLETE);
	return SEVERITY_INCOMPLETE;
    }

    STEPnode * item = 0;
    if(!assignVal)  // if not assigning values only need one node.  So only 
		    // one node is created. It is used to read the values
	item = (STEPnode*)NewNode();

	// ')' is the end of the aggregate
    while (in.good() && (c != ')') )
    {
	value_cnt++;
	if(assignVal) // create a new node each time through the loop
	    item = (STEPnode*)NewNode();

	errdesc.ClearErrorMsg();

	if(exchangeFileFormat)
	    item->STEPread(in, &errdesc);
	else
	    item->StrToVal(in, &errdesc);

	// read up to the next delimiter and set errors if garbage is
	// found before specified delims (i.e. comma and quote)
	CheckRemainingInput(in, &errdesc, elem_type->AttrTypeName(), ",)");

	if (errdesc.severity() < SEVERITY_INCOMPLETE)
	{
	    sprintf (errmsg, "  index:  %d\n", value_cnt );
	    errdesc.PrependToDetailMsg(errmsg);
	    err->AppendFromErrorArg(&errdesc);
	}
	if(assignVal) // pass the node to STEPaggregate
	    AddNode (item);

	in >> ws; // skip white space (although should already be skipped)
	in.get(c); // read delim

	// CheckRemainingInput should have left the input right at the delim
	// so that it would be read in in.get() above.  Since it did not find
	// the delim this does not know how to find it either!
	if( (c != ',') && (c != ')') )
	{
	    // cannot recover so give up and let STEPattribute recover
	    err->GreaterSeverity(SEVERITY_INPUT_ERROR);
	    return SEVERITY_INPUT_ERROR;
/*
	    // error read until you find a delimiter
	    SCLstring tmp;
	    while(in.good() && !strchr(",)", c) )
	    {
		in.get(c);
		tmp.Append(c);
	    }
	    // BUG could overwrite the error message buffer
	    sprintf(errmsg, "ERROR aggr. elem. followed by \'%s\' garbage.\n",
		    tmp.chars());
	    err->AppendToDetailMsg(errmsg);
	    err->AppendToUserMsg(errmsg);
	    err->GreaterSeverity(SEVERITY_WARNING);
	    if(!in.good())
		return err->severity();
*/
	}
    }
    return err->severity();
}

Severity 
STEPaggregate::StrToVal(const char *s, ErrorDescriptor *err, 
			const TypeDescriptor *elem_type, InstMgr *insts, 
			int addFileId)
{
    istrstream in((char *)s);
    return ReadValue(in, err, elem_type, insts, addFileId, 1, 0);
}

///////////////////////////////////////////////////////////////////////////////

Severity 
STEPaggregate::STEPread(istream& in, ErrorDescriptor *err, 
			const TypeDescriptor *elem_type, 
			InstMgr *insts, int addFileId)
{
    return ReadValue(in, err, elem_type, insts, addFileId, 1, 1);
}

const char * 
STEPaggregate::asStr(SCLstring & s) const
{
    s.set_null();

    if(!_null)
    {
	s = "(";
	STEPnode * n = (STEPnode *) head;
	SCLstring tmp;
	while (n)
	{
	    s.Append( n->STEPwrite(tmp) );
	    if (n = (STEPnode *) n -> NextNode ())
		s.Append(',');
	}
	s.Append(')');
    }
    return s.chars();
}

void
STEPaggregate::STEPwrite(ostream& out) const
{
    if(!_null)
    {
	out << '(';
	STEPnode * n = (STEPnode *)head;
	SCLstring s;
	while (n)  
	{
	    out << n->STEPwrite (s);
	    if ( n = (STEPnode *)(n -> NextNode ()) )
		out <<  ',';
	}
	out << ')';
    }
    else
	out << '$';
}

SingleLinkNode *
STEPaggregate::NewNode () 
{
  G4cerr << "Internal error:  " << __FILE__ << ": " <<  __LINE__ << "\n" ;
  G4cerr << "function:  STEPaggregate::NewNode \n" << _POC_ << "\n";
  return 0;
}

void
STEPaggregate::AddNode(SingleLinkNode *n)
{
    SingleLinkList::AppendNode(n);
    _null = 0;
}

void
STEPaggregate::Empty()
{
    SingleLinkList::Empty();
    _null = 1;
}


///////////////////////////////////////////////////////////////////////////////
// STEPnode
///////////////////////////////////////////////////////////////////////////////

Severity
STEPnode::StrToVal(const char *s, ErrorDescriptor *err)
{
    // defined in subtypes
  G4cerr << "Internal error:  " << __FILE__ << ": " <<  __LINE__ << "\n" ;
    err->AppendToDetailMsg(
	" function: STEPnode::StrToVal() called instead of virtual function.\n"
			   );
    err->AppendToDetailMsg("Aggr. attr value: '\n");
    err->AppendToDetailMsg("not assigned.\n");
    err->AppendToDetailMsg(_POC_);
    err->GreaterSeverity(SEVERITY_BUG);
    return SEVERITY_BUG;
}

Severity
STEPnode::StrToVal(istream &in, ErrorDescriptor *err)
{
    // defined in subtypes
  G4cerr << "Internal error:  " << __FILE__ << ": " <<  __LINE__ << "\n" ;
    err->AppendToDetailMsg(
	" function: STEPnode::StrToVal() called instead of virtual function.\n"
			   );
    err->AppendToDetailMsg("Aggr. attr value: '\n");
    err->AppendToDetailMsg("not assigned.\n");
    err->AppendToDetailMsg(_POC_);
    err->GreaterSeverity(SEVERITY_BUG);
    return SEVERITY_BUG;
}

Severity 
STEPnode::STEPread(const char *s, ErrorDescriptor *err)
{
    //  defined in subclasses
  G4cerr << "Internal error:  " << __FILE__ << ": " <<  __LINE__ << "\n" ;
  G4cerr << "function:  STEPnode::STEPread called instead of virtual function.\n"
       << _POC_ << "\n";

    err->AppendToDetailMsg(
	" function: STEPnode::STEPread() called instead of virtual function.\n"
			   );
    err->AppendToDetailMsg(_POC_);
    err->GreaterSeverity(SEVERITY_BUG);

    return SEVERITY_BUG;
}

Severity 
STEPnode::STEPread(istream &in, ErrorDescriptor *err)
{
  G4cerr << "Internal error:  " << __FILE__ << ": " <<  __LINE__ << "\n" ;
  G4cerr << "function:  STEPnode::STEPread called instead of virtual function.\n"
       << _POC_ << "\n";

    err->AppendToDetailMsg(
	" function: STEPnode::STEPread() called instead of virtual function.\n"
			   );
    err->AppendToDetailMsg(_POC_);
    err->GreaterSeverity(SEVERITY_BUG);
    return SEVERITY_BUG;
}

const char *
STEPnode::asStr(SCLstring &s)
{
    //  defined in subclasses
  G4cerr << "Internal error:  " << __FILE__ << ": " <<  __LINE__ << "\n" ;
  G4cerr << "function:  STEPnode::asStr called instead of virtual function.\n" 
       << _POC_ << "\n";

  return "";
}

const char *
STEPnode::STEPwrite (SCLstring &s)
{
  G4cerr << "Internal error:  " << __FILE__ << ": " <<  __LINE__ << "\n" ;
 G4cerr << "function:  STEPnode::STEPwrite called instead of virtual function.\n"
      << _POC_ << "\n";

 return "";
}

void 
STEPnode::STEPwrite (ostream& out)
{
  G4cerr << "Internal error:  " << __FILE__ << ": " <<  __LINE__ << "\n" ;
 G4cerr << "function:  STEPnode::STEPwrite called instead of virtual function.\n"
       << _POC_ << "\n";
}

///////////////////////////////////////////////////////////////////////////////
// GenericAggregate
///////////////////////////////////////////////////////////////////////////////

SingleLinkNode *
GenericAggregate::NewNode () 
{
    return new GenericAggrNode();
}

STEPaggregate& 
GenericAggregate::ShallowCopy (const STEPaggregate& a)
{
    Empty();
    
    SingleLinkNode* next = a.GetHead();
    SingleLinkNode* copy;

    while (next) 
    {
	copy = new GenericAggrNode (*(GenericAggrNode*)next);
	AddNode(copy);
	next = next->NextNode();
    }
    if(head)
	_null = 0;
    else
	_null = 1;
    return *this;
    
}

///////////////////////////////////////////////////////////////////////////////
// GenericAggrNode
///////////////////////////////////////////////////////////////////////////////

GenericAggrNode::GenericAggrNode (const char *str)
{  
    value = str;
}

GenericAggrNode::GenericAggrNode (GenericAggrNode& gan)
{  
    value = gan.value;
}

GenericAggrNode::GenericAggrNode()
{
}

GenericAggrNode::~GenericAggrNode()
{
}

SingleLinkNode *
GenericAggrNode::NewNode () 
{
    return new GenericAggrNode();
}

Severity 
GenericAggrNode::StrToVal(const char *s, ErrorDescriptor *err)
{
    return value.STEPread(s, err);
}

//TODO
Severity 
GenericAggrNode::StrToVal(istream &in, ErrorDescriptor *err)
{
    return value.STEPread(in, err);
}

Severity 
GenericAggrNode::STEPread(const char *s, ErrorDescriptor *err)
{
    istrstream in((char *) s);
    return value.STEPread(in, err);
}

Severity 
GenericAggrNode::STEPread(istream &in, ErrorDescriptor *err)
{
    return value.STEPread(in, err);
}

const char *
GenericAggrNode::asStr(SCLstring &s)  
{
    s.set_null();
    value.asStr(s);
    return s.chars();
}

const char *
GenericAggrNode::STEPwrite(SCLstring &s)
{
    return value.STEPwrite(s);
/*
// CHECK do we write dollar signs for nulls within an aggregate? DAS
    if(_value)
	return (const char *)_value;
    else
	return "$";
*/
}

void
GenericAggrNode::STEPwrite (ostream& out)
{
    value.STEPwrite(out);
}

///////////////////////////////////////////////////////////////////////////////
// EntityAggregate
///////////////////////////////////////////////////////////////////////////////
 
EntityAggregate::EntityAggregate () 
{
}

EntityAggregate::~EntityAggregate ()
{
//    delete v;
}


// if exchangeFileFormat == 1 then delims are required.

Severity 
EntityAggregate::ReadValue(istream &in, ErrorDescriptor *err, 
			 const TypeDescriptor *elem_type, 
			 InstMgr *insts, int addFileId, 
			 int assignVal, int exchangeFileFormat)
{
    ErrorDescriptor errdesc;
    char errmsg[BUFSIZ];
    int value_cnt = 0;

    if(assignVal)
	Empty ();  // read new values and discard existing ones

    char c;
    int validDelims = 1;

    in >> ws; // skip white space

    c = in.peek(); // does not advance input

    if(in.eof() || (c == '$') )
    {
	_null = 1;
	err->GreaterSeverity(SEVERITY_INCOMPLETE);
	return SEVERITY_INCOMPLETE;
    }

    if(c == '(')
    {
	in.get(c);
	validDelims = 0; // signal expectation for end delim
    }
    else if(exchangeFileFormat)
    {	// error did not find opening delim
	// give up because you do not know where to stop reading.
	err->GreaterSeverity(SEVERITY_INPUT_ERROR);
	return SEVERITY_INPUT_ERROR;
    }
    else if(!in.good())
    {// this should actually have been caught by skipping white space above
	err->GreaterSeverity(SEVERITY_INCOMPLETE);
	return SEVERITY_INCOMPLETE;
    }

    EntityNode * item = 0;
    if(!assignVal)  // if not assigning values only need one node.  So only 
		    // one node is created. It is used to read the values
	item = new EntityNode();

    while (in.good() && (c != ')') )
    {
	value_cnt++;
	if(assignVal) // create a new node each time through the loop
	    item = new EntityNode();

	errdesc.ClearErrorMsg();

	if(exchangeFileFormat)
	    item->STEPread(in, &errdesc, elem_type, insts, addFileId);
	else
	    item->StrToVal(in, &errdesc, elem_type, insts, addFileId);

	// read up to the next delimiter and set errors if garbage is
	// found before specified delims (i.e. comma and quote)
	CheckRemainingInput(in, &errdesc, elem_type->AttrTypeName(), ",)");

	if (errdesc.severity() < SEVERITY_INCOMPLETE)
	{
	    sprintf (errmsg, "  index:  %d\n", value_cnt );
	    errdesc.PrependToDetailMsg(errmsg);
	    err->AppendFromErrorArg(&errdesc);
	}
	if(assignVal)
	    AddNode (item);

	in >> ws; // skip white space (although should already be skipped)
	in.get(c); // read delim

	// CheckRemainingInput should have left the input right at the delim
	// so that it would be read in in.get() above.  Since it did not find
	// the delim this does not know how to find it either!
	if( (c != ',') && (c != ')') )
	{
	    // cannot recover so give up and let STEPattribute recover
	    err->GreaterSeverity(SEVERITY_INPUT_ERROR);
	    return SEVERITY_INPUT_ERROR;
/*
	    // error read until you find a delimiter
	    SCLstring tmp;
	    while(in.good() && !strchr(",)", c) )
	    {
		in.get(c);
		tmp.Append(c);
	    }
	    // BUG could overwrite the error message buffer
	    sprintf(errmsg, "ERROR aggr. elem. followed by \'%s\' garbage.\n",
		    tmp.chars());
	    err->AppendToDetailMsg(errmsg);
	    err->AppendToUserMsg(errmsg);
	    err->GreaterSeverity(SEVERITY_WARNING);
	    if(!in.good())
		return err->severity();
*/
	}
    }
    return err->severity();
}


STEPaggregate& 
EntityAggregate::ShallowCopy (const STEPaggregate& a)  
{
    const EntityNode * tmp = (const EntityNode*) a.GetHead ();
    while (tmp) 
    {
	AddNode (new EntityNode (tmp -> node));
	tmp = (const EntityNode*) tmp -> NextNode ();
    }
    if(head)
	_null = 0;
    else
	_null = 1;

    return *this;
}


SingleLinkNode *	
EntityAggregate::NewNode ()  
{
    return new EntityNode ();
}

///////////////////////////////////////////////////////////////////////////////
// EntityNode
///////////////////////////////////////////////////////////////////////////////

EntityNode::EntityNode  (STEPentity * e) : node (e) 
{
}

SingleLinkNode *	
EntityNode::NewNode ()  
{
    return new EntityNode ();
}
///////////////////////////////////////////////////////////////////////////////

Severity 
EntityNode::StrToVal(const char *s, ErrorDescriptor *err, 
		     const TypeDescriptor *elem_type,
		     InstMgr *insts, int addFileId)
{
    STEPentity *se = ReadEntityRef(s, err, ",)", insts, addFileId);
//    STEPentity *se = ReadEntityRef(s, err, 0, insts, addFileId);
    if( se != S_ENTITY_NULL )
    {
	ErrorDescriptor error;
	if(EntityValidLevel(se, elem_type, &error) == SEVERITY_NULL)
	    node = se;
	else
	{
	    node = S_ENTITY_NULL;
	    err->AppendToDetailMsg(error.DetailMsg());
	    err->AppendToUserMsg(error.UserMsg());
	    err->GreaterSeverity(error.severity());
	}
    }
    else
	node = S_ENTITY_NULL;
    return err->severity();
}

Severity 
EntityNode::StrToVal(istream &in, ErrorDescriptor *err, 
		     const TypeDescriptor *elem_type,
		     InstMgr *insts, int addFileId)
{
    return STEPread(in, err, elem_type, insts, addFileId);
}

Severity 
EntityNode::STEPread(const char *s, ErrorDescriptor *err, 
		     const TypeDescriptor *elem_type, 
		     InstMgr *insts, int addFileId)
{
    istrstream in((char *)s);
    return STEPread(in, err, elem_type, insts, addFileId);
}

Severity 
EntityNode::STEPread(istream &in, ErrorDescriptor *err, 
		     const TypeDescriptor *elem_type,
		     InstMgr *insts, int addFileId)
{
    STEPentity *se = ReadEntityRef(in, err, ",)", insts, addFileId);
    if( se != S_ENTITY_NULL )
    {
	ErrorDescriptor error;
	if( EntityValidLevel(se, elem_type, &error) == SEVERITY_NULL )
	    node = se;
	else
	{
	    node = S_ENTITY_NULL;
	    err->AppendToDetailMsg(error.DetailMsg());
	    err->AppendToUserMsg(error.UserMsg());
	    err->GreaterSeverity(error.severity());
	}
    }
    else
	node = S_ENTITY_NULL;
//    CheckRemainingInput(in, err, "enumeration", ",)");
    return err->severity();
}

const char *
EntityNode::asStr (SCLstring &s)  
{
    s.set_null();
    if (!node || (node == S_ENTITY_NULL))     //  nothing
	return "";
    else // otherwise return entity id
    {
	char tmp [64];
	sprintf(tmp, "#%d", node->STEPfile_id);
	s = tmp;
    }
    return s.chars();
}    

const char *
EntityNode::STEPwrite(SCLstring &s)
{
    if (!node || (node == S_ENTITY_NULL) )     //  nothing
    {
	s = "$";
	return s.chars();
    }
    asStr(s);
    return s.chars();
}

void 
EntityNode::STEPwrite(ostream& out)
{
    if (!node || (node == S_ENTITY_NULL))     //  nothing
      out << "$";
    SCLstring s;
    out << asStr(s);
}


///////////////////////////////////////////////////////////////////////////////
// SelectAggregate
///////////////////////////////////////////////////////////////////////////////
 
SelectAggregate::SelectAggregate () 
{
}

SelectAggregate::~SelectAggregate ()
{
//    delete v;
}


// if exchangeFileFormat == 1 then delims are required.

Severity 
SelectAggregate::ReadValue(istream &in, ErrorDescriptor *err, 
			   const TypeDescriptor *elem_type, 
			   InstMgr *insts, int addFileId, 
			   int assignVal, int exchangeFileFormat)
{
    ErrorDescriptor errdesc;
    char errmsg[BUFSIZ];
    int value_cnt = 0;

    if(assignVal)
	Empty ();  // read new values and discard existing ones

    char c;
    int validDelims = 1;

    in >> ws; // skip white space

    c = in.peek(); // does not advance input

    if(in.eof() || (c == '$') )
    {
	_null = 1;
	err->GreaterSeverity(SEVERITY_INCOMPLETE);
	return SEVERITY_INCOMPLETE;
    }

    if(c == '(')
    {
	in.get(c);
	validDelims = 0; // signal expectation for end delim
    }
    else if(exchangeFileFormat)
    {	// error did not find opening delim
	// give up because you do not know where to stop reading.
	err->GreaterSeverity(SEVERITY_INPUT_ERROR);
	return SEVERITY_INPUT_ERROR;
    }
    else if(!in.good())
    {// this should actually have been caught by skipping white space above
	err->GreaterSeverity(SEVERITY_INCOMPLETE);
	return SEVERITY_INCOMPLETE;
    }

    SelectNode * item = 0;
    if(!assignVal)  // if not assigning values only need one node.  So only 
		    // one node is created. It is used to read the values
    item = (SelectNode *) NewNode ();

    while (in.good() && (c != ')') )
    {
	value_cnt++;
	if(assignVal) // create a new node each time through the loop
	  item = (SelectNode *) NewNode ();

	errdesc.ClearErrorMsg();

	if(exchangeFileFormat)
	    item->STEPread(in, &errdesc, elem_type, insts, addFileId);
	else
	    item->StrToVal(in, &errdesc, elem_type, insts, addFileId);

	// read up to the next delimiter and set errors if garbage is
	// found before specified delims (i.e. comma and quote)
	CheckRemainingInput(in, &errdesc, elem_type->AttrTypeName(), ",)");

	if (errdesc.severity() < SEVERITY_INCOMPLETE)
	{
	    sprintf (errmsg, "  index:  %d\n", value_cnt );
	    errdesc.PrependToDetailMsg(errmsg);
	    err->AppendFromErrorArg(&errdesc);
//	    err->AppendToDetailMsg(errdesc.DetailMsg());
//	    err->AppendToUserMsg(errdesc.UserMsg());
	}
	if(assignVal)
	    AddNode (item);

	in >> ws; // skip white space (although should already be skipped)
	in.get(c); // read delim

	// CheckRemainingInput should have left the input right at the delim
	// so that it would be read in in.get() above.  Since it did not find
	// the delim this does not know how to find it either!
	if( (c != ',') && (c != ')') )
	{
	    // cannot recover so give up and let STEPattribute recover
	    err->GreaterSeverity(SEVERITY_INPUT_ERROR);
	    return SEVERITY_INPUT_ERROR;
	}
    }
    return err->severity();
}


STEPaggregate& 
SelectAggregate::ShallowCopy (const STEPaggregate& a)  
{
    const SelectNode * tmp = (const SelectNode*) a.GetHead ();
    while (tmp) 
    {
	AddNode (new SelectNode (tmp -> node));
	tmp = (const SelectNode*) tmp -> NextNode ();
    }
    if(head)
	_null = 0;
    else
	_null = 1;

    return *this;
}


SingleLinkNode *	
SelectAggregate::NewNode ()  
{
    return new SelectNode ();
}

///////////////////////////////////////////////////////////////////////////////
// SelectNode
///////////////////////////////////////////////////////////////////////////////

SingleLinkNode *	
SelectNode::NewNode ()  
{
    return new SelectNode ();
}
///////////////////////////////////////////////////////////////////////////////

Severity 
SelectNode::StrToVal(const char *s, ErrorDescriptor *err, 
		     const TypeDescriptor *elem_type,
		     InstMgr *insts, int addFileId)
{
/*
    STEPentity *se = ReadEntityRef(s, err, ",)", insts, addFileId);
    if( se != S_ENTITY_NULL )
    {
	ErrorDescriptor error;
	if(SelectValidLevel(se, elem_type, &error) == SEVERITY_NULL)
	    node = se;
	else
	{
	    node = S_ENTITY_NULL;
	    err->AppendToDetailMsg(error.DetailMsg());
	    err->AppendToUserMsg(error.UserMsg());
	    err->GreaterSeverity(error.severity());
	}
    }
    else
	node = S_ENTITY_NULL;
*/
	// KC you will have to decide what to do here
    istrstream in((char*)s);
    if (err->severity( node->STEPread(in, err, insts) ) != SEVERITY_NULL)
	err->AppendToDetailMsg (node ->Error ());
    return err->severity();
}

Severity 
SelectNode::StrToVal(istream &in, ErrorDescriptor *err, 
		     const TypeDescriptor *elem_type,
		     InstMgr *insts, int addFileId)
{
    return STEPread(in, err, elem_type, insts, addFileId);
}

Severity 
SelectNode::STEPread(const char *s, ErrorDescriptor *err, 
		     const TypeDescriptor *elem_type, 
		     InstMgr *insts, int addFileId)
{
    istrstream in((char *)s);
    return STEPread(in, err, elem_type, insts, addFileId);
}

Severity 
SelectNode::STEPread(istream &in, ErrorDescriptor *err, 
		     const TypeDescriptor *elem_type,
		     InstMgr *insts, int addFileId)
{
  if (!node)  {
    G4cerr << "Internal error:  " << __FILE__ << ": " <<  __LINE__ << "\n" 
	 << _POC_ "\n";
    G4cerr << "function:  SelectNode::STEPread \n" << "\n";
    return SEVERITY_BUG;
  }
  err->severity( node->STEPread(in, err, insts, addFileId));
  CheckRemainingInput(in, err, "select", ",)");
  return err->severity();
}

const char *
SelectNode::asStr (SCLstring &s)  
{
    s.set_null();
    if ( !node || (node->is_null()) )     //  nothing
	return "";
    else // otherwise return entity id
    {
      node -> STEPwrite (s);
      return s.chars ();
    }
}    

const char *
SelectNode::STEPwrite(SCLstring &s)
{
    s.set_null();
    if ( !node || (node->is_null()) )     //  nothing
    {
	s = "$";
	return "$";
    }
    node -> STEPwrite (s);
    return s.chars();
}

void 
SelectNode::STEPwrite(ostream& out)
{
    if ( !node || (node->is_null()) )     //  nothing
      out << "$";
    SCLstring s;
    out << asStr(s);
}

///////////////////////////////////////////////////////////////////////////////
// StringAggregate
///////////////////////////////////////////////////////////////////////////////

/******************************************************************
STEPaggregate& 
StringAggregate::ShallowCopy (const STEPaggregate&);
******************************************************************/

STEPaggregate& 
StringAggregate::ShallowCopy (const STEPaggregate& a)
{
    Empty();
    
    SingleLinkNode* next = a.GetHead();
    SingleLinkNode* copy;

    while (next) 
    {
	copy = new StringNode (*(StringNode*)next);
	AddNode(copy);
	next = next->NextNode();
    }
    if(head)
	_null = 0;
    else
	_null = 1;
    return *this;
    
}

SingleLinkNode *
StringAggregate::NewNode () 
{
    return new StringNode ();
}

///////////////////////////////////////////////////////////////////////////////
// StringNode
///////////////////////////////////////////////////////////////////////////////

StringNode::StringNode(StringNode& sn)
{
    value = sn.value.chars();
}

StringNode::StringNode(const char * sStr)
{
    // value is an SdaiString (the memory is copied)
    value = sStr;

/*
  // I do not think that you are expecting sStr in exchange file format
    ErrorDescriptor err;
    if(value.STEPread(sStr, &err) < SEVERITY_USERMSG)
	value.set_null();
*/
}

SingleLinkNode *
StringNode::NewNode () 
{
    return new StringNode ();
}

///////////////////////////////////////////////////////////////////////////////
// non-whitespace chars following s are considered garbage and is an error.
// a valid value will still be assigned if it exists before the garbage.
///////////////////////////////////////////////////////////////////////////////

Severity 
StringNode::StrToVal(const char *s, ErrorDescriptor *err)
{
    return STEPread(s, err);
}

// this function assumes you will check for garbage following input

Severity 
StringNode::StrToVal(istream &in, ErrorDescriptor *err)
{
    return value.STEPread(in, err);
}

// non-whitespace chars following s are considered garbage and is an error.
// a valid value will still be assigned if it exists before the garbage.
Severity 
StringNode::STEPread(const char *s, ErrorDescriptor *err)
{
    istrstream in((char *)s);

    value.STEPread(in, err);
    CheckRemainingInput(in, err, "string", ",)");
    return err->severity();
}

// this function assumes you will check for garbage following input

Severity 
StringNode::STEPread(istream &in, ErrorDescriptor *err)
{
    return value.STEPread(in, err);
}

const char *
StringNode::asStr(SCLstring &s)
{
//    return value.asStr(); // this does not put quotes around the value

    value.asStr(s);
    return s.chars();
}

const char *
StringNode::STEPwrite (SCLstring &s)
{
    value.STEPwrite(s);
    return s.chars();
}

void
StringNode::STEPwrite (ostream& out)
{
    value.STEPwrite(out);
}
///////////////////////////////////////////////////////////////////////////////
// BinaryAggregate
///////////////////////////////////////////////////////////////////////////////

/******************************************************************
STEPaggregate& 
BinaryAggregate::ShallowCopy (const STEPaggregate&);
******************************************************************/

STEPaggregate& 
BinaryAggregate::ShallowCopy (const STEPaggregate& a)
{
    Empty();
    
    SingleLinkNode* next = a.GetHead();
    SingleLinkNode* copy;

    while (next) 
    {
	copy = new BinaryNode (*(BinaryNode*)next);
	AddNode(copy);
	next = next->NextNode();
    }
    if(head)
	_null = 0;
    else
	_null = 1;
    return *this;
    
}

SingleLinkNode *
BinaryAggregate::NewNode () 
{
    return new BinaryNode ();
}

///////////////////////////////////////////////////////////////////////////////
// BinaryNode
///////////////////////////////////////////////////////////////////////////////

BinaryNode::BinaryNode(BinaryNode& bn)
{
    value = bn.value.chars();
}

BinaryNode::BinaryNode(const char *sStr)
{
    // value is an SdaiBinary (the memory is copied)
    value = sStr;
}

SingleLinkNode *
BinaryNode::NewNode () 
{
    return new BinaryNode ();
}

///////////////////////////////////////////////////////////////////////////////
// non-whitespace chars following s are considered garbage and is an error.
// a valid value will still be assigned if it exists before the garbage.
///////////////////////////////////////////////////////////////////////////////

Severity 
BinaryNode::StrToVal(const char *s, ErrorDescriptor *err)
{
    return STEPread(s, err);
}

// this function assumes you will check for garbage following input

Severity 
BinaryNode::StrToVal(istream &in, ErrorDescriptor *err)
{
    return value.STEPread(in, err);
}

// non-whitespace chars following s are considered garbage and is an error.
// a valid value will still be assigned if it exists before the garbage.

Severity 
BinaryNode::STEPread(const char *s, ErrorDescriptor *err)
{
    istrstream in((char *)s);

    value.STEPread(in, err);
    CheckRemainingInput(in, err, "binary", ",)");
    return err->severity();
}

// this function assumes you will check for garbage following input

Severity 
BinaryNode::STEPread(istream &in, ErrorDescriptor *err)
{
    return value.STEPread(in, err);
}

const char *
BinaryNode::asStr(SCLstring &s)
{
    s = value.chars();
    return s.chars();
}

const char *
BinaryNode::STEPwrite (SCLstring &s)
{
    value.STEPwrite(s);
    return s.chars();
}

void
BinaryNode::STEPwrite (ostream& out)
{
    value.STEPwrite(out);
}

///////////////////////////////////////////////////////////////////////////////
// EnumAggregate
///////////////////////////////////////////////////////////////////////////////

// COPY
STEPaggregate& 
EnumAggregate::ShallowCopy (const STEPaggregate& a)
{
    const EnumNode * tmp = (const EnumNode *) a.GetHead();
    EnumNode * to;
    
    while (tmp) 
    {
	to = (EnumNode *) NewNode ();
	to -> node -> put (tmp -> node ->asInt());
	AddNode (to);
	tmp = (const EnumNode *) tmp -> NextNode ();
    }
    if(head)
	_null = 0;
    else
	_null = 1;

    return *this;
}

EnumAggregate::EnumAggregate ()
{
    
}

EnumAggregate::~EnumAggregate ()
{
    
}

/******************************************************************
 ** Procedure:  EnumAggregate::NewNode
 ** Parameters:  
 ** Returns:  a new EnumNode which is of the correct derived type
 ** Description:  creates a node to put in an List of enumerated values
 **               function is virtual so that the right node will be 
 **               created
 ** Side Effects:  
 ** Status:  ok 2/91
 ******************************************************************/

/*EnumNode **/

SingleLinkNode *
EnumAggregate::NewNode ()
{
    //  defined in subclass
  G4cerr << "Internal error:  " << __FILE__ << ": " <<  __LINE__ << "\n" ;
  G4cerr << "function:  EnumAggregate::NewNode () called instead of virtual function. \n" 
       << _POC_ << "\n";
  return 0;
}

///////////////////////////////////////////////////////////////////////////////
// EnumNode
///////////////////////////////////////////////////////////////////////////////

SingleLinkNode *
EnumNode::NewNode ()
{
    //  defined in subclass
  G4cerr << "Internal error:  " << __FILE__ << ": " <<  __LINE__ << "\n" ;
  G4cerr << "function:  EnumNode::NewNode () called instead of virtual function. \n" 
       << _POC_ << "\n";
  return 0;
}

/*
// insts and addFileId are used for the EntityNode virtual definition of this
// function
Severity
EnumNode::StrToVal(const char *s, ErrorDescriptor *err)
{
    char messageBuf[BUFSIZ];
    messageBuf[0] = '\0';

    Severity sev = SEVERITY_NULL;
    int len = strlen (s);
    char *val = new char [len + 1];
    val[0] = '\0';
    char *saveForDelete = val;

    int numFound = sscanf((char *)s," %s", val);
    if(numFound != EOF)
    {
	if(val [0] == '.')  // strip the delims
	{
	    val++;
	    char * pos = strchr(val, '.');
	    if (pos) 
		*pos = '\0';
	    else 
	    {
		sev = SEVERITY_WARNING;
		err->GreaterSeverity(SEVERITY_WARNING);
		err->AppendToDetailMsg("Missing matching \'.\' delimiter.\n");
	    }
	}
	if(elem_type)
	{
	    if(elem_type->BaseType() == BOOLEAN_TYPE)
	    {
		switch(val[0])
		{
		  case 't':
		  case 'f':
		  case 'T':
		  case 'F':
		    break;
		  default:
		    sev = SEVERITY_WARNING;
		    err->GreaterSeverity(SEVERITY_WARNING);
		    sprintf(messageBuf, "Invalid boolean value: \'%s\'.\n",
			    val);
		    err->AppendToDetailMsg(messageBuf);
		    break;
		}
	    }
	}
	    // assign based on the result of this element (the error descriptor
	    // contains the error level for the whole aggregate).
//	if(assignVal && sev > SEVERITY_WARNING)
	if(1 && sev > SEVERITY_WARNING)
	    node -> put (val);
    }
//STEPenumeration::EnumValidLevel(const char *value, ErrorDescriptor *err,
//				int optional, char *tokenList, 
//				int needDelims, int clearError)

    node->EnumValidLevel((char *)val, err, 0, 0, 0, 0);
    delete [] saveForDelete;

    // an element being null shouldn't make the aggregate incomplete!!
    if(err->severity() == SEVERITY_INCOMPLETE)
	err->severity(SEVERITY_NULL);
    return err->severity();
}
*/

///////////////////////////////////////////////////////////////////////////////
// non-whitespace chars following s are considered garbage and is an error.
// a valid value will still be assigned if it exists before the garbage.
///////////////////////////////////////////////////////////////////////////////

Severity
EnumNode::StrToVal(const char *s, ErrorDescriptor *err)
{
    return STEPread(s, err);
}

// this function assumes you will check for garbage following input

Severity 
EnumNode::StrToVal(istream &in, ErrorDescriptor *err)
{
    return node->STEPread(in, err);
}

// non-whitespace chars following s are considered garbage and is an error.
// a valid value will still be assigned if it exists before the garbage.

Severity 
EnumNode::STEPread(const char *s, ErrorDescriptor *err)
{
    istrstream in((char *)s); // sz defaults to length of s

    int nullable = 0;
    node->STEPread (in, err,  nullable);
    CheckRemainingInput(in, err, "enumeration", ",)");
    return err->severity();
}

// this function assumes you will check for garbage following input

Severity 
EnumNode::STEPread(istream &in, ErrorDescriptor *err)
{
    int nullable = 0;
    node->STEPread (in, err,  nullable);
    return err->severity();
}

const char *
EnumNode::asStr (SCLstring &s)  
{
    node -> asStr(s);
    return s.chars();
}

const char *
EnumNode::STEPwrite (SCLstring &s)
{
    node->STEPwrite(s);
    return s.chars();

/*
    static char buf[BUFSIZ];
    buf[0] = '\0';
    
//    const char *str = asStr();
//    if( (strlen(str) > 0) && (str[0] != '$') )

    if(!(node->is_null()))
    {
	buf[0] = '.';
	buf[1] = '\0';
//	strcat(buf, str);
	strcat(buf, asStr());
	strcat(buf, ".");
    }
    return buf;
 */
}

void 
EnumNode::STEPwrite (ostream& out)
{
//    out << '.' << asStr() << '.';
    node->STEPwrite(out);
}

///////////////////////////////////////////////////////////////////////////////
// Logicals
///////////////////////////////////////////////////////////////////////////////

//EnumNode * 

SingleLinkNode *
Logicals::NewNode ()  
{
    return new EnumNode (new Logical);
}	

///////////////////////////////////////////////////////////////////////////////
// Booleans
///////////////////////////////////////////////////////////////////////////////

//EnumNode * 

SingleLinkNode *
Booleans::NewNode ()  
{
    return new EnumNode (new Boolean);
}	

///////////////////////////////////////////////////////////////////////////////
// RealAggregate
///////////////////////////////////////////////////////////////////////////////

SingleLinkNode *
RealAggregate::NewNode ()  
{
    return new RealNode();
}	

// COPY
STEPaggregate& 
RealAggregate::ShallowCopy (const STEPaggregate& a)
{
    const RealNode * tmp = (const RealNode *) a.GetHead();
    RealNode * to;
    
    while (tmp) 
    {
	to = (RealNode *) NewNode ();
	to -> value = tmp -> value;
	AddNode (to);
	tmp = (const RealNode *) tmp -> NextNode ();
    }
    if(head)
	_null = 0;
    else
	_null = 1;
    return *this;
}

///////////////////////////////////////////////////////////////////////////////
// IntAggregate
///////////////////////////////////////////////////////////////////////////////

SingleLinkNode *
IntAggregate::NewNode ()  
{
    return new IntNode();
}	

// COPY
STEPaggregate& 
IntAggregate::ShallowCopy (const STEPaggregate& a)
{
    const IntNode * tmp = (const IntNode *) a.GetHead();
    IntNode * to;
    
    while (tmp) 
    {
	to = (IntNode *) NewNode ();
	to -> value = tmp -> value;
	AddNode (to);
	tmp = (const IntNode *) tmp -> NextNode ();
    }
    if(head)
	_null = 0;
    else
	_null = 1;
    return *this;
}

///////////////////////////////////////////////////////////////////////////////
// RealNode
///////////////////////////////////////////////////////////////////////////////

SingleLinkNode *
RealNode::NewNode ()  
{
    return new RealNode();
}	

Severity 
RealNode::StrToVal(const char *s, ErrorDescriptor *err)
{
    if( ReadReal(value, s, err, ",)") ) // returns true if value is assigned
	_null = 0;
    else
    {
	set_null();
	value = S_REAL_NULL;
    }
    return err->severity ();
}

Severity 
RealNode::StrToVal(istream &in, ErrorDescriptor *err)
{
    if( ReadReal(value, in, err, ",)") ) // returns true if value is assigned
	_null = 0;
    else
    {
	set_null();
	value = S_REAL_NULL;
    }
    return err->severity ();
}


Severity 
RealNode::STEPread(const char *s, ErrorDescriptor *err)
{
    if( ReadReal(value, s, err, ",)") ) // returns true if value is assigned
	_null = 0;
    else
    {
	set_null();
	value = S_REAL_NULL;
    }
    return err->severity ();
}

Severity 
RealNode::STEPread(istream &in, ErrorDescriptor *err)
{
    if( ReadReal(value, in, err, ",)") ) // returns true if value is assigned
	_null = 0;
    else
    {
	set_null();
	value = S_REAL_NULL;
    }
    return err->severity ();
}

const char *
RealNode::asStr(SCLstring &s)
{
    STEPwrite(s);
    return s.chars();
}

const char *
RealNode::STEPwrite(SCLstring &s)
{
    char tmp[BUFSIZ];
    if(value != S_REAL_NULL)
    {
//	sprintf(tmp, "%.15g", value);
	sprintf(tmp, "%.*g", Real_Num_Precision, value);
	s = tmp;
    }
    else
	s.set_null();
    return s.chars();
}

void 
RealNode::STEPwrite(ostream& out)
{
    SCLstring s;
    out << STEPwrite(s);
}

///////////////////////////////////////////////////////////////////////////////
// IntNode
///////////////////////////////////////////////////////////////////////////////

SingleLinkNode *
IntNode::NewNode ()  
{
    return new IntNode();
}	

Severity 
IntNode::StrToVal(const char *s, ErrorDescriptor *err)
{
    if( ReadInteger(value, s, err, ",)") ) // returns true if value is assigned
	_null = 0;
    else
    {
	set_null();
	value = S_INT_NULL;
    }
    return err->severity ();
}

Severity 
IntNode::StrToVal(istream &in, ErrorDescriptor *err)
{
    if( ReadInteger(value, in, err, ",)") )// returns true if value is assigned
	_null = 0;
    else
    {
	set_null();
	value = S_INT_NULL;
    }
    return err->severity ();
}

Severity 
IntNode::STEPread(const char *s, ErrorDescriptor *err)
{
    if( ReadInteger(value, s, err, ",)") ) // returns true if value is assigned
	_null = 0;
    else
    {
	set_null();
	value = S_INT_NULL;
    }
    return err->severity ();
}

Severity 
IntNode::STEPread(istream &in, ErrorDescriptor *err)
{
    if( ReadInteger(value, in, err, ",)") ) // returns true if value is assigned
	_null = 0;
    else
    {
	set_null();
	value = S_INT_NULL;
    }
    return err->severity ();
}

const char *
IntNode::asStr(SCLstring &s)
{
    STEPwrite(s);
    return s.chars();
}

const char *
IntNode::STEPwrite(SCLstring &s)
{
    char tmp[BUFSIZ];
    if(value != S_INT_NULL)
    {
	sprintf(tmp, "%d", value);
	s = tmp;
    }
    else
	s.set_null();
    return s.chars();
}

void 
IntNode::STEPwrite(ostream& out)
{
    SCLstring s;
    out << STEPwrite(s);
}
