

//



//
// $Id: STEPcomplex.cc,v 1.3 1999-12-15 14:50:16 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include <ctype.h>

#include <STEPcomplex.h>

extern const char *
ReadStdKeyword(G4std::istream& in, SCLstring &buf, int skipInitWS);


STEPcomplex::STEPcomplex(Registry *registry, int fileid)
: STEPentity(fileid, 1),  sc(0), _registry(registry), visited(0)
{
    head = this;
/*
    _complex = 1;
    _registry = registry;
    sc = 0;
    visited = 0;
*/
}

STEPcomplex::STEPcomplex(Registry *registry, const SCLstring **names, 
			 int fileid) 
: STEPentity(fileid, 1),  sc(0), _registry(registry), visited(0)
{
    head = this;
/*
    _complex = 1;
    _registry = registry;
    sc = 0;
    visited = 0;
*/

    if(names[0])
    {
	BuildAttrs( names[0]->chars() );
    }
    int i = 1;
    while(!eDesc && names[i])
    {
	// at least 1 entity part did not have a valid name
	_error.GreaterSeverity(SEVERITY_INCOMPLETE);
	BuildAttrs( names[i]->chars() );
	i++;
    }
    if(!eDesc) // no entity part had a valid name
	_error.GreaterSeverity(SEVERITY_WARNING);

    while(names[i])
    {
	AddEntityPart( names[i]->chars() );
	i++;
    }
    AssignDerives();
}

STEPcomplex::STEPcomplex(Registry *registry, const char **names, 
			 int fileid)
: STEPentity(fileid, 1),  sc(0), _registry(registry), visited(0)
{
    head = this;
/*
    _complex = 1;
    _registry = registry;
    sc = 0;
    visited = 0;
*/

    if(names[0])
    {
	BuildAttrs(names[0]);
    }
    int i = 1;
    while(!eDesc && names[i])
    {
	// at least 1 entity part did not have a valid name
	_error.GreaterSeverity(SEVERITY_INCOMPLETE);
	BuildAttrs( names[i] );
	i++;
    }
    if(!eDesc) // no entity part had a valid name
	_error.GreaterSeverity(SEVERITY_WARNING);

    while(names[i])
    {
	AddEntityPart( names[i] );
	i++;
    }
    AssignDerives();
}

STEPcomplex::~STEPcomplex()
{
    if(sc)
	delete sc;
}

void 
STEPcomplex::AssignDerives()
{
    const char *nm;
    STEPattribute * a = 0;
    STEPcomplex *scomp1 = head;
    STEPcomplex *scomp2;

    const AttrDescriptorList *attrList;
    AttrDescLinkNode *attrPtr;
    const AttrDescriptor *ad;

    // find out how many attrs there are
//	int attrCount = attrList->EntryCount();

    while(scomp1 && scomp1->eDesc)
    {
	a = 0;
	attrList = &( scomp1->eDesc->ExplicitAttr() );
	attrPtr = (AttrDescLinkNode *)attrList->GetHead();

	// assign nm to be derived attr
	// while( more derived attr for entity part )
	while( attrPtr != 0 )
	{
	    ad = attrPtr->AttrDesc();
	    if( (LOGICAL)( ad->Derived() ) == sdaiTRUE)
	    {
		const char *nm = ad->Name();
		const char *attrNm = 0;
		if(strrchr(nm,'.'))
		{
		    attrNm = strrchr(nm,'.');
		    attrNm++;
		}
		else
		    attrNm = nm;
		scomp2 = head;
	        while(scomp2 && !a)
		{
		    if(scomp1 != scomp2)
		    {
//			scomp2->MakeDerived ( ad->Name() );
			scomp2->MakeDerived ( attrNm );
//			a = scomp2->GetSTEPattribute( ad->Name() );
			a = scomp2->GetSTEPattribute( attrNm );
		    }
		    scomp2 = scomp2->sc;
		}
//		if (a)  a ->Derive ();
	    }
	    // Increment attr
	    attrPtr = (AttrDescLinkNode *)attrPtr->NextNode();
	}
	scomp1 = scomp1->sc;
    }
}

// this function should only be called for the head entity
// in the List of entity parts.

void 
STEPcomplex::AddEntityPart(const char *name)
{
    STEPcomplex *scomplex;

    if(name)
    {
	scomplex = new STEPcomplex(_registry, STEPfile_id);
	scomplex->BuildAttrs(name);
	if(scomplex->eDesc)
	{
	    scomplex->head = this;
//	    scomplex->STEPfile_id = STEPfile_id;
	    AppendEntity(scomplex);
	}
	else
	{
	    G4cout << scomplex->_error.DetailMsg() << G4endl;
	    delete scomplex;
	}
    }
}

STEPcomplex *
STEPcomplex::EntityPart(const char *name)
{
    STEPcomplex *scomp = head;
    SCLstring s1, s2;
    while(scomp)
    {	
	if(scomp->eDesc)
	{
	    if( !strcmp( StrToUpper(name, s1), 
			 StrToUpper(scomp->eDesc->Name(), s2) ) )
		return scomp;
	}
	else
	    G4cout << "Bug in STEPcomplex::EntityPart(): entity part has "
		 << "no EntityDescriptor\n";
	scomp = scomp->sc;
    }
    return 0;
}

int 
STEPcomplex::EntityExists(const char *name)
{
    return (EntityPart(name) ? 1 : 0);
}


Severity 
STEPcomplex::ValidLevel(ErrorDescriptor *error, InstMgr *im, 
			int clearError)
{
    G4cout << "STEPcomplex::ValidLevel() not implemented.\n";
    return SEVERITY_NULL;
}

void 
STEPcomplex::AppendEntity(STEPcomplex *stepc)
{
    if(sc)
	sc->AppendEntity(stepc);
    else
	sc = stepc;
}

// READ
Severity 
STEPcomplex::STEPread(int id, int addFileId, class InstMgr * instance_set,
		 G4std::istream& in)
{
    char c;
    SCLstring typeNm;
    STEPcomplex *stepc = 0;

    ClearError(1);
    STEPfile_id = id;
    
    stepc = head;
    while(stepc)
    {
	stepc->visited = 0;
	stepc = stepc->sc;
    }

    in >> G4std::ws;
    in.get(c);
    if(c == '(') // opening paren for subsuperRecord
    {
	in >> G4std::ws;
	c = in.peek();
	while(c != ')')
	{
	    typeNm.set_null();
	    in >> G4std::ws;
	    ReadStdKeyword(in, typeNm, 1); // read the type name
	    in >> G4std::ws;
	    c = in.peek();
	    if(c != '(')
	    {
		_error.AppendToDetailMsg("Missing open paren before entity attr values.\n");
		G4cout << "ERROR: missing open paren\n";
		_error.GreaterSeverity(SEVERITY_INPUT_ERROR);
		STEPread_error(c,0,in);
		return _error.severity();
	    }

	    stepc = EntityPart(typeNm.chars());
	    if(stepc)
		stepc->STEPentity::STEPread(id, addFileId, instance_set, in);
	    else
	    {
		G4cout << "ERROR: complex entity part does not exist.\n";
		_error.AppendToDetailMsg("Complex entity part of instance does not exist.\n");
		G4cout << "ERROR: missing open paren\n";
		_error.GreaterSeverity(SEVERITY_INPUT_ERROR);
		STEPread_error(c,0,in);
		return _error.severity();
	    }
	    in >> G4std::ws;
	    c = in.peek();
	}
	if(c != ')')
	    G4cout << "ERROR: missing ending paren for complex entity instance.\n";
	else
	    in.get(c); // read the closing paren
    }
    return _error.severity();
}

#ifdef buildwhileread
// READ
Severity 
STEPcomplex::STEPread(int id, int addFileId, class InstMgr * instance_set,
		 G4std::istream& in)
{
    ClearError(1);
    STEPfile_id = id;
    
    STEPcomplex stepc = head;
    while(stepc)
    {
	stepc->visited = 0;
	stepc = stepc->sc;
    }

    char c;
    in >> G4std::ws;
    in.get(c);
    if(c == '(')
    {
	SCLstring s;
	in >> G4std::ws;
	in.get(c);
	while( in && (c != '(') && !isspace(c) ) // get the entity name
	{
	    s.Append(c);
	    in.get(c);
	}
	if(isspace(c))
	{
	    in >> G4std::ws;
	    in.get(c);
	}
//    STEPcomplex *EntityPart(const char *name);

	if(c != '(')
	{
	    _error.AppendToDetailMsg(
				     "Missing open paren before entity attr values.\n");
	    G4cout << "ERROR: missing open paren\n";
	    _error.GreaterSeverity(SEVERITY_INPUT_ERROR);
	    STEPread_error(c,0,in);
	    return _error.severity();
	}
	else // c == '('
	    in.putback(c);
	
	G4cout << s.chars() << G4endl;
	BuildAttrs( s.chars() );
	STEPentity::STEPread(id, addFileId, instance_set, in);
	
	in >> G4std::ws;
	in.get(c);
	while(c != ')')
	{
	    s.set_null();
	    while( in && (c != '(') && !isspace(c) ) // get the entity name
	    {
		s.Append(c);
		in.get(c);
	    }
	    if(isspace(c))
	    {
		in >> G4std::ws;
		in.get(c);
	    }
	    if(c != '(')
	    {
		_error.AppendToDetailMsg(
					 "Missing open paren before entity attr values.\n");
		G4cout << "ERROR: missing open paren\n";
		_error.GreaterSeverity(SEVERITY_INPUT_ERROR);
		STEPread_error(c,0,in);
		return _error.severity();
	    }
	    else // c == '('
		in.putback(c);

	    G4cout << s.chars() << G4endl; // diagnostics DAS
	    
	    STEPcomplex *stepc = new STEPcomplex( _registry );
	    AppendEntity(stepc);
	    stepc->BuildAttrs( s.chars() );
	    stepc->STEPentity::STEPread(id, addFileId, instance_set, in);
	    in >> G4std::ws;
	    in.get(c);
	}
    }
    return _error.severity();
}

#endif

void 
STEPcomplex::BuildAttrs(const char *s )
{
    // assign inherited member variable
    eDesc = (class EntityDescriptor *)_registry->FindEntity(s);

    if(eDesc)
    {
	const AttrDescriptorList *attrList;
	attrList = &( eDesc->ExplicitAttr() );

      //////////////////////////////////////////////
      // find out how many attrs there are
      //////////////////////////////////////////////
	int attrCount = attrList->EntryCount();

	STEPattribute * a = 0;

	AttrDescLinkNode *attrPtr = (AttrDescLinkNode *)attrList->GetHead();
	while( attrPtr != 0)
	{
	    const AttrDescriptor *ad = attrPtr->AttrDesc();

	    if( (LOGICAL)( ad->Derived() ) != sdaiTRUE)
	    {

		switch ( ad->NonRefType() )
		{
		  case INTEGER_TYPE:
		    a = new STEPattribute (*ad,  new SdaiInteger);
		    break;

		  case STRING_TYPE:
		    a = new STEPattribute (*ad,  new SdaiString);
		    break;

		  case BINARY_TYPE:
		    a = new STEPattribute (*ad,  new SdaiBinary);
		    break;

		  case REAL_TYPE:
		    a = new STEPattribute (*ad,  new SdaiReal);
		    break;

		  case BOOLEAN_TYPE:
		    a = new STEPattribute (*ad,  new Boolean);
		    break;

		  case LOGICAL_TYPE:
		    a = new STEPattribute (*ad,  new Logical);
		    break;

		  case ENTITY_TYPE:
		    a = new STEPattribute (*ad,  new (STEPentity *) );
		    break;

		  case ENUM_TYPE:
		  {
		    EnumTypeDescriptor * enumD = 
				(EnumTypeDescriptor *)ad->ReferentType();
		    a = new STEPattribute (*ad,  enumD->CreateEnum() );
		    break;
		  }
		  case SELECT_TYPE:
		  {
		    SelectTypeDescriptor * selectD = 
				(SelectTypeDescriptor *)ad->ReferentType();
		    a = new STEPattribute (*ad,  selectD->CreateSelect() );
		    break;
		  }
		  case AGGREGATE_TYPE:
		  case ARRAY_TYPE:		// DAS
		  case BAG_TYPE:		// DAS
		  case SET_TYPE:		// DAS
		  case LIST_TYPE:		// DAS
		  {
		    AggrTypeDescriptor * aggrD = 
				(AggrTypeDescriptor *)ad->ReferentType();
		    a = new STEPattribute (*ad,  aggrD->CreateAggregate() );
		    break;
		  }
		}

		a -> set_null ();
		attributes.push (a);
	    }

/*	    // for when inverse information is included
	    else if( (LOGICAL)( ad->Inverse() ) == sdaiTRUE)
	    {
		str.Append('(');
		endchar = ')';
	    }
*/
	    attrPtr = (AttrDescLinkNode *)attrPtr->NextNode();
	}
    }
    else
    {
	_error.AppendToDetailMsg("Entity does not exist.\n");
	_error.GreaterSeverity(SEVERITY_INPUT_ERROR);
    }
}

void 
STEPcomplex::STEPread_error(char c, int index, G4std::istream& in)
{
    G4cout << "STEPcomplex::STEPread_error() \n";
}

// WRITE
void 
STEPcomplex::STEPwrite(G4std::ostream& out, int writeComment)
{
    if(writeComment && p21Comment && !p21Comment->is_null() )
	out << p21Comment->chars();
    out << "#" << STEPfile_id << "=(";
    WriteExtMapEntities(out);
    out << ");\n";
}

const char * 
STEPcomplex::STEPwrite(SCLstring &buf)
{
    buf.set_null();

    buf.Append('#');
    buf.Append(STEPfile_id);
    buf.Append('=');
/*
    char instanceInfo[BUFSIZ];
    sprintf(instanceInfo, "#%d=", STEPfile_id );
*/
/*
    G4std::strstream ss;
    ss << "#" << STEPfile_id << "=(";
    WriteExtMapEntities(ss);
    ss << ");";
    ss << G4std::ends;

    char *tmpstr = ss.str();
    buf.Append(tmpstr);
    delete tmpstr;
*/
    WriteExtMapEntities(buf);
    buf.Append(");");

    return buf.chars();
}

void 
STEPcomplex::WriteExtMapEntities(G4std::ostream& out)
{
    SCLstring tmp;
    out << StrToUpper (EntityName(), tmp);
    out << "(";
    int n = attributes.list_length();

    for (int i = 0 ; i < n; i++) {
	(attributes[i]).STEPwrite(out);
	if (i < n-1) out << ",";
    }
    out << ")";
    if(sc)
    {
	sc->WriteExtMapEntities(out);
    }
}

const char * 
STEPcomplex::WriteExtMapEntities(SCLstring &buf)
{
    char instanceInfo[BUFSIZ];
    
    SCLstring tmp;
//    sprintf(instanceInfo, "%s(", (char *)StrToUpper( EntityName(), tmp ) );
//    buf.Append(instanceInfo);

    buf.Append( (char *)StrToUpper(EntityName(),tmp) );
    buf.Append( '(' );

    int n = attributes.list_length();

    for (int i = 0 ; i < n; i++) {
	attributes[i].asStr(tmp) ;
	buf.Append (tmp);
	if (i < n-1) {
	    buf.Append( ',' );
	}
    }    
    buf.Append( ")" );

    if(sc)
    {
	sc->WriteExtMapEntities(buf);
    }

    return buf.chars();
}

void 
STEPcomplex::CopyAs (STEPentity *)
{
    G4cout << "ERROR: STEPcomplex::CopyAs() not implemented.\n";
}
