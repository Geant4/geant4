
/*
* NIST Utils Class Library
* clutils/errordesc.cc
* February, 1993
* David Sauder
* K. C. Morris

* Development of this software was funded by the United States Government,
* and is not subject to copyright.
*/

/* $Id: errordesc.cc,v 1.4 2000-01-21 13:43:10 gcosmo Exp $  */ 

#include <errordesc.h>
#include <Str.h>

DebugLevel ErrorDescriptor::_debug_level = DEBUG_OFF;
G4std::ostream *  ErrorDescriptor::_out = 0;

void 
ErrorDescriptor::PrintContents(G4std::ostream &out) const
{
    SCLstring s;
    out << "Severity: " << severity(s) << G4endl;
    if(_userMsg)
    {
	out << "User message in parens:" << G4endl << "(";
	out << UserMsg() << ")" << G4endl;
    }
    if(_detailMsg)
    {
	out << "Detailed message in parens:" << G4endl << "(";
	out << DetailMsg() << ")" << G4endl;
    }
}

const char *
ErrorDescriptor::severity(SCLstring &s) const
{
    s.set_null();
    switch(severity())
    {
      case SEVERITY_NULL :
      {
	  s = "SEVERITY_NULL";
	  break;
      }
      case SEVERITY_USERMSG :
      {
	  s = "SEVERITY_USERMSG";
	  break;
      }
      case SEVERITY_INCOMPLETE :
      {
	  s = "SEVERITY_INCOMPLETE";
	  break;
      }
      case SEVERITY_WARNING :
      {
	  s = "SEVERITY_WARNING";
	  break;
      }
      case SEVERITY_INPUT_ERROR :
      {
	  s = "SEVERITY_INPUT_ERROR";
	  break;
      }
      case SEVERITY_BUG :
      {
	  s = "SEVERITY_BUG";
	  break;
      }
      case SEVERITY_EXIT :
      {
	  s = "SEVERITY_EXIT";
	  break;
      }
      case SEVERITY_DUMP :
      {
	  s = "SEVERITY_DUMP";
	  break;
      }
      case SEVERITY_MAX :
      {
	  s = "SEVERITY_MAX";
	  break;
      }
    }
    return s.chars();
}


Severity 
ErrorDescriptor::GetCorrSeverity(const char *s)
{
//    G4cout << "s is (" << s << ") \n";
    if(s && s[0] != 0)
    {
	SCLstring s2;
	StrToUpper(s,s2);
//	G4cout << "s after if is (" << s << ") \n" << "s2 is (" << s2.chars() << ")\n";
	if(!strcmp(s2.chars(),"SEVERITY_NULL"))
	{
//	    G4cout << "SEVERITY_NULL" << G4endl;
	    return SEVERITY_NULL;
	}
	if(!strcmp(s2.chars(),"SEVERITY_USERMSG"))
	{
//	    G4cout << "SEVERITY_USERMSG" << G4endl;
	    return SEVERITY_USERMSG;
	}
	if(!strcmp(s2.chars(),"SEVERITY_INCOMPLETE"))
	{
//	    G4cout << "SEVERITY_INCOMPLETE" << G4endl;
	    return SEVERITY_INCOMPLETE;
	}
	if(!strcmp(s2.chars(),"SEVERITY_WARNING"))
	{
//	    G4cout << "SEVERITY_WARNING" << G4endl;
	    return SEVERITY_WARNING;
	}
	if(!strcmp(s2.chars(),"SEVERITY_INPUT_ERROR"))
	{
//	    G4cout << "SEVERITY_INPUT_ERROR" << G4endl;
	    return SEVERITY_INPUT_ERROR;
	}
	if(!strcmp(s2.chars(),"SEVERITY_BUG"))
	{
//	    G4cout << "SEVERITY_BUG" << G4endl;
	    return SEVERITY_BUG;
	}
	if(!strcmp(s2.chars(),"SEVERITY_EXIT"))
	{
//	    G4cout << "SEVERITY_EXIT" << G4endl;
	    return SEVERITY_EXIT;
	}
	if(!strcmp(s2.chars(),"SEVERITY_DUMP"))
	{
//	    G4cout << "SEVERITY_DUMP" << G4endl;
	    return SEVERITY_DUMP;
	}
	if(!strcmp(s2.chars(),"SEVERITY_MAX"))
	{
//	    G4cout << "SEVERITY_MAX" << G4endl;
	    return SEVERITY_MAX;
	}
    }
    G4cerr << "Internal error:  " << __FILE__ <<  __LINE__
      << "\n" << _POC_ "\n";
    G4cerr << "Calling ErrorDescriptor::GetCorrSeverity() with null string\n";
    return SEVERITY_BUG;
}

ErrorDescriptor::ErrorDescriptor ( Severity s,  DebugLevel d)
:     _severity (s), 
      _userMsg (0),
      _detailMsg (0)
{
  if (d  != DEBUG_OFF) 
    _debug_level = d;
}

/*
ErrorDescriptor::ErrorDescriptor ( Severity s,  Enforcement e,  DebugLevel d)
:     _severity (s), 
      _enforcement_level (e), 
      _userMsg (0),
      _detailMsg (0)
{
  if (d  != DEBUG_OFF) 
    _debug_level = d;
}
*/
 
const char *
ErrorDescriptor::UserMsg () const
{
    if(_userMsg)
	return _userMsg->chars();
    else
	return "";
}

void
ErrorDescriptor::UserMsg ( const char * msg)  
{
    if(!_userMsg)
#ifdef __OSTORE__
	_userMsg = new (os_database::of(this), 
			SCLstring::get_os_typespec()) SCLstring;
#else
	_userMsg = new SCLstring;
#endif
    *_userMsg = msg;
}

void  
ErrorDescriptor::PrependToUserMsg ( const char * msg)  
{
    if(!_userMsg)
#ifdef __OSTORE__
	_userMsg = new (os_database::of(this), 
			SCLstring::get_os_typespec()) SCLstring;
#else
	_userMsg = new SCLstring;
#endif
    _userMsg -> Prepend (msg);
}

void 
ErrorDescriptor::AppendToUserMsg ( const char c)
{
    if(!_userMsg) 
#ifdef __OSTORE__
	_userMsg = new (os_database::of(this), 
			SCLstring::get_os_typespec()) SCLstring;
#else
	_userMsg = new SCLstring;
#endif
    _userMsg -> Append(c);
}    

void  
ErrorDescriptor::AppendToUserMsg ( const char * msg)  
{
    if(!_userMsg)
#ifdef __OSTORE__
	_userMsg = new (os_database::of(this), 
			SCLstring::get_os_typespec()) SCLstring;
#else
	_userMsg = new SCLstring;
#endif
    _userMsg -> Append (msg);
}

 
const char *
ErrorDescriptor::DetailMsg ()  const
{
    if(_detailMsg)
	return _detailMsg->chars();
    else
	return "";
}

void
ErrorDescriptor::DetailMsg ( const char * msg)  
{
    if(!_detailMsg)
#ifdef __OSTORE__
	_detailMsg = new (os_database::of(this), 
			SCLstring::get_os_typespec()) SCLstring;
#else
	_detailMsg = new SCLstring;
#endif
    *_detailMsg = msg;
    // G4cerr << "D " << _detailMsg->chars() << '\n';
}

void
ErrorDescriptor::PrependToDetailMsg (const char * msg)  
{
    if(!_detailMsg)
#ifdef __OSTORE__
	_detailMsg = new (os_database::of(this), 
			SCLstring::get_os_typespec()) SCLstring;
#else
	_detailMsg = new SCLstring;
#endif
    _detailMsg -> Prepend (msg);
}    

void 
ErrorDescriptor::AppendToDetailMsg ( const char c)
{
    if(!_detailMsg)
#ifdef __OSTORE__
	_detailMsg = new (os_database::of(this), 
			SCLstring::get_os_typespec()) SCLstring;
#else
	_detailMsg = new SCLstring;
#endif
    _detailMsg -> Append(c);
}    

void
ErrorDescriptor::AppendToDetailMsg (const char * msg)  
{
    if(!_detailMsg)
#ifdef __OSTORE__
	_detailMsg = new (os_database::of(this), 
			SCLstring::get_os_typespec()) SCLstring;
#else
	_detailMsg = new SCLstring;
#endif
    _detailMsg -> Append (msg);
}    
