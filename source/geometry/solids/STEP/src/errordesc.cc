// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: errordesc.cc,v 1.1 1999-01-07 16:08:18 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

/*
* NIST Utils Class Library
* clutils/errordesc.cc
* February, 1993
* David Sauder
* K. C. Morris

* Development of this software was funded by the United States Government,
* and is not subject to copyright.
*/

/*   */ 

#include <errordesc.h>

DebugLevel ErrorDescriptor::_debug_level = DEBUG_OFF;
ostream *  ErrorDescriptor::_out = 0;

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
	_userMsg = new SCLstring;
    *_userMsg = msg;
}

void  
ErrorDescriptor::PrependToUserMsg ( const char * msg)  
{
    if(!_userMsg)
	_userMsg = new SCLstring;
    _userMsg -> Prepend (msg);
}

void  
ErrorDescriptor::AppendToUserMsg ( const char * msg)  
{
    if(!_userMsg)
	_userMsg = new SCLstring;
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
	_detailMsg = new SCLstring;
    *_detailMsg = msg;
    // G4cerr << "D " << _detailMsg->chars() << '\n';
}

void
ErrorDescriptor::PrependToDetailMsg (const char * msg)  
{
    if(!_detailMsg)
	_detailMsg = new SCLstring;
    _detailMsg -> Prepend (msg);
}    

void
ErrorDescriptor::AppendToDetailMsg (const char * msg)  
{
    if(!_detailMsg)
	_detailMsg = new SCLstring;
    _detailMsg -> Append (msg);
}    
