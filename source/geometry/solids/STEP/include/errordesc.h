

//



//
// $Id: errordesc.h,v 1.2 1999-05-21 20:20:40 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef ERRORDESC_H
#define	ERRORDESC_H

/*
* NIST Utils Class Library
* clutils/errordesc.h
* May 1995
* David Sauder
* K. C. Morris

* Development of this software was funded by the United States Government,
* and is not subject to copyright.
*/

/*   */ 

#ifdef __O3DB__
#include <OpenOODB.h>
#endif

#include <scl_string.h>
#include "G4ios.hh"

typedef enum Severity {
    SEVERITY_MAX	= -5,
    SEVERITY_DUMP	= -4,
    SEVERITY_EXIT	= -3,	// fatal
    SEVERITY_BUG	= -2,	// non-recoverable error -- probably bug
    SEVERITY_INPUT_ERROR = -1,	// non-recoverable error
    SEVERITY_WARNING	= 0,	// recoverable error
    SEVERITY_INCOMPLETE	= 1,	// incomplete data
    SEVERITY_USERMSG	= 2,	// possibly an error
    SEVERITY_NULL	= 3	// no error or message
  } Severity;

/* DAS not used
enum Enforcement {
    ENFORCE_OFF,
    ENFORCE_OPTIONALITY,
    ENFORCE_ALL
    
    } ;
*/    

enum  DebugLevel  {
    DEBUG_OFF	=0,
    DEBUG_USR	=1,
    DEBUG_ALL	=2
    
};

/******************************************************************
 ** Class:  ErrorDescriptor
 ** Data Members:  
 **	severity level of error
 **	user message
 **	detailed message
 ** Description:  
 **	the error is a detailed error message + a severity level
 **	also keeps a user message separately
 **	detailed message gets sent to ostream
 **	uses SCLstring class to keep the user messages
 **	keeps severity of error
 **	created with or without error
 ** Status:  
 ******************************************************************/

class ErrorDescriptor {
//  friend     istream &operator<< ( istream&, ErrorDescriptor& );
  protected:
  
    Severity	_severity;

    static DebugLevel	_debug_level;
    static ostream* _out;
    
    SCLstring *_userMsg;
    SCLstring *_detailMsg;
  public:

    ErrorDescriptor (Severity s    = SEVERITY_NULL, 
			  DebugLevel d  = DEBUG_OFF);
    ~ErrorDescriptor () { delete _userMsg; delete _detailMsg; }

    void ClearErrorMsg() {
	_severity = SEVERITY_NULL;
	delete _userMsg;   _userMsg = 0;
	delete _detailMsg; _detailMsg = 0;
    }

    Severity severity() const        { return _severity; }
    Severity severity(Severity s) {  return (_severity = s); }
    Severity GreaterSeverity(Severity s)
	{ return ((s < _severity) ?  _severity = s : _severity); }

    const char * UserMsg () const;  //  UserMsg is from String
    void UserMsg ( const char *);
    void AppendToUserMsg ( const char *);
    void PrependToUserMsg ( const char *);
    void AppendToUserMsg ( const char c)
    {
      if(!_userMsg) _userMsg = new SCLstring;
      _userMsg -> Append(c);
    }    

    const char * DetailMsg () const;
    void DetailMsg ( const char *);
    void AppendToDetailMsg ( const char *);
    void PrependToDetailMsg ( const char *);
    void AppendToDetailMsg ( const char c)
    {
      if(!_detailMsg) _detailMsg = new SCLstring;
      _detailMsg -> Append(c);
    }    

    Severity AppendFromErrorArg(ErrorDescriptor *err)
    {
	GreaterSeverity( err->severity() );
	AppendToDetailMsg( err->DetailMsg() );
	AppendToUserMsg( err->UserMsg() );
	return severity();
    }

    DebugLevel debug_level() const        { return _debug_level; }
    void debug_level(DebugLevel d)   { _debug_level = d; }
    void SetOutput(ostream *o)      { _out = o; }
    
    // void operator =  (const char* y)  {  SCLstring::operator = (y);  }
    
/*
//    Enforcement	_enforcement_level;	
    ErrorDescriptor (Severity s    = SEVERITY_NULL, 
			  Enforcement e = ENFORCE_OFF,
			  DebugLevel d  = DEBUG_OFF);
    Enforcement enforcement() const       { return _enforcement_level; }
    void enforcement(Enforcement e) { _enforcement_level = e; }
*/
} ;

#endif
