#ifndef SAX_ERROR_EVENTS_H
#define SAX_ERROR_EVENTS_H 1

#include <string>

#include "SAXEvent.hh"

class SAXErrorEventBase
{
public:
  SAXErrorEventBase( const std::string& msg )
  {
    fMessage = msg;
  }
  SAXErrorEventBase( const char* msg ) : fMessage( msg )
  {
  }
  ~SAXErrorEventBase()
  {
  }
  
  const std::string& Message() const
  {
    return fMessage;
  }
  
private:
  std::string fMessage;
};

class SAXEventWarning : virtual public SAXEvent, public SAXErrorEventBase
{
public:
  SAXEventWarning( const std::string& msg ) : SAXErrorEventBase( msg )
  {
  }
  SAXEventWarning( const char* msg ) : SAXErrorEventBase( msg )
  {
  }
  ~SAXEventWarning()
  {
  }
  
  virtual const SAXEvent::EventType Type() const
  {
    return SAXEvent::eWarning;
  }
};

class SAXEventError : virtual public SAXEvent, public SAXErrorEventBase
{
public:
  SAXEventError( const std::string& msg ) : SAXErrorEventBase( msg )
  {
  }
  SAXEventError( const char* msg ) : SAXErrorEventBase( msg )
  {
  }
  ~SAXEventError()
  {
  }
  
  virtual const SAXEvent::EventType Type() const
  {
    return SAXEvent::eError;
  }
};

class SAXEventFatalError : virtual public SAXEvent, public SAXErrorEventBase
{
public:
  SAXEventFatalError( const std::string& msg ) : SAXErrorEventBase( msg )
  {
  }
  SAXEventFatalError( const char* msg ) : SAXErrorEventBase( msg )
  {
  }
  ~SAXEventFatalError()
  {
  }
  
  virtual const SAXEvent::EventType Type() const
  {
    return SAXEvent::eFatalError;
  }
};

#endif // SAX_ERROR_EVENTS_H

