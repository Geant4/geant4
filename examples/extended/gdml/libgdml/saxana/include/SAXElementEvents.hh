#ifndef SAX_ELEMENT_EVENTS_H
#define SAX_ELEMENT_EVENTS_H 1

#include <string>

#include "SAXEvent.hh"
#include "ASCIIAttributeList.hh"

class SAXEventStartElement : virtual public SAXEvent
{
public:
  SAXEventStartElement( const std::string& name, const ASCIIAttributeList& attrs )
  {
    fName       = name;
    fAttributes = attrs;
  }
  SAXEventStartElement( const char* name, const ASCIIAttributeList& attrs )
  : fName( name ), fAttributes( attrs )
  {
  }
  ~SAXEventStartElement()
  {
  }
  
  virtual const SAXEvent::EventType Type() const
  {
    return SAXEvent::eStartElement;
  }
  
  const std::string& Name() const
  {
    return fName;
  }
  
  const ASCIIAttributeList& Attributes() const
  {
    return fAttributes;
  }
  
private:
  std::string        fName;
  ASCIIAttributeList fAttributes;
};

class SAXEventEndElement : virtual public SAXEvent
{
public:
  SAXEventEndElement( const std::string& name )
  {
    fName = name;
  }
  SAXEventEndElement( const char* name )
  : fName( name )
  {
  }
  ~SAXEventEndElement()
  {
  }
  
  virtual const SAXEvent::EventType Type() const
  {
    return SAXEvent::eEndElement;
  }
  
  const std::string& Name() const
  {
    return fName;
  }

private:
  std::string fName;
};

class SAXEventCharacters : virtual public SAXEvent
{
public:
  SAXEventCharacters( const std::string& data )
  {
    fData = data;
  }
  SAXEventCharacters( const char* data )
  : fData( data )
  {
  }
  ~SAXEventCharacters()
  {
  }
  
  virtual const SAXEvent::EventType Type() const
  {
    return SAXEvent::eCharacters;
  }
  
  const std::string& Data() const
  {
    return fData;
  }

private:
  std::string fData;
};

class SAXEventPI : virtual public SAXEvent
{
public:
  SAXEventPI( const std::string& target, const std::string& data )
  {
    fTarget = target; fData = data;
  }
  SAXEventPI( const char* target, const char* data )
  : fTarget( target), fData( data )
  {
  }
  ~SAXEventPI()
  {
  }
  
  virtual const SAXEvent::EventType Type() const
  {
    return SAXEvent::ePI;
  }

  const std::string& Target() const
  {
    return fTarget;
  }
  const std::string& Data() const
  {
    return fData;
  }

private:
  std::string fTarget;
  std::string fData;
};

#endif // SAX_ELEMENT_EVENTS_H

