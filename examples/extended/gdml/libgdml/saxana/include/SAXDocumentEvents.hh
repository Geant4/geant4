#ifndef SAX_DOCUMENT_EVENTS_H
#define SAX_DOCUMENT_EVENTS_H 1

#include "SAXEvent.hh"

class SAXEventStartDocument : virtual public SAXEvent
{
public:
  SAXEventStartDocument() {}
  ~SAXEventStartDocument() {}
  
  virtual const SAXEvent::EventType Type() const
  {
    return SAXEvent::eStartDocument;
  }
};

class SAXEventEndDocument : virtual public SAXEvent
{
public:
  SAXEventEndDocument() {}
  ~SAXEventEndDocument() {}
  
  virtual const SAXEvent::EventType Type() const
  {
    return SAXEvent::eEndDocument;
  }
};

#endif // SAX_DOCUMENT_EVENTS_H

