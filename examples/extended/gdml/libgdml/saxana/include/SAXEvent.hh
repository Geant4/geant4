#ifndef SAX_EVENT_H
#define SAX_EVENT_H 1

class SAXEvent
{
public:
  typedef enum
  {
    eStartDocument,
    eEndDocument,
    eStartElement,
    eEndElement,
    eCharacters,
    ePI,
    eComment,
    eWarning,
    eError,
    eFatalError
  } EventType;
public:
  SAXEvent() {}
  virtual ~SAXEvent() {}
  
  virtual const SAXEvent::EventType Type() const = 0;
};

#endif // SAX_EVENT_H

