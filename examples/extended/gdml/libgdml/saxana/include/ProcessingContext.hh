#ifndef PROCESSING_CONTEXT_H
#define PROCESSING_CONTEXT_H 1

class StateStack;
class SAXEvent;
class SAXEventGun;
class ProcessingConfigurator;
class SAXObject;

class ProcessingContext
{
public:
  virtual const StateStack*             GetStack() const = 0;
  virtual const SAXEvent*               GetLastEvent() const = 0;
  virtual const SAXEventGun*            GetSAXEventGun() const = 0; 
  virtual void                          SetSAXEventGun( const SAXEventGun* gun ) = 0;
  virtual const ProcessingConfigurator* GetConfig() const = 0;
  virtual SAXObject**                   GetTopObject() const = 0;
};

#endif // PROCESSING_CONTEXT_H

