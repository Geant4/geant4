#ifndef SAX_PROCESSOR_H
#define SAX_PROCESSOR_H 1

#include "ProcessingContext.hh"
#include "StatusCode.hh"

class SAXEvent;
class SAXEventGun;
class StateStack;
class StateProcessMap;
class SAXSubscriberPool;
class ProcessingConfigurator;
class SAXObject;

class SAXProcessor : virtual public ProcessingContext
{
public:
  SAXProcessor();
  ~SAXProcessor();
  
  // Mandatory interface from ProcessingContext
  virtual const StateStack*             GetStack() const;
  virtual const SAXEvent*               GetLastEvent() const;
  virtual const SAXEventGun*            GetSAXEventGun() const; 
  virtual void                          SetSAXEventGun( const SAXEventGun* gun );
  virtual const ProcessingConfigurator* GetConfig() const;
  virtual SAXObject**                   GetTopObject() const;
  
  StatusCode Initialize();
  StatusCode Configure( ProcessingConfigurator* config );
  StatusCode Run();
  StatusCode Finalize();
  
  void ProcessEvent( const SAXEvent* const event );
  
private:
  StateProcessMap*        fMap;
  SAXSubscriberPool*      fPool;
  StateStack*             fStack;
  StateStack*             fNotifyStack;
  ProcessingConfigurator* fConfig;
  const SAXEvent*         fCurrentEvent;
  const SAXEventGun*      fCurrentGun;
  bool                    fIgnoring;
};

#endif // SAX_PROCESSOR_H

