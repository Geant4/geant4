#ifndef ACTION_PROCESING_STATE_H
#define ACTION_PROCESING_STATE_H 1

#include "ProcessingState.hh"
#include "SAXEventGun.hh"

class ActionProcessingState : virtual public ProcessingState
{
public:
  //typedef RCObjectHandle<SAXEventGun> EventGun;
  typedef SAXEventGun* EventGun;
  
  ActionProcessingState();
  ~ActionProcessingState();

  virtual ProcessingState::EType Type() const
  {
    return ProcessingState::eAction;
  }

private:
  EventGun fPrevious;
};

#endif // ACTION_PROCESING_STATE_H

