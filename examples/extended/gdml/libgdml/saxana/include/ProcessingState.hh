#ifndef PROCESSING_STATE_H
#define PROCESSING_STATE_H 1

class ProcessingState
{
public:
  typedef enum
  {
    eState,
    eAction
  } EType;
  
  virtual EType Type() const = 0;
};

#endif // PROCESSING_STATE_H

