#ifndef ACTION_H
#define ACTION_H 1

#include "StatusCode.hh"
#include "SXComponentObject.hh"

class ProcessingContext;

class Action : virtual public SAXComponentObject
{
public:
  virtual const SAXComponentObject* Build() const
  {
    return 0;
  }
  virtual const SAXComponentObject::EType Type() const
  {
    return SAXComponentObject::eAction;
  }

public:
  virtual StatusCode Run( const ProcessingContext* const context ) = 0;
};

#endif // ACTION_H

