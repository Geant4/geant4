#ifndef SAX_COMPONENT_FACTORY_BASE_H
#define SAX_COMPONENT_FACTORY_BASE_H 1

#include "SAXComponentObject.hh"

class SAXComponentFactoryBase
{
public:
  virtual void                Load() const {}
  virtual SAXComponentObject* Build() const = 0;
  virtual SAXComponentObject::EType Type() const = 0;
};

#endif // SAX_COMPONENT_FACTORY_BASE_H

