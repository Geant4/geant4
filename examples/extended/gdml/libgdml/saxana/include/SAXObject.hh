#ifndef SAX_OBJECT_H
#define SAX_OBJECT_H 1

#include "SAXObjectBase.hh"

#include <string>

class SAXObject : virtual public SAXObjectBase
{
public:
  typedef enum {
    element,             // element
    contentGroup,        // choice, sequence, all
    attributeGroup,      // reusable group of attributes
    elementGroup         // reusable group of elements
  } Type;
    
public:
  virtual ~SAXObject() {
  }

  virtual bool IsCacheable() const
  {
    return false;
  }
  
  virtual Type type() = 0;
};

#endif // SAX_OBJECT_H

