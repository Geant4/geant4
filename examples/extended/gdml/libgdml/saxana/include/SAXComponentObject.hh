#ifndef SAX_COMPONENT_OBJECT_H
#define SAX_COMPONENT_OBJECT_H 1

class SAXComponentObject
{
public:
  typedef enum
  {
    eProcess,
    eAction,
    eSubscriber
  } EType;
  
  SAXComponentObject() {}
  virtual ~SAXComponentObject() {}
  
  virtual const SAXComponentObject* Build() const = 0;
  virtual EType Type() const = 0;
};

#endif // SAX_COMPONENT_OBJECT_H

