#ifndef SAX_OBJECT_BASE_H
#define SAX_OBJECT_BASE_H 1

class SAXObjectBase
{
public:
  virtual bool IsCacheable() const = 0;
};

#endif // SAX_OBJECT_BASE_H

