#ifndef COMPOSITE_H
#define COMPOSITE_H 1

#include "ReferenceType.hh"

class composite : virtual public SAXObject, public ReferenceType
{
public:
  composite() {
  }
  virtual ~composite() {
  }
  virtual SAXObject::Type type() {
    return SAXObject::element;
  }
  const std::string& get_n() const {
    return m_n;
  }
  void set_n( const std::string& n ) {
    m_n = n;
  }
  
private:
  // xs:positiveInteger
  std::string m_n;
};

#endif // COMPOSITE_H
