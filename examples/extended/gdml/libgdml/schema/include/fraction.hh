#ifndef FRACTION_H
#define FRACTION_H 1

#include "ReferenceType.hh"

class fraction : virtual public SAXObject, public ReferenceType
{
public:
  fraction() {
  }
  virtual ~fraction() {
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
  // xs:double
  std::string m_n;
};

#endif // FRACTION_H
