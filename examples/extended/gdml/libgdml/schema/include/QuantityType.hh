#ifndef QUANTITYTYPE_H
#define QUANTITYTYPE_H 1

#include "ConstantType.hh"

class QuantityType : public ConstantType {
public:
  QuantityType() {
  }
  ~QuantityType() {
  }
  void set_unit( std::string& s ) {
    m_unit = s;
  }
  void set_unit( char* s ) {
    m_unit = s;
  }
  std::string get_unit() const {
    return m_unit;
  }
  void set_type( std::string& s ) {
    m_type = s;
  }
  void set_type( char* s ) {
    m_type = s;
  }
  std::string get_type() const {
    return m_type;
  }
private:
  // string
  std::string m_unit;
  // string
  std::string m_type;
};



#endif // QUANTITYTYPE_H
