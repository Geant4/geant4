#ifndef IDENTIFIABLEQUANTITYTYPE_H
#define IDENTIFIABLEQUANTITYTYPE_H 1

#include "QuantityType.hh"


class IdentifiableQuantityType : public QuantityType {
public:
  IdentifiableQuantityType() {
  }
  ~IdentifiableQuantityType() {
  }
  void set_name( std::string& s ) {
    m_name = s;
  }
  void set_name( char* s ) {
    m_name = s;
  }
  std::string get_name() const {
    return m_name;
  }
  
private:
  std::string m_name;
};



#endif // IDENTIFIABLEQUANTITYTYPE_H
