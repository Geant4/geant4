#ifndef IDENTIFIABLEQUANTITYVECTORTYPE_H
#define IDENTIFIABLEQUANTITYVECTORTYPE_H 1


#include "QuantityVectorType.hh"


class IdentifiableQuantityVectorType : public QuantityVectorType {
public:
  IdentifiableQuantityVectorType() {
  }
  ~IdentifiableQuantityVectorType() {
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
  // ID required
  std::string m_name;    
};


#endif // IDENTIFIABLEQUANTITYVECTORTYPE_H
