#ifndef IDENTIFIABLECONSTANTTYPE_H
#define IDENTIFIABLECONSTANTTYPE_H 1

#include "ConstantType.hh"

class IdentifiableConstantType : public ConstantType
{
public:
  IdentifiableConstantType() {
  }
  IdentifiableConstantType( std::string& s ) : m_name( s ) {
  }
  IdentifiableConstantType( char* s ) {
    m_name = s;
  }
  ~IdentifiableConstantType() {
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

#endif // IDENTIFIABLECONSTANTTYPE_H
