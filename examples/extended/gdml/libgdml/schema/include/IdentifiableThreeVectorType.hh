#ifndef IDENTIFIABLETHREEVECTORTYPE_H
#define IDENTIFIABLETHREEVECTORTYPE_H 1


#include "ThreeVectorType.hh"


class IdentifiableThreeVectorType : public ThreeVectorType {
public:
  IdentifiableThreeVectorType() {
  }
  ~IdentifiableThreeVectorType() {
  }
  void set_name( std::string& s ) {
    m_name = s;
  }
  void set_name( char* s ) {
    m_name = s;
  }
  std::string& name() const {
    return m_name;
  }
private:
  // ID required
  std::string m_name;    
};


#endif // IDENTIFIABLETHREEVECTORTYPE_H
