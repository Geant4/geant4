#ifndef IDENTIFIABLEEXPRESSIONTYPE_H
#define IDENTIFIABLEEXPRESSIONTYPE_H 1


#include "ExpressionType.hh"


class IdentifiableExpressionType : public ExpressionType {
public:
  IdentifiableExpressionType() {
  }
  ~IdentifiableExpressionType() {
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


#endif // IDENTIFIABLEEXPRESSIONTYPE_H
