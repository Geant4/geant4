#ifndef REFERENCETYPE_H
#define REFERENCETYPE_H 1

#include <string>

class ReferenceType {
public:
  ReferenceType() {
  }
  ~ReferenceType() {
  }
  void set_ref( std::string& s ) {
    m_ref = s;
  }
  void set_ref( char* s ) {
    m_ref = s;
  }
  const std::string& get_ref() const {
    return m_ref;
  }
 
private:
  // IDREF
  std::string m_ref;
};


#endif // REFERENCETYPE_H
