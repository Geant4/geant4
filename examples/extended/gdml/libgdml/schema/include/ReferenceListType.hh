#ifndef REFERENCELISTTYPE_H
#define REFERENCELISTTYPE_H 1

#include <string>


class ReferenceListType {
public:
  ReferenceListType() {
  }
  ~ReferenceListType() {
  }
  void set_refs( std::string& s ) {
    m_refs = s;
  }
  void set_refs( char* s ) {
    m_refs = s;
  }
  std::string& refs() const {
    return m_refs;
  }
 
private:
  // IDREFS
  std::string m_refs;
};


#endif // REFERENCELISTTYPE_H
