#ifndef SOLIDTYPE_H
#define SOLIDTYPE_H 1

#include <string>

class SolidType {
public:
  SolidType() {
  }
  ~SolidType() {
  }
  
  const std::string& get_lunit() const {
    return m_lunit;
  }
  const std::string& get_aunit() const {
    return m_aunit;
  }
  const std::string& get_name() const {
    return m_name;
  }
  
  void set_lunit( const std::string& lu ) {
    m_lunit = lu;
  }
  void set_aunit( const std::string& au ) {
    m_aunit = au;
  }
  void set_name( const std::string& name ) {
    m_name = name;
  }
  
private:
  std::string m_lunit;
  std::string m_aunit;
  std::string m_name; 
};



#endif // SOLIDTYPE_H
