#ifndef CONSTANTTYPE_H
#define CONSTANTTYPE_H 1

#include <string>

class ConstantType {
  // An anonymous, local scope, value
  public:
    ConstantType() {
    }
    ConstantType( std::string& s ) : m_value( s ) {
    }
    ConstantType( char* s ) {
      m_value = s;
    }
    ~ConstantType() {
    }
    void set_value( std::string& s ) {
      m_value = s;
    } 
    std::string get_value() const {
      return m_value;
    }
  private:
    // double    
    std::string m_value;
};

#endif // CONSTANTTYPE_H
