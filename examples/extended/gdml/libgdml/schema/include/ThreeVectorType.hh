#ifndef THREEVECTORTYPE_H
#define THREEVECTORTYPE_H 1

#include <string>

class ThreeVectorType {
public:
  ThreeVectorType() {
  }
  ~ThreeVectorType() {
  }
  void set_x( std::string& s ) {
    m_x = s;
  }
  void set_x( char* s ) {
    m_x = s;
  }
  std::string get_x() const {
    return m_x;
  }
  void set_y( std::string& s ) {
    m_y = s;
  }
  void set_y( char* s ) {
    m_y = s;
  }
  std::string get_y() const {
    return m_y;
  }
  void set_z( std::string& s ) {
    m_z = s;
  }
  void set_z( char* s ) {
    m_z = s;
  }
  std::string get_z() const {
    return m_z;
  }
  
private:
  std::string m_x;
  std::string m_y;
  std::string m_z;
};


#endif // THREEVECTORTYPE_H
