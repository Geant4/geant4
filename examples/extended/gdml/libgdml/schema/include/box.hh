#ifndef BOX_H
#define BOX_H 1

#include "SAXObject.hh"
#include "SolidType.hh"

class box : public SAXObject, public SolidType {
public:
  box() {
  }
  virtual ~box() {
  }
  virtual SAXObject::Type type() {
    return SAXObject::element;
  }
  
  const std::string& get_x() const {
    return m_x;
  }
  const std::string& get_y() const {
    return m_y;
  }
  const std::string& get_z() const {
    return m_z;
  }
  
  void set_x( const std::string& x ) {
    m_x = x;
  }
  void set_y( const std::string& y ) {
    m_y = y;
  }
  void set_z( const std::string& z ) {
    m_z = z;
  }
  
private:
  std::string m_x;
  std::string m_y;
  std::string m_z;
};



#endif // BOX_H
