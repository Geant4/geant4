#ifndef TRD_H
#define TRD_H 1

#include "SAXObject.hh"
#include "SolidType.hh"

class trd : public SAXObject, public SolidType {
public:
  trd() {
  }
  virtual ~trd() {
  }
  virtual SAXObject::Type type() {
    return SAXObject::element;
  }
  
  const std::string& get_x1() const {
    return m_x1;
  }
  const std::string& get_x2() const {
    return m_x2;
  }
  const std::string& get_y1() const {
    return m_y1;
  }
  const std::string& get_y2() const {
    return m_y2;
  }
  const std::string& get_z() const {
    return m_z;
  }
  
  void set_x1( const std::string& x1 ) {
    m_x1 = x1;
  }
  void set_x2( const std::string& x2 ) {
    m_x2 = x2;
  }
  void set_y1( const std::string& y1 ) {
    m_y1 = y1;
  }
  void set_y2( const std::string& y2 ) {
    m_y2 = y2;
  }
  void set_z( const std::string& z ) {
    m_z = z;
  }
  
private:
  std::string m_x1;
  std::string m_x2;
  std::string m_y1;
  std::string m_y2;
  std::string m_z;
};

#endif // TRD_H
