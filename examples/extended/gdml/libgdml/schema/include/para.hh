#ifndef PARA_H
#define PARA_H 1

#include "SAXObject.hh"
#include "SolidType.hh"

class para : public SAXObject, public SolidType {
public:
  para() {
  }
  virtual ~para() {
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
  const std::string& get_alpha() const {
    return m_alpha;
  }
  const std::string& get_theta() const {
    return m_theta;
  }
  const std::string& get_phi() const {
    return m_phi;
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
  void set_alpha( const std::string& alpha ) {
    m_alpha = alpha;
  }
  void set_theta( const std::string& theta ) {
    m_theta = theta;
  }
  void set_phi( const std::string& phi ) {
    m_phi = phi;
  }
  
private:
  std::string m_x;
  std::string m_y;
  std::string m_z;
  std::string m_alpha;
  std::string m_theta;
  std::string m_phi;
};

#endif // PARA_H
