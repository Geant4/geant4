#ifndef TRAP_H
#define TRAP_H 1

#include "SAXObject.hh"
#include "SolidType.hh"

class trap : public SAXObject, public SolidType {
public:
  trap() {
  }
  virtual ~trap() {
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
  const std::string& get_x3() const {
    return m_x3;
  }
  const std::string& get_x4() const {
    return m_x4;
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
  const std::string& get_theta() const {
    return m_theta;
  }
  const std::string& get_phi() const {
    return m_phi;
  }
  const std::string& get_alpha1() const {
    return m_theta;
  }
  const std::string& get_alpha2() const {
    return m_alpha2;
  }
  
  void set_x1( const std::string& x1 ) {
    m_x1 = x1;
  }
  void set_x2( const std::string& x2 ) {
    m_x2 = x2;
  }
  void set_x3( const std::string& x3 ) {
    m_x3 = x3;
  }
  void set_x4( const std::string& x4 ) {
    m_x4 = x4;
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
  void set_theta( const std::string& theta ) {
    m_theta = theta;
  }
  void set_phi( const std::string& phi ) {
    m_phi = phi;
  }
  void set_alpha1( const std::string& alpha1 ) {
    m_alpha1 = alpha1;
  }
  void set_alpha2( const std::string& alpha2 ) {
    m_alpha2 = alpha2;
  }
  
private:
  std::string m_x1;
  std::string m_x2;
  std::string m_x3;
  std::string m_x4;
  std::string m_y1;
  std::string m_y2;
  std::string m_z;
  std::string m_theta;
  std::string m_phi;
  std::string m_alpha1;
  std::string m_alpha2;
};

#endif // TRAP_H
