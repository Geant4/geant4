#ifndef SPHERE_H
#define SPHERE_H 1

#include "SAXObject.hh"
#include "SolidType.hh"

class sphere : public SAXObject, public SolidType {
public:
  sphere() {
  }
  virtual ~sphere() {
  }
  virtual SAXObject::Type type() {
    return SAXObject::element;
  }
  
  const std::string& get_rmin() const {
    return m_rmin;
  }
  const std::string& get_rmax() const {
    return m_rmax;
  }
  const std::string& get_startphi() const {
    return m_startphi;
  }
  const std::string& get_deltaphi() const {
    return m_deltaphi;
  }
  const std::string& get_starttheta() const {
    return m_starttheta;
  }
  const std::string& get_deltatheta() const {
    return m_deltatheta;
  }
  
  void set_rmin( const std::string& rmin ) {
    m_rmin = rmin;
  }
  void set_rmax( const std::string& rmax ) {
    m_rmax = rmax;
  }
  void set_startphi( const std::string& startphi ) {
    m_startphi = startphi;
  }
  void set_deltaphi( const std::string& deltaphi ) {
    m_deltaphi = deltaphi;
  }
  void set_starttheta( const std::string& starttheta ) {
    m_starttheta = starttheta;
  }
  void set_deltatheta( const std::string& deltatheta ) {
    m_deltatheta = deltatheta;
  }
  
private:
  std::string m_rmin;
  std::string m_rmax;
  std::string m_startphi;
  std::string m_deltaphi;
  std::string m_starttheta;
  std::string m_deltatheta;
};



#endif // SPHERE_H
