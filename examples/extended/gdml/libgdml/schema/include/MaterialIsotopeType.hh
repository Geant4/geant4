#ifndef MATERIALISOTOPETYPE_H
#define MATERIALISOTOPETYPE_H 1


#include "MaterialType.hh"
#include "DorDref.hh"
#include "atom.hh"

class MaterialIsotopeType : public MaterialType {
public:
  MaterialIsotopeType() : m_density( 0 ) {
  }
  ~MaterialIsotopeType() {
    if( m_density != 0 )
      delete m_density;
  }
  
  const std::string& get_N() const {
    return m_N;
  }
  const std::string& get_Z() const {
    return m_Z;
  }
  const atom& get_atom() const {
    return m_atom;
  }
  const SAXObject* get_DorDref() const {
    const SAXObject* so = 0;
    if( m_density != 0 ) {
      so = m_density->get_content();
    }
    return so;
  }
  
  void set_N( std::string& n ) {
    m_N = n;
  }
  void set_Z( std::string& z ) {
    m_Z = z;
  }
  void set_atom( const atom& a ) {
    m_atom = a;
  }
  void set_DorDref( const std::string& tag, SAXObject* so ) {
    if( !m_density ) {
      m_density = new DorDref;
    }
    m_density->set_content( tag, so);
  }
  
  
private:
  DorDref*    m_density;
  atom        m_atom;
  std::string m_N;
  std::string m_Z;
};

#endif // MATERIALISOTOPETYPE_H
