#ifndef MATERIALMIXTURETYPE_H
#define MATERIALMIXTURETYPE_H 1

#include "MaterialType.hh"
#include "DorDref.hh"
#include "TorTref.hh"
#include "PorPref.hh"
#include "atom.hh"

class MaterialMixtureType : public MaterialType {
public:
  MaterialMixtureType() {
  }
  
  ~MaterialMixtureType() {
  }
  
  const std::string& get_Z() const {
    return m_Z;
  }
  
  const SAXObject* get_DorDref() const {
    return( m_density.get_content() );
  }
  
  ContentChoice* get_choice() {
    return &m_choice;
  }
  
  const ContentChoice* get_choice() const {
    return &m_choice;
  }
  
  void set_Z( const std::string& z ) {
    m_Z = z;
  }
  
  void set_DorDref( const std::string& tag, SAXObject* so ) {
    m_density.set_content( tag, so);
  }
  
  void set_choice( const std::string& tag, SAXObject* so ) {
    ContentGroup::ContentItem ci = { tag, so };
    m_choice.set_content( ci );
  }
  
private:
  DorDref         m_density;
  ContentChoice   m_choice;
  std::string     m_Z;
};

#endif // MATERIALMIXTURETYPE_H
