#ifndef SETUP_H
#define SETUP_H 1

#include "SAXObject.hh"
#include "ReferenceType.hh"

class setup : virtual public SAXObject
{
public:
  class world : public SAXObject, public ReferenceType {
    public:
      world() {
      }
      virtual ~world() {
      }
      virtual SAXObject::Type type() {
        return SAXObject::element;
      }
  };

public:
  setup() : m_match( false ) {
  }
  virtual ~setup() {
  }
  virtual SAXObject::Type type() {
    return SAXObject::element;
  }
  
  const setup::world& get_world() const {
    return m_world;
  }
  const std::string& get_name() const {
    return m_ID;
  }
  const std::string& get_version() const {
    return m_version;
  }
  bool get_match() const {
    return m_match;
  }
  
  void set_world( const setup::world& w ) {
    m_world = w;
  }
  void set_name( const std::string& n ) {
    m_ID = n;
  }
  void set_version( const std::string& v ) {
    m_version = v;
  }
  void set_match( bool m ) {
    m_match = m;
  }
  
private:
  setup::world m_world;
  std::string  m_ID;
  std::string  m_version;
  bool         m_match;
};

#endif // SETUP_H
