#ifndef DEFINETYPE_H
#define DEFINETYPE_H 1

#include "ContentGroup.hh"

#include <string>
#include <vector>

class defineType {
   
public:
  defineType() : m_choices( 0 ) {
  }
  
  ~defineType() {
    
    std::cout << "defineType desctructor called..." << std::endl;
    
    if( m_choices != 0 ) {
      if( m_choices->size() > 0 ) {
        for( unsigned int idx = 0; idx < m_choices->size(); idx++ ) {
          if( (*m_choices)[idx] != 0 ) {
            delete (*m_choices)[idx];
          }
        }
        delete m_choices;
        m_choices = 0;
      }
    }
  }
  
  void add_to_content( const std::string& tag, SAXObject* so ) {
    ContentGroup::ContentItem ci = {tag,so};
    if( m_choices == 0 ) {
      m_choices = new std::vector<ContentChoice*>;
    }
    ContentChoice* choice = new ContentChoice;
    choice->set_content( ci );
    m_choices->push_back( choice );
  }

  SAXObject* get_content( unsigned int idx ) const {
    return( (m_choices->size()>0) ? (*m_choices)[idx] : 0  );
  }

  unsigned int size() {
    unsigned int ret = 0;
    if( m_choices != 0 ) {
      ret = m_choices->size();
    }
    return ret;
  }

  const unsigned int size() const {
    unsigned int ret = 0;
    if( m_choices != 0 ) {
      ret = m_choices->size();
    }
    return ret;
  }

private:
  // 0..* choice
  std::vector<ContentChoice*>* m_choices;
};


#endif // DEFINETYPE_H
