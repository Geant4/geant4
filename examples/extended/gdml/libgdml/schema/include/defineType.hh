//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: defineType.hh,v 1.3 2002-06-03 12:09:34 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------
//
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
    
    //std::cout << "defineType desctructor called..." << std::endl;
    
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
