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
// $Id: ContentGroup.hh,v 1.3 2002-06-03 12:09:34 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------
//
#ifndef CONTENTGROUP_H
#define CONTENTGROUP_H 1

#include "SAXObject.hh"

#include <vector>

#include <iostream>

class ContentGroup : virtual public SAXObject {
public:
  typedef enum {
    sequence,
    choice,
    all
  } ContentType;
  struct ContentItem {
    std::string tag;
    SAXObject*  object;
  };
public:
  ContentGroup() {
  }
  virtual ~ContentGroup() {
  }
  virtual SAXObject::Type type() {
    return SAXObject::contentGroup;
  }
  virtual ContentGroup::ContentType content_type() const = 0;
  virtual bool               empty() = 0;
};

class ContentChoice : virtual public ContentGroup {
public:
  ContentChoice() {
    m_choice.tag    = "";
    m_choice.object = 0;
  }
  virtual ~ContentChoice() {
    if( m_choice.object != 0 ) {
      //std::cout << "ContentChoice destructor deleting " << m_choice.tag << std::endl;
      delete m_choice.object;
    }
  }
  virtual ContentGroup::ContentType content_type() const {
    return ContentGroup::choice;
  }
  virtual bool empty() {
    return( m_choice.object == 0 );
  }
  void set_content( ContentItem& ci ) {
    m_choice = ci;
  }
  ContentItem& content() {
    return m_choice;
  }
  const ContentItem& content() const {
    return m_choice;
  }

private:
  ContentItem m_choice;
};

class ContentSequence : virtual public ContentGroup {
public:
  ContentSequence() {
  }
  virtual ~ContentSequence() {
    for( unsigned int idx= 0; idx < m_sequence.size(); idx++ ) {
      if( m_sequence[idx].object != 0 ) {
        //std::cout << "ContentSequence destructor deleting " << m_sequence[idx].tag << std::endl;
        delete m_sequence[idx].object;
      }
    }
  }
  virtual ContentGroup::ContentType content_type() const {
    return ContentGroup::sequence;
  }
  virtual bool empty() {
    return( m_sequence.size() == 0 );
  }
  size_t size() const {
    return( m_sequence.size() );
  }
  void add_content( ContentItem& ci ) {
    m_sequence.push_back( ci );
  }
  ContentItem content( unsigned int idx ) const {
    return m_sequence[idx];
  }
  void clear() {
    for( unsigned int idx= 0; idx < m_sequence.size(); idx++ ) {
      if( m_sequence[idx].object != 0 ) {
        //std::cout << "ContentSequence destructor clearing " << m_sequence[idx].tag << std::endl;
        delete m_sequence[idx].object;
      }
    }
    m_sequence.clear();
  }
  
private:
  std::vector<ContentItem> m_sequence;
};

class ContentAll : virtual public ContentGroup {
public:
  ContentAll() : m_all( 0 ) {
  }
  virtual ~ContentAll() {
    for( unsigned int idx= 0; idx < m_all.size(); idx++ ) {
      if( m_all[idx].object != 0 ) {
        //std::cout << "ContentAll destructor deleting " << m_all[idx].tag << std::endl;
        delete m_all[idx].object;
      }
    }
  }
  virtual ContentGroup::ContentType content_type() const {
    return ContentGroup::all;
  }
  virtual bool empty() {
    return( m_all.size() == 0 );
  }
  size_t size() const {
    return( m_all.size() );
  }
  void add_content( ContentItem& ci ) {
    m_all.push_back( ci );
  }
  ContentItem content( unsigned int idx ) const {
    return m_all[idx];
  }
  void clear() {
    for( unsigned int idx= 0; idx < m_all.size(); idx++ ) {
      if( m_all[idx].object != 0 ) {
        //std::cout << "ContentAll clearing " << m_all[idx].tag << std::endl;
        delete m_all[idx].object;
      }
    }
    m_all.clear();
  }
  
private:
  std::vector<ContentItem> m_all;
};

#endif // CONTENTGROUP_H
