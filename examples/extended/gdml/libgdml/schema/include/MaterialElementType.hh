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
// $Id: MaterialElementType.hh,v 1.2 2002-06-03 12:09:34 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------
//
#ifndef MATERIALELEMENTTYPE_H
#define MATERIALELEMENTTYPE_H 1

#include "ContentGroup.hh"
#include "MaterialType.hh"
#include "DorDref.hh"
#include "atom.hh"

class AtomOrFraction : public ContentChoice
{
public:
  AtomOrFraction() : isAtom( false ) {
  }
  ~AtomOrFraction() {
  }
  const SAXObject* get_content() const {
    const SAXObject* so;
    if( isAtom )
      so = this;
    else
      so = &m_fractions;
    return so;
  }

  void set_content( const std::string& tag, SAXObject* so ) {
    if( !isAtom ) {
      ContentGroup::ContentItem ci = {tag,so};
      if( tag == "fraction" ) {
        // Switch to sequence content model
        m_fractions.add_content( ci );
      } else if( tag == "atom" ) {
        ContentChoice::set_content( ci );
        isAtom = true;
        // No more updates should happen from now on....
      }
    }
  }

private:
  bool isAtom;
  ContentSequence m_fractions;
};

class MaterialElementType : public MaterialType {
public:
  MaterialElementType() : m_density( 0 ) {
  }
  ~MaterialElementType() {
    if( m_density != 0 )
      delete m_density;
  }
  
  const std::string& get_N() const {
    return m_N;
  }
  const std::string& get_Z() const {
    return m_Z;
  }
  const SAXObject* get_AtomOrFraction() const {
    return m_aorf.get_content();
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
  void set_AtomOrFraction( const std::string& tag, SAXObject* so ) {
    m_aorf.set_content( tag, so );
  }
  void set_DorDref( const std::string& tag, SAXObject* so ) {
    if( !m_density ) {
      m_density = new DorDref;
    }
    m_density->set_content( tag, so);
  }
  
  
private:
  DorDref*        m_density;
  AtomOrFraction  m_aorf;
  std::string     m_N;
  std::string     m_Z;
};

#endif // MATERIALELEMENTTYPE_H
