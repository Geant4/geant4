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
// $Id: MaterialMixtureType.hh,v 1.2 2002-06-03 12:09:34 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------
//
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
