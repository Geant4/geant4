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
// $Id: trd.hh,v 1.2 2002-06-03 12:09:34 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------
//
#ifndef TRD_H
#define TRD_H 1

#include "SAXObject.hh"
#include "SolidType.hh"

class trd : public SAXObject, public SolidType {
public:
  trd() {
  }
  virtual ~trd() {
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
  const std::string& get_y1() const {
    return m_y1;
  }
  const std::string& get_y2() const {
    return m_y2;
  }
  const std::string& get_z() const {
    return m_z;
  }
  
  void set_x1( const std::string& x1 ) {
    m_x1 = x1;
  }
  void set_x2( const std::string& x2 ) {
    m_x2 = x2;
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
  
private:
  std::string m_x1;
  std::string m_x2;
  std::string m_y1;
  std::string m_y2;
  std::string m_z;
};

#endif // TRD_H
