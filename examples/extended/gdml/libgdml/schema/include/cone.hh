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
// $Id: cone.hh,v 1.2 2002-06-03 12:09:34 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------
//
#ifndef CONE_H
#define CONE_H 1

#include "SAXObject.hh"
#include "SolidType.hh"

class cone : public SAXObject, public SolidType {
public:
  cone() {
  }
  virtual ~cone() {
  }
  virtual SAXObject::Type type() {
    return SAXObject::element;
  }
  
  const std::string& get_rmin1() const {
    return m_rmin1;
  }
  const std::string& get_rmin2() const {
    return m_rmin2;
  }
  const std::string& get_rmax1() const {
    return m_rmax1;
  }
  const std::string& get_rmax2() const {
    return m_rmax2;
  }
  const std::string& get_z() const {
    return m_z;
  }
  const std::string& get_startphi() const {
    return m_startphi;
  }
  const std::string& get_deltaphi() const {
    return m_deltaphi;
  }
  
  void set_rmin1( const std::string& rmin1 ) {
    m_rmin1 = rmin1;
  }
  void set_rmin2( const std::string& rmin2 ) {
    m_rmin2 = rmin2;
  }
  void set_rmax1( const std::string& rmax1 ) {
    m_rmax1 = rmax1;
  }
  void set_rmax2( const std::string& rmax2 ) {
    m_rmax2 = rmax2;
  }
  void set_z( const std::string& z ) {
    m_z = z;
  }
  void set_startphi( const std::string& startphi ) {
    m_startphi = startphi;
  }
  void set_deltaphi( const std::string& deltaphi ) {
    m_deltaphi = deltaphi;
  }
  
private:
  std::string m_rmin1;
  std::string m_rmin2;
  std::string m_rmax1;
  std::string m_rmax2;
  std::string m_z;
  std::string m_startphi;
  std::string m_deltaphi;
};


#endif // CONE_H
