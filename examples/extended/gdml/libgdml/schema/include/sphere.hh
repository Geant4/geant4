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
// $Id: sphere.hh,v 1.2 2002-06-03 12:09:34 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------
//
#ifndef SPHERE_H
#define SPHERE_H 1

#include "SAXObject.hh"
#include "SolidType.hh"

class sphere : public SAXObject, public SolidType {
public:
  sphere() {
  }
  virtual ~sphere() {
  }
  virtual SAXObject::Type type() {
    return SAXObject::element;
  }
  
  const std::string& get_rmin() const {
    return m_rmin;
  }
  const std::string& get_rmax() const {
    return m_rmax;
  }
  const std::string& get_startphi() const {
    return m_startphi;
  }
  const std::string& get_deltaphi() const {
    return m_deltaphi;
  }
  const std::string& get_starttheta() const {
    return m_starttheta;
  }
  const std::string& get_deltatheta() const {
    return m_deltatheta;
  }
  
  void set_rmin( const std::string& rmin ) {
    m_rmin = rmin;
  }
  void set_rmax( const std::string& rmax ) {
    m_rmax = rmax;
  }
  void set_startphi( const std::string& startphi ) {
    m_startphi = startphi;
  }
  void set_deltaphi( const std::string& deltaphi ) {
    m_deltaphi = deltaphi;
  }
  void set_starttheta( const std::string& starttheta ) {
    m_starttheta = starttheta;
  }
  void set_deltatheta( const std::string& deltatheta ) {
    m_deltatheta = deltatheta;
  }
  
private:
  std::string m_rmin;
  std::string m_rmax;
  std::string m_startphi;
  std::string m_deltaphi;
  std::string m_starttheta;
  std::string m_deltatheta;
};



#endif // SPHERE_H
