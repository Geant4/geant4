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
// $Id: ThreeVectorType.hh,v 1.2 2002-06-03 12:09:34 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------
//
#ifndef THREEVECTORTYPE_H
#define THREEVECTORTYPE_H 1

#include <string>

class ThreeVectorType {
public:
  ThreeVectorType() {
  }
  ~ThreeVectorType() {
  }
  void set_x( std::string& s ) {
    m_x = s;
  }
  void set_x( char* s ) {
    m_x = s;
  }
  std::string get_x() const {
    return m_x;
  }
  void set_y( std::string& s ) {
    m_y = s;
  }
  void set_y( char* s ) {
    m_y = s;
  }
  std::string get_y() const {
    return m_y;
  }
  void set_z( std::string& s ) {
    m_z = s;
  }
  void set_z( char* s ) {
    m_z = s;
  }
  std::string get_z() const {
    return m_z;
  }
  
private:
  std::string m_x;
  std::string m_y;
  std::string m_z;
};


#endif // THREEVECTORTYPE_H
