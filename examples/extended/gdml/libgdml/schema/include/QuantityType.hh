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
// $Id: QuantityType.hh,v 1.3 2002-07-25 10:41:21 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------
//
#ifndef QUANTITYTYPE_H
#define QUANTITYTYPE_H 1

#include "ConstantType.hh"

class QuantityType : public ConstantType {
public:
  QuantityType() {
  }
  ~QuantityType() {
  }
  void set_unit( const std::string& s ) {
    m_unit = s;
  }
  void set_unit( const char* s ) {
    m_unit = s;
  }
  std::string get_unit() const {
    return m_unit;
  }
  void set_type( const std::string& s ) {
    m_type = s;
  }
  void set_type( const char* s ) {
    m_type = s;
  }
  std::string get_type() const {
    return m_type;
  }
private:
  // string
  std::string m_unit;
  // string
  std::string m_type;
};



#endif // QUANTITYTYPE_H
