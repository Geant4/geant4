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
// $Id: ConstantType.hh,v 1.2 2002-06-03 12:09:33 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------
//
#ifndef CONSTANTTYPE_H
#define CONSTANTTYPE_H 1

#include <string>

class ConstantType {
  // An anonymous, local scope, value
  public:
    ConstantType() {
    }
    ConstantType( std::string& s ) : m_value( s ) {
    }
    ConstantType( char* s ) {
      m_value = s;
    }
    ~ConstantType() {
    }
    void set_value( std::string& s ) {
      m_value = s;
    } 
    std::string get_value() const {
      return m_value;
    }
  private:
    // double    
    std::string m_value;
};

#endif // CONSTANTTYPE_H
