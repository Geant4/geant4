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
// $Id: composite.hh,v 1.2 2002-06-03 12:09:34 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------
//
#ifndef COMPOSITE_H
#define COMPOSITE_H 1

#include "ReferenceType.hh"

class composite : virtual public SAXObject, public ReferenceType
{
public:
  composite() {
  }
  virtual ~composite() {
  }
  virtual SAXObject::Type type() {
    return SAXObject::element;
  }
  const std::string& get_n() const {
    return m_n;
  }
  void set_n( const std::string& n ) {
    m_n = n;
  }
  
private:
  // xs:positiveInteger
  std::string m_n;
};

#endif // COMPOSITE_H
