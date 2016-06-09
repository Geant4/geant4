//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: XBase.hh,v 1.4 2006-06-29 15:39:00 gunter Exp $
// ====================================================================
//   XBase.hh
//
//                                         2005 Q
// ====================================================================
#ifndef XBASE_H
#define XBASE_H

#include <string>

// ====================================================================
//
// class definition
//
// ====================================================================

class XBase {
protected:
  int ival;

public:
  XBase();
  virtual ~XBase();

  inline void SetIVal(int i);
  inline int GetIVal() const;
  
  virtual std::string PVMethod()=0;

};

// ====================================================================
// inline functions
// ====================================================================

inline void XBase::SetIVal(int i) { ival= i; }
inline int XBase::GetIVal() const { return ival; }

#endif
