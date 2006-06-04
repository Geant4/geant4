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
// $Id: XBase.hh,v 1.3 2006-06-04 21:35:59 kmura Exp $
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
