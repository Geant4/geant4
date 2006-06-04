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
//                                              2004 Q
// ====================================================================
#ifndef XBASE_H
#define XBASE_H

// ====================================================================
//
// class definition
//
// ====================================================================

class XBase {
protected:
  int ival;
  double dval;

public:
  XBase();
  virtual ~XBase();

  void SetIVal(int aval);
  int GetIVal() const;

  void SetDVal(double aval);
  double GetDVal() const;

  void AMethod();
  virtual int VMethod(const XBase* abase) const=0;

};

// ====================================================================
// inline functions
// ====================================================================
inline void XBase::SetIVal(int i) { ival= i; }
inline int XBase::GetIVal() const { return ival; }

inline void XBase::SetDVal(double d) { dval= d; }
inline double XBase::GetDVal() const  { return dval; }

#endif
