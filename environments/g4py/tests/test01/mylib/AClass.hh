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
// $Id: AClass.hh,v 1.3 2006-06-04 21:35:59 kmura Exp $
// ====================================================================
//   AClass.hh
//
//                                         2005 Q
// ====================================================================
#ifndef ACLASS_H
#define ACLASS_H

// ====================================================================
//
// class definition
//
// ====================================================================

class AClass {
private:
  int ival;
  double dval;
  
public:
  AClass();
  AClass(int i, double d=0.);

  ~AClass();

  inline void SetIVal(int i);
  inline int GetIVal() const;

  inline void SetDVal(double d);
  inline double GetDVal() const;

  void AMethod();

};

// ====================================================================
// inline functions
// ====================================================================

inline void AClass::SetIVal(int i) { ival= i; }
inline int AClass::GetIVal() const { return ival; }

inline void AClass::SetDVal(double d) { dval= d; }
inline double AClass::GetDVal() const { return dval; }

#endif
