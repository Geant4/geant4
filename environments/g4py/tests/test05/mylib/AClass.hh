// $Id: AClass.hh,v 1.2 2006-04-25 07:54:28 kmura Exp $
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

public:
  AClass();
  AClass(int i, double d=0.);
  ~AClass();

  void SetIVal(int i);
  int GetIVal() const;

  int AMethod();
  int AMethod(int i);
  int AMethod(int i, double d);

  int  BMethod();
  double BMethod(double d);
  
  double CMethod(int i, double d1=1., double d2=2.);

};

// ====================================================================
// inline functions
// ====================================================================

inline void AClass::SetIVal(int i) { ival= i; }
inline int AClass::GetIVal() const { return ival; }

#endif

