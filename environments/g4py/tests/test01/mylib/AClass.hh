// $Id: AClass.hh,v 1.1 2006-02-27 10:05:24 kmura Exp $
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
