// $Id: XBase.hh,v 1.2 2006-04-25 07:54:28 kmura Exp $
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
