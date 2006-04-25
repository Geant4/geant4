// $Id: XBase.hh,v 1.2 2006-04-25 07:54:29 kmura Exp $
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
