// $Id: AClass.hh,v 1.2 2006-04-25 07:54:28 kmura Exp $
// ====================================================================
//   AClass.hh
//
//                                              2004 Q
// ====================================================================
#ifndef ACLASS_H
#define ACLASS_H

#include "XBase.hh"

// ====================================================================
//
// class definition
//
// ====================================================================

class AClass : public XBase {
public:
  AClass();
  ~AClass();

  void AMethod(); // overrided
  virtual int VMethod(const XBase* abase) const;

};

#endif
