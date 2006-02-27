// $Id: AClass.hh,v 1.1 2006-02-27 10:08:30 kmura Exp $
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
