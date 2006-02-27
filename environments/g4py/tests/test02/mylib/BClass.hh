// $Id: BClass.hh,v 1.1 2006-02-27 10:08:30 kmura Exp $
// ====================================================================
//   BClass.hh
//
//                                              2004 Q
// ====================================================================
#ifndef BCLASS_H
#define BCLASS_H

#include "XBase.hh"

// ====================================================================
//
// class definition
//
// ====================================================================

class BClass : public XBase {
public:
  BClass();
  ~BClass();

  void AMethod(); // overrided
  virtual int VMethod(const XBase* abase) const;

};

#endif
