// $Id: ZBase.hh,v 1.1 2006-02-27 10:05:25 kmura Exp $
// ====================================================================
//   ZBase.hh
//
//                                         2005 Q
// ====================================================================
#ifndef ZBASE_H
#define ZBASE_H

#include <string>

// ====================================================================
//
// class definition
//
// ====================================================================

class ZBase {
public:
  ZBase();
  virtual ~ZBase();

  void AMethod();
  virtual void VMethod(std::string message);

};

#endif
