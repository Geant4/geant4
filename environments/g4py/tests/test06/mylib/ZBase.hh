// $Id: ZBase.hh,v 1.2 2006-04-25 07:54:29 kmura Exp $
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
