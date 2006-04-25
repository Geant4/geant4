// $Id: AClass.hh,v 1.2 2006-04-25 07:54:28 kmura Exp $
// ====================================================================
//   AClass.hh
//
//                                              2004 Q
// ====================================================================
#ifndef ACLASS_H
#define ACLASS_H

// ====================================================================
//
// class definition
//
// ====================================================================

class AClass {
protected:
  static AClass* thePointer;
  AClass();

public:
  ~AClass();

  static AClass* GetPointer();

};

#endif
