//
// $RCSfile: Constants.hh,v $
//
// $Revision: 1.1 $
// $Date: 1999-12-09 11:50:35 $
// $Author: sscherer $
// $Locker:  $
// $State: Exp $
// DOSfile: constant.h
// UNIXfile: Constants.H
//
// $Log: not supported by cvs2svn $
// Revision 1.1.1.1  1998/09/22 16:31:05  mhofmann
// U++ V0.9
//
// Revision 1.1.1.1  1997/07/17 12:54:54  mhofmann
// Initial import
//
// Revision 1.1.1.1  1996/10/04 14:37:50  mhofmann
// Phase Transition Project
//
//
 
#ifndef _Constants_H
#define _Constants_H

#ifdef __MSDOS__
  #include "definiti.h"
#else
  #include "Definitions.hh"
#endif

#define Real double

class PhysConstants {
  public:
    static const Real PionRange;
    static const Real hc;
    static const Real hcSquared;
    static const Real alpha;
    static const Real e;
};

class MathConstants {
  public:
    static const Real pi;
    static const Real fourPi;
};

class NumConstants {
  public:
    static const Real hugeFloat;
    static const Real tinyFloat;
    static const Real hugeReal;
    static const Real tinyReal;
    static const int hugeInt;
};

#endif // _Constants_H

