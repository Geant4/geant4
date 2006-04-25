// $Id: XFunc.cc,v 1.2 2006-04-25 07:54:29 kmura Exp $
// ====================================================================
//   XFunc.cc
//
//                                         2005 Q
// ====================================================================
#include "XFunc.hh"

// ====================================================================
//
// class description
//
// ====================================================================

int func1(int i)
{
  return i;
}

double func1(int i, double d)
{
  return i*d;
}


int func2(int i)
{
  return i;
}

int func2(int i0, int i1)
{
  return i0*i1;
}

int func3(int i)
{
  return 0;
}

int func3(double d)
{
  return 1;
}

