#ifdef IS_GCC
#pragma implementation
#include "double.hh"
//REAL double_a;
//const REAL double_b=double_a,double_c=double_a;
//REAL double_x=G4std::max(double_b,double_c);
#else
#include "double.hh"
#endif

#include "Error.hh"
#include <math.h>
#include "g4std/iostream"
#include <algo.h>

const double Double::Epsilon = 1e-12;
const double Double::EpsWarning = 1e-8;
const Double Double::Null = 0.0;
const Double Double::One = 1.0;
const Double Double::Infinity = 1e+30;
const Double Double::NotDefined = 1e+60;

const Double mathConstants::Pi = 4*atan(1); 
const Double mathConstants::hc = 0.1973; // [GeV*fm] 
const Double mathConstants::gamma = 0.5277; // [GeV*fm] 

//const Double NotDefined(IsNotDefined);

G4std::ostream& operator<<(G4std::ostream& o,const Double& a)
{
  if (fabs(a) == Double::Infinity) 
    o << "Double::Infinity";
  else
    o << a.x;
  return o;
}

Double operator/(const Double& a,const Double& b)
{ 
  double y;
  if (b!=0.0)
    if (a!=0)
      y = a.x/b.x;
    else
      y = 0;
  else 
    if (a != 0)
      y = sign(a)*Double::Infinity;
    else {
      //      G4cerr << "WARNING: Division by Zero" << G4endl;
      //exit(1);
      throw DividedByZero(a.x); 
    }
  return y;
}

Double operator/(const Double& a,double b)
{ 
  double y;
  if (Double(b)!=0.0)
    if (a!=0) 
      y = a.x/b;
    else
      y = 0;
  else {
    //    G4cerr << "WARNING: Division by Zero" << G4endl;
    //exit(1);
    throw DividedByZero(a.x); 
  }
  return y;
}

Double operator/(double a,const Double& b)
{ 
  Double a1(a);
  double y;
  if (b!=0.0)
    if (Double(a)!=0) 
      y = a/b.x;
    else
      y = 0;
  else 
    if (a1 != 0)
      y = sign(a1)*Double::Infinity;
    else {
      //  G4cerr << "WARNING: Division by Zero" << G4endl;
      // exit(1);
      throw DividedByZero(a); 
    }
  return y;
}

Double operator/(const Double& a,int b)
{ 
  double y;
  if (b!=0)
    if (a!=0)
      y = a.x/double(b);
    else
      y = 0;
  else 
    if (a != 0)
      y = sign(a)*Double::Infinity;
    else {
      G4cerr << "WARNING: Division by Zero" << G4endl;
      throw DividedByZero(a.x); 
  }
  return y;
}

Double operator/(int a,const Double& b)
{ 
  double y;
  if (b!=0.0)
    if (a!=0)
      y = double(a)/b.x;
    else
      y = 0;
  else 
    if (a) 
      y = a*Double::Infinity;
    else {
      //      G4cerr << "WARNING: Division by Zero" << G4endl;
      //exit(1);
      throw DividedByZero(a); 
    }
  return y;
}

Double sqrt(const Double& y)
{
  if (y >= 0)
    return sqrt(G4std::max(y.x,0.0));
  else {
    //    G4cerr << "WARNING: SQRT: x = " << y.x << G4endl;
    //exit(1);
    //return 0;
    if (fabs(y.x)>Double::EpsWarning)
      throw WrongArgument("sqrt",y.x);
    else
      return 0;
  }
}

Double log(const Double& y)
{
  if (y > 0)
    return log(y.x);
  else {
    //    G4cerr << "LOG: x = " << y.x << G4endl;
    //exit(1);
    throw WrongArgument("log (x>0)",y.x);
  }
}

Double log(const Double& b,const Double& y)
{
  if (y > 0)
    return log(y.x)/log(b);
  else {
    //    G4cerr << "LOG_B: x = " << y.x << G4endl;
    //exit(1);
    throw WrongArgument("log_b (x>0)",y.x);
  }
}

Double tan(const Double& a) 
{
  if (cos(a) != 0)
    return tan(a.x);
  else {
    G4cerr << "TAN: x = " << a.x << G4endl;
//    exit(1);
//    return 0;
    throw WrongArgument("tan",a.x);
  }
}

Double cot(const Double& a) 
{
  if (sin(a) != 0)
    return 1/tan(a.x);
  else {
    //    G4cerr << "COTANGENS: x = " << a.x << G4endl;
    //exit(1);
    throw WrongArgument("cot",a.x);
  }
}

Double arcsin(const Double& y) 
{
  if (fabs(y) <= 1) {
    if (fabs(y) < 1) 
      return asin(y.x);
    else
      return sign(y.x)*mathConstants::Pi/2;
  }
  else {
    if (fabs(fabs(y.x)-1)>Double::EpsWarning) {
      //      G4cerr << "WARNING: ARCSIN: x = " << y.x << ":" << fabs(y.x) - 1 << G4endl;
      //exit(1);
      throw WrongArgument("arcsin (-1<=x<=1)",y.x);
    }
    else
      return sign(y.x)*mathConstants::Pi/2;
  }
}

Double arccos(const Double& y) 
{
  if (fabs(y) <= 1) {
    if (fabs(y) < 1) 
      return acos(y.x);
    else
      return (1-sign(y.x))/2*mathConstants::Pi;
  }
  else {
    if (fabs(fabs(y.x)-1)>Double::EpsWarning) {
      //      G4cerr << "WARNING: ARCCOS: x = " << y.x << G4endl;
      //exit(1);
      throw WrongArgument("arccos (-1<=x<=1)",y.x);
    }
    else
      return (1-sign(y.x))/2*mathConstants::Pi;
  }
}

inline Double arccot(const Double& y) 
{ 
  return mathConstants::Pi/2 - atan(y.x); 
}
  
Double coth(const Double& a) 
{
  if (tanh(a) != 0)
    return 1/tanh(a.x);
  else {
    //    G4cerr << "COTH: x = " << a.x << G4endl;
    //exit(1);
    throw WrongArgument("coth",a.x);
  }
}

Double arcosh(const Double& y) 
{
  if (y >= 1)
    return acosh(G4std::max(y.x,1.0));
  else {
    if (fabs(y.x-1)>Double::EpsWarning) {
      //    G4cerr << "ARCOSH: x = " << y.x << G4endl;
      //exit(1);
      throw WrongArgument("arcosh (x>=1)",y.x);
    }
    else
      return 0;
  }
}

Double artanh(const Double& y) 
{
  if (fabs(y) < 1)
    return atanh(y.x);
  else {
    //    G4cerr << "ARTANH: x = " << y.x << G4endl;
    //exit(1);
    throw WrongArgument("artanh (|x|<1)",y.x);
  }
}

Double arcoth(const Double& y) 
{
  if (fabs(y) > 1)
    return atanh(y.x);
  else {
    //    G4cerr << "ARCOTH: x = " << y.x << G4endl;
    //exit(1);
    throw WrongArgument("arcoth (|x|>1)",y.x);
  }
}

G4std::ostream& operator<<(G4std::ostream& o,const Xdouble& x) 
{
  if (Valid(x))
    o << x.x;
  else
    o << "NotDef";
  return o;
}

BinOStream& operator<<(BinOStream& o,const Xdouble& x) 
{
  o << double(x.x) << char(x.status);
  return o;
}

BinIStream& operator>>(BinIStream& o,Xdouble& x) 
{
  char c;
  double y;
  o >> y >> c;
  x = y;
  x.status = c;
  return o;
}




