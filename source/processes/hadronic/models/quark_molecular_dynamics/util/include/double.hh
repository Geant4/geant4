#ifndef __DOUBLE__
#define __DOUBLE__

#ifdef IS_GCC
#pragma interface
#endif

#include "g4std/iostream"
#include <math.h>
#include "BinIOStream.hh"
#include "Error.hh"
/*
template<class T>
T max(const T& a,const T& b) { return (a>b) ? a : b; }

template<class T>
T min(const T& a,const T& b) { return (a<b) ? a : b; }
*/
class Xdouble;

class Double
{
  friend Xdouble;
  friend G4std::ostream& operator<<(G4std::ostream& o,const Double& a);
  friend G4std::istream& operator>>(G4std::istream& o,Double& a) { return o >> a.x; }
  friend BinOStream& operator<<(BinOStream& o,const Double& y) { return o << y.x; };
  friend BinIStream& operator>>(BinIStream& in,Double& y) { return in >> y.x; }
friend int operator==(const Double& a,const Double& b) { return a.IsEqual(b.x); }
friend int operator==(double a,const Double& b) { return b.IsEqual(a); }
friend int operator==(const Double& a,double b) { return a.IsEqual(b); }
friend int operator==(const Double& a,int b) { return a.IsEqual((double)b); }
friend int operator==(int a,const Double& b) { return b.IsEqual((double)a); }
#ifndef GCC27
friend int operator!=(const Double& a,const Double& b) { return !a.IsEqual(b.x); }
#endif
friend int operator!=(double a,const Double& b) { return !b.IsEqual(a); }
friend int operator!=(const Double& a,double b) { return !a.IsEqual(b); }
friend int operator!=(int a,const Double& b) { return !b.IsEqual(a); }
friend int operator!=(const Double& a,int b) { return !a.IsEqual(b); }
#ifndef GCC27
friend int operator<=(const Double& a,const Double& b) { return (a==b) ? 1 : (a.x<b.x); }
#endif
friend int operator<=(double a,const Double& b) { return (Double(a)==b) ? 1 : (a<b.x); }
friend int operator<=(const Double& a,double b) { return (a==Double(b)) ? 1 : (a.x<b); }
friend int operator<=(int a,const Double& b) { return (Double(a)==b) ? 1 : (a<b.x); }
friend int operator<=(const Double& a,int b) { return (a==Double(b)) ? 1 : (a.x<b); }
#ifndef GCC27
  friend int operator>=(const Double& a,const Double& b) { return (a==b) ? 1 : (a.x>b.x); }
#endif
friend int operator>=(double a,const Double& b) { return (Double(a)==b) ? 1 : (a>b.x); }
friend int operator>=(const Double& a,double b) { return (a==Double(b)) ? 1 : (a.x>b); }
friend int operator>=(int a,const Double& b) { return (Double(a)==b) ? 1 : (a>b.x); }
friend int operator>=(const Double& a,int b) { return (a==Double(b)) ? 1 : (a.x>b); }
friend int operator<(const Double& a,const Double& b) { return (a==b) ? 0 : (a.x<b.x); }
friend int operator<(double a,const Double& b) { return (Double(a)==b) ? 0 : (a<b.x); }
friend int operator<(const Double& a,double b) { return (a==Double(b)) ? 0 : (a.x<b); }
friend int operator<(int a,const Double& b) { return (Double(a)==b) ? 0 : (a<b.x); }
friend int operator<(const Double& a,int b) { return (a==Double(b)) ? 0 : (a.x<b); }
#ifndef GCC27
  friend int operator>(const Double& a,const Double& b) { return (a==b) ? 0 : (a.x>b.x); }
#endif
friend int operator>(double a,const Double& b) { return (Double(a)==b) ? 0 : (a>b.x); }
friend int operator>(const Double& a,double b) { return (a==Double(b)) ? 0 : (a.x>b); }
friend int operator>(int a,const Double& b) { return (Double(a)==b) ? 0 : (a>b.x); }
friend int operator>(const Double& a,int b) { return (a==Double(b)) ? 0 : (a.x>b); }
inline friend Double operator-(const Double& a) { return -a.x; }
friend Double operator+(const Double& a,const Double& b) { return a.x+b.x; }
friend Double operator+(double a,const Double& b) { return a+b.x; }
friend Double operator+(const Double& a,double b) { return a.x+b; }
friend Double operator+(int a,const Double& b) { return a+b.x; }
friend Double operator+(const Double& a,int b) { return a.x+b; }
friend Double operator-(const Double& a,const Double& b) { return a.x-b.x; }
friend Double operator-(const Double& a,double b) { return a.x-b; }
friend Double operator-(double a,const Double& b) { return a-b.x; }
friend Double operator-(const Double& a,int b) { return a.x-b; }
friend Double operator-(int a,const Double& b) { return a-b.x; }
friend Double operator*(const Double& a,const Double& b) { return a.x*b.x; }
friend Double operator*(const Double& a,double b) { return a.x*b; }
friend Double operator*(double a,const Double& b) { return a*b.x; }
friend Double operator*(const Double& a,int b) { return a.x*b; }
friend Double operator*(int a,const Double& b) { return a*b.x; }
friend Double operator/(const Double& a,const Double& b);
friend Double operator/(double a,const Double& b);
friend Double operator/(const Double& a,double b);
friend Double operator/(int a,const Double& b);
friend Double operator/(const Double& a,int b);

friend Double sign(Double x) { return (x==0.0) ? 0.0 : ((x < 0.0) ? -1.0 : 1.0); }
friend int round(Double y) { return int(sign(y.x)*(int)(fabs(y.x)+Epsilon)); }
friend Double fabs(const Double& y) { return (y == 0.0) ? 0.0 : fabs(y.x); }
friend Double sqr(const Double& y) { return y.x*y.x; }
  //friend Double pow(const Double& y,int n);

friend Double sqrt(const Double& y);
friend Double exp(const Double& y) { return exp(y.x); }
friend Double log(const Double& y);
friend Double log(const Double& b,const Double& y);

friend Double sin(const Double& y) { return sin(y.x); }
friend Double cos(const Double& y) { return cos(y.x); }
friend Double tan(const Double& y);
friend Double cot(const Double& y);

friend Double arcsin(const Double& y);
friend Double arccos(const Double& y);
friend Double arctan(const Double& y) { return atan(y.x); }
friend Double arccot(const Double& y);

friend Double sinh(const Double& y) { return sinh(y.x); }
friend Double cosh(const Double& y) { return cosh(y.x); }
friend Double tanh(const Double& y) { return tanh(y.x); }
friend Double coth(const Double& y);

friend Double arsinh(const Double& y) { return asinh(y.x); }
friend Double arcosh(const Double& y);
friend Double artanh(const Double& y);
friend Double arcoth(const Double& y);

public:
  Double() : x(0) {}
  Double(double y) : x(y) {}
  Double(const Double& y) : x(y.x) {}
  operator double() { return x; }
  operator double() const { return x; }
  Double& operator=(const Double& y)  { x  = y.x; return *this; }
  Double& operator+=(const Double& y) { x += y.x; return *this; }
  Double& operator-=(const Double& y) { x -= y.x; return *this; }
  Double& operator*=(const Double& y) { x *= y.x; return *this; }
  Double& operator/=(const Double& y) { x /= y.x; return *this; }
  static const double Epsilon;
  static const double EpsWarning;
  static const Double Infinity;
  static const Double NotDefined;
  static const Double Null;
  static const Double One;
protected:
  double x;
  int IsEqual(double y = 0.0) const { 
#ifndef DOUBLE_NO_EQUAL_CHECK
    return (fabs(x-y)<Epsilon); 
#else
    return x==y; 
#endif
  }
};

class Xdouble : public Double
{
  friend int Valid(const Xdouble& y) { return  (y.status ==0); }
  friend G4std::ostream& operator<<(G4std::ostream& o,const Xdouble& x);
  friend BinOStream& operator<<(BinOStream&,const Xdouble&);
  friend BinIStream& operator>>(BinIStream&,Xdouble&);
public:
  enum { OK = 0, NotDef = 1, NotSet = 2 };
  //  Xdouble(int err = 0) : Double(0),status(err) {}
  Xdouble() : Double(0),status(NotSet) {}
  Xdouble(const Double& y,int err=0) : Double(y),status(err) {}
  Xdouble(const double y,int err=0) : Double(y),status(err) {}
  Xdouble(const Xdouble& y) : Double(y.x),status(y.status) {}
  void SetStatus(int err) { x = 0.0; status = err; }
  int GetStatus() { return status; }
  Xdouble& operator=(const Xdouble& y) { x = y.x; status = y.status; return *this; }
  //  Xdouble& operator=(const Double& y)  {  status = 0; x  = y.x; return *this; }
  //  Xdouble& operator=(int err)  { x  = 0; status = err; return *this; }
  Xdouble& operator+=(const Xdouble& y) { if (!status) x += y.x; return *this; }
  Xdouble& operator-=(const Xdouble& y) { if (!status) x -= y.x; return *this; }
  Xdouble& operator*=(const Xdouble& y) { if (!status) x *= y.x; return *this; }
  Xdouble& operator/=(const Xdouble& y) { if (!status) x /= y.x; return *this; }
private:
  int IsEqual(double y = 0.0) const { return status ? 0 : (fabs(x-y)<Epsilon); }
  int status;
};

class mathConstants
{
public:
  static const Double Pi;
  static const Double hc;
  static const Double gamma;
};


typedef Double REAL;
#endif





