#ifndef __MATHERR__
#define __MATHERR__

#include "globals.hh"
#include "g4std/iostream"

class Error
{
public:
  virtual void writeMessage(G4std::ostream& o) const { o << "ERROR occured...\n"; }
  virtual ~Error() {}
};

class IndexOutOfRange : public Error
{
public:
//
// sps-0801: use of G4std::max(max_ind) produces error "max cannot be initialized in a constructor."
//
  IndexOutOfRange(int i,int max_ind) : index(i),max(max_ind) {} 
  void writeMessage(G4std::ostream& o) const { o << "Index out of Range: " << index << " > " << max << " !" << G4endl; } 
private:
  int index,max;
};

class WrongNumberOfArguments : public Error
{
public:
  WrongNumberOfArguments(int n,int m) : need(n),got(m) {}
  void writeMessage(G4std::ostream& o) const { o << "Wrong number of arguments: Need " << need << " and got only " << got << G4endl; }
private:
  int need,got;
};

class MathErr : public Error 
{
public:
  MathErr(char* text,double arg) : CallFun(text),Argument(arg) {}
  MathErr(double arg) : CallFun("Unknown"),Argument(arg) {}
  MathErr() : CallFun("Unknown"),Argument(0) {}
  void writeMessage(G4std::ostream& o) const { 
    _print(o); 
    o << CallFun << " called with argument x = " << Argument << G4endl; 
  } 
private:
  char* CallFun;
  double Argument;
protected:
  virtual void _print(G4std::ostream& o) const { o << "MATHERR: "; }
};

class DividedByZero : public MathErr
{
public:
  DividedByZero(double arg) : MathErr("operator/()",arg) {}
protected:
  void _print(G4std::ostream& o) const { o << "DIVIDED BY ZERO: "; }
};

class WrongArgument : public MathErr
{
public:
  WrongArgument(char* text,double arg) : MathErr(text,arg) {}
protected:
  void _print(G4std::ostream& o) const { o << "WRONG ARGUMENT: "; }
};

class Zero_Not_Found : public Error
{
public:
  Zero_Not_Found(double x0,double x1,double f0,double f1) : a(x0),b(x1),fa(f0),fb(f1) {}
  void writeMessage(G4std::ostream& o) const { o << "Zero not found!\n"; _print(o); }
  double a,b,fa,fb;
protected:
  virtual void _print(G4std::ostream& o) const { 
    o << "a = " << a << ", b = " << b << ": f(a) = " << fa << ", f(b) = " << fb << G4endl; 
  }
};

template <class t>
class AmbiguousResult : public Error
{
  friend int dim(const AmbiguousResult<t>& m) { return m.N; }
public:
  AmbiguousResult(int n,t* a) : N(n),x(new t[n]) 
  { for (int i=0; i<n; i++) x[i] = a[i]; }
  ~AmbiguousResult() { delete [] x; };
  t& operator[](int i) { return x[i]; }
  void writeMessage(G4std::ostream& o) const { o << "Ambiguos Result: " << N << " solutions detected.\nx = "; for (int i=0; i<N; i++) o << x[i] << "  "; o << G4endl;} 
private:
  int N;
  t* x;
};

class TooMuchIter : public Zero_Not_Found
{
public:
  TooMuchIter(double x0,double x1,double f0,double f1) : Zero_Not_Found(x0,x1,f0,f1) {}
  void _print(G4std::ostream& o) const { o << "Too many iterations: " << a << "  " << b << " : " << fa << "  " << fb << G4endl; }
};

class NoZeroDetected : public Zero_Not_Found
{
public:
  NoZeroDetected(double x0,double x1,double f0,double f1) 
    : Zero_Not_Found(x0,x1,f0,f1) {}
  void _print(G4std::ostream& o) const { 
    o << "Interval has NO or MORE THAN ONE zeros. Cannot treat this...\n"; 
  }
};

class VariableNowhereDefined : public Error
{
public:
  VariableNowhereDefined() {}
  void writeMessage(G4std::ostream& o) const { o << "Variable nowhere defined!" << G4endl; }
};

#endif
