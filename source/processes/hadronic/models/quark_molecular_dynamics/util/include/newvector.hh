#ifndef __NEWVECTOR__
#define __NEWVECTOR__

#include "g4std/iostream"
#include <math.h>
#include "String.hh"

class Vektor
{
  friend int dim(const Vektor& v) { return v.N; }
  friend double length(const Vektor& v) { return sqrt(v*v); }
  friend double square(const Vektor& v) { return v*v; }
  friend G4std::ostream& operator<<(G4std::ostream& o,const Vektor& v) { o << "("; for (int k=1; k<v.N; k++) o << v[k] << ","; o << v[v.N] << ")"; return o; }
  friend G4std::istream& operator>>(G4std::istream& in,Vektor& v);
  friend Vektor operator-(const Vektor& v) { Vektor x(v.N); for (int k=0; k<v.N; k++) x.a[k] = -v.a[k]; return x; }
  friend Vektor operator+(const Vektor& v,const Vektor& w) { Vektor x(v.N); for (int k=0; k<v.N; k++) x.a[k] = v.a[k]+w.a[k]; return x; }
  friend Vektor operator-(const Vektor& v,const Vektor& w) { Vektor x(v.N); for (int k=0; k<v.N; k++) x.a[k] = v.a[k]-w.a[k]; return x; }
  friend double operator*(const Vektor& v,const Vektor& w) { double s = 0.0; for (int k=0; k<v.N; k++) s += v.a[k]*w.a[k]; return s; }
  friend Vektor operator*(const Vektor& v,const double& y) { Vektor x(v.N); for (int k=0; k<v.N; k++) x.a[k] = v.a[k]*y; return x; }
  friend Vektor operator*(const double& y,const Vektor& v) { Vektor x(v.N); for (int k=0; k<v.N; k++) x.a[k] = v.a[k]*y; return x; }
  friend Vektor operator/(const Vektor& v,const double& y) { Vektor x(v.N); for (int k=0; k<v.N; k++) x.a[k] = v.a[k]/y; return x; }
  friend Vektor product(const Vektor& v,const Vektor& w) { Vektor x(v.N); for (int k=0; k<v.N; k++) x.a[k] = v.a[k]*w.a[k]; return x; }
  
  int N;
protected:
  double* a;
public:
  Vektor() : N(0),a(0) {}
  Vektor(int n) : N(n),a(new double[N]) { for (int k=0; k<N; k++) a[k] = 0.0; }
  Vektor(int n,double* b) : N(n),a(new double[N]) { memcpy(a,b,N*sizeof(double)); }
  Vektor(const Vektor& v) : N(v.N),a(new double[N]) { memcpy(a,v.a,N*sizeof(double)); }
  Vektor(int n,const String&);
  virtual ~Vektor() { delete [] a; }
  void reset(int dN,int where = -1);
  void remove(int pos,int dN);
  Vektor& operator=(double* b) { memcpy(a,b,N*sizeof(double)); return *this; }
  Vektor& operator=(const Vektor& v) { if ( N != v.N ) reset(v.N); memcpy(a,v.a,N*sizeof(double)); return *this; }
  Vektor& operator+=(const Vektor& v) { for (int k=0; k<N; k++) a[k] += v.a[k]; return *this; }
  Vektor& operator-=(const Vektor& v) { for (int k=0; k<N; k++) a[k] -= v.a[k]; return *this; }
  Vektor& operator+=(const double& y) { for (int k=0; k<N; k++) a[k] += y; return *this; }
  Vektor& operator-=(const double& y) { for (int k=0; k<N; k++) a[k] -= y; return *this; }
  Vektor operator+(const double& y) { Vektor x(N); for (int k=0; k<N; k++) x.a[k] = a[k]+y; return x; }
  Vektor operator-(const double& y) { Vektor x(N); for (int k=0; k<N; k++) x.a[k] = a[k]-y; return x; }
  Vektor& operator*=(const double& y) { for (int k=0; k<N; k++) a[k] *= y; return *this; }
  Vektor& operator/=(const double& y) { for (int k=0; k<N; k++) a[k] /= y; return *this; }
  operator double*() { return a; }
  operator void*() const { return a; }
  virtual double operator[](int k) const { return a[k-1]; }
  virtual double& operator[](int k) { return a[k-1]; }
};

class Vektor3 : public Vektor
{
public:
  Vektor3() : Vektor(3) {}
  Vektor3(const Vektor& v) : Vektor(v) {}
  Vektor3(double* b) : Vektor(3,b) {}
  Vektor3(double x,double y,double z) : Vektor(3) { a[0] = x; a[1] = y; a[2] = z; }
  Vektor3(const String& s) : Vektor(3,s) {}
  static Vektor3 isotropy(double R);
  //  Vektor3(const three_vector& x) : Vektor(3,
  //  operator three_vector() { return three_vector(a[0],a[1],a[2]); }
};

class Vektor4 : public Vektor
{
  friend Vektor3 spatial(const Vektor4& x) { return Vektor3(x.a[1],x.a[2],x.a[3]); }
  friend double operator*(const Vektor4& v,const Vektor4& w) { double s = v.a[0]*w.a[0]; for (int k=1; k<4; k++) s -= v.a[k]*w.a[k]; return s; }
 public:
  Vektor4() : Vektor(4) {}
  Vektor4(const Vektor& v) : Vektor(v) {}
  Vektor4(double* b) : Vektor(4,b) {}
  Vektor4(const Vektor3& x,double t) : Vektor(4) { a[1] = x[1]; a[2] = x[2]; a[3] = x[3]; a[0] = t; }
  Vektor4(double x,double y,double z,double t) : Vektor(4) { a[1] = x; a[2] = y; a[3] = z; a[0] = t;}
  Vektor4& operator=(const Vektor4& x) { for (int i=0; i<4; i++) a[i] = x[i]; return *this; }
  void setSpatial(const Vektor3& x) { for (int i=1; i<4; i++) a[i] = x[i]; }
  virtual double operator[](int k) const { return a[k]; }
  virtual double& operator[](int k) { return a[k]; }
};

class Matrize
{
  friend G4std::ostream& operator<<(G4std::ostream& o,const Matrize& A);
  friend int Rows(const Matrize& A) { return A.N; }
  friend int Cols(const Matrize& A) { return A.M; }
  friend Matrize operator*(const Matrize& A,const Matrize& B);
  friend Vektor operator*(const Matrize& A,const Vektor& B);
  friend Matrize transpose(const Matrize& A);

  int N,M;
  Vektor* a;
public:
  Matrize() : N(0),M(0),a(0) {}
  Matrize(int n);
  Matrize(int n,int m);
  void set(int n,int m);
  double operator()(int i,int j) const { return (a[j-1])[i]; }
  double& operator()(int i,int j) { return (a[j-1])[i]; }
  Vektor column(int k) const { return a[k-1]; }
  Vektor& column(int k) { return a[k-1]; }
};

#endif
