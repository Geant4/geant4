#ifndef __GEOMETRY__
#define __GEOMETRY__

#ifdef IS_GCC
#pragma interface
#endif

#ifndef IS_SGI
#include "globals.hh"
#endif

#include "newvector.hh"
#include "double.hh"

class Geometry
{
public:
  Geometry() : V(0) {}
  Geometry(REAL v) : V(v),V1(0) {}
  Geometry(REAL v,REAL v1) : V(v),V1(v1) {}
  REAL getVolume() { return V; }
  void newSize(REAL V_) { V = V_; V1 = V_; dSize(); }
  void newSize(REAL V_,REAL V1_) { V = V_; V1 = V1_; dSize(); }
  virtual void adjustSize(REAL) = 0;
  virtual REAL totalSurface() { Vektor3 x; return surface(x); }
  virtual REAL size() const = 0;
  virtual REAL surface(const Vektor3&) const = 0;
  virtual REAL expansion(REAL velocity) { throw "Geometry: expansion not defined!"; };
  virtual bool isInSurface(const Vektor3&) const { return true; }
  virtual bool isInside(const Vektor3&) { return true; }
  virtual REAL hitSurface(const Vektor3& v,const Vektor3& x0,Vektor3&) = 0;
  virtual bool homogeneous(Vektor4&) = 0;
  virtual Matrize rotateFrame(const Vektor3&) = 0;
  virtual void whatAmI(ostream& o) { o << "Error..."; }
  virtual void reflect(Vektor3& x) const { }
  bool moveToSurface(Vektor4& x,const Vektor4& p);
  bool moveToSurface(Vektor4& x,Vektor4& p,Matrize&,Vektor3&);
  virtual Vektor3 dr(const Vektor3& x1,const Vektor3& x2) const { return x1-x2; }
protected:
  virtual REAL dSize() = 0;
  REAL shrinkVolume(REAL dV) { REAL d = pow(1+dV/V,1.0/3.0); V += dV; return d; }
  REAL V,V1;
};

class halfSpace : public Geometry
{
public:
  halfSpace(REAL l,REAL dl = 0.0) : Geometry(l*frontArea,(l-dl)*frontArea),L(l),dL(dl) {}
  REAL hitSurface(const Vektor3& v,const Vektor3& x0,Vektor3& r0); 
  REAL expansion(REAL v) { return v/L; }
  bool isInside(const Vektor3& x) { return (x[3] <= 0); }
  //  bool isInSurface(const Vektor3& x) const { return (fabs(x[1])<=0.5 && fabs(x[2])<=0.5*L); }
  bool homogeneous(Vektor4&);
  Matrize rotateFrame(const Vektor3&);
  REAL size() const { return L; }
  REAL surface(const Vektor3&) const { return frontArea; }
  void whatAmI(ostream& o) { o << "half space (L = " << L << ")"; }
  void newSize(REAL V_) { L = V_/frontArea; V = V_; }
  void adjustSize(REAL dV) { L *= shrinkVolume(dV); }

  static REAL frontArea;
protected:
  REAL dSize() { return dL = (V-V1)/frontArea; }
  REAL shrinkVolume(REAL dV) { REAL d = dV/V; V += dV; return d; }
private:
  REAL L,dL;
};

class Box : public Geometry
{
public:
  Box(REAL l,REAL dl = 0.0) : Geometry(l*l*l,pow(l-dl,3)),L(l),dL(dl) {}
  bool isInside(const Vektor3& x) { 
    bool y = true; 
    for (int i=1; i<=3; i++) 
      if ( fabs(x[i]) > 0.5*L ) { 
	y = false; break; 
      }
    return y; 
  }
  REAL hitSurface(const Vektor3& v,const Vektor3& x0,Vektor3& r0);
  bool homogeneous(Vektor4&);
  Matrize rotateFrame(const Vektor3&);
  REAL size() const { return L; }
  REAL surface(const Vektor3&) const { return 6*L*L; }
  void whatAmI(ostream& o) { o << "Box (L = " << L << ")"; }
  void newSize(REAL V_) { L = pow(V,1.0/3.0); V = V_; }
  void adjustSize(REAL dV) { L *= shrinkVolume(dV); }
  virtual void reflect(Vektor3& x) const;
protected:
  REAL dSize() { return dL = (V-V1)/(3*sqr(L)); }

  REAL L,dL;
};

class InfiniteBox : public Box
{
public:
  InfiniteBox(REAL l) : Box(l) {}
  Vektor3 dr(const Vektor3& x1,const Vektor3& x2) const;
};

class Sphere : public Geometry
{
public:
  Sphere(REAL r,REAL dr = 0.0) 
    : Geometry(4.0/3.0*mathConstants::Pi*r*r*r,4.0/3.0*mathConstants::Pi*pow((r-dr),3)),R(r),dR(dr) {}
  REAL expansion(REAL v) { return 3*v/R; }
  bool isInside(const Vektor3& x) { return square(x) < sqr(R); }
  REAL hitSurface(const Vektor3& v,const Vektor3& x0,Vektor3& r0);
  bool homogeneous(Vektor4&);
  REAL size() const { return R; }
  Matrize rotateFrame(const Vektor3&);
  REAL totalSurface() { return 4.0*mathConstants::Pi*R*R; }
  REAL surface(const Vektor3& x) const { return 4.0*mathConstants::Pi*square(x); }
  void whatAmI(ostream& o) { o << "Sphere (R = " << R << ")"; }
  void adjustSize(REAL dV) { R *= shrinkVolume(dV); }
  void newSize(REAL V_) { R = pow(3.0/4.0/mathConstants::Pi*V_,1.0/3.0); V = V_; }
protected:
  REAL dSize() { return dR = 1.0/3.0*R*(1-V1/V); }
  REAL R,dR;
};

class Tube : public Geometry
{
public:
  Tube(REAL r,REAL l) : Geometry(mathConstants::Pi*r*r*l),R(r),L(l) {}
  bool isInside(const Vektor3& x) { return double(sqr(x[1])+sqr(x[2]))<=R*R && fabs(x[3]) <= double(L/2); }
  REAL hitSurface(const Vektor3& v,const Vektor3& x0,Vektor3& r0);
  bool homogeneous(Vektor4&);
  REAL size() const { return R; }
  Matrize rotateFrame(const Vektor3&);
  REAL surface(const Vektor3&) const { return 2.0*mathConstants::Pi*R*(R+L); }
  void whatAmI(ostream& o) { o << "Tube (R = " << R << ",L = " << L << ")"; }
  void adjustSize(REAL dV) { REAL d = shrinkVolume(dV); R*=d; L*=d; }
public:
  REAL dSize() { return 0.0; }
  REAL R,L;
};

#endif
