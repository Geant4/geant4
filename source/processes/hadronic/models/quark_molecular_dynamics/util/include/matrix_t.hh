#ifndef __MATRIX_T__
#define __MATRIX_T__

#ifdef IS_GCC
#pragma interface
void terminate();
#endif

#include "g4std/iostream"
#include <stdlib.h>
#include <stdarg.h>
#include "Subscript.hh"
#include "heapwr.hh"
#include "globals.hh"

template <class t>
class mathType
{
public:
    static t null;
    static t one;
};


template<class t>
class arrayBase 
{
  friend G4std::ostream& operator<<(G4std::ostream& o,const arrayBase<t>& a);
  friend G4std::istream& operator>>(G4std::istream& in,const arrayBase<t>& a);
  friend int operator==(const arrayBase<t>& a,const arrayBase<t>& b);
  friend int operator!=(const arrayBase<t>& a,const arrayBase<t>& b);
public:
  virtual Subscript dim() const { return size; }
  virtual Subscript rows() const { return size; }
  virtual Subscript cols() const { return 1; }
  //These operators may vary the starting index (0 or 1) as defined in the
  //deried classes! Per default they start at i=0!
  virtual t& operator[](Subscript i)       { return coeff(i); }
  virtual t  operator[](Subscript i) const { return coeff(i); }
  // coeff ALWAYS start at i=0!!! These function are called from classes 
  // arrayBase, mathArrayBase and vector. Whatever you do in derived classes,
  // be careful and aware: Don't overload them unless you are sure what you
  // are doing!!!!
  virtual t& coeff(Subscript i) ;
  virtual t coeff(Subscript i) const;
  virtual int isEqual(const arrayBase<t>& b) const;
  virtual ~arrayBase();
  void set(Subscript n);
  void set(Subscript n,const t& def);
  void set(Subscript n,const t* defArray);
  void set(const arrayBase<t>& a);
  virtual void writeToStream(G4std::ostream& o) const;
  virtual void readFromStream(G4std::istream& in);
  operator t*() { return x; }
protected:
  arrayBase() : size(0),x(0) {}
  arrayBase(Subscript n);
  arrayBase(Subscript n,const t& def);
  arrayBase(Subscript n,const t* defArray);
  arrayBase(const arrayBase<t>& a);
protected:
  virtual void writeOut(G4std::ostream& o) const;
  virtual void readIn(G4std::istream& i) const;
  Subscript size;
  t* x;
};

template <class t>
class mathArrayBase : public arrayBase<t>
{
  friend t square(const mathArrayBase<t>& a) { return a.norm(); }
public:
    mathArrayBase<t> add(const mathArrayBase<t>& b) const;
    mathArrayBase<t> sub(const mathArrayBase<t>& b) const;
    mathArrayBase<t> mul_t(const t& b) const;
    mathArrayBase<t>& addSelf(const mathArrayBase<t>& b)  ;
    mathArrayBase<t>& subSelf(const mathArrayBase<t>& b)  ;
    mathArrayBase<t>& mulSelf_t(const t& y) ;
  mathArrayBase<t>& flipSign() ;
  virtual t norm() const;
protected:
    mathArrayBase() : arrayBase<t>() {}
    mathArrayBase(Subscript n) : arrayBase<t>(n) {}
    mathArrayBase(Subscript n,const t& def) : arrayBase<t>(n,def) {}
    mathArrayBase(Subscript n,const t* defArray) : arrayBase<t>(n,defArray) {}
    mathArrayBase(const arrayBase<t>& a) : arrayBase<t>(a) {}
};

template <class t> class Matrix;

template <class t>
class Vector : public mathArrayBase<t>
{
  friend Vector<t> operator-(const Vector<t>& a);
  friend Vector<t> operator+(const Vector<t>& a,const Vector<t>& b);
  friend Vector<t> operator-(const Vector<t>& a,const Vector<t>& b);
  friend Vector<t> operator*(const Vector<t>& a,const t& y);
  friend Vector<t> operator*(const t& y,const Vector<t>& a);
  friend t scalarP(const Vector<t>& a,const Vector<t>& b);
public:
  Vector() : mathArrayBase<t>() {}
  Vector(Subscript n) : mathArrayBase<t>(n) {}
  Vector(Subscript n,const t& def) : mathArrayBase<t>(n,def) {}
  Vector(Subscript n,const t* defArray) : mathArrayBase<t>(n,defArray) {}
  Vector(const arrayBase<t>& a) : mathArrayBase<t>(a) {}
  ~Vector() {}
  /*  Vector(Subscript n,REAL x1 ...) :  mathArrayBase<t>(n)
  {
    va_list ap;
    va_start(ap,n);

    this->coeff(0) = x1;
    for (Subscript i=1; i<n; i++)
      coeff(i) = va_arg(ap,t);
  }
  */
  virtual t& operator[](Subscript j) { return coeff(j-1); }
  virtual t operator[](Subscript j) const { return coeff(j-1); }
  Vector<t>& operator=(const Vector<t>& y);
  Vector<t>& operator+=(const Vector<t>& y);
  Matrix<t> concat(const Matrix<t>& B) const;
};

template <class t>
class Matrix : public mathArrayBase<t>
{
  friend Matrix<t> operator-(const Matrix<t>& a);
  friend Matrix<t> operator+(const Matrix<t>& a,const Matrix<t>& b);
  friend Matrix<t> operator-(const Matrix<t>& a,const Matrix<t>& b);
  friend Matrix<t> operator*(const Matrix<t>& a,const t& y);
  friend Matrix<t> operator*(const t& y,const Matrix<t>& a);
  friend Matrix<t> operator*(const Matrix<t>& a,const Matrix<t>& b);
  friend Vector<t> operator*(const Matrix<t>& a,const Vector<t>& b);
  //  friend Matrix<t> diag(int n,const t* y);
  //  friend Matrix<t> diag(int n,const t& y);
public:
  Matrix() : mathArrayBase<t>(),N(0),M(0) {}
  Matrix(Subscript n,Subscript m) : mathArrayBase<t>(n*m),N(n),M(m) {}
  Matrix(Subscript n,Subscript m,const t& def) 
      : mathArrayBase<t>(n*m,def),N(n),M(m) {}
  Matrix(Subscript n,Subscript m,const t* defArray) 
      : mathArrayBase<t>(n*m,defArray),N(n),M(m) {}
  Matrix(const Matrix<t>& a);
  Matrix(const arrayBase<t>& a);
  void set(Subscript n,Subscript m) ;
  void set(Subscript n,Subscript m,const t* defArray);
  void set(const Matrix<t>& a);
  void set(const arrayBase<t>& a);
  virtual Subscript rows() const { return N; }
  virtual Subscript cols() const { return M; }
  virtual t& operator()(Subscript j,Subscript k) { return coeff(j-1+(k-1)*N); }
  virtual t  operator()(Subscript j,Subscript k) const { return coeff(j-1+(k-1)*N); }
    
  Matrix<t>& addMatrix(const Matrix<t>& b) ;
  Matrix<t>& operator=(const Matrix<t>& B) { set(B); return *this; }
  Matrix<t>& operator*=(const Matrix<t>& B);
  virtual Matrix<t> concat(const Matrix<t>& B) const;
  virtual Matrix<t> concat(const Vector<t>& B) const;
  virtual Matrix<t> pow(int n) const;
  Matrix<t> transpose() const;
  virtual Matrix<t> inverse() const;
  Matrix<t> deleteRC(int row,int col) const;
  Vector<t> colVector(Subscript i);
  Matrix<t>& gaussAlgorithm(const Matrix<t>& B,Matrix<t>& C) const;
  Matrix<t> eraseZeros(Vector<t>& match) const;
  Vector<t>& gaussAlgorithm(const Vector<t>& B,Vector<t>& D) const;
protected:
  Matrix<t>& multiply(const Matrix<t>& B,Matrix<t>& C) const;
private:
  void writeOut(G4std::ostream& o ) const;
  void addLine(Subscript j,Subscript k,const t& y);
  void swapLines(Subscript j,Subscript k);
  Subscript N,M;
};

template <class t>
class diagonalMatrix : public Matrix<t>
{
  friend diagonalMatrix<t> operator-(const diagonalMatrix<t>& a);
  friend diagonalMatrix<t> operator+(const diagonalMatrix<t>& a,const diagonalMatrix<t>& b);
  friend diagonalMatrix<t> operator-(const diagonalMatrix<t>& a,const diagonalMatrix<t>& b);
  friend diagonalMatrix<t> operator*(const diagonalMatrix<t>& a,const t& y);
  friend diagonalMatrix<t> operator*(const t& y,const diagonalMatrix<t>& a);
public:
  diagonalMatrix() : Matrix<t>() {}
  diagonalMatrix(Subscript n);
  diagonalMatrix(Subscript n,const t& def);
  diagonalMatrix(Subscript n,const t* defArray);
  diagonalMatrix(const Vector<t>& a);
  diagonalMatrix(const Matrix<t>& a);
  Subscript rows() const { return N; }
  Subscript cols() const { return N; }
  t& operator()(Subscript j,Subscript k) ;
  t operator()(Subscript j,Subscript k) const ;
  diagonalMatrix<t> transpose() const;
  Matrix<t> inverse() const;
private:
  Subscript N;
};


template <class t>
class tridiagonalMatrix : public Matrix<t>
{
  friend tridiagonalMatrix<t> operator-(const tridiagonalMatrix<t>& a);
  friend tridiagonalMatrix<t> operator+(const tridiagonalMatrix<t>& a,const tridiagonalMatrix<t>& b);
  friend tridiagonalMatrix<t> operator+(const tridiagonalMatrix<t>& a,const diagonalMatrix<t>& b);
  friend tridiagonalMatrix<t> operator-(const tridiagonalMatrix<t>& a,const tridiagonalMatrix<t>& b);
  friend tridiagonalMatrix<t> operator-(const tridiagonalMatrix<t>& a,const diagonalMatrix<t>& b);
  friend tridiagonalMatrix<t> operator*(const tridiagonalMatrix<t>& a,const t& y);
  friend tridiagonalMatrix<t> operator*(const t& y,const tridiagonalMatrix<t>& a);
public:
  tridiagonalMatrix() : Matrix<t>() {}
  tridiagonalMatrix(Subscript n) : Matrix<t>(),N(n) {arrayBase<t>::set(3*n);}
  tridiagonalMatrix(Subscript n,const t& def);
  tridiagonalMatrix(Subscript n,const t* defArray);
  tridiagonalMatrix(const Matrix<t>& a);
  Subscript rows() const { return N; }
  Subscript cols() const { return N; }
  t& operator()(Subscript j,Subscript k);
  t operator()(Subscript j,Subscript k) const;
  
  tridiagonalMatrix<t> transpose() const;
  Matrix<t> inverse() const;
private:
  Subscript N;
};

#include "matrix_t.tcc"
#endif




























