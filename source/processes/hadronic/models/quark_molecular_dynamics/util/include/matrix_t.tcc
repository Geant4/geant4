#include <iostream.h>
#include <stdlib.h>
#include <stdarg.h>
#ifndef GCC
#include "globals.hh"
#endif 

template<class t>
ostream& operator<<(ostream& o,const arrayBase<t>& a) 
{ 
  a.writeOut(o); 
  return o; 
}

template<class t>
istream& operator>>(istream& in,const arrayBase<t>& a) 
{ 
  a.readIn(in); 
  return in; 
}
template<class t>
int operator==(const arrayBase<t>& a,const arrayBase<t>& b) 
{ 
  return a.isEqual(b); 
}

template<class t>
int operator!=(const arrayBase<t>& a,const arrayBase<t>& b) 
{ 
  return !a.isEqual(b); 
}

template<class t>
inline t& arrayBase<t>::coeff(int i) 
{ 
  if (i>=size) {
    cerr << i << endl;
    throw "arrayBase::coeff: Index out of range"; 
  }
  return x[i]; 
}

template<class t>
inline t arrayBase<t>::coeff(int i) const { 
  if (i>=size) {
    cerr << i << endl;
    throw "arrayBase::coeff: Index out of range"; 
  }
  return x[i]; 
}
  
template<class t>
int arrayBase<t>::isEqual(const arrayBase<t>& b) const
{
  int y = ( rows() == b.rows() && cols() == b.cols() && dim() == b.dim() );
  for (int i=0; i<dim() && y; i++)
    y = y & (coeff(i) == b.coeff(i));
  return y;
}

template<class t>
arrayBase<t>::~arrayBase() 
{ 
  if (x) delete [] x; 
}

template<class t>
void arrayBase<t>::set(int n)
{
  if (size && n!=size)
    throw "arrayBase<t>::set: Incompatible array dimensions";
  size = n;
  if (!x)
    x = new t[size];
  for (int i=0; i<size; i++)
    x[i] = mathType<t>::null;
}

template<class t>
void arrayBase<t>::set(int n,const t& def)
{
  if (size != 0 && n != size) 
    throw "arrayBase<t>::set: Incompatible array dimensions";
  size = n;
  if (!x)
    x = new t[size];
  for (int i=0; i<size; i++)
    x[i] = def;
}

template<class t>
void arrayBase<t>::set(int n,const t* defArray)
{
  if (size != 0 && n != size) 
    throw "arrayBase<t>::set: Incompatible array dimensions";
  size = n;
  if (!x)
    x = new t[size];
  for (int i=0; i<size; i++)
    x[i] = defArray[i];
}

template<class t>
void arrayBase<t>::set(const arrayBase<t>& a)
    {
	if (size != 0 && size != a.size)
	    throw "arrayBase<t>::set: Incompatible array dimensions";
	size = a.dim();
	if (!x)
	    x = new t[size];
	for (int i=0; i<size; i++)
	    x[i] = a.coeff(i);
    }

template<class t>
void arrayBase<t>::writeToStream(ostream& o) const
{
  o << size << " ";
  for (int i=0; i<size; i++)
    o << x[i] << " ";
  o.flush();
}

template<class t>
void arrayBase<t>::readFromStream(istream& in)
{
  in >> size;
  set(size);
  for (int i=0; i<size; i++)
    in >> x[i];
}

template<class t>
arrayBase<t>::arrayBase(int n) : size(n),x(new t[size])
{
  for (int i=0; i<size; i++)
    x[i] = mathType<t>::null;
}
template<class t>
arrayBase<t>::arrayBase(int n,const t& def): size(n),x(new t[size])
{
  for (int i=0; i<size; i++)
    x[i] = def;
}

template<class t>
arrayBase<t>::arrayBase(int n,const t* defArray) : size(n),x(new t[size])
{
  for (int i=0; i<size; i++)
    x[i] = defArray[i];
}

template<class t>
arrayBase<t>::arrayBase(const arrayBase<t>& a) : size(a.dim()),x(new t[size])
{
  for (int i=0; i<size; i++)
    x[i] = a.coeff(i);
}
template<class t>
void arrayBase<t>::writeOut(ostream& o) const
{
  o << "(";
  for (int i=0; i<dim()-1; i++)
    o << coeff(i) << ",";
  o << coeff(dim()-1) << ")";
}

template<class t>
void arrayBase<t>::readIn(istream& i) const
{
  char c;
  t y;
  do {
    i >> c;
  }
  while (i.good() && c == ' ');
  if (i.good())
    if (c=='(') {
      for (int j=0; j<rows(); j++) {
	i >> y >> c;
	if ( ! i.good() ) return; 
	x[j] = y;
      }
      if ( c != ')' ) {
      }
    }
    else 
      i.putback(c);
}

template<class t>
mathArrayBase<t> mathArrayBase<t>::add(const mathArrayBase<t>& b) const
{ 
  mathArrayBase<t> c(*this); 
  c.addSelf(b); 
  return c; 
}

template<class t>
mathArrayBase<t> mathArrayBase<t>::sub(const mathArrayBase<t>& b) const
{ 
  mathArrayBase<t> c(*this); 
  c.subSelf(b); 
  return c; 
}

template<class t>
mathArrayBase<t> mathArrayBase<t>::mul_t(const t& b) const
{ 
  mathArrayBase<t> c(*this); 
  c.mulSelf_t(b); 
  return c; 
}

template<class t>
mathArrayBase<t>& mathArrayBase<t>::addSelf(const mathArrayBase<t>& b)  
{
  for (int i=0; i<dim(); i++) 
    coeff(i) += b.coeff(i);
  return *this;
}

template<class t>
mathArrayBase<t>& mathArrayBase<t>::subSelf(const mathArrayBase<t>& b)  
{
  for (int i=0; i<dim(); i++) 
    coeff(i) -= b.coeff(i);
  return *this;
}

template<class t>
mathArrayBase<t>& mathArrayBase<t>::mulSelf_t(const t& y) 
{
  for (int i=0; i<dim(); i++) 
    coeff(i) *= y;
  return *this;
}

template<class t>
mathArrayBase<t>& mathArrayBase<t>::flipSign() 
{
  for (int i=0; i<dim(); i++) 
    coeff(i) = -coeff(i);
  return *this;
}

template<class t>
t mathArrayBase<t>::norm() const
{
  t s = coeff(0)*coeff(0);   // coeff always starts at i=0!
  for (int i=1; i<dim(); i++)
    s += coeff(i)*coeff(i);
  return s;
}

template<class t>
Vector<t> operator-(const Vector<t>& a)
{
  Vector<t> c(a);
  c.flipSign();
  return c;
}

template<class t>
Vector<t> operator+(const Vector<t>& a,const Vector<t>& b)
{
  if (a.dim() != b.dim())
    throw "Vector::operator+: Incompatible array dimensions";
  Vector<t> c(a);
  c.addSelf(b);
  return c;
}

template<class t>
Vector<t> operator-(const Vector<t>& a,const Vector<t>& b)
{
  if (a.dim() != b.dim())
    throw "Vector::operator-: Incompatible array dimensions";
  Vector<t> c(a);
  c.subSelf(b);
  return c;
}

template<class t>
Vector<t> operator*(const Vector<t>& a,const t& y)
{
  Vector<t> c(a);
  c.mulSelf_t(y);
  return c;
}

template<class t>
Vector<t> operator*(const t& y,const Vector<t>& a)
{
  Vector<t> c(a);
  c.mulSelf_t(y);
  return c;
}

template<class t>
t scalarP(const Vector<t>& a,const Vector<t>& b)
{
  t y = a.coeff(0)*b.coeff(0);
  for (int i=1; i<a.dim(); i++)
    y += a.coeff(i)*b.coeff(i);
  return y;
}

template<class t>
Vector<t>& Vector<t>::operator=(const Vector<t>& y) {
  for (int k=1; k<=y.rows(); k++)
    (*this)[k] = y[k];
  return *this;
}

template<class t>
Vector<t>& Vector<t>::operator+=(const Vector<t>& y) {
  if (dim() != y.dim() )
    throw "Vector::operator+=: Incompatible array dimensions";
  addSelf(y);
  return *this;
}
template<class t>
Matrix<t> Vector<t>::concat(const Matrix<t>& B) const
{
  if (rows() != B.rows())
    throw "Vector::concat: Incompatible array dimensions";
  Matrix<t> C(rows(),cols()+B.cols());
  for (int i=1; i<=rows(); i++) {
    C(i,1) = (*this)[i];
    for (int k=1; k<=B.cols(); k++)
      C(i,k+cols()) = B(i,k);
  }
  return C;
}

template<class t>
Matrix<t> operator-(const Matrix<t>& a)
    {
	Matrix<t> c(a);
	c.flipSign();
	return c;
    }
template<class t>
Matrix<t> operator+(const Matrix<t>& a,const Matrix<t>& b)
{
  if (a.dim() != b.dim())
    throw "Matrix::operator+: Incompatible array dimensions";
  
  Matrix<t> c(a);
  c.addSelf(b);
  return c;
}

template<class t>
Matrix<t> operator-(const Matrix<t>& a,const Matrix<t>& b)
{
  if (a.dim() != b.dim())
    throw "Matrix::operator-: Incompatible array dimensions";
  
  Matrix<t> c(a);
  c.subSelf(b);
  return c;
}

template<class t>
Matrix<t> operator*(const Matrix<t>& a,const t& y)
{
  Matrix<t> c(a);
  c.mulSelf_t(y);
  return c;
}

template<class t>
Matrix<t> operator*(const t& y,const Matrix<t>& a)
{
  Matrix<t> c(a);
  c.mulSelf_t(y);
  return c;
}

template<class t>
Matrix<t> operator*(const Matrix<t>& a,const Matrix<t>& b)
{
  if (a.cols() != b.rows())
    throw "Matrix::operator*: Incompatible array dimensions";
  
  Matrix<t> c(a.rows(),b.cols());
  return a.multiply(b,c);
}
template<class t>
Vector<t> operator*(const Matrix<t>& a,const Vector<t>& b)
{
  if (a.cols() != b.dim())
    throw "Matrix::operator*: Incompatible array dimensions";
  
  Matrix<t> c(a.rows(),1);
  Matrix<t> B(b);
  return (Vector<t>&)a.multiply(B,c);
}

template<class t>
Matrix<t> diag(int n,const t* y)
{
  Matrix<t> C(n,n,mathType<t>::null);
  for (int i=1; i<=n; i++)
    C(i,i) = y[i];
  return C;
}

template<class t>
Matrix<t> diag(int n,const t& y)
{
  Matrix<t> C(n,n,mathType<t>::null);
  for (int i=1; i<=n; i++)
    C(i,i) = y;
  return C;
}

template<class t>
Matrix<t>::Matrix(const Matrix<t>& a) 
  : mathArrayBase<t>(a.rows()*a.cols()),N(a.rows()),M(a.cols()) 
{
  for (int i=1; i<=N; i++) 
    for (int j=1; j<=M; j++)
      (*this)(i,j) = a(i,j);
}

template<class t>
Matrix<t>::Matrix(const arrayBase<t>& a) 
  : mathArrayBase<t>(a),N(a.dim()),M(1) {}

template<class t>
void Matrix<t>::set(int n,int m) 
{ 
  if (dim() && (n!=N || m!=M))
    throw "Matrix already initialized!";
  N = n; M = m; 
  arrayBase<t>::set(n*m); 
}

template<class t>
void Matrix<t>::set(int n,int m,const t* defArray)
{
  if (dim() && (n!=N || m!=M))
    throw "Matrix already initialized!";
  N = n; M = m; 
  arrayBase<t>::set(n*m,defArray);
}

template<class t>
void Matrix<t>::set(const Matrix<t>& a)
{
  if (dim() && (a.rows()!=N || a.cols()!=M))
    throw "Matrix already initialized!";
  N = a.rows(); M = a.cols();
  arrayBase<t>::set(a);
}
template<class t>
void Matrix<t>::set(const arrayBase<t>& a)
{
  if (dim() && (a.dim()!=N || M!=1))
    throw "Matrix already initialized!";
  N = a.dim(); M = 1;
  arrayBase<t>::set(a);
}

template<class t>
Matrix<t>& Matrix<t>::addMatrix(const Matrix<t>& b) 
  {
    for (int i=1; i<=cols(); i++)
      for (int j=1; j<=rows(); j++)
	(*this)(i,j) += b(i,j);
    return *this;
  }
template<class t>
Matrix<t>& Matrix<t>::operator*=(const Matrix<t>& B)
{
  if (!(cols() == rows() && B.cols() == B.rows() && cols() == B.cols()))
    throw "Matrix::operator*=: Incompatible array dimensions";
  Matrix<t> C(*this);
  set(multiply(B,C));
  return *this;
}
template<class t>
Matrix<t> Matrix<t>::concat(const Matrix<t>& B) const
{
  if (rows() != B.rows())
    throw "Matrix::concat: Incompatible array dimensions";
  
  Matrix<t> C(rows(),cols()+B.cols());
  for (int i=1; i<=rows(); i++) {
    for (int j=1; j<=cols(); j++)
      C(i,j) = (*this)(i,j);
    for (int k=1; k<=B.cols(); k++)
      C(i,k+cols()) = B(i,k);
  }
  return C;
}
template<class t>
Matrix<t> Matrix<t>::concat(const Vector<t>& B) const
{
  if (rows() != B.rows())
    throw "Matrix::concat: Incompatible array dimensions";
  
  Matrix<t> C(rows(),cols()+1);
  for (int i=1; i<=rows(); i++) {
    for (int j=1; j<=cols(); j++)
      C(i,j) = (*this)(i,j);
    C(i,1+cols()) = B[i];
  }
  return C;
}
template<class t>
Matrix<t> Matrix<t>::pow(int n) const
{
  if (rows() != cols())
    throw "Matrix::pow: Incompatible array dimensions";
  
  Matrix<t> C;
  if (n==0)
    C = diag(rows(),mathType<t>::one);
  else
    if (n > 0) {
      C.set(*this);
      for (int i=1; i<n; i++)
	C *= *this;
    }
  return C;
}
template<class t>
Matrix<t> Matrix<t>::transpose() const
{
  Matrix<t> C(cols(),rows());
  for (int i=1; i<=rows(); i++)
    for (int j=1; j<=cols(); j++)
      C(j,i) = (*this)(i,j);
  return C;
}
template<class t>
Matrix<t> Matrix<t>::inverse() const
{
  Matrix<t> C(cols(),rows());
  diagonalMatrix<t> B(cols(),1);
  gaussAlgorithm(B,C);
  return C;
}
template<class t>
Matrix<t> Matrix<t>::deleteRC(int row,int col) const
{
  Matrix<t> C(rows()-1,cols()-1);
  int di=0,dj;
  for (int i=1; i<=rows(); i++) {
    dj = 0;
    if (i==row) 
      di = 1;
    else
      for (int j=1; j<=cols(); j++) {
	if (j==col) 
	  dj = 1;
	else
	  C(i-di,j-dj) = (*this)(i,j);
      }
  }
  return C;
}
template<class t>
Vector<t> Matrix<t>::colVector(int i) 
{ 
  Vector<t> c(rows(),&x[(i-1)*rows()]); return c; 
}
template<class t>
Matrix<t>& Matrix<t>::gaussAlgorithm(const Matrix<t>& B,Matrix<t>& C) const
{
  if ( rows() != B.rows() )
    throw "Matrix::gauss: Incompatible array dimensions";
  
  Matrix<t> A = concat(B);
  for (int i=1; i<=rows(); i++) {
    int j = i;
    while ( j<=A.rows() && A(j,i) == 0.0 )
      ++j;
    if ( j>i && j<=A.rows() ) {
      A.swapLines(i,j);
    }
    else
      if (j>A.rows())
	throw "Matrix::gauss: No Pivot";
    
    for (int k=i+1; k<=rows(); k++) 
      A.addLine(k,i,-A(k,i)/A(i,i));
  }
  for (int ii=rows(); ii>=1; ii--)
    for (int j=1; j<=C.cols(); j++) {
      C(ii,j) = A(ii,cols()+j);
      for (int k=ii+1; k<=rows(); k++)
	C(ii,j) -= A(ii,k)*C(k,j);
      C(ii,j) /= A(ii,ii);
    }
  return C;
}
template<class t>
Matrix<t> Matrix<t>::eraseZeros(Vector<t>& match) const
{
  int* z1 = new int[rows()];
  int no = 0;
  for (int l1=1; l1<=rows(); l1++) {
    int l2=1;
    match[l1] = mathType<t>::one;
    while ( l2<=cols() && (*this)(l1,l2) == 0.0 ) ++l2;
    if ( l2>cols() ) {
      int l3=1;
      while ( l3<=rows() && (*this)(l3,l1) == 0.0 ) ++l3;
      if ( l3>rows() ) {
	z1[++no] = l1;
	match[l1] = mathType<t>::null;
      }
      else
	throw "Matrix::gauss::Unterbestimmt!!";
    }
  }
  Matrix<t> C(rows()-no,cols()-no);
  if ( no ) {
    for (int i=1; i<=rows(); i++) {
      for (int j=1; j<=cols(); j++) {
	int di=0,dj=0;
	bool cross = false;
	for (int k=1; k<=no; k++) {
	  if (i>=z1[k]) ++di;
	  if (j>=z1[k]) ++dj;
	  cross = cross || (i==z1[k] || j==z1[k]);
	}
	if (!cross)
	  C(i-di,j-dj) = (*this)(i,j);
      }
    }
  }
  else
    C = *this;
  delete [] z1;
  return C;
}
template<class t>
Vector<t>& Matrix<t>::gaussAlgorithm(const Vector<t>& B,Vector<t>& D) const
{
  if ( rows() != B.rows() )
    throw "Matrix::gauss: Incompatible array dimensions";
  
  Vector<t> match(rows());
  Matrix<t> A = concat(B).eraseZeros(match);
  for (int i=1; i<=A.rows(); i++) {
    int j = i;
    while ( j<=A.rows() && A(j,i) == 0.0 )
      ++j;
    if ( j>i && j<=A.rows() ) {
      A.swapLines(i,j);
    }
    else
      if (j>A.rows()) {
	throw "Matrix::gauss: No Pivot";
      }
    
    for (int k=i+1; k<=A.rows(); k++) 
      A.addLine(k,i,-A(k,i)/A(i,i));
  }
  Vector<t> C(A.rows());
  for (int ii=A.rows(); ii>=1; ii--) {
    C[ii] = A(ii,A.rows()+1);
    for (int k=ii+1; k<=A.rows(); k++)
      C[ii] -= A(ii,k)*C[k];
    C[ii] /= A(ii,ii);
  }
  int djj = 0;
  for (int jj=1; jj<=rows(); jj++) {
    if (match[jj] == 0) {
      ++djj;
      D[jj] = 0;
    }
    else
      D[jj] = C[jj-djj];
  }
  return D;
}
template<class t>
Matrix<t>& Matrix<t>::multiply(const Matrix<t>& B,Matrix<t>& C) const
{
  if (cols() != B.rows())
    throw "Matrix::multiply: Incompatible array dimensions";
  
  for (int i=1; i<=rows(); i++)
    for (int j=1; j<=B.cols(); j++) {
      t s = mathType<t>::null;
      for (int k=1; k<=cols(); k++)
	s += (*this)(i,k)*B(k,j);
      C(i,j) = s;
    }
  return C;
}
template<class t>
void Matrix<t>::writeOut(ostream& o ) const
{
  for (int i=1; i<=rows(); i++) {
    for (int j=1; j<=cols(); j++) {
      o.width(12);
      o << (*this)(i,j) << "  ";
    }
    o << endl;
  }
}
template<class t>
void Matrix<t>::addLine(int j,int k,const t& y)
{
  for (int i=1; i<cols(); i++)
    (*this)(j,i) += (*this)(k,i)*y;
}
template<class t>
void Matrix<t>::swapLines(int j,int k)
{
  for (int i=1; i<cols(); i++) {
    t y = (*this)(j,i);
    (*this)(j,i) = (*this)(k,i);
    (*this)(k,i) = y;
  }
}
template<class t>
diagonalMatrix<t> operator-(const diagonalMatrix<t>& a)
{
  diagonalMatrix<t> c(a);
  c.flipSign();
  return c;
}

template<class t>
diagonalMatrix<t> operator+(const diagonalMatrix<t>& a,const diagonalMatrix<t>& b)
{
  if (a.dim() != b.dim())
    throw "diagonal::operator+: Incompatible array dimensions";
  
  diagonalMatrix<t> c(a);
  c.addSelf(b);
  return c;
}
template<class t>
diagonalMatrix<t> operator-(const diagonalMatrix<t>& a,const diagonalMatrix<t>& b)
{
  if (a.dim() != b.dim())
    throw "diagonal::operator-: Incompatible array dimensions";
  
  diagonalMatrix<t> c(a);
  c.subSelf(b);
  return c;
}
template<class t>
diagonalMatrix<t> operator*(const diagonalMatrix<t>& a,const t& y)
{
  diagonalMatrix<t> c(a);
  c.mulSelf_t(y);
  return c;
}
template<class t>
diagonalMatrix<t> operator*(const t& y,const diagonalMatrix<t>& a)
{
  diagonalMatrix<t> c(a);
  c.mulSelf_t(y);
  return c;
}
template<class t>
diagonalMatrix<t>::diagonalMatrix(int n) : Matrix<t>(),N(n) 
{
  arrayBase<t>::set(n);
}
template<class t>
diagonalMatrix<t>::diagonalMatrix(int n,const t& def): Matrix<t>(),N(n)
{ 
  arrayBase<t>::set(n,def); 
}
template<class t>
diagonalMatrix<t>::diagonalMatrix(int n,const t* defArray): Matrix<t>(),N(n)
{
	arrayBase<t>::set(n,defArray);
}
template<class t>
diagonalMatrix<t>::diagonalMatrix(const Vector<t>& a) : Matrix<t>()
{
  N = a.dim();
  arrayBase<t>::set(N);
  for (int i=1; i<=N; i++)
    (*this)(i,i) = a[i];
}
template<class t>
diagonalMatrix<t>::diagonalMatrix(const Matrix<t>& a)
{
  if (a.rows()!=a.cols())
    throw "diagonal::constructor: Incompatible array dimensions";
  
  N = a.rows();
  arrayBase<t>::set(N);
  for (int i=1; i<=N; i++)
    (*this)(i,i) = a(i,i);
}
template<class t>
t& diagonalMatrix<t>::operator()(int j,int k) 
{ 
  return (j==k) ? this->coeff(j-1) : mathType<t>::null; 
}
template<class t>
t diagonalMatrix<t>::operator()(int j,int k) const 
{ 
  return (j==k) ? this->coeff(j-1) : mathType<t>::null; 
}
template<class t>
diagonalMatrix<t> diagonalMatrix<t>::transpose() const
{
  return *this;
}
template<class t>
Matrix<t> diagonalMatrix<t>::inverse() const
{
  diagonalMatrix<t> I(rows());
  for (int i=0; i<rows(); i++)
    I.coeff(i) = 1.0/coeff(i);
  return I;
}

template <class t>
tridiagonalMatrix<t> operator-(const tridiagonalMatrix<t>& a)
{
  tridiagonalMatrix<t> c(a);
  c.flipSign();
  return c;
}
template<class t>
tridiagonalMatrix<t> operator+(const tridiagonalMatrix<t>& a,const tridiagonalMatrix<t>& b)
{
  if (a.dim() != b.dim())
    throw "tridiag::operator+: Incompatible array dimensions";
  
  tridiagonalMatrix<t> c(a);
  c.addSelf(b);
  return c;
}
template<class t>
tridiagonalMatrix<t> operator+(const tridiagonalMatrix<t>& a,const diagonalMatrix<t>& b)
{
  if (a.rows() != b.rows() || a.cols() != b.cols())
    throw "tridiag::operator+: Incompatible array dimensions";
  
  tridiagonalMatrix<t> c(a);
  c.addMatrix(b);
  return c;
}
template<class t>
tridiagonalMatrix<t> operator-(const tridiagonalMatrix<t>& a,const tridiagonalMatrix<t>& b)
{
  if (a.dim() != b.dim())
    throw "tridiag::operator-: Incompatible array dimensions";
  
  tridiagonalMatrix<t> c(a);
  c.subSelf(b);
  return c;
}
template<class t>
tridiagonalMatrix<t> operator-(const tridiagonalMatrix<t>& a,const diagonalMatrix<t>& b)
{
  if (a.rows() != b.rows() || a.cols() != b.cols())
    throw "tridiag::operator-: Incompatible array dimensions";
  
  tridiagonalMatrix<t> c(a);
  c.subSelf(b);
  return c;
}
template<class t>
tridiagonalMatrix<t> operator*(const tridiagonalMatrix<t>& a,const t& y)
{
  tridiagonalMatrix<t> c(a);
  c.mulSelf_t(y);
  return c;
}
template<class t>
tridiagonalMatrix<t> operator*(const t& y,const tridiagonalMatrix<t>& a)
{
  tridiagonalMatrix<t> c(a);
  c.mulSelf_t(y);
  return c;
}
template<class t>
tridiagonalMatrix<t>::tridiagonalMatrix(int n,const t& def): Matrix<t>(),N(n){ arrayBase<t>::set(3*n,def); }

template<class t>
tridiagonalMatrix<t>::tridiagonalMatrix(int n,const t* defArray): Matrix<t>(),N(n)
{
  arrayBase<t>::set(3*n,defArray);
}
template<class t>
tridiagonalMatrix<t>::tridiagonalMatrix(const Matrix<t>& a)
{
  if (a.rows()!=a.cols())
    throw "tridiag::operator+: Incompatible array dimensions";
  
  N = a.rows();
  arrayBase<t>::set(3*N);
  for (int i=1; i<=N; i++)
    for (int j=i-1; j<=i+1; j++)
      (*this)(i,j) = a(i,j);
}
template<class t>
t&  tridiagonalMatrix<t>::operator()(int j,int k) { return abs((int)j-(int)k)<=1 ? (*this)[j-1+N*(j-k+1)] : mathType<t>::null; }

template<class t>
t  tridiagonalMatrix<t>::operator()(int j,int k) const { return abs((int)j-(int)k)<=1 ? (*this)[j-1+N*(j-k+1)] : mathType<t>::null; }

template<class t>
tridiagonalMatrix<t>  tridiagonalMatrix<t>::transpose() const
{
  tridiagonalMatrix<t> C(*this);
  t* a = &C[0];
  t* b = &C[2*N+1];
  for (int j=0; j<rows(); j++) {
    t y = a[j];
    a[j] = b[j];
    b[j] = y;
  }
  return C;
}
template<class t>
Matrix<t>  tridiagonalMatrix<t>::inverse() const
{
  t* a = &x[N];
  t* b = &x[0];
  t* c = &x[2*N+1];
  Matrix<t> E = diag(N,mathType<t>::one);
  Matrix<t> I(rows(),rows());
  t* d = new t[rows()],* p = new t[rows()],* g = new t[rows()];
  d[0] = a[0];
  for (int j=1; j<rows(); j++) {
    p[j] = -c[j-1]/d[j-1];
    d[j] = a[j]+p[j]*b[j-1];
  }
  for (int i=1; i<=rows(); i++) {
    g[0] = E(i,1);
    for (int j=1; j<rows(); j++)
      g[j] = E(i,j+1)+p[j]*g[j-1];
    I(N,i) = g[N-1]/d[N-1];
    for (int k=N-1; k>=1; k--)
      I(k,i) = (-b[k-1]*I(k+1,i)+g[k-1])/d[k-1];
  }
  return I;
}



























