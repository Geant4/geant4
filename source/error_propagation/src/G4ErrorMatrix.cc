//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file
// ------------------------------------------------------------

#ifdef GNUPRAGMA
#pragma implementation
#endif

#include "globals.hh"
#include <string.h>
#include <float.h>        // for DBL_EPSILON
#include <cmath>
#include <stdlib.h>
#include <iostream>

#include "G4ErrorMatrix.hh"
#include "G4ErrorSymMatrix.hh"

// Simple operation for all elements

#define SIMPLE_UOP(OPER)                            \
    mIter a=m.begin();                      \
    mIter e=m.end();                        \
   for(;a!=e; a++) (*a) OPER t;

#define SIMPLE_BOP(OPER)                            \
    mIter a=m.begin();                      \
    mcIter b=m2.m.begin();                  \
    mIter e=m.end();                        \
   for(;a!=e; a++, b++) (*a) OPER (*b);

#define SIMPLE_TOP(OPER)                            \
    mcIter a=m1.m.begin();       \
    mcIter b=m2.m.begin();       \
    mIter t=mret.m.begin();      \
    mcIter e=m1.m.end();         \
   for(;a!=e; a++, b++, t++) (*t) = (*a) OPER (*b);

// Static functions.

#define CHK_DIM_2(r1,r2,c1,c2,fun) \
   if (r1!=r2 || c1!=c2)  { \
    G4ErrorMatrix::error("Range error in Matrix function " #fun "(1)."); \
   }

#define CHK_DIM_1(c1,r2,fun) \
   if (c1!=r2) { \
     G4ErrorMatrix::error("Range error in Matrix function " #fun "(2)."); \
   }

// Constructors. (Default constructors are inlined and in .icc file)

G4ErrorMatrix::G4ErrorMatrix(int p,int q)
   : m(p*q), nrow(p), ncol(q)
{
  size = nrow * ncol;
}

G4ErrorMatrix::G4ErrorMatrix(int p,int q,int init)
   : m(p*q), nrow(p), ncol(q)
{
   size = nrow * ncol;

   if (size > 0) {
      switch(init)
      {
      case 0:
	 break;

      case 1:
	 {
	    if ( ncol == nrow ) {
 	       mIter a = m.begin();
 	       mIter b = m.end();
	       for( ; a<b; a+=(ncol+1)) *a = 1.0;
	    } else {
	       error("Invalid dimension in G4ErrorMatrix(int,int,1).");
	    }
	    break;
	 }
      default:
	 error("G4ErrorMatrix: initialization must be either 0 or 1.");
      }
   }
}

//
// Destructor
//
G4ErrorMatrix::~G4ErrorMatrix() {
}

G4ErrorMatrix::G4ErrorMatrix(const G4ErrorMatrix &m1)
   : m(m1.size), nrow(m1.nrow), ncol(m1.ncol), size(m1.size)
{
   m = m1.m;

}

G4ErrorMatrix::G4ErrorMatrix(const G4ErrorSymMatrix &m1)
   : m(m1.nrow*m1.nrow), nrow(m1.nrow), ncol(m1.nrow)
{
   size = nrow * ncol;

   int n = ncol;
   mcIter sjk = m1.m.begin();
   mIter m1j = m.begin();
   mIter mj  = m.begin();
   // j >= k
   for(int j=1;j<=nrow;j++) {
      mIter mjk = mj;
      mIter mkj = m1j;
      for(int k=1;k<=j;k++) {
	 *(mjk++) = *sjk;
	 if(j!=k) *mkj = *sjk;
	 sjk++;
	 mkj += n;
      }
      mj += n;
      m1j++;
   }
}

//
//
// Sub matrix
//
//

G4ErrorMatrix G4ErrorMatrix::sub(int min_row, int max_row,
			 int min_col,int max_col) const
#ifdef HEP_GNU_OPTIMIZED_RETURN
return mret(max_row-min_row+1,max_col-min_col+1);
{
#else
{
  G4ErrorMatrix mret(max_row-min_row+1,max_col-min_col+1);
#endif
  if(max_row > num_row() || max_col >num_col())
    error("G4ErrorMatrix::sub: Index out of range");
  mIter a = mret.m.begin();
  int nc = num_col();
  mcIter b1 = m.begin() + (min_row - 1) * nc + min_col - 1;
  
  for(int irow=1; irow<=mret.num_row(); irow++) {
    mcIter brc = b1;
    for(int icol=1; icol<=mret.num_col(); icol++) {
      *(a++) = *(brc++);
    }
    b1 += nc;
  }
  return mret;
}

void G4ErrorMatrix::sub(int row,int col,const G4ErrorMatrix &m1)
{
  if(row <1 || row+m1.num_row()-1 > num_row() || 
     col <1 || col+m1.num_col()-1 > num_col()   )
    error("G4ErrorMatrix::sub: Index out of range");
  mcIter a = m1.m.begin();
  int nc = num_col();
  mIter b1 = m.begin() + (row - 1) * nc + col - 1;
  
  for(int irow=1; irow<=m1.num_row(); irow++) {
    mIter brc = b1;
    for(int icol=1; icol<=m1.num_col(); icol++) {
      *(brc++) = *(a++);
    }
    b1 += nc;
  }
}

//
// Direct sum of two matricies
//

G4ErrorMatrix dsum(const G4ErrorMatrix &m1, const G4ErrorMatrix &m2)
#ifdef HEP_GNU_OPTIMIZED_RETURN
  return mret(m1.num_row() + m2.num_row(), m1.num_col() + m2.num_col(),
	      0);
{
#else
{
  G4ErrorMatrix mret(m1.num_row() + m2.num_row(), m1.num_col() + m2.num_col(),
		 0);
#endif
  mret.sub(1,1,m1);
  mret.sub(m1.num_row()+1,m1.num_col()+1,m2);
  return mret;
}

/* -----------------------------------------------------------------------
   This section contains support routines for matrix.h. This section contains
   The two argument functions +,-. They call the copy constructor and +=,-=.
   ----------------------------------------------------------------------- */
G4ErrorMatrix G4ErrorMatrix::operator- () const 
#ifdef HEP_GNU_OPTIMIZED_RETURN
      return m2(nrow, ncol);
{
#else
{
   G4ErrorMatrix m2(nrow, ncol);
#endif
   mcIter a=m.begin();
   mIter b=m2.m.begin();
   mcIter e=m.end();
   for(;a<e; a++, b++) (*b) = -(*a);
   return m2;
}

   

G4ErrorMatrix operator+(const G4ErrorMatrix &m1,const G4ErrorMatrix &m2)
#ifdef HEP_GNU_OPTIMIZED_RETURN
     return mret(m1.nrow, m1.ncol);
{
#else
{
  G4ErrorMatrix mret(m1.nrow, m1.ncol);
#endif
  CHK_DIM_2(m1.num_row(),m2.num_row(), m1.num_col(),m2.num_col(),+);
  SIMPLE_TOP(+)
    return mret;
}

//
// operator -
//

G4ErrorMatrix operator-(const G4ErrorMatrix &m1,const G4ErrorMatrix &m2)
#ifdef HEP_GNU_OPTIMIZED_RETURN
     return mret(m1.num_row(), m1.num_col());
{
#else
{
  G4ErrorMatrix mret(m1.num_row(), m1.num_col());
#endif
  CHK_DIM_2(m1.num_row(),m2.num_row(),
			 m1.num_col(),m2.num_col(),-);
  SIMPLE_TOP(-)
  return mret;
}

/* -----------------------------------------------------------------------
   This section contains support routines for matrix.h. This file contains
   The two argument functions *,/. They call copy constructor and then /=,*=.
   ----------------------------------------------------------------------- */

G4ErrorMatrix operator/(
const G4ErrorMatrix &m1,double t)
#ifdef HEP_GNU_OPTIMIZED_RETURN
     return mret(m1);
{
#else
{
  G4ErrorMatrix mret(m1);
#endif
  mret /= t;
  return mret;
}

G4ErrorMatrix operator*(const G4ErrorMatrix &m1,double t)
#ifdef HEP_GNU_OPTIMIZED_RETURN
     return mret(m1);
{
#else
{
  G4ErrorMatrix mret(m1);
#endif
  mret *= t;
  return mret;
}

G4ErrorMatrix operator*(double t,const G4ErrorMatrix &m1)
#ifdef HEP_GNU_OPTIMIZED_RETURN
     return mret(m1);
{
#else
{
  G4ErrorMatrix mret(m1);
#endif
  mret *= t;
  return mret;
}

G4ErrorMatrix operator*(const G4ErrorMatrix &m1,const G4ErrorMatrix &m2)
#ifdef HEP_GNU_OPTIMIZED_RETURN
     return mret(m1.nrow,m2.ncol,0);
{
#else
{
  // initialize matrix to 0.0
  G4ErrorMatrix mret(m1.nrow,m2.ncol,0);
#endif
  CHK_DIM_1(m1.ncol,m2.nrow,*);

  int m1cols = m1.ncol;
  int m2cols = m2.ncol;

  for (int i=0; i<m1.nrow; i++)
  {
     for (int j=0; j<m1cols; j++) 
     {
	 double temp = m1.m[i*m1cols+j];
	 mIter pt = mret.m.begin() + i*m2cols;
	
	// Loop over k (the column index in matrix m2)
	mcIter pb = m2.m.begin() + m2cols*j;
	const mcIter pblast = pb + m2cols;
	while (pb < pblast)
	{
	   (*pt) += temp * (*pb);
	   pb++;
	   pt++;
	}
     }
  }

  return mret;
}

/* -----------------------------------------------------------------------
   This section contains the assignment and inplace operators =,+=,-=,*=,/=.
   ----------------------------------------------------------------------- */

G4ErrorMatrix & G4ErrorMatrix::operator+=(const G4ErrorMatrix &m2)
{
  CHK_DIM_2(num_row(),m2.num_row(),num_col(),m2.num_col(),+=);
  SIMPLE_BOP(+=)
  return (*this);
}

G4ErrorMatrix & G4ErrorMatrix::operator-=(const G4ErrorMatrix &m2)
{
  CHK_DIM_2(num_row(),m2.num_row(),num_col(),m2.num_col(),-=);
  SIMPLE_BOP(-=)
  return (*this);
}

G4ErrorMatrix & G4ErrorMatrix::operator/=(double t)
{
  SIMPLE_UOP(/=)
  return (*this);
}

G4ErrorMatrix & G4ErrorMatrix::operator*=(double t)
{
  SIMPLE_UOP(*=)
  return (*this);
}

G4ErrorMatrix & G4ErrorMatrix::operator=(const G4ErrorMatrix &m1)
{
   if(m1.nrow*m1.ncol != size) //??fixme?? m1.size != size
   {
      size = m1.nrow * m1.ncol;
      m.resize(size); //??fixme?? if (size < m1.size) m.resize(m1.size);
   }
   nrow = m1.nrow;
   ncol = m1.ncol;
   m = m1.m;
   return (*this);
}


// Print the G4ErrorMatrix.

std::ostream& operator<<(std::ostream &s, const G4ErrorMatrix &q)
{
  s << "\n";
  /* Fixed format needs 3 extra characters for field, while scientific needs 7 */
  int width;
  if(s.flags() & std::ios::fixed)
    width = s.precision()+3;
  else
    width = s.precision()+7;
  for(int irow = 1; irow<= q.num_row(); irow++)
    {
      for(int icol = 1; icol <= q.num_col(); icol++)
	{
	  s.width(width);
	  s << q(irow,icol) << " ";
	}
      s << G4endl;
    }
  return s;
}

G4ErrorMatrix G4ErrorMatrix::T() const
#ifdef HEP_GNU_OPTIMIZED_RETURN
return mret(ncol,nrow);
{
#else
{
   G4ErrorMatrix mret(ncol,nrow);
#endif
    mcIter pl = m.end();
    mcIter pme = m.begin();
    mIter pt = mret.m.begin();
    mIter ptl = mret.m.end();
   for (; pme < pl; pme++, pt+=nrow)
   {
      if (pt >= ptl) pt -= (size-1);
      (*pt) = (*pme);
   }
   return mret;
}

G4ErrorMatrix G4ErrorMatrix::apply(double (*f)(double, int, int)) const
#ifdef HEP_GNU_OPTIMIZED_RETURN
return mret(num_row(),num_col());
{
#else
{
  G4ErrorMatrix mret(num_row(),num_col());
#endif
  mcIter a = m.begin();
  mIter b = mret.m.begin();
  for(int ir=1;ir<=num_row();ir++) {
    for(int ic=1;ic<=num_col();ic++) {
      *(b++) = (*f)(*(a++), ir, ic);
    }
  }
  return mret;
}

int G4ErrorMatrix::dfinv_matrix(int *ir) {
  if (num_col()!=num_row())
    error("dfinv_matrix: G4ErrorMatrix is not NxN");
   int n = num_col();
  if (n==1) return 0;

  double s31, s32;
   double s33, s34;

  mIter m11 = m.begin();
  mIter m12 = m11 + 1;
  mIter m21 = m11 + n;
  mIter m22 = m12 + n;
  *m21 = -(*m22) * (*m11) * (*m21);
  *m12 = -(*m12);
  if (n>2) {
    mIter mi = m.begin() + 2 * n;
    mIter mii= m.begin() + 2 * n + 2;
    mIter mimim = m.begin() + n + 1;
    for (int i=3;i<=n;i++) {
      int im2 = i - 2;
      mIter mj = m.begin();
      mIter mji = mj + i - 1;
      mIter mij = mi;
      for (int j=1;j<=im2;j++) { 
	s31 = 0.0;
	s32 = *mji;
	mIter mkj = mj + j - 1;
	mIter mik = mi + j - 1;
	mIter mjkp = mj + j;
	mIter mkpi = mj + n + i - 1;
	for (int k=j;k<=im2;k++) {
	  s31 += (*mkj) * (*(mik++));
	  s32 += (*(mjkp++)) * (*mkpi);
	  mkj += n;
	  mkpi += n;
	}
	*mij = -(*mii) * (((*(mij-n)))*( (*(mii-1)))+(s31));
	*mji = -s32;
	mj += n;
	mji += n;
	mij++;
      }
      *(mii-1) = -(*mii) * (*mimim) * (*(mii-1));
      *(mimim+1) = -(*(mimim+1));
      mi += n;
      mimim += (n+1);
      mii += (n+1);
    }
  }
  mIter mi = m.begin();
  mIter mii = m.begin();
  for (int i=1;i<n;i++) {
    int ni = n - i;
    mIter mij = mi;
    int j;
    for (j=1; j<=i;j++) {
      s33 = *mij;
       mIter mikj = mi + n + j - 1;
       mIter miik = mii + 1;
      mIter min_end = mi + n;
      for (;miik<min_end;) {
	s33 += (*mikj) * (*(miik++));
	mikj += n;
      }
      *(mij++) = s33;
    }
    for (j=1;j<=ni;j++) {
      s34 = 0.0;
      mIter miik = mii + j;
      mIter mikij = mii + j * n + j;
      for (int k=j;k<=ni;k++) {
	s34 += *mikij * (*(miik++));
	mikij += n;
      }
      *(mii+j) = s34;
    }
    mi += n;
    mii += (n+1);
  }
  int nxch = ir[n];
  if (nxch==0) return 0;
  for (int mm=1;mm<=nxch;mm++) {
    int k = nxch - mm + 1;
    int ij = ir[k];
    int i = ij >> 12;
    int j = ij%4096;
    mIter mki = m.begin() + i - 1;
    mIter mkj = m.begin() + j - 1;
    for (k=1; k<=n;k++) {
      // 2/24/05 David Sachs fix of improper swap bug that was present
      // for many years:
      double ti = *mki; // 2/24/05
      *mki = *mkj;
      *mkj = ti;	// 2/24/05
      mki += n;
      mkj += n;
    }
  }
  return 0;
}

int G4ErrorMatrix::dfact_matrix(double &det, int *ir) {
  if (ncol!=nrow)
     error("dfact_matrix: G4ErrorMatrix is not NxN");

  int ifail, jfail;
   int n = ncol;

  double tf;
  double g1 = 1.0e-19, g2 = 1.0e19;

  double p, q, t;
  double s11, s12;

  double epsilon = 8*DBL_EPSILON;
  // could be set to zero (like it was before)
  // but then the algorithm often doesn't detect
  // that a matrix is singular

  int normal = 0, imposs = -1;
  int jrange = 0, jover = 1, junder = -1;
  ifail = normal;
  jfail = jrange;
  int nxch = 0;
  det = 1.0;
  mIter mj = m.begin();
  mIter mjj = mj;
  for (int j=1;j<=n;j++) {
    int k = j;
    p = (fabs(*mjj));
    if (j!=n) {
      mIter mij = mj + n + j - 1; 
      for (int i=j+1;i<=n;i++) {
	q = (fabs(*(mij)));
	if (q > p) {
	  k = i;
	  p = q;
	}
	mij += n;
      }
      if (k==j) {
	if (p <= epsilon) {
	  det = 0;
	  ifail = imposs;
	  jfail = jrange;
	  return ifail;
	}
	det = -det; // in this case the sign of the determinant
	            // must not change. So I change it twice. 
      }
      mIter mjl = mj;
      mIter mkl = m.begin() + (k-1)*n;
      for (int l=1;l<=n;l++) {
        tf = *mjl;
        *(mjl++) = *mkl;
        *(mkl++) = tf;
      }
      nxch = nxch + 1;  // this makes the determinant change its sign
      ir[nxch] = (((j)<<12)+(k));
    } else {
      if (p <= epsilon) {
	det = 0.0;
	ifail = imposs;
	jfail = jrange;
	return ifail;
      }
    }
    det *= *mjj;
    *mjj = 1.0 / *mjj;
    t = (fabs(det));
    if (t < g1) {
      det = 0.0;
      if (jfail == jrange) jfail = junder;
    } else if (t > g2) {
      det = 1.0;
      if (jfail==jrange) jfail = jover;
    }
    if (j!=n) {
      mIter mk = mj + n;
      mIter mkjp = mk + j;
      mIter mjk = mj + j;
      for (k=j+1;k<=n;k++) {
	s11 = - (*mjk);
	s12 = - (*mkjp);
	if (j!=1) {
	  mIter mik = m.begin() + k - 1;
	  mIter mijp = m.begin() + j;
	  mIter mki = mk;
	  mIter mji = mj;
	  for (int i=1;i<j;i++) {
	    s11 += (*mik) * (*(mji++));
	    s12 += (*mijp) * (*(mki++));
	    mik += n;
	    mijp += n;
	  }
	}
	*(mjk++) = -s11 * (*mjj);
	*(mkjp) = -(((*(mjj+1)))*((*(mkjp-1)))+(s12));
	mk += n;
	mkjp += n;
      }
    }
    mj += n;
    mjj += (n+1);
  }
  if (nxch%2==1) det = -det;
  if (jfail !=jrange) det = 0.0;
  ir[n] = nxch;
  return 0;
}

void G4ErrorMatrix::invert(int &ierr) {
  if(ncol != nrow)
     error("G4ErrorMatrix::invert: G4ErrorMatrix is not NxN");

  static int max_array = 20;
  static int *ir = new int [max_array+1];

  if (ncol > max_array) {
    delete [] ir;
    max_array = nrow;
    ir = new int [max_array+1];
  }
  double t1, t2, t3;
  double det, temp, s;
  int ifail;
  switch(nrow) {
  case 3:
    double c11,c12,c13,c21,c22,c23,c31,c32,c33;
    ifail = 0;
    c11 = (*(m.begin()+4)) * (*(m.begin()+8)) - (*(m.begin()+5)) * (*(m.begin()+7));
    c12 = (*(m.begin()+5)) * (*(m.begin()+6)) - (*(m.begin()+3)) * (*(m.begin()+8));
    c13 = (*(m.begin()+3)) * (*(m.begin()+7)) - (*(m.begin()+4)) * (*(m.begin()+6));
    c21 = (*(m.begin()+7)) * (*(m.begin()+2)) - (*(m.begin()+8)) * (*(m.begin()+1));
    c22 = (*(m.begin()+8)) * (*m.begin()) - (*(m.begin()+6)) * (*(m.begin()+2));
    c23 = (*(m.begin()+6)) * (*(m.begin()+1)) - (*(m.begin()+7)) * (*m.begin());
    c31 = (*(m.begin()+1)) * (*(m.begin()+5)) - (*(m.begin()+2)) * (*(m.begin()+4));
    c32 = (*(m.begin()+2)) * (*(m.begin()+3)) - (*m.begin()) * (*(m.begin()+5));
    c33 = (*m.begin()) * (*(m.begin()+4)) - (*(m.begin()+1)) * (*(m.begin()+3));
    t1 = fabs(*m.begin());
    t2 = fabs(*(m.begin()+3));
    t3 = fabs(*(m.begin()+6));
    if (t1 >= t2) {
      if (t3 >= t1) {
      temp = *(m.begin()+6);
      det = c23*c12-c22*c13;
      } else {
	temp = *(m.begin());
	det = c22*c33-c23*c32;
      }
    } else if (t3 >= t2) {
      temp = *(m.begin()+6);
      det = c23*c12-c22*c13;
    } else {
      temp = *(m.begin()+3);
      det = c13*c32-c12*c33;
    }
    if (det==0) {
      ierr = 1;
      return;
    }
    {
      double s = temp/det;
      mIter mm = m.begin();
      *(mm++) = s*c11;
      *(mm++) = s*c21;
      *(mm++) = s*c31;
      *(mm++) = s*c12;
      *(mm++) = s*c22;
      *(mm++) = s*c32;
      *(mm++) = s*c13;
      *(mm++) = s*c23;
      *(mm) = s*c33;
    }
    break;
  case 2:
    ifail = 0;
    det = (*m.begin())*(*(m.begin()+3)) - (*(m.begin()+1))*(*(m.begin()+2));
    if (det==0) {
      ierr = 1;
      return;
    }
    s = 1.0/det;
    temp = s*(*(m.begin()+3));
    *(m.begin()+1) *= -s;
    *(m.begin()+2) *= -s;
    *(m.begin()+3) = s*(*m.begin());
    *(m.begin()) = temp;
    break;
  case 1:
    ifail = 0;
    if ((*(m.begin()))==0) {
      ierr = 1;
      return;
    }
    *(m.begin()) = 1.0/(*(m.begin()));
    break;
  case 4:
    invertHaywood4(ierr);
    return;
  case 5:
    invertHaywood5(ierr);
    return;
  case 6:
    invertHaywood6(ierr);
    return;
  default:
    ifail = dfact_matrix(det, ir);
    if(ifail) {
      ierr = 1;
      return;
    }
    dfinv_matrix(ir);
    break;
  }
  ierr = 0;
  return;
}

double G4ErrorMatrix::determinant() const {
  static int max_array = 20;
  static int *ir = new int [max_array+1];
  if(ncol != nrow)
    error("G4ErrorMatrix::determinant: G4ErrorMatrix is not NxN");
  if (ncol > max_array) {
    delete [] ir;
    max_array = nrow;
    ir = new int [max_array+1];
  }
  double det;
  G4ErrorMatrix mt(*this);
  int i = mt.dfact_matrix(det, ir);
  if(i==0) return det;
  return 0;
}

double G4ErrorMatrix::trace() const {
   double t = 0.0;
   for (mcIter d = m.begin(); d < m.end(); d += (ncol+1) )
      t += *d;
   return t;
}

void G4ErrorMatrix::error(const char *s)
{
  G4cerr << s << G4endl;
  G4cerr << "---Exiting to System." << G4endl;
  abort();
}


// Aij are indices for a 6x6 matrix.
// Mij are indices for a 5x5 matrix.
// Fij are indices for a 4x4 matrix.

#define A00 0
#define A01 1
#define A02 2
#define A03 3
#define A04 4
#define A05 5

#define A10 6
#define A11 7
#define A12 8
#define A13 9
#define A14 10
#define A15 11

#define A20 12
#define A21 13
#define A22 14
#define A23 15
#define A24 16
#define A25 17

#define A30 18
#define A31 19
#define A32 20
#define A33 21
#define A34 22
#define A35 23

#define A40 24
#define A41 25
#define A42 26
#define A43 27
#define A44 28
#define A45 29

#define A50 30
#define A51 31
#define A52 32
#define A53 33
#define A54 34
#define A55 35

#define M00 0
#define M01 1
#define M02 2
#define M03 3
#define M04 4

#define M10 5
#define M11 6
#define M12 7
#define M13 8
#define M14 9

#define M20 10
#define M21 11
#define M22 12
#define M23 13
#define M24 14

#define M30 15
#define M31 16
#define M32 17
#define M33 18
#define M34 19

#define M40 20
#define M41 21
#define M42 22
#define M43 23
#define M44 24

#define F00 0
#define F01 1
#define F02 2
#define F03 3

#define F10 4
#define F11 5
#define F12 6
#define F13 7

#define F20 8
#define F21 9
#define F22 10
#define F23 11

#define F30 12
#define F31 13
#define F32 14
#define F33 15


void G4ErrorMatrix::invertHaywood4  (int & ifail) {

  ifail = 0;

  // Find all NECESSARY 2x2 dets:  (18 of them)

  double Det2_12_01 = m[F10]*m[F21] - m[F11]*m[F20];
  double Det2_12_02 = m[F10]*m[F22] - m[F12]*m[F20];
  double Det2_12_03 = m[F10]*m[F23] - m[F13]*m[F20];			//
  double Det2_12_13 = m[F11]*m[F23] - m[F13]*m[F21];			//
  double Det2_12_23 = m[F12]*m[F23] - m[F13]*m[F22];			//
  double Det2_12_12 = m[F11]*m[F22] - m[F12]*m[F21];
  double Det2_13_01 = m[F10]*m[F31] - m[F11]*m[F30];
  double Det2_13_02 = m[F10]*m[F32] - m[F12]*m[F30];
  double Det2_13_03 = m[F10]*m[F33] - m[F13]*m[F30];
  double Det2_13_12 = m[F11]*m[F32] - m[F12]*m[F31];
  double Det2_13_13 = m[F11]*m[F33] - m[F13]*m[F31];
  double Det2_13_23 = m[F12]*m[F33] - m[F13]*m[F32];			//
  double Det2_23_01 = m[F20]*m[F31] - m[F21]*m[F30];
  double Det2_23_02 = m[F20]*m[F32] - m[F22]*m[F30];
  double Det2_23_03 = m[F20]*m[F33] - m[F23]*m[F30];
  double Det2_23_12 = m[F21]*m[F32] - m[F22]*m[F31];
  double Det2_23_13 = m[F21]*m[F33] - m[F23]*m[F31];
  double Det2_23_23 = m[F22]*m[F33] - m[F23]*m[F32];

  // Find all NECESSARY 3x3 dets:   (16 of them)

  double Det3_012_012 = m[F00]*Det2_12_12 - m[F01]*Det2_12_02 
				+ m[F02]*Det2_12_01;
  double Det3_012_013 = m[F00]*Det2_12_13 - m[F01]*Det2_12_03 
				+ m[F03]*Det2_12_01;			//
  double Det3_012_023 = m[F00]*Det2_12_23 - m[F02]*Det2_12_03 
				+ m[F03]*Det2_12_02;			//
  double Det3_012_123 = m[F01]*Det2_12_23 - m[F02]*Det2_12_13 
				+ m[F03]*Det2_12_12;			//
  double Det3_013_012 = m[F00]*Det2_13_12 - m[F01]*Det2_13_02 
				+ m[F02]*Det2_13_01;
  double Det3_013_013 = m[F00]*Det2_13_13 - m[F01]*Det2_13_03
				+ m[F03]*Det2_13_01;
  double Det3_013_023 = m[F00]*Det2_13_23 - m[F02]*Det2_13_03
				+ m[F03]*Det2_13_02;			//
  double Det3_013_123 = m[F01]*Det2_13_23 - m[F02]*Det2_13_13
				+ m[F03]*Det2_13_12;			//
  double Det3_023_012 = m[F00]*Det2_23_12 - m[F01]*Det2_23_02 
				+ m[F02]*Det2_23_01;
  double Det3_023_013 = m[F00]*Det2_23_13 - m[F01]*Det2_23_03
				+ m[F03]*Det2_23_01;
  double Det3_023_023 = m[F00]*Det2_23_23 - m[F02]*Det2_23_03
				+ m[F03]*Det2_23_02;
  double Det3_023_123 = m[F01]*Det2_23_23 - m[F02]*Det2_23_13
				+ m[F03]*Det2_23_12;			//
  double Det3_123_012 = m[F10]*Det2_23_12 - m[F11]*Det2_23_02 
				+ m[F12]*Det2_23_01;
  double Det3_123_013 = m[F10]*Det2_23_13 - m[F11]*Det2_23_03 
				+ m[F13]*Det2_23_01;
  double Det3_123_023 = m[F10]*Det2_23_23 - m[F12]*Det2_23_03 
				+ m[F13]*Det2_23_02;
  double Det3_123_123 = m[F11]*Det2_23_23 - m[F12]*Det2_23_13 
				+ m[F13]*Det2_23_12;

  // Find the 4x4 det:

  double det =    m[F00]*Det3_123_123 
		- m[F01]*Det3_123_023 
		+ m[F02]*Det3_123_013 
		- m[F03]*Det3_123_012;

  if ( det == 0 ) {  
#ifdef SINGULAR_DIAGNOSTICS
    std::cerr << "Kramer's rule inversion of a singular 4x4 matrix: "
	<< *this << "\n";
#endif
    ifail = 1;
    return;
  } 

  double oneOverDet = 1.0/det;
  double mn1OverDet = - oneOverDet;

  m[F00] =  Det3_123_123 * oneOverDet;
  m[F01] =  Det3_023_123 * mn1OverDet;
  m[F02] =  Det3_013_123 * oneOverDet;
  m[F03] =  Det3_012_123 * mn1OverDet;

  m[F10] =  Det3_123_023 * mn1OverDet;
  m[F11] =  Det3_023_023 * oneOverDet;
  m[F12] =  Det3_013_023 * mn1OverDet;
  m[F13] =  Det3_012_023 * oneOverDet;

  m[F20] =  Det3_123_013 * oneOverDet;
  m[F21] =  Det3_023_013 * mn1OverDet;
  m[F22] =  Det3_013_013 * oneOverDet;
  m[F23] =  Det3_012_013 * mn1OverDet;

  m[F30] =  Det3_123_012 * mn1OverDet;
  m[F31] =  Det3_023_012 * oneOverDet;
  m[F32] =  Det3_013_012 * mn1OverDet;
  m[F33] =  Det3_012_012 * oneOverDet;

  return;
}



void G4ErrorMatrix::invertHaywood5  (int & ifail) {

  ifail = 0;

  // Find all NECESSARY 2x2 dets:  (30 of them)

  double Det2_23_01 = m[M20]*m[M31] - m[M21]*m[M30];
  double Det2_23_02 = m[M20]*m[M32] - m[M22]*m[M30];
  double Det2_23_03 = m[M20]*m[M33] - m[M23]*m[M30];
  double Det2_23_04 = m[M20]*m[M34] - m[M24]*m[M30];
  double Det2_23_12 = m[M21]*m[M32] - m[M22]*m[M31];	//
  double Det2_23_13 = m[M21]*m[M33] - m[M23]*m[M31];
  double Det2_23_14 = m[M21]*m[M34] - m[M24]*m[M31];	//
  double Det2_23_23 = m[M22]*m[M33] - m[M23]*m[M32];
  double Det2_23_24 = m[M22]*m[M34] - m[M24]*m[M32];	//
  double Det2_23_34 = m[M23]*m[M34] - m[M24]*m[M33];	//
  double Det2_24_01 = m[M20]*m[M41] - m[M21]*m[M40];
  double Det2_24_02 = m[M20]*m[M42] - m[M22]*m[M40];
  double Det2_24_03 = m[M20]*m[M43] - m[M23]*m[M40];
  double Det2_24_04 = m[M20]*m[M44] - m[M24]*m[M40];
  double Det2_24_12 = m[M21]*m[M42] - m[M22]*m[M41];
  double Det2_24_13 = m[M21]*m[M43] - m[M23]*m[M41];
  double Det2_24_14 = m[M21]*m[M44] - m[M24]*m[M41];
  double Det2_24_23 = m[M22]*m[M43] - m[M23]*m[M42];
  double Det2_24_24 = m[M22]*m[M44] - m[M24]*m[M42];
  double Det2_24_34 = m[M23]*m[M44] - m[M24]*m[M43];	//
  double Det2_34_01 = m[M30]*m[M41] - m[M31]*m[M40];
  double Det2_34_02 = m[M30]*m[M42] - m[M32]*m[M40];
  double Det2_34_03 = m[M30]*m[M43] - m[M33]*m[M40];
  double Det2_34_04 = m[M30]*m[M44] - m[M34]*m[M40];
  double Det2_34_12 = m[M31]*m[M42] - m[M32]*m[M41];
  double Det2_34_13 = m[M31]*m[M43] - m[M33]*m[M41];
  double Det2_34_14 = m[M31]*m[M44] - m[M34]*m[M41];
  double Det2_34_23 = m[M32]*m[M43] - m[M33]*m[M42];
  double Det2_34_24 = m[M32]*m[M44] - m[M34]*m[M42];
  double Det2_34_34 = m[M33]*m[M44] - m[M34]*m[M43];

  // Find all NECESSARY 3x3 dets:   (40 of them)

  double Det3_123_012 = m[M10]*Det2_23_12 - m[M11]*Det2_23_02 
				+ m[M12]*Det2_23_01;
  double Det3_123_013 = m[M10]*Det2_23_13 - m[M11]*Det2_23_03 
				+ m[M13]*Det2_23_01;
  double Det3_123_014 = m[M10]*Det2_23_14 - m[M11]*Det2_23_04 
				+ m[M14]*Det2_23_01;		//
  double Det3_123_023 = m[M10]*Det2_23_23 - m[M12]*Det2_23_03 
				+ m[M13]*Det2_23_02;
  double Det3_123_024 = m[M10]*Det2_23_24 - m[M12]*Det2_23_04 
				+ m[M14]*Det2_23_02;		//
  double Det3_123_034 = m[M10]*Det2_23_34 - m[M13]*Det2_23_04 
				+ m[M14]*Det2_23_03;		//
  double Det3_123_123 = m[M11]*Det2_23_23 - m[M12]*Det2_23_13 
				+ m[M13]*Det2_23_12;
  double Det3_123_124 = m[M11]*Det2_23_24 - m[M12]*Det2_23_14 
				+ m[M14]*Det2_23_12;		//
  double Det3_123_134 = m[M11]*Det2_23_34 - m[M13]*Det2_23_14 
				+ m[M14]*Det2_23_13;		//
  double Det3_123_234 = m[M12]*Det2_23_34 - m[M13]*Det2_23_24 
				+ m[M14]*Det2_23_23;		//
  double Det3_124_012 = m[M10]*Det2_24_12 - m[M11]*Det2_24_02 
				+ m[M12]*Det2_24_01;
  double Det3_124_013 = m[M10]*Det2_24_13 - m[M11]*Det2_24_03 
				+ m[M13]*Det2_24_01;
  double Det3_124_014 = m[M10]*Det2_24_14 - m[M11]*Det2_24_04 
				+ m[M14]*Det2_24_01;
  double Det3_124_023 = m[M10]*Det2_24_23 - m[M12]*Det2_24_03 
				+ m[M13]*Det2_24_02;
  double Det3_124_024 = m[M10]*Det2_24_24 - m[M12]*Det2_24_04 
				+ m[M14]*Det2_24_02;
  double Det3_124_034 = m[M10]*Det2_24_34 - m[M13]*Det2_24_04 
				+ m[M14]*Det2_24_03;		//
  double Det3_124_123 = m[M11]*Det2_24_23 - m[M12]*Det2_24_13 
				+ m[M13]*Det2_24_12;
  double Det3_124_124 = m[M11]*Det2_24_24 - m[M12]*Det2_24_14 
				+ m[M14]*Det2_24_12;
  double Det3_124_134 = m[M11]*Det2_24_34 - m[M13]*Det2_24_14 
				+ m[M14]*Det2_24_13;		//
  double Det3_124_234 = m[M12]*Det2_24_34 - m[M13]*Det2_24_24 
				+ m[M14]*Det2_24_23;		//
  double Det3_134_012 = m[M10]*Det2_34_12 - m[M11]*Det2_34_02 
				+ m[M12]*Det2_34_01;
  double Det3_134_013 = m[M10]*Det2_34_13 - m[M11]*Det2_34_03 
				+ m[M13]*Det2_34_01;
  double Det3_134_014 = m[M10]*Det2_34_14 - m[M11]*Det2_34_04 
				+ m[M14]*Det2_34_01;
  double Det3_134_023 = m[M10]*Det2_34_23 - m[M12]*Det2_34_03 
				+ m[M13]*Det2_34_02;
  double Det3_134_024 = m[M10]*Det2_34_24 - m[M12]*Det2_34_04 
				+ m[M14]*Det2_34_02;
  double Det3_134_034 = m[M10]*Det2_34_34 - m[M13]*Det2_34_04 
				+ m[M14]*Det2_34_03;
  double Det3_134_123 = m[M11]*Det2_34_23 - m[M12]*Det2_34_13 
				+ m[M13]*Det2_34_12;
  double Det3_134_124 = m[M11]*Det2_34_24 - m[M12]*Det2_34_14 
				+ m[M14]*Det2_34_12;
  double Det3_134_134 = m[M11]*Det2_34_34 - m[M13]*Det2_34_14 
				+ m[M14]*Det2_34_13;
  double Det3_134_234 = m[M12]*Det2_34_34 - m[M13]*Det2_34_24 
				+ m[M14]*Det2_34_23;		//
  double Det3_234_012 = m[M20]*Det2_34_12 - m[M21]*Det2_34_02 
				+ m[M22]*Det2_34_01;
  double Det3_234_013 = m[M20]*Det2_34_13 - m[M21]*Det2_34_03 
				+ m[M23]*Det2_34_01;
  double Det3_234_014 = m[M20]*Det2_34_14 - m[M21]*Det2_34_04 
				+ m[M24]*Det2_34_01;
  double Det3_234_023 = m[M20]*Det2_34_23 - m[M22]*Det2_34_03 
				+ m[M23]*Det2_34_02;
  double Det3_234_024 = m[M20]*Det2_34_24 - m[M22]*Det2_34_04 
				+ m[M24]*Det2_34_02;
  double Det3_234_034 = m[M20]*Det2_34_34 - m[M23]*Det2_34_04 
				+ m[M24]*Det2_34_03;
  double Det3_234_123 = m[M21]*Det2_34_23 - m[M22]*Det2_34_13 
				+ m[M23]*Det2_34_12;
  double Det3_234_124 = m[M21]*Det2_34_24 - m[M22]*Det2_34_14 
				+ m[M24]*Det2_34_12;
  double Det3_234_134 = m[M21]*Det2_34_34 - m[M23]*Det2_34_14 
				+ m[M24]*Det2_34_13;
  double Det3_234_234 = m[M22]*Det2_34_34 - m[M23]*Det2_34_24 
				+ m[M24]*Det2_34_23;

  // Find all NECESSARY 4x4 dets:   (25 of them)

  double Det4_0123_0123 = m[M00]*Det3_123_123 - m[M01]*Det3_123_023 
				+ m[M02]*Det3_123_013 - m[M03]*Det3_123_012;
  double Det4_0123_0124 = m[M00]*Det3_123_124 - m[M01]*Det3_123_024 
				+ m[M02]*Det3_123_014 - m[M04]*Det3_123_012; //
  double Det4_0123_0134 = m[M00]*Det3_123_134 - m[M01]*Det3_123_034 
				+ m[M03]*Det3_123_014 - m[M04]*Det3_123_013; //
  double Det4_0123_0234 = m[M00]*Det3_123_234 - m[M02]*Det3_123_034 
				+ m[M03]*Det3_123_024 - m[M04]*Det3_123_023; //
  double Det4_0123_1234 = m[M01]*Det3_123_234 - m[M02]*Det3_123_134 
				+ m[M03]*Det3_123_124 - m[M04]*Det3_123_123; //
  double Det4_0124_0123 = m[M00]*Det3_124_123 - m[M01]*Det3_124_023 
				+ m[M02]*Det3_124_013 - m[M03]*Det3_124_012;
  double Det4_0124_0124 = m[M00]*Det3_124_124 - m[M01]*Det3_124_024 
				+ m[M02]*Det3_124_014 - m[M04]*Det3_124_012;
  double Det4_0124_0134 = m[M00]*Det3_124_134 - m[M01]*Det3_124_034 
				+ m[M03]*Det3_124_014 - m[M04]*Det3_124_013; //
  double Det4_0124_0234 = m[M00]*Det3_124_234 - m[M02]*Det3_124_034 
				+ m[M03]*Det3_124_024 - m[M04]*Det3_124_023; //
  double Det4_0124_1234 = m[M01]*Det3_124_234 - m[M02]*Det3_124_134 
				+ m[M03]*Det3_124_124 - m[M04]*Det3_124_123; //
  double Det4_0134_0123 = m[M00]*Det3_134_123 - m[M01]*Det3_134_023 
				+ m[M02]*Det3_134_013 - m[M03]*Det3_134_012;
  double Det4_0134_0124 = m[M00]*Det3_134_124 - m[M01]*Det3_134_024 
				+ m[M02]*Det3_134_014 - m[M04]*Det3_134_012;
  double Det4_0134_0134 = m[M00]*Det3_134_134 - m[M01]*Det3_134_034 
				+ m[M03]*Det3_134_014 - m[M04]*Det3_134_013;
  double Det4_0134_0234 = m[M00]*Det3_134_234 - m[M02]*Det3_134_034 
				+ m[M03]*Det3_134_024 - m[M04]*Det3_134_023; //
  double Det4_0134_1234 = m[M01]*Det3_134_234 - m[M02]*Det3_134_134 
				+ m[M03]*Det3_134_124 - m[M04]*Det3_134_123; //
  double Det4_0234_0123 = m[M00]*Det3_234_123 - m[M01]*Det3_234_023 
				+ m[M02]*Det3_234_013 - m[M03]*Det3_234_012;
  double Det4_0234_0124 = m[M00]*Det3_234_124 - m[M01]*Det3_234_024 
				+ m[M02]*Det3_234_014 - m[M04]*Det3_234_012;
  double Det4_0234_0134 = m[M00]*Det3_234_134 - m[M01]*Det3_234_034 
				+ m[M03]*Det3_234_014 - m[M04]*Det3_234_013;
  double Det4_0234_0234 = m[M00]*Det3_234_234 - m[M02]*Det3_234_034 
				+ m[M03]*Det3_234_024 - m[M04]*Det3_234_023;
  double Det4_0234_1234 = m[M01]*Det3_234_234 - m[M02]*Det3_234_134 
				+ m[M03]*Det3_234_124 - m[M04]*Det3_234_123; //
  double Det4_1234_0123 = m[M10]*Det3_234_123 - m[M11]*Det3_234_023 
				+ m[M12]*Det3_234_013 - m[M13]*Det3_234_012;
  double Det4_1234_0124 = m[M10]*Det3_234_124 - m[M11]*Det3_234_024 
				+ m[M12]*Det3_234_014 - m[M14]*Det3_234_012;
  double Det4_1234_0134 = m[M10]*Det3_234_134 - m[M11]*Det3_234_034 
				+ m[M13]*Det3_234_014 - m[M14]*Det3_234_013;
  double Det4_1234_0234 = m[M10]*Det3_234_234 - m[M12]*Det3_234_034 
				+ m[M13]*Det3_234_024 - m[M14]*Det3_234_023;
  double Det4_1234_1234 = m[M11]*Det3_234_234 - m[M12]*Det3_234_134 
				+ m[M13]*Det3_234_124 - m[M14]*Det3_234_123;

  // Find the 5x5 det:

  double det =    m[M00]*Det4_1234_1234 
	 	- m[M01]*Det4_1234_0234 
		+ m[M02]*Det4_1234_0134 
		- m[M03]*Det4_1234_0124 
		+ m[M04]*Det4_1234_0123;

  if ( det == 0 ) {  
#ifdef SINGULAR_DIAGNOSTICS
    std::cerr << "Kramer's rule inversion of a singular 5x5 matrix: "
	<< *this << "\n";
#endif
    ifail = 1;
    return;
  } 

  double oneOverDet = 1.0/det;
  double mn1OverDet = - oneOverDet;

  m[M00] =  Det4_1234_1234 * oneOverDet;
  m[M01] =  Det4_0234_1234 * mn1OverDet;
  m[M02] =  Det4_0134_1234 * oneOverDet;
  m[M03] =  Det4_0124_1234 * mn1OverDet;
  m[M04] =  Det4_0123_1234 * oneOverDet;

  m[M10] =  Det4_1234_0234 * mn1OverDet;
  m[M11] =  Det4_0234_0234 * oneOverDet;
  m[M12] =  Det4_0134_0234 * mn1OverDet;
  m[M13] =  Det4_0124_0234 * oneOverDet;
  m[M14] =  Det4_0123_0234 * mn1OverDet;

  m[M20] =  Det4_1234_0134 * oneOverDet;
  m[M21] =  Det4_0234_0134 * mn1OverDet;
  m[M22] =  Det4_0134_0134 * oneOverDet;
  m[M23] =  Det4_0124_0134 * mn1OverDet;
  m[M24] =  Det4_0123_0134 * oneOverDet;

  m[M30] =  Det4_1234_0124 * mn1OverDet;
  m[M31] =  Det4_0234_0124 * oneOverDet;
  m[M32] =  Det4_0134_0124 * mn1OverDet;
  m[M33] =  Det4_0124_0124 * oneOverDet;
  m[M34] =  Det4_0123_0124 * mn1OverDet;

  m[M40] =  Det4_1234_0123 * oneOverDet;
  m[M41] =  Det4_0234_0123 * mn1OverDet;
  m[M42] =  Det4_0134_0123 * oneOverDet;
  m[M43] =  Det4_0124_0123 * mn1OverDet;
  m[M44] =  Det4_0123_0123 * oneOverDet;

  return;
}

void G4ErrorMatrix::invertHaywood6  (int & ifail) {

  ifail = 0;

  // Find all NECESSARY 2x2 dets:  (45 of them)

  double Det2_34_01 = m[A30]*m[A41] - m[A31]*m[A40];
  double Det2_34_02 = m[A30]*m[A42] - m[A32]*m[A40];
  double Det2_34_03 = m[A30]*m[A43] - m[A33]*m[A40];
  double Det2_34_04 = m[A30]*m[A44] - m[A34]*m[A40];
  double Det2_34_05 = m[A30]*m[A45] - m[A35]*m[A40];		//
  double Det2_34_12 = m[A31]*m[A42] - m[A32]*m[A41];
  double Det2_34_13 = m[A31]*m[A43] - m[A33]*m[A41];
  double Det2_34_14 = m[A31]*m[A44] - m[A34]*m[A41];
  double Det2_34_15 = m[A31]*m[A45] - m[A35]*m[A41];		//
  double Det2_34_23 = m[A32]*m[A43] - m[A33]*m[A42];
  double Det2_34_24 = m[A32]*m[A44] - m[A34]*m[A42];
  double Det2_34_25 = m[A32]*m[A45] - m[A35]*m[A42];		//
  double Det2_34_34 = m[A33]*m[A44] - m[A34]*m[A43];
  double Det2_34_35 = m[A33]*m[A45] - m[A35]*m[A43];		//
  double Det2_34_45 = m[A34]*m[A45] - m[A35]*m[A44];		//
  double Det2_35_01 = m[A30]*m[A51] - m[A31]*m[A50];
  double Det2_35_02 = m[A30]*m[A52] - m[A32]*m[A50];
  double Det2_35_03 = m[A30]*m[A53] - m[A33]*m[A50];
  double Det2_35_04 = m[A30]*m[A54] - m[A34]*m[A50];
  double Det2_35_05 = m[A30]*m[A55] - m[A35]*m[A50];
  double Det2_35_12 = m[A31]*m[A52] - m[A32]*m[A51];
  double Det2_35_13 = m[A31]*m[A53] - m[A33]*m[A51];
  double Det2_35_14 = m[A31]*m[A54] - m[A34]*m[A51];
  double Det2_35_15 = m[A31]*m[A55] - m[A35]*m[A51];
  double Det2_35_23 = m[A32]*m[A53] - m[A33]*m[A52];
  double Det2_35_24 = m[A32]*m[A54] - m[A34]*m[A52];
  double Det2_35_25 = m[A32]*m[A55] - m[A35]*m[A52];
  double Det2_35_34 = m[A33]*m[A54] - m[A34]*m[A53];
  double Det2_35_35 = m[A33]*m[A55] - m[A35]*m[A53];
  double Det2_35_45 = m[A34]*m[A55] - m[A35]*m[A54];		//
  double Det2_45_01 = m[A40]*m[A51] - m[A41]*m[A50];
  double Det2_45_02 = m[A40]*m[A52] - m[A42]*m[A50];
  double Det2_45_03 = m[A40]*m[A53] - m[A43]*m[A50];
  double Det2_45_04 = m[A40]*m[A54] - m[A44]*m[A50];
  double Det2_45_05 = m[A40]*m[A55] - m[A45]*m[A50];
  double Det2_45_12 = m[A41]*m[A52] - m[A42]*m[A51];
  double Det2_45_13 = m[A41]*m[A53] - m[A43]*m[A51];
  double Det2_45_14 = m[A41]*m[A54] - m[A44]*m[A51];
  double Det2_45_15 = m[A41]*m[A55] - m[A45]*m[A51];
  double Det2_45_23 = m[A42]*m[A53] - m[A43]*m[A52];
  double Det2_45_24 = m[A42]*m[A54] - m[A44]*m[A52];
  double Det2_45_25 = m[A42]*m[A55] - m[A45]*m[A52];
  double Det2_45_34 = m[A43]*m[A54] - m[A44]*m[A53];
  double Det2_45_35 = m[A43]*m[A55] - m[A45]*m[A53];
  double Det2_45_45 = m[A44]*m[A55] - m[A45]*m[A54];

  // Find all NECESSARY 3x3 dets:  (80 of them)

  double Det3_234_012 = m[A20]*Det2_34_12 - m[A21]*Det2_34_02 
						+ m[A22]*Det2_34_01;
  double Det3_234_013 = m[A20]*Det2_34_13 - m[A21]*Det2_34_03 
						+ m[A23]*Det2_34_01;
  double Det3_234_014 = m[A20]*Det2_34_14 - m[A21]*Det2_34_04 
						+ m[A24]*Det2_34_01;
  double Det3_234_015 = m[A20]*Det2_34_15 - m[A21]*Det2_34_05
						+ m[A25]*Det2_34_01;	//
  double Det3_234_023 = m[A20]*Det2_34_23 - m[A22]*Det2_34_03 
						+ m[A23]*Det2_34_02;
  double Det3_234_024 = m[A20]*Det2_34_24 - m[A22]*Det2_34_04 
						+ m[A24]*Det2_34_02;
  double Det3_234_025 = m[A20]*Det2_34_25 - m[A22]*Det2_34_05 
						+ m[A25]*Det2_34_02;	//
  double Det3_234_034 = m[A20]*Det2_34_34 - m[A23]*Det2_34_04 
						+ m[A24]*Det2_34_03;
  double Det3_234_035 = m[A20]*Det2_34_35 - m[A23]*Det2_34_05 		
						+ m[A25]*Det2_34_03;	//
  double Det3_234_045 = m[A20]*Det2_34_45 - m[A24]*Det2_34_05 
						+ m[A25]*Det2_34_04;	//
  double Det3_234_123 = m[A21]*Det2_34_23 - m[A22]*Det2_34_13 
						+ m[A23]*Det2_34_12;
  double Det3_234_124 = m[A21]*Det2_34_24 - m[A22]*Det2_34_14 
						+ m[A24]*Det2_34_12;
  double Det3_234_125 = m[A21]*Det2_34_25 - m[A22]*Det2_34_15 
						+ m[A25]*Det2_34_12;	//
  double Det3_234_134 = m[A21]*Det2_34_34 - m[A23]*Det2_34_14 
						+ m[A24]*Det2_34_13;
  double Det3_234_135 = m[A21]*Det2_34_35 - m[A23]*Det2_34_15 
						+ m[A25]*Det2_34_13;	//
  double Det3_234_145 = m[A21]*Det2_34_45 - m[A24]*Det2_34_15 
						+ m[A25]*Det2_34_14;	//
  double Det3_234_234 = m[A22]*Det2_34_34 - m[A23]*Det2_34_24 
						+ m[A24]*Det2_34_23;
  double Det3_234_235 = m[A22]*Det2_34_35 - m[A23]*Det2_34_25
						+ m[A25]*Det2_34_23;	//
  double Det3_234_245 = m[A22]*Det2_34_45 - m[A24]*Det2_34_25
						+ m[A25]*Det2_34_24;	//
  double Det3_234_345 = m[A23]*Det2_34_45 - m[A24]*Det2_34_35
						+ m[A25]*Det2_34_34;	//
  double Det3_235_012 = m[A20]*Det2_35_12 - m[A21]*Det2_35_02 
						+ m[A22]*Det2_35_01;
  double Det3_235_013 = m[A20]*Det2_35_13 - m[A21]*Det2_35_03 
						+ m[A23]*Det2_35_01;
  double Det3_235_014 = m[A20]*Det2_35_14 - m[A21]*Det2_35_04 
						+ m[A24]*Det2_35_01;
  double Det3_235_015 = m[A20]*Det2_35_15 - m[A21]*Det2_35_05 
						+ m[A25]*Det2_35_01;
  double Det3_235_023 = m[A20]*Det2_35_23 - m[A22]*Det2_35_03 
						+ m[A23]*Det2_35_02;
  double Det3_235_024 = m[A20]*Det2_35_24 - m[A22]*Det2_35_04 
						+ m[A24]*Det2_35_02;
  double Det3_235_025 = m[A20]*Det2_35_25 - m[A22]*Det2_35_05 
						+ m[A25]*Det2_35_02;
  double Det3_235_034 = m[A20]*Det2_35_34 - m[A23]*Det2_35_04 
						+ m[A24]*Det2_35_03;
  double Det3_235_035 = m[A20]*Det2_35_35 - m[A23]*Det2_35_05 
						+ m[A25]*Det2_35_03;
  double Det3_235_045 = m[A20]*Det2_35_45 - m[A24]*Det2_35_05 
						+ m[A25]*Det2_35_04;	//
  double Det3_235_123 = m[A21]*Det2_35_23 - m[A22]*Det2_35_13 
						+ m[A23]*Det2_35_12;
  double Det3_235_124 = m[A21]*Det2_35_24 - m[A22]*Det2_35_14 
						+ m[A24]*Det2_35_12;
  double Det3_235_125 = m[A21]*Det2_35_25 - m[A22]*Det2_35_15 
						+ m[A25]*Det2_35_12;
  double Det3_235_134 = m[A21]*Det2_35_34 - m[A23]*Det2_35_14 
						+ m[A24]*Det2_35_13;
  double Det3_235_135 = m[A21]*Det2_35_35 - m[A23]*Det2_35_15 
						+ m[A25]*Det2_35_13;
  double Det3_235_145 = m[A21]*Det2_35_45 - m[A24]*Det2_35_15 
						+ m[A25]*Det2_35_14;	//
  double Det3_235_234 = m[A22]*Det2_35_34 - m[A23]*Det2_35_24 
						+ m[A24]*Det2_35_23;
  double Det3_235_235 = m[A22]*Det2_35_35 - m[A23]*Det2_35_25 
						+ m[A25]*Det2_35_23;
  double Det3_235_245 = m[A22]*Det2_35_45 - m[A24]*Det2_35_25 
						+ m[A25]*Det2_35_24;	//
  double Det3_235_345 = m[A23]*Det2_35_45 - m[A24]*Det2_35_35 
						+ m[A25]*Det2_35_34;	//
  double Det3_245_012 = m[A20]*Det2_45_12 - m[A21]*Det2_45_02 
						+ m[A22]*Det2_45_01;
  double Det3_245_013 = m[A20]*Det2_45_13 - m[A21]*Det2_45_03 
						+ m[A23]*Det2_45_01;
  double Det3_245_014 = m[A20]*Det2_45_14 - m[A21]*Det2_45_04 
						+ m[A24]*Det2_45_01;
  double Det3_245_015 = m[A20]*Det2_45_15 - m[A21]*Det2_45_05 
						+ m[A25]*Det2_45_01;
  double Det3_245_023 = m[A20]*Det2_45_23 - m[A22]*Det2_45_03 
						+ m[A23]*Det2_45_02;
  double Det3_245_024 = m[A20]*Det2_45_24 - m[A22]*Det2_45_04 
						+ m[A24]*Det2_45_02;
  double Det3_245_025 = m[A20]*Det2_45_25 - m[A22]*Det2_45_05 
						+ m[A25]*Det2_45_02;
  double Det3_245_034 = m[A20]*Det2_45_34 - m[A23]*Det2_45_04 
						+ m[A24]*Det2_45_03;
  double Det3_245_035 = m[A20]*Det2_45_35 - m[A23]*Det2_45_05 
						+ m[A25]*Det2_45_03;
  double Det3_245_045 = m[A20]*Det2_45_45 - m[A24]*Det2_45_05 
						+ m[A25]*Det2_45_04;
  double Det3_245_123 = m[A21]*Det2_45_23 - m[A22]*Det2_45_13 
						+ m[A23]*Det2_45_12;
  double Det3_245_124 = m[A21]*Det2_45_24 - m[A22]*Det2_45_14 
						+ m[A24]*Det2_45_12;
  double Det3_245_125 = m[A21]*Det2_45_25 - m[A22]*Det2_45_15 
						+ m[A25]*Det2_45_12;
  double Det3_245_134 = m[A21]*Det2_45_34 - m[A23]*Det2_45_14 
						+ m[A24]*Det2_45_13;
  double Det3_245_135 = m[A21]*Det2_45_35 - m[A23]*Det2_45_15 
						+ m[A25]*Det2_45_13;
  double Det3_245_145 = m[A21]*Det2_45_45 - m[A24]*Det2_45_15 
						+ m[A25]*Det2_45_14;
  double Det3_245_234 = m[A22]*Det2_45_34 - m[A23]*Det2_45_24 
						+ m[A24]*Det2_45_23;
  double Det3_245_235 = m[A22]*Det2_45_35 - m[A23]*Det2_45_25 
						+ m[A25]*Det2_45_23;
  double Det3_245_245 = m[A22]*Det2_45_45 - m[A24]*Det2_45_25 
						+ m[A25]*Det2_45_24;
  double Det3_245_345 = m[A23]*Det2_45_45 - m[A24]*Det2_45_35 
						+ m[A25]*Det2_45_34;	//
  double Det3_345_012 = m[A30]*Det2_45_12 - m[A31]*Det2_45_02 
						+ m[A32]*Det2_45_01;
  double Det3_345_013 = m[A30]*Det2_45_13 - m[A31]*Det2_45_03 
						+ m[A33]*Det2_45_01;
  double Det3_345_014 = m[A30]*Det2_45_14 - m[A31]*Det2_45_04 
						+ m[A34]*Det2_45_01;
  double Det3_345_015 = m[A30]*Det2_45_15 - m[A31]*Det2_45_05 
						+ m[A35]*Det2_45_01;
  double Det3_345_023 = m[A30]*Det2_45_23 - m[A32]*Det2_45_03 
						+ m[A33]*Det2_45_02;
  double Det3_345_024 = m[A30]*Det2_45_24 - m[A32]*Det2_45_04 
						+ m[A34]*Det2_45_02;
  double Det3_345_025 = m[A30]*Det2_45_25 - m[A32]*Det2_45_05 
						+ m[A35]*Det2_45_02;
  double Det3_345_034 = m[A30]*Det2_45_34 - m[A33]*Det2_45_04 
						+ m[A34]*Det2_45_03;
  double Det3_345_035 = m[A30]*Det2_45_35 - m[A33]*Det2_45_05 
						+ m[A35]*Det2_45_03;
  double Det3_345_045 = m[A30]*Det2_45_45 - m[A34]*Det2_45_05 
						+ m[A35]*Det2_45_04;
  double Det3_345_123 = m[A31]*Det2_45_23 - m[A32]*Det2_45_13 
						+ m[A33]*Det2_45_12;
  double Det3_345_124 = m[A31]*Det2_45_24 - m[A32]*Det2_45_14 
						+ m[A34]*Det2_45_12;
  double Det3_345_125 = m[A31]*Det2_45_25 - m[A32]*Det2_45_15 
						+ m[A35]*Det2_45_12;
  double Det3_345_134 = m[A31]*Det2_45_34 - m[A33]*Det2_45_14 
						+ m[A34]*Det2_45_13;
  double Det3_345_135 = m[A31]*Det2_45_35 - m[A33]*Det2_45_15 
						+ m[A35]*Det2_45_13;
  double Det3_345_145 = m[A31]*Det2_45_45 - m[A34]*Det2_45_15 
						+ m[A35]*Det2_45_14;
  double Det3_345_234 = m[A32]*Det2_45_34 - m[A33]*Det2_45_24 
						+ m[A34]*Det2_45_23;
  double Det3_345_235 = m[A32]*Det2_45_35 - m[A33]*Det2_45_25 
						+ m[A35]*Det2_45_23;
  double Det3_345_245 = m[A32]*Det2_45_45 - m[A34]*Det2_45_25 
						+ m[A35]*Det2_45_24;
  double Det3_345_345 = m[A33]*Det2_45_45 - m[A34]*Det2_45_35 
						+ m[A35]*Det2_45_34;

  // Find all NECESSARY 4x4 dets:  (75 of them)

  double Det4_1234_0123 = m[A10]*Det3_234_123 - m[A11]*Det3_234_023 
			+ m[A12]*Det3_234_013 - m[A13]*Det3_234_012;
  double Det4_1234_0124 = m[A10]*Det3_234_124 - m[A11]*Det3_234_024 
			+ m[A12]*Det3_234_014 - m[A14]*Det3_234_012;
  double Det4_1234_0125 = m[A10]*Det3_234_125 - m[A11]*Det3_234_025 
			+ m[A12]*Det3_234_015 - m[A15]*Det3_234_012;	//
  double Det4_1234_0134 = m[A10]*Det3_234_134 - m[A11]*Det3_234_034 
			+ m[A13]*Det3_234_014 - m[A14]*Det3_234_013;
  double Det4_1234_0135 = m[A10]*Det3_234_135 - m[A11]*Det3_234_035
			+ m[A13]*Det3_234_015 - m[A15]*Det3_234_013;	//
  double Det4_1234_0145 = m[A10]*Det3_234_145 - m[A11]*Det3_234_045
			+ m[A14]*Det3_234_015 - m[A15]*Det3_234_014;	//
  double Det4_1234_0234 = m[A10]*Det3_234_234 - m[A12]*Det3_234_034 
			+ m[A13]*Det3_234_024 - m[A14]*Det3_234_023;
  double Det4_1234_0235 = m[A10]*Det3_234_235 - m[A12]*Det3_234_035 
			+ m[A13]*Det3_234_025 - m[A15]*Det3_234_023;	//
  double Det4_1234_0245 = m[A10]*Det3_234_245 - m[A12]*Det3_234_045 
			+ m[A14]*Det3_234_025 - m[A15]*Det3_234_024;	//
  double Det4_1234_0345 = m[A10]*Det3_234_345 - m[A13]*Det3_234_045 
			+ m[A14]*Det3_234_035 - m[A15]*Det3_234_034;	//
  double Det4_1234_1234 = m[A11]*Det3_234_234 - m[A12]*Det3_234_134 
			+ m[A13]*Det3_234_124 - m[A14]*Det3_234_123;
  double Det4_1234_1235 = m[A11]*Det3_234_235 - m[A12]*Det3_234_135 
			+ m[A13]*Det3_234_125 - m[A15]*Det3_234_123;	//
  double Det4_1234_1245 = m[A11]*Det3_234_245 - m[A12]*Det3_234_145 
			+ m[A14]*Det3_234_125 - m[A15]*Det3_234_124;	//
  double Det4_1234_1345 = m[A11]*Det3_234_345 - m[A13]*Det3_234_145 
			+ m[A14]*Det3_234_135 - m[A15]*Det3_234_134;	//
  double Det4_1234_2345 = m[A12]*Det3_234_345 - m[A13]*Det3_234_245 
			+ m[A14]*Det3_234_235 - m[A15]*Det3_234_234;	//
  double Det4_1235_0123 = m[A10]*Det3_235_123 - m[A11]*Det3_235_023 
			+ m[A12]*Det3_235_013 - m[A13]*Det3_235_012;
  double Det4_1235_0124 = m[A10]*Det3_235_124 - m[A11]*Det3_235_024 
			+ m[A12]*Det3_235_014 - m[A14]*Det3_235_012;
  double Det4_1235_0125 = m[A10]*Det3_235_125 - m[A11]*Det3_235_025 
			+ m[A12]*Det3_235_015 - m[A15]*Det3_235_012;
  double Det4_1235_0134 = m[A10]*Det3_235_134 - m[A11]*Det3_235_034 
			+ m[A13]*Det3_235_014 - m[A14]*Det3_235_013;
  double Det4_1235_0135 = m[A10]*Det3_235_135 - m[A11]*Det3_235_035 
			+ m[A13]*Det3_235_015 - m[A15]*Det3_235_013;
  double Det4_1235_0145 = m[A10]*Det3_235_145 - m[A11]*Det3_235_045 
			+ m[A14]*Det3_235_015 - m[A15]*Det3_235_014;	//
  double Det4_1235_0234 = m[A10]*Det3_235_234 - m[A12]*Det3_235_034 
			+ m[A13]*Det3_235_024 - m[A14]*Det3_235_023;
  double Det4_1235_0235 = m[A10]*Det3_235_235 - m[A12]*Det3_235_035 
			+ m[A13]*Det3_235_025 - m[A15]*Det3_235_023;
  double Det4_1235_0245 = m[A10]*Det3_235_245 - m[A12]*Det3_235_045 
			+ m[A14]*Det3_235_025 - m[A15]*Det3_235_024;	//
  double Det4_1235_0345 = m[A10]*Det3_235_345 - m[A13]*Det3_235_045 
			+ m[A14]*Det3_235_035 - m[A15]*Det3_235_034;	//
  double Det4_1235_1234 = m[A11]*Det3_235_234 - m[A12]*Det3_235_134 
			+ m[A13]*Det3_235_124 - m[A14]*Det3_235_123;
  double Det4_1235_1235 = m[A11]*Det3_235_235 - m[A12]*Det3_235_135 
			+ m[A13]*Det3_235_125 - m[A15]*Det3_235_123;
  double Det4_1235_1245 = m[A11]*Det3_235_245 - m[A12]*Det3_235_145 
			+ m[A14]*Det3_235_125 - m[A15]*Det3_235_124;	//
  double Det4_1235_1345 = m[A11]*Det3_235_345 - m[A13]*Det3_235_145 
			+ m[A14]*Det3_235_135 - m[A15]*Det3_235_134;	//
  double Det4_1235_2345 = m[A12]*Det3_235_345 - m[A13]*Det3_235_245 
			+ m[A14]*Det3_235_235 - m[A15]*Det3_235_234;	//
  double Det4_1245_0123 = m[A10]*Det3_245_123 - m[A11]*Det3_245_023 
			+ m[A12]*Det3_245_013 - m[A13]*Det3_245_012;
  double Det4_1245_0124 = m[A10]*Det3_245_124 - m[A11]*Det3_245_024 
			+ m[A12]*Det3_245_014 - m[A14]*Det3_245_012;
  double Det4_1245_0125 = m[A10]*Det3_245_125 - m[A11]*Det3_245_025 
			+ m[A12]*Det3_245_015 - m[A15]*Det3_245_012;
  double Det4_1245_0134 = m[A10]*Det3_245_134 - m[A11]*Det3_245_034 
			+ m[A13]*Det3_245_014 - m[A14]*Det3_245_013;
  double Det4_1245_0135 = m[A10]*Det3_245_135 - m[A11]*Det3_245_035 
			+ m[A13]*Det3_245_015 - m[A15]*Det3_245_013;
  double Det4_1245_0145 = m[A10]*Det3_245_145 - m[A11]*Det3_245_045 
			+ m[A14]*Det3_245_015 - m[A15]*Det3_245_014;
  double Det4_1245_0234 = m[A10]*Det3_245_234 - m[A12]*Det3_245_034 
			+ m[A13]*Det3_245_024 - m[A14]*Det3_245_023;
  double Det4_1245_0235 = m[A10]*Det3_245_235 - m[A12]*Det3_245_035 
			+ m[A13]*Det3_245_025 - m[A15]*Det3_245_023;
  double Det4_1245_0245 = m[A10]*Det3_245_245 - m[A12]*Det3_245_045 
			+ m[A14]*Det3_245_025 - m[A15]*Det3_245_024;
  double Det4_1245_0345 = m[A10]*Det3_245_345 - m[A13]*Det3_245_045 
			+ m[A14]*Det3_245_035 - m[A15]*Det3_245_034;	//
  double Det4_1245_1234 = m[A11]*Det3_245_234 - m[A12]*Det3_245_134 
			+ m[A13]*Det3_245_124 - m[A14]*Det3_245_123;
  double Det4_1245_1235 = m[A11]*Det3_245_235 - m[A12]*Det3_245_135 
			+ m[A13]*Det3_245_125 - m[A15]*Det3_245_123;
  double Det4_1245_1245 = m[A11]*Det3_245_245 - m[A12]*Det3_245_145 
			+ m[A14]*Det3_245_125 - m[A15]*Det3_245_124;
  double Det4_1245_1345 = m[A11]*Det3_245_345 - m[A13]*Det3_245_145 
			+ m[A14]*Det3_245_135 - m[A15]*Det3_245_134;	//
  double Det4_1245_2345 = m[A12]*Det3_245_345 - m[A13]*Det3_245_245 
			+ m[A14]*Det3_245_235 - m[A15]*Det3_245_234;	//
  double Det4_1345_0123 = m[A10]*Det3_345_123 - m[A11]*Det3_345_023 
			+ m[A12]*Det3_345_013 - m[A13]*Det3_345_012;
  double Det4_1345_0124 = m[A10]*Det3_345_124 - m[A11]*Det3_345_024 
			+ m[A12]*Det3_345_014 - m[A14]*Det3_345_012;
  double Det4_1345_0125 = m[A10]*Det3_345_125 - m[A11]*Det3_345_025 
			+ m[A12]*Det3_345_015 - m[A15]*Det3_345_012;
  double Det4_1345_0134 = m[A10]*Det3_345_134 - m[A11]*Det3_345_034 
			+ m[A13]*Det3_345_014 - m[A14]*Det3_345_013;
  double Det4_1345_0135 = m[A10]*Det3_345_135 - m[A11]*Det3_345_035 
			+ m[A13]*Det3_345_015 - m[A15]*Det3_345_013;
  double Det4_1345_0145 = m[A10]*Det3_345_145 - m[A11]*Det3_345_045 
			+ m[A14]*Det3_345_015 - m[A15]*Det3_345_014;
  double Det4_1345_0234 = m[A10]*Det3_345_234 - m[A12]*Det3_345_034 
			+ m[A13]*Det3_345_024 - m[A14]*Det3_345_023;
  double Det4_1345_0235 = m[A10]*Det3_345_235 - m[A12]*Det3_345_035 
			+ m[A13]*Det3_345_025 - m[A15]*Det3_345_023;
  double Det4_1345_0245 = m[A10]*Det3_345_245 - m[A12]*Det3_345_045 
			+ m[A14]*Det3_345_025 - m[A15]*Det3_345_024;
  double Det4_1345_0345 = m[A10]*Det3_345_345 - m[A13]*Det3_345_045 
			+ m[A14]*Det3_345_035 - m[A15]*Det3_345_034;
  double Det4_1345_1234 = m[A11]*Det3_345_234 - m[A12]*Det3_345_134 
			+ m[A13]*Det3_345_124 - m[A14]*Det3_345_123;
  double Det4_1345_1235 = m[A11]*Det3_345_235 - m[A12]*Det3_345_135 
			+ m[A13]*Det3_345_125 - m[A15]*Det3_345_123;
  double Det4_1345_1245 = m[A11]*Det3_345_245 - m[A12]*Det3_345_145 
			+ m[A14]*Det3_345_125 - m[A15]*Det3_345_124;
  double Det4_1345_1345 = m[A11]*Det3_345_345 - m[A13]*Det3_345_145 
			+ m[A14]*Det3_345_135 - m[A15]*Det3_345_134;
  double Det4_1345_2345 = m[A12]*Det3_345_345 - m[A13]*Det3_345_245 
			+ m[A14]*Det3_345_235 - m[A15]*Det3_345_234;	//
  double Det4_2345_0123 = m[A20]*Det3_345_123 - m[A21]*Det3_345_023 
			+ m[A22]*Det3_345_013 - m[A23]*Det3_345_012;
  double Det4_2345_0124 = m[A20]*Det3_345_124 - m[A21]*Det3_345_024 
			+ m[A22]*Det3_345_014 - m[A24]*Det3_345_012;
  double Det4_2345_0125 = m[A20]*Det3_345_125 - m[A21]*Det3_345_025 
			+ m[A22]*Det3_345_015 - m[A25]*Det3_345_012;
  double Det4_2345_0134 = m[A20]*Det3_345_134 - m[A21]*Det3_345_034 
			+ m[A23]*Det3_345_014 - m[A24]*Det3_345_013;
  double Det4_2345_0135 = m[A20]*Det3_345_135 - m[A21]*Det3_345_035 
			+ m[A23]*Det3_345_015 - m[A25]*Det3_345_013;
  double Det4_2345_0145 = m[A20]*Det3_345_145 - m[A21]*Det3_345_045 
			+ m[A24]*Det3_345_015 - m[A25]*Det3_345_014;
  double Det4_2345_0234 = m[A20]*Det3_345_234 - m[A22]*Det3_345_034 
			+ m[A23]*Det3_345_024 - m[A24]*Det3_345_023;
  double Det4_2345_0235 = m[A20]*Det3_345_235 - m[A22]*Det3_345_035 
			+ m[A23]*Det3_345_025 - m[A25]*Det3_345_023;
  double Det4_2345_0245 = m[A20]*Det3_345_245 - m[A22]*Det3_345_045 
			+ m[A24]*Det3_345_025 - m[A25]*Det3_345_024;
  double Det4_2345_0345 = m[A20]*Det3_345_345 - m[A23]*Det3_345_045 
			+ m[A24]*Det3_345_035 - m[A25]*Det3_345_034;
  double Det4_2345_1234 = m[A21]*Det3_345_234 - m[A22]*Det3_345_134 
			+ m[A23]*Det3_345_124 - m[A24]*Det3_345_123;
  double Det4_2345_1235 = m[A21]*Det3_345_235 - m[A22]*Det3_345_135 
			+ m[A23]*Det3_345_125 - m[A25]*Det3_345_123;
  double Det4_2345_1245 = m[A21]*Det3_345_245 - m[A22]*Det3_345_145 
			+ m[A24]*Det3_345_125 - m[A25]*Det3_345_124;
  double Det4_2345_1345 = m[A21]*Det3_345_345 - m[A23]*Det3_345_145 
			+ m[A24]*Det3_345_135 - m[A25]*Det3_345_134;
  double Det4_2345_2345 = m[A22]*Det3_345_345 - m[A23]*Det3_345_245 
			+ m[A24]*Det3_345_235 - m[A25]*Det3_345_234;

  // Find all NECESSARY 5x5 dets:  (36 of them)

  double Det5_01234_01234 = m[A00]*Det4_1234_1234 - m[A01]*Det4_1234_0234 
    + m[A02]*Det4_1234_0134 - m[A03]*Det4_1234_0124 + m[A04]*Det4_1234_0123;
  double Det5_01234_01235 = m[A00]*Det4_1234_1235 - m[A01]*Det4_1234_0235
    + m[A02]*Det4_1234_0135 - m[A03]*Det4_1234_0125 + m[A05]*Det4_1234_0123; //
  double Det5_01234_01245 = m[A00]*Det4_1234_1245 - m[A01]*Det4_1234_0245
    + m[A02]*Det4_1234_0145 - m[A04]*Det4_1234_0125 + m[A05]*Det4_1234_0124; //
  double Det5_01234_01345 = m[A00]*Det4_1234_1345 - m[A01]*Det4_1234_0345
    + m[A03]*Det4_1234_0145 - m[A04]*Det4_1234_0135 + m[A05]*Det4_1234_0134; //
  double Det5_01234_02345 = m[A00]*Det4_1234_2345 - m[A02]*Det4_1234_0345
    + m[A03]*Det4_1234_0245 - m[A04]*Det4_1234_0235 + m[A05]*Det4_1234_0234; //
  double Det5_01234_12345 = m[A01]*Det4_1234_2345 - m[A02]*Det4_1234_1345
    + m[A03]*Det4_1234_1245 - m[A04]*Det4_1234_1235 + m[A05]*Det4_1234_1234; //
  double Det5_01235_01234 = m[A00]*Det4_1235_1234 - m[A01]*Det4_1235_0234 
    + m[A02]*Det4_1235_0134 - m[A03]*Det4_1235_0124 + m[A04]*Det4_1235_0123;
  double Det5_01235_01235 = m[A00]*Det4_1235_1235 - m[A01]*Det4_1235_0235 
    + m[A02]*Det4_1235_0135 - m[A03]*Det4_1235_0125 + m[A05]*Det4_1235_0123;
  double Det5_01235_01245 = m[A00]*Det4_1235_1245 - m[A01]*Det4_1235_0245 
    + m[A02]*Det4_1235_0145 - m[A04]*Det4_1235_0125 + m[A05]*Det4_1235_0124; //
  double Det5_01235_01345 = m[A00]*Det4_1235_1345 - m[A01]*Det4_1235_0345 
    + m[A03]*Det4_1235_0145 - m[A04]*Det4_1235_0135 + m[A05]*Det4_1235_0134; //
  double Det5_01235_02345 = m[A00]*Det4_1235_2345 - m[A02]*Det4_1235_0345 
    + m[A03]*Det4_1235_0245 - m[A04]*Det4_1235_0235 + m[A05]*Det4_1235_0234; //
  double Det5_01235_12345 = m[A01]*Det4_1235_2345 - m[A02]*Det4_1235_1345 
    + m[A03]*Det4_1235_1245 - m[A04]*Det4_1235_1235 + m[A05]*Det4_1235_1234; //
  double Det5_01245_01234 = m[A00]*Det4_1245_1234 - m[A01]*Det4_1245_0234 
    + m[A02]*Det4_1245_0134 - m[A03]*Det4_1245_0124 + m[A04]*Det4_1245_0123;
  double Det5_01245_01235 = m[A00]*Det4_1245_1235 - m[A01]*Det4_1245_0235 
    + m[A02]*Det4_1245_0135 - m[A03]*Det4_1245_0125 + m[A05]*Det4_1245_0123;
  double Det5_01245_01245 = m[A00]*Det4_1245_1245 - m[A01]*Det4_1245_0245 
    + m[A02]*Det4_1245_0145 - m[A04]*Det4_1245_0125 + m[A05]*Det4_1245_0124;
  double Det5_01245_01345 = m[A00]*Det4_1245_1345 - m[A01]*Det4_1245_0345 
    + m[A03]*Det4_1245_0145 - m[A04]*Det4_1245_0135 + m[A05]*Det4_1245_0134; //
  double Det5_01245_02345 = m[A00]*Det4_1245_2345 - m[A02]*Det4_1245_0345 
    + m[A03]*Det4_1245_0245 - m[A04]*Det4_1245_0235 + m[A05]*Det4_1245_0234; //
  double Det5_01245_12345 = m[A01]*Det4_1245_2345 - m[A02]*Det4_1245_1345 
    + m[A03]*Det4_1245_1245 - m[A04]*Det4_1245_1235 + m[A05]*Det4_1245_1234; //
  double Det5_01345_01234 = m[A00]*Det4_1345_1234 - m[A01]*Det4_1345_0234 
    + m[A02]*Det4_1345_0134 - m[A03]*Det4_1345_0124 + m[A04]*Det4_1345_0123;
  double Det5_01345_01235 = m[A00]*Det4_1345_1235 - m[A01]*Det4_1345_0235 
    + m[A02]*Det4_1345_0135 - m[A03]*Det4_1345_0125 + m[A05]*Det4_1345_0123;
  double Det5_01345_01245 = m[A00]*Det4_1345_1245 - m[A01]*Det4_1345_0245 
    + m[A02]*Det4_1345_0145 - m[A04]*Det4_1345_0125 + m[A05]*Det4_1345_0124;
  double Det5_01345_01345 = m[A00]*Det4_1345_1345 - m[A01]*Det4_1345_0345 
    + m[A03]*Det4_1345_0145 - m[A04]*Det4_1345_0135 + m[A05]*Det4_1345_0134;
  double Det5_01345_02345 = m[A00]*Det4_1345_2345 - m[A02]*Det4_1345_0345 
    + m[A03]*Det4_1345_0245 - m[A04]*Det4_1345_0235 + m[A05]*Det4_1345_0234; //
  double Det5_01345_12345 = m[A01]*Det4_1345_2345 - m[A02]*Det4_1345_1345 
    + m[A03]*Det4_1345_1245 - m[A04]*Det4_1345_1235 + m[A05]*Det4_1345_1234; //
  double Det5_02345_01234 = m[A00]*Det4_2345_1234 - m[A01]*Det4_2345_0234 
    + m[A02]*Det4_2345_0134 - m[A03]*Det4_2345_0124 + m[A04]*Det4_2345_0123;
  double Det5_02345_01235 = m[A00]*Det4_2345_1235 - m[A01]*Det4_2345_0235 
    + m[A02]*Det4_2345_0135 - m[A03]*Det4_2345_0125 + m[A05]*Det4_2345_0123;
  double Det5_02345_01245 = m[A00]*Det4_2345_1245 - m[A01]*Det4_2345_0245 
    + m[A02]*Det4_2345_0145 - m[A04]*Det4_2345_0125 + m[A05]*Det4_2345_0124;
  double Det5_02345_01345 = m[A00]*Det4_2345_1345 - m[A01]*Det4_2345_0345 
    + m[A03]*Det4_2345_0145 - m[A04]*Det4_2345_0135 + m[A05]*Det4_2345_0134;
  double Det5_02345_02345 = m[A00]*Det4_2345_2345 - m[A02]*Det4_2345_0345 
    + m[A03]*Det4_2345_0245 - m[A04]*Det4_2345_0235 + m[A05]*Det4_2345_0234;
  double Det5_02345_12345 = m[A01]*Det4_2345_2345 - m[A02]*Det4_2345_1345 
    + m[A03]*Det4_2345_1245 - m[A04]*Det4_2345_1235 + m[A05]*Det4_2345_1234; //
  double Det5_12345_01234 = m[A10]*Det4_2345_1234 - m[A11]*Det4_2345_0234 
    + m[A12]*Det4_2345_0134 - m[A13]*Det4_2345_0124 + m[A14]*Det4_2345_0123;
  double Det5_12345_01235 = m[A10]*Det4_2345_1235 - m[A11]*Det4_2345_0235 
    + m[A12]*Det4_2345_0135 - m[A13]*Det4_2345_0125 + m[A15]*Det4_2345_0123;
  double Det5_12345_01245 = m[A10]*Det4_2345_1245 - m[A11]*Det4_2345_0245 
    + m[A12]*Det4_2345_0145 - m[A14]*Det4_2345_0125 + m[A15]*Det4_2345_0124;
  double Det5_12345_01345 = m[A10]*Det4_2345_1345 - m[A11]*Det4_2345_0345 
    + m[A13]*Det4_2345_0145 - m[A14]*Det4_2345_0135 + m[A15]*Det4_2345_0134;
  double Det5_12345_02345 = m[A10]*Det4_2345_2345 - m[A12]*Det4_2345_0345 
    + m[A13]*Det4_2345_0245 - m[A14]*Det4_2345_0235 + m[A15]*Det4_2345_0234;
  double Det5_12345_12345 = m[A11]*Det4_2345_2345 - m[A12]*Det4_2345_1345 
    + m[A13]*Det4_2345_1245 - m[A14]*Det4_2345_1235 + m[A15]*Det4_2345_1234;

  // Find the determinant 

  double det =    m[A00]*Det5_12345_12345 
	     	- m[A01]*Det5_12345_02345 
	     	+ m[A02]*Det5_12345_01345 
		- m[A03]*Det5_12345_01245 
		+ m[A04]*Det5_12345_01235 
		- m[A05]*Det5_12345_01234;

  if ( det == 0 ) {  
#ifdef SINGULAR_DIAGNOSTICS
    std::cerr << "Kramer's rule inversion of a singular 6x6 matrix: "
	<< *this << "\n";
#endif
    ifail = 1;
    return;
  } 

  double oneOverDet = 1.0/det;
  double mn1OverDet = - oneOverDet;

  m[A00] =  Det5_12345_12345*oneOverDet;
  m[A01] =  Det5_02345_12345*mn1OverDet;
  m[A02] =  Det5_01345_12345*oneOverDet;
  m[A03] =  Det5_01245_12345*mn1OverDet;
  m[A04] =  Det5_01235_12345*oneOverDet;
  m[A05] =  Det5_01234_12345*mn1OverDet;

  m[A10] =  Det5_12345_02345*mn1OverDet;
  m[A11] =  Det5_02345_02345*oneOverDet;
  m[A12] =  Det5_01345_02345*mn1OverDet;
  m[A13] =  Det5_01245_02345*oneOverDet;
  m[A14] =  Det5_01235_02345*mn1OverDet;
  m[A15] =  Det5_01234_02345*oneOverDet;

  m[A20] =  Det5_12345_01345*oneOverDet;
  m[A21] =  Det5_02345_01345*mn1OverDet;
  m[A22] =  Det5_01345_01345*oneOverDet;
  m[A23] =  Det5_01245_01345*mn1OverDet;
  m[A24] =  Det5_01235_01345*oneOverDet;
  m[A25] =  Det5_01234_01345*mn1OverDet;

  m[A30] =  Det5_12345_01245*mn1OverDet;
  m[A31] =  Det5_02345_01245*oneOverDet;
  m[A32] =  Det5_01345_01245*mn1OverDet;
  m[A33] =  Det5_01245_01245*oneOverDet;
  m[A34] =  Det5_01235_01245*mn1OverDet;
  m[A35] =  Det5_01234_01245*oneOverDet;

  m[A40] =  Det5_12345_01235*oneOverDet;
  m[A41] =  Det5_02345_01235*mn1OverDet;
  m[A42] =  Det5_01345_01235*oneOverDet;
  m[A43] =  Det5_01245_01235*mn1OverDet;
  m[A44] =  Det5_01235_01235*oneOverDet;
  m[A45] =  Det5_01234_01235*mn1OverDet;

  m[A50] =  Det5_12345_01234*mn1OverDet;
  m[A51] =  Det5_02345_01234*oneOverDet;
  m[A52] =  Det5_01345_01234*mn1OverDet;
  m[A53] =  Det5_01245_01234*oneOverDet;
  m[A54] =  Det5_01235_01234*mn1OverDet;
  m[A55] =  Det5_01234_01234*oneOverDet;

  return;
}


