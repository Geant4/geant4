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
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file
// ------------------------------------------------------------

#include "globals.hh"
#include <iostream>
#include <cmath>

#include "G4ErrorSymMatrix.hh"
#include "G4ErrorMatrix.hh"

// Simple operation for all elements

#define SIMPLE_UOP(OPER)                                                       \
  G4ErrorMatrixIter a = m.begin();                                             \
  G4ErrorMatrixIter e = m.begin() + num_size();                                \
  for(; a < e; a++)                                                            \
    (*a) OPER t;

#define SIMPLE_BOP(OPER)                                                       \
  G4ErrorMatrixIter a      = m.begin();                                        \
  G4ErrorMatrixConstIter b = mat2.m.begin();                                   \
  G4ErrorMatrixConstIter e = m.begin() + num_size();                           \
  for(; a < e; a++, b++)                                                       \
    (*a) OPER(*b);

#define SIMPLE_TOP(OPER)                                                       \
  G4ErrorMatrixConstIter a = mat1.m.begin();                                   \
  G4ErrorMatrixConstIter b = mat2.m.begin();                                   \
  G4ErrorMatrixIter t      = mret.m.begin();                                   \
  G4ErrorMatrixConstIter e = mat1.m.begin() + mat1.num_size();                 \
  for(; a < e; a++, b++, t++)                                                  \
    (*t) = (*a) OPER(*b);

#define CHK_DIM_2(r1, r2, c1, c2, fun)                                         \
  if(r1 != r2 || c1 != c2)                                                     \
  {                                                                            \
    G4ErrorMatrix::error("Range error in Matrix function " #fun "(1).");       \
  }

#define CHK_DIM_1(c1, r2, fun)                                                 \
  if(c1 != r2)                                                                 \
  {                                                                            \
    G4ErrorMatrix::error("Range error in Matrix function " #fun "(2).");       \
  }

// Constructors. (Default constructors are inlined and in .icc file)

G4ErrorSymMatrix::G4ErrorSymMatrix(G4int p)
  : m(p * (p + 1) / 2)
  , nrow(p)
{
  size = nrow * (nrow + 1) / 2;
  m.assign(size, 0);
}

G4ErrorSymMatrix::G4ErrorSymMatrix(G4int p, G4int init)
  : m(p * (p + 1) / 2)
  , nrow(p)
{
  size = nrow * (nrow + 1) / 2;

  m.assign(size, 0);
  switch(init)
  {
    case 0:
      break;

    case 1:
    {
      G4ErrorMatrixIter a = m.begin();
      for(G4int i = 1; i <= nrow; i++)
      {
        *a = 1.0;
        a += (i + 1);
      }
      break;
    }
    default:
      G4ErrorMatrix::error("G4ErrorSymMatrix: initialization must be 0 or 1.");
  }
}

//
// Destructor
//

G4ErrorSymMatrix::~G4ErrorSymMatrix() {}

G4ErrorSymMatrix::G4ErrorSymMatrix(const G4ErrorSymMatrix& mat1)
  : m(mat1.size)
  , nrow(mat1.nrow)
  , size(mat1.size)
{
  m = mat1.m;
}

//
//
// Sub matrix
//
//

G4ErrorSymMatrix G4ErrorSymMatrix::sub(G4int min_row, G4int max_row) const
{
  G4ErrorSymMatrix mret(max_row - min_row + 1);
  if(max_row > num_row())
  {
    G4ErrorMatrix::error("G4ErrorSymMatrix::sub: Index out of range");
  }
  G4ErrorMatrixIter a       = mret.m.begin();
  G4ErrorMatrixConstIter b1 = m.begin() + (min_row + 2) * (min_row - 1) / 2;
  for(G4int irow = 1; irow <= mret.num_row(); irow++)
  {
    G4ErrorMatrixConstIter b = b1;
    for(G4int icol = 1; icol <= irow; icol++)
    {
      *(a++) = *(b++);
    }
    b1 += irow + min_row - 1;
  }
  return mret;
}

G4ErrorSymMatrix G4ErrorSymMatrix::sub(G4int min_row, G4int max_row)
{
  G4ErrorSymMatrix mret(max_row - min_row + 1);
  if(max_row > num_row())
  {
    G4ErrorMatrix::error("G4ErrorSymMatrix::sub: Index out of range");
  }
  G4ErrorMatrixIter a  = mret.m.begin();
  G4ErrorMatrixIter b1 = m.begin() + (min_row + 2) * (min_row - 1) / 2;
  for(G4int irow = 1; irow <= mret.num_row(); irow++)
  {
    G4ErrorMatrixIter b = b1;
    for(G4int icol = 1; icol <= irow; icol++)
    {
      *(a++) = *(b++);
    }
    b1 += irow + min_row - 1;
  }
  return mret;
}

void G4ErrorSymMatrix::sub(G4int row, const G4ErrorSymMatrix& mat1)
{
  if(row < 1 || row + mat1.num_row() - 1 > num_row())
  {
    G4ErrorMatrix::error("G4ErrorSymMatrix::sub: Index out of range");
  }
  G4ErrorMatrixConstIter a = mat1.m.begin();
  G4ErrorMatrixIter b1     = m.begin() + (row + 2) * (row - 1) / 2;
  for(G4int irow = 1; irow <= mat1.num_row(); irow++)
  {
    G4ErrorMatrixIter b = b1;
    for(G4int icol = 1; icol <= irow; icol++)
    {
      *(b++) = *(a++);
    }
    b1 += irow + row - 1;
  }
}

//
// Direct sum of two matricies
//

G4ErrorSymMatrix dsum(const G4ErrorSymMatrix& mat1,
                      const G4ErrorSymMatrix& mat2)
{
  G4ErrorSymMatrix mret(mat1.num_row() + mat2.num_row(), 0);
  mret.sub(1, mat1);
  mret.sub(mat1.num_row() + 1, mat2);
  return mret;
}

/* -----------------------------------------------------------------------
   This section contains support routines for matrix.h. This section contains
   The two argument functions +,-. They call the copy constructor and +=,-=.
   ----------------------------------------------------------------------- */

G4ErrorSymMatrix G4ErrorSymMatrix::operator-() const
{
  G4ErrorSymMatrix mat2(nrow);
  G4ErrorMatrixConstIter a = m.begin();
  G4ErrorMatrixIter b      = mat2.m.begin();
  G4ErrorMatrixConstIter e = m.begin() + num_size();
  for(; a < e; a++, b++)
  {
    (*b) = -(*a);
  }
  return mat2;
}

G4ErrorMatrix operator+(const G4ErrorMatrix& mat1, const G4ErrorSymMatrix& mat2)
{
  G4ErrorMatrix mret(mat1);
  CHK_DIM_2(mat1.num_row(), mat2.num_row(), mat1.num_col(), mat2.num_col(), +);
  mret += mat2;
  return mret;
}

G4ErrorMatrix operator+(const G4ErrorSymMatrix& mat1, const G4ErrorMatrix& mat2)
{
  G4ErrorMatrix mret(mat2);
  CHK_DIM_2(mat1.num_row(), mat2.num_row(), mat1.num_col(), mat2.num_col(), +);
  mret += mat1;
  return mret;
}

G4ErrorSymMatrix operator+(const G4ErrorSymMatrix& mat1,
                           const G4ErrorSymMatrix& mat2)
{
  G4ErrorSymMatrix mret(mat1.nrow);
  CHK_DIM_1(mat1.nrow, mat2.nrow, +);
  SIMPLE_TOP(+)
  return mret;
}

//
// operator -
//

G4ErrorMatrix operator-(const G4ErrorMatrix& mat1, const G4ErrorSymMatrix& mat2)
{
  G4ErrorMatrix mret(mat1);
  CHK_DIM_2(mat1.num_row(), mat2.num_row(), mat1.num_col(), mat2.num_col(), -);
  mret -= mat2;
  return mret;
}

G4ErrorMatrix operator-(const G4ErrorSymMatrix& mat1, const G4ErrorMatrix& mat2)
{
  G4ErrorMatrix mret(mat1);
  CHK_DIM_2(mat1.num_row(), mat2.num_row(), mat1.num_col(), mat2.num_col(), -);
  mret -= mat2;
  return mret;
}

G4ErrorSymMatrix operator-(const G4ErrorSymMatrix& mat1,
                           const G4ErrorSymMatrix& mat2)
{
  G4ErrorSymMatrix mret(mat1.num_row());
  CHK_DIM_1(mat1.num_row(), mat2.num_row(), -);
  SIMPLE_TOP(-)
  return mret;
}

/* -----------------------------------------------------------------------
   This section contains support routines for matrix.h. This file contains
   The two argument functions *,/. They call copy constructor and then /=,*=.
   ----------------------------------------------------------------------- */

G4ErrorSymMatrix operator/(const G4ErrorSymMatrix& mat1, G4double t)
{
  G4ErrorSymMatrix mret(mat1);
  mret /= t;
  return mret;
}

G4ErrorSymMatrix operator*(const G4ErrorSymMatrix& mat1, G4double t)
{
  G4ErrorSymMatrix mret(mat1);
  mret *= t;
  return mret;
}

G4ErrorSymMatrix operator*(G4double t, const G4ErrorSymMatrix& mat1)
{
  G4ErrorSymMatrix mret(mat1);
  mret *= t;
  return mret;
}

G4ErrorMatrix operator*(const G4ErrorMatrix& mat1, const G4ErrorSymMatrix& mat2)
{
  G4ErrorMatrix mret(mat1.num_row(), mat2.num_col());
  CHK_DIM_1(mat1.num_col(), mat2.num_row(), *);
  G4ErrorMatrixConstIter mit1, mit2, sp, snp;  // mit2=0
  G4double temp;
  G4ErrorMatrixIter mir = mret.m.begin();
  for(mit1 = mat1.m.begin();
      mit1 < mat1.m.begin() + mat1.num_row() * mat1.num_col(); mit1 = mit2)
  {
    snp = mat2.m.begin();
    for(int step = 1; step <= mat2.num_row(); ++step)
    {
      mit2 = mit1;
      sp   = snp;
      snp += step;
      temp = 0;
      while(sp < snp)  // Loop checking, 06.08.2015, G.Cosmo
      {
        temp += *(sp++) * (*(mit2++));
      }
      if(step < mat2.num_row())
      {  // only if we aren't on the last row
        sp += step - 1;
        for(int stept = step + 1; stept <= mat2.num_row(); stept++)
        {
          temp += *sp * (*(mit2++));
          if(stept < mat2.num_row())
            sp += stept;
        }
      }  // if(step
      *(mir++) = temp;
    }  // for(step
  }    // for(mit1
  return mret;
}

G4ErrorMatrix operator*(const G4ErrorSymMatrix& mat1, const G4ErrorMatrix& mat2)
{
  G4ErrorMatrix mret(mat1.num_row(), mat2.num_col());
  CHK_DIM_1(mat1.num_col(), mat2.num_row(), *);
  G4int step, stept;
  G4ErrorMatrixConstIter mit1, mit2, sp, snp;
  G4double temp;
  G4ErrorMatrixIter mir = mret.m.begin();
  for(step = 1, snp = mat1.m.begin(); step <= mat1.num_row(); snp += step++)
  {
    for(mit1 = mat2.m.begin(); mit1 < mat2.m.begin() + mat2.num_col(); mit1++)
    {
      mit2 = mit1;
      sp   = snp;
      temp = 0;
      while(sp < snp + step)  // Loop checking, 06.08.2015, G.Cosmo
      {
        temp += *mit2 * (*(sp++));
        mit2 += mat2.num_col();
      }
      sp += step - 1;
      for(stept = step + 1; stept <= mat1.num_row(); stept++)
      {
        temp += *mit2 * (*sp);
        mit2 += mat2.num_col();
        sp += stept;
      }
      *(mir++) = temp;
    }
  }
  return mret;
}

G4ErrorMatrix operator*(const G4ErrorSymMatrix& mat1,
                        const G4ErrorSymMatrix& mat2)
{
  G4ErrorMatrix mret(mat1.num_row(), mat1.num_row());
  CHK_DIM_1(mat1.num_col(), mat2.num_row(), *);
  G4int step1, stept1, step2, stept2;
  G4ErrorMatrixConstIter snp1, sp1, snp2, sp2;
  G4double temp;
  G4ErrorMatrixIter mr = mret.m.begin();
  for(step1 = 1, snp1 = mat1.m.begin(); step1 <= mat1.num_row();
      snp1 += step1++)
  {
    for(step2 = 1, snp2 = mat2.m.begin(); step2 <= mat2.num_row();)
    {
      sp1 = snp1;
      sp2 = snp2;
      snp2 += step2;
      temp = 0;
      if(step1 < step2)
      {
        while(sp1 < snp1 + step1)  // Loop checking, 06.08.2015, G.Cosmo
        {
          temp += (*(sp1++)) * (*(sp2++));
        }
        sp1 += step1 - 1;
        for(stept1 = step1 + 1; stept1 != step2 + 1; sp1 += stept1++)
        {
          temp += (*sp1) * (*(sp2++));
        }
        sp2 += step2 - 1;
        for(stept2 = ++step2; stept2 <= mat2.num_row();
            sp1 += stept1++, sp2 += stept2++)
        {
          temp += (*sp1) * (*sp2);
        }
      }
      else
      {
        while(sp2 < snp2)  // Loop checking, 06.08.2015, G.Cosmo
        {
          temp += (*(sp1++)) * (*(sp2++));
        }
        sp2 += step2 - 1;
        for(stept2 = ++step2; stept2 != step1 + 1; sp2 += stept2++)
        {
          temp += (*(sp1++)) * (*sp2);
        }
        sp1 += step1 - 1;
        for(stept1 = step1 + 1; stept1 <= mat1.num_row();
            sp1 += stept1++, sp2 += stept2++)
        {
          temp += (*sp1) * (*sp2);
        }
      }
      *(mr++) = temp;
    }
  }
  return mret;
}

/* -----------------------------------------------------------------------
   This section contains the assignment and inplace operators =,+=,-=,*=,/=.
   ----------------------------------------------------------------------- */

G4ErrorMatrix& G4ErrorMatrix::operator+=(const G4ErrorSymMatrix& mat2)
{
  CHK_DIM_2(num_row(), mat2.num_row(), num_col(), mat2.num_col(), +=);
  G4int n                    = num_col();
  G4ErrorMatrixConstIter sjk = mat2.m.begin();
  G4ErrorMatrixIter m1j      = m.begin();
  G4ErrorMatrixIter mj       = m.begin();
  // j >= k
  for(G4int j = 1; j <= num_row(); j++)
  {
    G4ErrorMatrixIter mjk = mj;
    G4ErrorMatrixIter mkj = m1j;
    for(G4int k = 1; k <= j; k++)
    {
      *(mjk++) += *sjk;
      if(j != k)
        *mkj += *sjk;
      sjk++;
      mkj += n;
    }
    mj += n;
    m1j++;
  }
  return (*this);
}

G4ErrorSymMatrix& G4ErrorSymMatrix::operator+=(const G4ErrorSymMatrix& mat2)
{
  CHK_DIM_2(num_row(), mat2.num_row(), num_col(), mat2.num_col(), +=);
  SIMPLE_BOP(+=)
  return (*this);
}

G4ErrorMatrix& G4ErrorMatrix::operator-=(const G4ErrorSymMatrix& mat2)
{
  CHK_DIM_2(num_row(), mat2.num_row(), num_col(), mat2.num_col(), -=);
  G4int n                    = num_col();
  G4ErrorMatrixConstIter sjk = mat2.m.begin();
  G4ErrorMatrixIter m1j      = m.begin();
  G4ErrorMatrixIter mj       = m.begin();
  // j >= k
  for(G4int j = 1; j <= num_row(); j++)
  {
    G4ErrorMatrixIter mjk = mj;
    G4ErrorMatrixIter mkj = m1j;
    for(G4int k = 1; k <= j; k++)
    {
      *(mjk++) -= *sjk;
      if(j != k)
        *mkj -= *sjk;
      sjk++;
      mkj += n;
    }
    mj += n;
    m1j++;
  }
  return (*this);
}

G4ErrorSymMatrix& G4ErrorSymMatrix::operator-=(const G4ErrorSymMatrix& mat2)
{
  CHK_DIM_2(num_row(), mat2.num_row(), num_col(), mat2.num_col(), -=);
  SIMPLE_BOP(-=)
  return (*this);
}

G4ErrorSymMatrix& G4ErrorSymMatrix::operator/=(G4double t)
{
  SIMPLE_UOP(/=)
  return (*this);
}

G4ErrorSymMatrix& G4ErrorSymMatrix::operator*=(G4double t)
{
  SIMPLE_UOP(*=)
  return (*this);
}

G4ErrorMatrix& G4ErrorMatrix::operator=(const G4ErrorSymMatrix& mat1)
{
  if(mat1.nrow * mat1.nrow != size)
  {
    size = mat1.nrow * mat1.nrow;
    m.resize(size);
  }
  nrow                       = mat1.nrow;
  ncol                       = mat1.nrow;
  G4int n                    = ncol;
  G4ErrorMatrixConstIter sjk = mat1.m.begin();
  G4ErrorMatrixIter m1j      = m.begin();
  G4ErrorMatrixIter mj       = m.begin();
  // j >= k
  for(G4int j = 1; j <= num_row(); j++)
  {
    G4ErrorMatrixIter mjk = mj;
    G4ErrorMatrixIter mkj = m1j;
    for(G4int k = 1; k <= j; k++)
    {
      *(mjk++) = *sjk;
      if(j != k)
        *mkj = *sjk;
      sjk++;
      mkj += n;
    }
    mj += n;
    m1j++;
  }
  return (*this);
}

G4ErrorSymMatrix& G4ErrorSymMatrix::operator=(const G4ErrorSymMatrix& mat1)
{
  if(&mat1 == this)
  {
    return *this;
  }
  if(mat1.nrow != nrow)
  {
    nrow = mat1.nrow;
    size = mat1.size;
    m.resize(size);
  }
  m = mat1.m;
  return (*this);
}

// Print the Matrix.

std::ostream& operator<<(std::ostream& os, const G4ErrorSymMatrix& q)
{
  os << G4endl;

  // Fixed format needs 3 extra characters for field,
  // while scientific needs 7

  std::size_t width;
  if(os.flags() & std::ios::fixed)
  {
    width = os.precision() + 3;
  }
  else
  {
    width = os.precision() + 7;
  }
  for(G4int irow = 1; irow <= q.num_row(); ++irow)
  {
    for(G4int icol = 1; icol <= q.num_col(); ++icol)
    {
      os.width(width);
      os << q(irow, icol) << " ";
    }
    os << G4endl;
  }
  return os;
}

G4ErrorSymMatrix G4ErrorSymMatrix::apply(G4double (*f)(G4double, G4int,
                                                       G4int)) const
{
  G4ErrorSymMatrix mret(num_row());
  G4ErrorMatrixConstIter a = m.cbegin();
  G4ErrorMatrixIter b      = mret.m.begin();
  for(G4int ir = 1; ir <= num_row(); ++ir)
  {
    for(G4int ic = 1; ic <= ir; ++ic)
    {
      *(b++) = (*f)(*(a++), ir, ic);
    }
  }
  return mret;
}

void G4ErrorSymMatrix::assign(const G4ErrorMatrix& mat1)
{
  if(mat1.nrow != nrow)
  {
    nrow = mat1.nrow;
    size = nrow * (nrow + 1) / 2;
    m.resize(size);
  }
  G4ErrorMatrixConstIter a = mat1.m.begin();
  G4ErrorMatrixIter b      = m.begin();
  for(G4int r = 1; r <= nrow; r++)
  {
    G4ErrorMatrixConstIter d = a;
    for(G4int c = 1; c <= r; c++)
    {
      *(b++) = *(d++);
    }
    a += nrow;
  }
}

G4ErrorSymMatrix G4ErrorSymMatrix::similarity(const G4ErrorMatrix& mat1) const
{
  G4ErrorSymMatrix mret(mat1.num_row());
  G4ErrorMatrix temp = mat1 * (*this);

  // If mat1*(*this) has correct dimensions, then so will the mat1.T
  // multiplication. So there is no need to check dimensions again.

  G4int n                  = mat1.num_col();
  G4ErrorMatrixIter mr     = mret.m.begin();
  G4ErrorMatrixIter tempr1 = temp.m.begin();
  for(G4int r = 1; r <= mret.num_row(); r++)
  {
    G4ErrorMatrixConstIter m1c1 = mat1.m.begin();
    for(G4int c = 1; c <= r; c++)
    {
      G4double tmp                = 0.0;
      G4ErrorMatrixIter tempri    = tempr1;
      G4ErrorMatrixConstIter m1ci = m1c1;
      for(G4int i = 1; i <= mat1.num_col(); i++)
      {
        tmp += (*(tempri++)) * (*(m1ci++));
      }
      *(mr++) = tmp;
      m1c1 += n;
    }
    tempr1 += n;
  }
  return mret;
}

G4ErrorSymMatrix G4ErrorSymMatrix::similarity(
  const G4ErrorSymMatrix& mat1) const
{
  G4ErrorSymMatrix mret(mat1.num_row());
  G4ErrorMatrix temp       = mat1 * (*this);
  G4int n                  = mat1.num_col();
  G4ErrorMatrixIter mr     = mret.m.begin();
  G4ErrorMatrixIter tempr1 = temp.m.begin();
  for(G4int r = 1; r <= mret.num_row(); r++)
  {
    G4ErrorMatrixConstIter m1c1 = mat1.m.begin();
    G4int c;
    for(c = 1; c <= r; c++)
    {
      G4double tmp                = 0.0;
      G4ErrorMatrixIter tempri    = tempr1;
      G4ErrorMatrixConstIter m1ci = m1c1;
      G4int i;
      for(i = 1; i < c; i++)
      {
        tmp += (*(tempri++)) * (*(m1ci++));
      }
      for(i = c; i <= mat1.num_col(); i++)
      {
        tmp += (*(tempri++)) * (*(m1ci));
        m1ci += i;
      }
      *(mr++) = tmp;
      m1c1 += c;
    }
    tempr1 += n;
  }
  return mret;
}

G4ErrorSymMatrix G4ErrorSymMatrix::similarityT(const G4ErrorMatrix& mat1) const
{
  G4ErrorSymMatrix mret(mat1.num_col());
  G4ErrorMatrix temp       = (*this) * mat1;
  G4int n                  = mat1.num_col();
  G4ErrorMatrixIter mrc    = mret.m.begin();
  G4ErrorMatrixIter temp1r = temp.m.begin();
  for(G4int r = 1; r <= mret.num_row(); r++)
  {
    G4ErrorMatrixConstIter m11c = mat1.m.begin();
    for(G4int c = 1; c <= r; c++)
    {
      G4double tmp                = 0.0;
      G4ErrorMatrixIter tempir    = temp1r;
      G4ErrorMatrixConstIter m1ic = m11c;
      for(G4int i = 1; i <= mat1.num_row(); i++)
      {
        tmp += (*(tempir)) * (*(m1ic));
        tempir += n;
        m1ic += n;
      }
      *(mrc++) = tmp;
      m11c++;
    }
    temp1r++;
  }
  return mret;
}

void G4ErrorSymMatrix::invert(G4int& ifail)
{
  ifail = 0;

  switch(nrow)
  {
    case 3:
    {
      G4double det, temp;
      G4double t1, t2, t3;
      G4double c11, c12, c13, c22, c23, c33;
      c11 = (*(m.begin() + 2)) * (*(m.begin() + 5)) -
            (*(m.begin() + 4)) * (*(m.begin() + 4));
      c12 = (*(m.begin() + 4)) * (*(m.begin() + 3)) -
            (*(m.begin() + 1)) * (*(m.begin() + 5));
      c13 = (*(m.begin() + 1)) * (*(m.begin() + 4)) -
            (*(m.begin() + 2)) * (*(m.begin() + 3));
      c22 = (*(m.begin() + 5)) * (*m.begin()) -
            (*(m.begin() + 3)) * (*(m.begin() + 3));
      c23 = (*(m.begin() + 3)) * (*(m.begin() + 1)) -
            (*(m.begin() + 4)) * (*m.begin());
      c33 = (*m.begin()) * (*(m.begin() + 2)) -
            (*(m.begin() + 1)) * (*(m.begin() + 1));
      t1 = std::fabs(*m.begin());
      t2 = std::fabs(*(m.begin() + 1));
      t3 = std::fabs(*(m.begin() + 3));
      if(t1 >= t2)
      {
        if(t3 >= t1)
        {
          temp = *(m.begin() + 3);
          det  = c23 * c12 - c22 * c13;
        }
        else
        {
          temp = *m.begin();
          det  = c22 * c33 - c23 * c23;
        }
      }
      else if(t3 >= t2)
      {
        temp = *(m.begin() + 3);
        det  = c23 * c12 - c22 * c13;
      }
      else
      {
        temp = *(m.begin() + 1);
        det  = c13 * c23 - c12 * c33;
      }
      if(det == 0)
      {
        ifail = 1;
        return;
      }
      {
        G4double ss          = temp / det;
        G4ErrorMatrixIter mq = m.begin();
        *(mq++)              = ss * c11;
        *(mq++)              = ss * c12;
        *(mq++)              = ss * c22;
        *(mq++)              = ss * c13;
        *(mq++)              = ss * c23;
        *(mq)                = ss * c33;
      }
    }
    break;
    case 2:
    {
      G4double det, temp, ss;
      det = (*m.begin()) * (*(m.begin() + 2)) -
            (*(m.begin() + 1)) * (*(m.begin() + 1));
      if(det == 0)
      {
        ifail = 1;
        return;
      }
      ss = 1.0 / det;
      *(m.begin() + 1) *= -ss;
      temp             = ss * (*(m.begin() + 2));
      *(m.begin() + 2) = ss * (*m.begin());
      *m.begin()       = temp;
      break;
    }
    case 1:
    {
      if((*m.begin()) == 0)
      {
        ifail = 1;
        return;
      }
      *m.begin() = 1.0 / (*m.begin());
      break;
    }
    case 5:
    {
      invert5(ifail);
      return;
    }
    case 6:
    {
      invert6(ifail);
      return;
    }
    case 4:
    {
      invert4(ifail);
      return;
    }
    default:
    {
      invertBunchKaufman(ifail);
      return;
    }
  }
  return;  // inversion successful
}

G4double G4ErrorSymMatrix::determinant() const
{
  static const G4int max_array = 20;

  // ir must point to an array which is ***1 longer than*** nrow

  static std::vector<G4int> ir_vec(max_array + 1);
  if(ir_vec.size() <= static_cast<unsigned int>(nrow))
  {
    ir_vec.resize(nrow + 1);
  }
  G4int* ir = &ir_vec[0];

  G4double det;
  G4ErrorMatrix mt(*this);
  G4int i = mt.dfact_matrix(det, ir);
  if(i == 0)
  {
    return det;
  }
  return 0.0;
}

G4double G4ErrorSymMatrix::trace() const
{
  G4double t = 0.0;
  for(G4int i = 0; i < nrow; i++)
  {
    t += *(m.begin() + (i + 3) * i / 2);
  }
  return t;
}

void G4ErrorSymMatrix::invertBunchKaufman(G4int& ifail)
{
  // Bunch-Kaufman diagonal pivoting method
  // It is decribed in J.R. Bunch, L. Kaufman (1977).
  // "Some Stable Methods for Calculating Inertia and Solving Symmetric
  // Linear Systems", Math. Comp. 31, p. 162-179. or in Gene H. Golub,
  // Charles F. van Loan, "Matrix Computations" (the second edition
  // has a bug.) and implemented in "lapack"
  // Mario Stanke, 09/97

  G4int i, j, k, ss;
  G4int pivrow;

  // Establish the two working-space arrays needed:  x and piv are
  // used as pointers to arrays of doubles and ints respectively, each
  // of length nrow.  We do not want to reallocate each time through
  // unless the size needs to grow.  We do not want to leak memory, even
  // by having a new without a delete that is only done once.

  static const G4int max_array                     = 25;
  static G4ThreadLocal std::vector<G4double>* xvec = 0;
  if(!xvec)
    xvec = new std::vector<G4double>(max_array);
  static G4ThreadLocal std::vector<G4int>* pivv = 0;
  if(!pivv)
    pivv = new std::vector<G4int>(max_array);
  typedef std::vector<G4int>::iterator pivIter;
  if(xvec->size() < static_cast<unsigned int>(nrow))
    xvec->resize(nrow);
  if(pivv->size() < static_cast<unsigned int>(nrow))
    pivv->resize(nrow);
  // Note - resize should do  nothing if the size is already larger than nrow,
  //        but on VC++ there are indications that it does so we check.
  // Note - the data elements in a vector are guaranteed to be contiguous,
  //        so x[i] and piv[i] are optimally fast.
  G4ErrorMatrixIter x = xvec->begin();
  // x[i] is used as helper storage, needs to have at least size nrow.
  pivIter piv = pivv->begin();
  // piv[i] is used to store details of exchanges

  G4double temp1, temp2;
  G4ErrorMatrixIter ip, mjj, iq;
  G4double lambda, sigma;
  const G4double alpha   = .6404;  // = (1+sqrt(17))/8
  const G4double epsilon = 32 * DBL_EPSILON;
  // whenever a sum of two doubles is below or equal to epsilon
  // it is set to zero.
  // this constant could be set to zero but then the algorithm
  // doesn't neccessarily detect that a matrix is singular

  for(i = 0; i < nrow; i++)
  {
    piv[i] = i + 1;
  }

  ifail = 0;

  // compute the factorization P*A*P^T = L * D * L^T
  // L is unit lower triangular, D is direct sum of 1x1 and 2x2 matrices
  // L and D^-1 are stored in A = *this, P is stored in piv[]

  for(j = 1; j < nrow; j += ss)  // main loop over columns
  {
    mjj    = m.begin() + j * (j - 1) / 2 + j - 1;
    lambda = 0;  // compute lambda = max of A(j+1:n,j)
    pivrow = j + 1;
    ip     = m.begin() + (j + 1) * j / 2 + j - 1;
    for(i = j + 1; i <= nrow; ip += i++)
    {
      if(std::fabs(*ip) > lambda)
      {
        lambda = std::fabs(*ip);
        pivrow = i;
      }
    }
    if(lambda == 0)
    {
      if(*mjj == 0)
      {
        ifail = 1;
        return;
      }
      ss   = 1;
      *mjj = 1. / *mjj;
    }
    else
    {
      if(std::fabs(*mjj) >= lambda * alpha)
      {
        ss     = 1;
        pivrow = j;
      }
      else
      {
        sigma = 0;  // compute sigma = max A(pivrow, j:pivrow-1)
        ip    = m.begin() + pivrow * (pivrow - 1) / 2 + j - 1;
        for(k = j; k < pivrow; k++)
        {
          if(std::fabs(*ip) > sigma)
            sigma = std::fabs(*ip);
          ip++;
        }
        if(sigma * std::fabs(*mjj) >= alpha * lambda * lambda)
        {
          ss     = 1;
          pivrow = j;
        }
        else if(std::fabs(*(m.begin() + pivrow * (pivrow - 1) / 2 + pivrow -
                            1)) >= alpha * sigma)
        {
          ss = 1;
        }
        else
        {
          ss = 2;
        }
      }
      if(pivrow == j)  // no permutation neccessary
      {
        piv[j - 1] = pivrow;
        if(*mjj == 0)
        {
          ifail = 1;
          return;
        }
        temp2 = *mjj = 1. / *mjj;  // invert D(j,j)

        // update A(j+1:n, j+1,n)
        for(i = j + 1; i <= nrow; i++)
        {
          temp1 = *(m.begin() + i * (i - 1) / 2 + j - 1) * temp2;
          ip    = m.begin() + i * (i - 1) / 2 + j;
          for(k = j + 1; k <= i; k++)
          {
            *ip -= temp1 * *(m.begin() + k * (k - 1) / 2 + j - 1);
            if(std::fabs(*ip) <= epsilon)
            {
              *ip = 0;
            }
            ip++;
          }
        }
        // update L
        ip = m.begin() + (j + 1) * j / 2 + j - 1;
        for(i = j + 1; i <= nrow; ip += i++)
        {
          *ip *= temp2;
        }
      }
      else if(ss == 1)  // 1x1 pivot
      {
        piv[j - 1] = pivrow;

        // interchange rows and columns j and pivrow in
        // submatrix (j:n,j:n)
        ip = m.begin() + pivrow * (pivrow - 1) / 2 + j;
        for(i = j + 1; i < pivrow; i++, ip++)
        {
          temp1 = *(m.begin() + i * (i - 1) / 2 + j - 1);
          *(m.begin() + i * (i - 1) / 2 + j - 1) = *ip;
          *ip                                    = temp1;
        }
        temp1 = *mjj;
        *mjj  = *(m.begin() + pivrow * (pivrow - 1) / 2 + pivrow - 1);
        *(m.begin() + pivrow * (pivrow - 1) / 2 + pivrow - 1) = temp1;
        ip = m.begin() + (pivrow + 1) * pivrow / 2 + j - 1;
        iq = ip + pivrow - j;
        for(i = pivrow + 1; i <= nrow; ip += i, iq += i++)
        {
          temp1 = *iq;
          *iq   = *ip;
          *ip   = temp1;
        }

        if(*mjj == 0)
        {
          ifail = 1;
          return;
        }
        temp2 = *mjj = 1. / *mjj;  // invert D(j,j)

        // update A(j+1:n, j+1:n)
        for(i = j + 1; i <= nrow; i++)
        {
          temp1 = *(m.begin() + i * (i - 1) / 2 + j - 1) * temp2;
          ip    = m.begin() + i * (i - 1) / 2 + j;
          for(k = j + 1; k <= i; k++)
          {
            *ip -= temp1 * *(m.begin() + k * (k - 1) / 2 + j - 1);
            if(std::fabs(*ip) <= epsilon)
            {
              *ip = 0;
            }
            ip++;
          }
        }
        // update L
        ip = m.begin() + (j + 1) * j / 2 + j - 1;
        for(i = j + 1; i <= nrow; ip += i++)
        {
          *ip *= temp2;
        }
      }
      else  // ss=2, ie use a 2x2 pivot
      {
        piv[j - 1] = -pivrow;
        piv[j]     = 0;  // that means this is the second row of a 2x2 pivot

        if(j + 1 != pivrow)
        {
          // interchange rows and columns j+1 and pivrow in
          // submatrix (j:n,j:n)
          ip = m.begin() + pivrow * (pivrow - 1) / 2 + j + 1;
          for(i = j + 2; i < pivrow; i++, ip++)
          {
            temp1 = *(m.begin() + i * (i - 1) / 2 + j);
            *(m.begin() + i * (i - 1) / 2 + j) = *ip;
            *ip                                = temp1;
          }
          temp1 = *(mjj + j + 1);
          *(mjj + j + 1) =
            *(m.begin() + pivrow * (pivrow - 1) / 2 + pivrow - 1);
          *(m.begin() + pivrow * (pivrow - 1) / 2 + pivrow - 1) = temp1;
          temp1                                                 = *(mjj + j);
          *(mjj + j) = *(m.begin() + pivrow * (pivrow - 1) / 2 + j - 1);
          *(m.begin() + pivrow * (pivrow - 1) / 2 + j - 1) = temp1;
          ip = m.begin() + (pivrow + 1) * pivrow / 2 + j;
          iq = ip + pivrow - (j + 1);
          for(i = pivrow + 1; i <= nrow; ip += i, iq += i++)
          {
            temp1 = *iq;
            *iq   = *ip;
            *ip   = temp1;
          }
        }
        // invert D(j:j+1,j:j+1)
        temp2 = *mjj * *(mjj + j + 1) - *(mjj + j) * *(mjj + j);
        if(temp2 == 0)
        {
          G4Exception("G4ErrorSymMatrix::bunch_invert()",
                      "GEANT4e-Notification", JustWarning,
                      "Error in pivot choice!");
        }
        temp2 = 1. / temp2;

        // this quotient is guaranteed to exist by the choice
        // of the pivot

        temp1          = *mjj;
        *mjj           = *(mjj + j + 1) * temp2;
        *(mjj + j + 1) = temp1 * temp2;
        *(mjj + j)     = -*(mjj + j) * temp2;

        if(j < nrow - 1)  // otherwise do nothing
        {
          // update A(j+2:n, j+2:n)
          for(i = j + 2; i <= nrow; i++)
          {
            ip    = m.begin() + i * (i - 1) / 2 + j - 1;
            temp1 = *ip * *mjj + *(ip + 1) * *(mjj + j);
            if(std::fabs(temp1) <= epsilon)
            {
              temp1 = 0;
            }
            temp2 = *ip * *(mjj + j) + *(ip + 1) * *(mjj + j + 1);
            if(std::fabs(temp2) <= epsilon)
            {
              temp2 = 0;
            }
            for(k = j + 2; k <= i; k++)
            {
              ip = m.begin() + i * (i - 1) / 2 + k - 1;
              iq = m.begin() + k * (k - 1) / 2 + j - 1;
              *ip -= temp1 * *iq + temp2 * *(iq + 1);
              if(std::fabs(*ip) <= epsilon)
              {
                *ip = 0;
              }
            }
          }
          // update L
          for(i = j + 2; i <= nrow; i++)
          {
            ip    = m.begin() + i * (i - 1) / 2 + j - 1;
            temp1 = *ip * *mjj + *(ip + 1) * *(mjj + j);
            if(std::fabs(temp1) <= epsilon)
            {
              temp1 = 0;
            }
            *(ip + 1) = *ip * *(mjj + j) + *(ip + 1) * *(mjj + j + 1);
            if(std::fabs(*(ip + 1)) <= epsilon)
            {
              *(ip + 1) = 0;
            }
            *ip = temp1;
          }
        }
      }
    }
  }  // end of main loop over columns

  if(j == nrow)  // the the last pivot is 1x1
  {
    mjj = m.begin() + j * (j - 1) / 2 + j - 1;
    if(*mjj == 0)
    {
      ifail = 1;
      return;
    }
    else
    {
      *mjj = 1. / *mjj;
    }
  }  // end of last pivot code

  // computing the inverse from the factorization

  for(j = nrow; j >= 1; j -= ss)  // loop over columns
  {
    mjj = m.begin() + j * (j - 1) / 2 + j - 1;
    if(piv[j - 1] > 0)  // 1x1 pivot, compute column j of inverse
    {
      ss = 1;
      if(j < nrow)
      {
        ip = m.begin() + (j + 1) * j / 2 + j - 1;
        for(i = 0; i < nrow - j; ip += 1 + j + i++)
        {
          x[i] = *ip;
        }
        for(i = j + 1; i <= nrow; i++)
        {
          temp2 = 0;
          ip    = m.begin() + i * (i - 1) / 2 + j;
          for(k = 0; k <= i - j - 1; k++)
          {
            temp2 += *ip++ * x[k];
          }
          for(ip += i - 1; k < nrow - j; ip += 1 + j + k++)
          {
            temp2 += *ip * x[k];
          }
          *(m.begin() + i * (i - 1) / 2 + j - 1) = -temp2;
        }
        temp2 = 0;
        ip    = m.begin() + (j + 1) * j / 2 + j - 1;
        for(k = 0; k < nrow - j; ip += 1 + j + k++)
        {
          temp2 += x[k] * *ip;
        }
        *mjj -= temp2;
      }
    }
    else  // 2x2 pivot, compute columns j and j-1 of the inverse
    {
      if(piv[j - 1] != 0)
      {
        std::ostringstream message;
        message << "Error in pivot: " << piv[j - 1];
        G4Exception("G4ErrorSymMatrix::invertBunchKaufman()",
                    "GEANT4e-Notification", JustWarning, message);
      }
      ss = 2;
      if(j < nrow)
      {
        ip = m.begin() + (j + 1) * j / 2 + j - 1;
        for(i = 0; i < nrow - j; ip += 1 + j + i++)
        {
          x[i] = *ip;
        }
        for(i = j + 1; i <= nrow; i++)
        {
          temp2 = 0;
          ip    = m.begin() + i * (i - 1) / 2 + j;
          for(k = 0; k <= i - j - 1; k++)
          {
            temp2 += *ip++ * x[k];
          }
          for(ip += i - 1; k < nrow - j; ip += 1 + j + k++)
          {
            temp2 += *ip * x[k];
          }
          *(m.begin() + i * (i - 1) / 2 + j - 1) = -temp2;
        }
        temp2 = 0;
        ip    = m.begin() + (j + 1) * j / 2 + j - 1;
        for(k = 0; k < nrow - j; ip += 1 + j + k++)
        {
          temp2 += x[k] * *ip;
        }
        *mjj -= temp2;
        temp2 = 0;
        ip    = m.begin() + (j + 1) * j / 2 + j - 2;
        for(i = j + 1; i <= nrow; ip += i++)
        {
          temp2 += *ip * *(ip + 1);
        }
        *(mjj - 1) -= temp2;
        ip = m.begin() + (j + 1) * j / 2 + j - 2;
        for(i = 0; i < nrow - j; ip += 1 + j + i++)
        {
          x[i] = *ip;
        }
        for(i = j + 1; i <= nrow; i++)
        {
          temp2 = 0;
          ip    = m.begin() + i * (i - 1) / 2 + j;
          for(k = 0; k <= i - j - 1; k++)
          {
            temp2 += *ip++ * x[k];
          }
          for(ip += i - 1; k < nrow - j; ip += 1 + j + k++)
          {
            temp2 += *ip * x[k];
          }
          *(m.begin() + i * (i - 1) / 2 + j - 2) = -temp2;
        }
        temp2 = 0;
        ip    = m.begin() + (j + 1) * j / 2 + j - 2;
        for(k = 0; k < nrow - j; ip += 1 + j + k++)
        {
          temp2 += x[k] * *ip;
        }
        *(mjj - j) -= temp2;
      }
    }

    // interchange rows and columns j and piv[j-1]
    // or rows and columns j and -piv[j-2]

    pivrow = (piv[j - 1] == 0) ? -piv[j - 2] : piv[j - 1];
    ip     = m.begin() + pivrow * (pivrow - 1) / 2 + j;
    for(i = j + 1; i < pivrow; i++, ip++)
    {
      temp1 = *(m.begin() + i * (i - 1) / 2 + j - 1);
      *(m.begin() + i * (i - 1) / 2 + j - 1) = *ip;
      *ip                                    = temp1;
    }
    temp1 = *mjj;
    *mjj  = *(m.begin() + pivrow * (pivrow - 1) / 2 + pivrow - 1);
    *(m.begin() + pivrow * (pivrow - 1) / 2 + pivrow - 1) = temp1;
    if(ss == 2)
    {
      temp1      = *(mjj - 1);
      *(mjj - 1) = *(m.begin() + pivrow * (pivrow - 1) / 2 + j - 2);
      *(m.begin() + pivrow * (pivrow - 1) / 2 + j - 2) = temp1;
    }

    ip = m.begin() + (pivrow + 1) * pivrow / 2 + j - 1;  // &A(i,j)
    iq = ip + pivrow - j;
    for(i = pivrow + 1; i <= nrow; ip += i, iq += i++)
    {
      temp1 = *iq;
      *iq   = *ip;
      *ip   = temp1;
    }
  }  // end of loop over columns (in computing inverse from factorization)

  return;  // inversion successful
}

G4ThreadLocal G4double G4ErrorSymMatrix::posDefFraction5x5 = 1.0;
G4ThreadLocal G4double G4ErrorSymMatrix::posDefFraction6x6 = 1.0;
G4ThreadLocal G4double G4ErrorSymMatrix::adjustment5x5     = 0.0;
G4ThreadLocal G4double G4ErrorSymMatrix::adjustment6x6     = 0.0;
const G4double G4ErrorSymMatrix::CHOLESKY_THRESHOLD_5x5    = .5;
const G4double G4ErrorSymMatrix::CHOLESKY_THRESHOLD_6x6    = .2;
const G4double G4ErrorSymMatrix::CHOLESKY_CREEP_5x5        = .005;
const G4double G4ErrorSymMatrix::CHOLESKY_CREEP_6x6        = .002;

// Aij are indices for a 6x6 symmetric matrix.
//     The indices for 5x5 or 4x4 symmetric matrices are the same,
//     ignoring all combinations with an index which is inapplicable.

#define A00 0
#define A01 1
#define A02 3
#define A03 6
#define A04 10
#define A05 15

#define A10 1
#define A11 2
#define A12 4
#define A13 7
#define A14 11
#define A15 16

#define A20 3
#define A21 4
#define A22 5
#define A23 8
#define A24 12
#define A25 17

#define A30 6
#define A31 7
#define A32 8
#define A33 9
#define A34 13
#define A35 18

#define A40 10
#define A41 11
#define A42 12
#define A43 13
#define A44 14
#define A45 19

#define A50 15
#define A51 16
#define A52 17
#define A53 18
#define A54 19
#define A55 20

void G4ErrorSymMatrix::invert5(G4int& ifail)
{
  if(posDefFraction5x5 >= CHOLESKY_THRESHOLD_5x5)
  {
    invertCholesky5(ifail);
    posDefFraction5x5 = .9 * posDefFraction5x5 + .1 * (1 - ifail);
    if(ifail != 0)  // Cholesky failed -- invert using Haywood
    {
      invertHaywood5(ifail);
    }
  }
  else
  {
    if(posDefFraction5x5 + adjustment5x5 >= CHOLESKY_THRESHOLD_5x5)
    {
      invertCholesky5(ifail);
      posDefFraction5x5 = .9 * posDefFraction5x5 + .1 * (1 - ifail);
      if(ifail != 0)  // Cholesky failed -- invert using Haywood
      {
        invertHaywood5(ifail);
        adjustment5x5 = 0;
      }
    }
    else
    {
      invertHaywood5(ifail);
      adjustment5x5 += CHOLESKY_CREEP_5x5;
    }
  }
  return;
}

void G4ErrorSymMatrix::invert6(G4int& ifail)
{
  if(posDefFraction6x6 >= CHOLESKY_THRESHOLD_6x6)
  {
    invertCholesky6(ifail);
    posDefFraction6x6 = .9 * posDefFraction6x6 + .1 * (1 - ifail);
    if(ifail != 0)  // Cholesky failed -- invert using Haywood
    {
      invertHaywood6(ifail);
    }
  }
  else
  {
    if(posDefFraction6x6 + adjustment6x6 >= CHOLESKY_THRESHOLD_6x6)
    {
      invertCholesky6(ifail);
      posDefFraction6x6 = .9 * posDefFraction6x6 + .1 * (1 - ifail);
      if(ifail != 0)  // Cholesky failed -- invert using Haywood
      {
        invertHaywood6(ifail);
        adjustment6x6 = 0;
      }
    }
    else
    {
      invertHaywood6(ifail);
      adjustment6x6 += CHOLESKY_CREEP_6x6;
    }
  }
  return;
}

void G4ErrorSymMatrix::invertHaywood5(G4int& ifail)
{
  ifail = 0;

  // Find all NECESSARY 2x2 dets:  (25 of them)

  G4double Det2_23_01 = m[A20] * m[A31] - m[A21] * m[A30];
  G4double Det2_23_02 = m[A20] * m[A32] - m[A22] * m[A30];
  G4double Det2_23_03 = m[A20] * m[A33] - m[A23] * m[A30];
  G4double Det2_23_12 = m[A21] * m[A32] - m[A22] * m[A31];
  G4double Det2_23_13 = m[A21] * m[A33] - m[A23] * m[A31];
  G4double Det2_23_23 = m[A22] * m[A33] - m[A23] * m[A32];
  G4double Det2_24_01 = m[A20] * m[A41] - m[A21] * m[A40];
  G4double Det2_24_02 = m[A20] * m[A42] - m[A22] * m[A40];
  G4double Det2_24_03 = m[A20] * m[A43] - m[A23] * m[A40];
  G4double Det2_24_04 = m[A20] * m[A44] - m[A24] * m[A40];
  G4double Det2_24_12 = m[A21] * m[A42] - m[A22] * m[A41];
  G4double Det2_24_13 = m[A21] * m[A43] - m[A23] * m[A41];
  G4double Det2_24_14 = m[A21] * m[A44] - m[A24] * m[A41];
  G4double Det2_24_23 = m[A22] * m[A43] - m[A23] * m[A42];
  G4double Det2_24_24 = m[A22] * m[A44] - m[A24] * m[A42];
  G4double Det2_34_01 = m[A30] * m[A41] - m[A31] * m[A40];
  G4double Det2_34_02 = m[A30] * m[A42] - m[A32] * m[A40];
  G4double Det2_34_03 = m[A30] * m[A43] - m[A33] * m[A40];
  G4double Det2_34_04 = m[A30] * m[A44] - m[A34] * m[A40];
  G4double Det2_34_12 = m[A31] * m[A42] - m[A32] * m[A41];
  G4double Det2_34_13 = m[A31] * m[A43] - m[A33] * m[A41];
  G4double Det2_34_14 = m[A31] * m[A44] - m[A34] * m[A41];
  G4double Det2_34_23 = m[A32] * m[A43] - m[A33] * m[A42];
  G4double Det2_34_24 = m[A32] * m[A44] - m[A34] * m[A42];
  G4double Det2_34_34 = m[A33] * m[A44] - m[A34] * m[A43];

  // Find all NECESSARY 3x3 dets:   (30 of them)

  G4double Det3_123_012 =
    m[A10] * Det2_23_12 - m[A11] * Det2_23_02 + m[A12] * Det2_23_01;
  G4double Det3_123_013 =
    m[A10] * Det2_23_13 - m[A11] * Det2_23_03 + m[A13] * Det2_23_01;
  G4double Det3_123_023 =
    m[A10] * Det2_23_23 - m[A12] * Det2_23_03 + m[A13] * Det2_23_02;
  G4double Det3_123_123 =
    m[A11] * Det2_23_23 - m[A12] * Det2_23_13 + m[A13] * Det2_23_12;
  G4double Det3_124_012 =
    m[A10] * Det2_24_12 - m[A11] * Det2_24_02 + m[A12] * Det2_24_01;
  G4double Det3_124_013 =
    m[A10] * Det2_24_13 - m[A11] * Det2_24_03 + m[A13] * Det2_24_01;
  G4double Det3_124_014 =
    m[A10] * Det2_24_14 - m[A11] * Det2_24_04 + m[A14] * Det2_24_01;
  G4double Det3_124_023 =
    m[A10] * Det2_24_23 - m[A12] * Det2_24_03 + m[A13] * Det2_24_02;
  G4double Det3_124_024 =
    m[A10] * Det2_24_24 - m[A12] * Det2_24_04 + m[A14] * Det2_24_02;
  G4double Det3_124_123 =
    m[A11] * Det2_24_23 - m[A12] * Det2_24_13 + m[A13] * Det2_24_12;
  G4double Det3_124_124 =
    m[A11] * Det2_24_24 - m[A12] * Det2_24_14 + m[A14] * Det2_24_12;
  G4double Det3_134_012 =
    m[A10] * Det2_34_12 - m[A11] * Det2_34_02 + m[A12] * Det2_34_01;
  G4double Det3_134_013 =
    m[A10] * Det2_34_13 - m[A11] * Det2_34_03 + m[A13] * Det2_34_01;
  G4double Det3_134_014 =
    m[A10] * Det2_34_14 - m[A11] * Det2_34_04 + m[A14] * Det2_34_01;
  G4double Det3_134_023 =
    m[A10] * Det2_34_23 - m[A12] * Det2_34_03 + m[A13] * Det2_34_02;
  G4double Det3_134_024 =
    m[A10] * Det2_34_24 - m[A12] * Det2_34_04 + m[A14] * Det2_34_02;
  G4double Det3_134_034 =
    m[A10] * Det2_34_34 - m[A13] * Det2_34_04 + m[A14] * Det2_34_03;
  G4double Det3_134_123 =
    m[A11] * Det2_34_23 - m[A12] * Det2_34_13 + m[A13] * Det2_34_12;
  G4double Det3_134_124 =
    m[A11] * Det2_34_24 - m[A12] * Det2_34_14 + m[A14] * Det2_34_12;
  G4double Det3_134_134 =
    m[A11] * Det2_34_34 - m[A13] * Det2_34_14 + m[A14] * Det2_34_13;
  G4double Det3_234_012 =
    m[A20] * Det2_34_12 - m[A21] * Det2_34_02 + m[A22] * Det2_34_01;
  G4double Det3_234_013 =
    m[A20] * Det2_34_13 - m[A21] * Det2_34_03 + m[A23] * Det2_34_01;
  G4double Det3_234_014 =
    m[A20] * Det2_34_14 - m[A21] * Det2_34_04 + m[A24] * Det2_34_01;
  G4double Det3_234_023 =
    m[A20] * Det2_34_23 - m[A22] * Det2_34_03 + m[A23] * Det2_34_02;
  G4double Det3_234_024 =
    m[A20] * Det2_34_24 - m[A22] * Det2_34_04 + m[A24] * Det2_34_02;
  G4double Det3_234_034 =
    m[A20] * Det2_34_34 - m[A23] * Det2_34_04 + m[A24] * Det2_34_03;
  G4double Det3_234_123 =
    m[A21] * Det2_34_23 - m[A22] * Det2_34_13 + m[A23] * Det2_34_12;
  G4double Det3_234_124 =
    m[A21] * Det2_34_24 - m[A22] * Det2_34_14 + m[A24] * Det2_34_12;
  G4double Det3_234_134 =
    m[A21] * Det2_34_34 - m[A23] * Det2_34_14 + m[A24] * Det2_34_13;
  G4double Det3_234_234 =
    m[A22] * Det2_34_34 - m[A23] * Det2_34_24 + m[A24] * Det2_34_23;

  // Find all NECESSARY 4x4 dets:   (15 of them)

  G4double Det4_0123_0123 = m[A00] * Det3_123_123 - m[A01] * Det3_123_023 +
                            m[A02] * Det3_123_013 - m[A03] * Det3_123_012;
  G4double Det4_0124_0123 = m[A00] * Det3_124_123 - m[A01] * Det3_124_023 +
                            m[A02] * Det3_124_013 - m[A03] * Det3_124_012;
  G4double Det4_0124_0124 = m[A00] * Det3_124_124 - m[A01] * Det3_124_024 +
                            m[A02] * Det3_124_014 - m[A04] * Det3_124_012;
  G4double Det4_0134_0123 = m[A00] * Det3_134_123 - m[A01] * Det3_134_023 +
                            m[A02] * Det3_134_013 - m[A03] * Det3_134_012;
  G4double Det4_0134_0124 = m[A00] * Det3_134_124 - m[A01] * Det3_134_024 +
                            m[A02] * Det3_134_014 - m[A04] * Det3_134_012;
  G4double Det4_0134_0134 = m[A00] * Det3_134_134 - m[A01] * Det3_134_034 +
                            m[A03] * Det3_134_014 - m[A04] * Det3_134_013;
  G4double Det4_0234_0123 = m[A00] * Det3_234_123 - m[A01] * Det3_234_023 +
                            m[A02] * Det3_234_013 - m[A03] * Det3_234_012;
  G4double Det4_0234_0124 = m[A00] * Det3_234_124 - m[A01] * Det3_234_024 +
                            m[A02] * Det3_234_014 - m[A04] * Det3_234_012;
  G4double Det4_0234_0134 = m[A00] * Det3_234_134 - m[A01] * Det3_234_034 +
                            m[A03] * Det3_234_014 - m[A04] * Det3_234_013;
  G4double Det4_0234_0234 = m[A00] * Det3_234_234 - m[A02] * Det3_234_034 +
                            m[A03] * Det3_234_024 - m[A04] * Det3_234_023;
  G4double Det4_1234_0123 = m[A10] * Det3_234_123 - m[A11] * Det3_234_023 +
                            m[A12] * Det3_234_013 - m[A13] * Det3_234_012;
  G4double Det4_1234_0124 = m[A10] * Det3_234_124 - m[A11] * Det3_234_024 +
                            m[A12] * Det3_234_014 - m[A14] * Det3_234_012;
  G4double Det4_1234_0134 = m[A10] * Det3_234_134 - m[A11] * Det3_234_034 +
                            m[A13] * Det3_234_014 - m[A14] * Det3_234_013;
  G4double Det4_1234_0234 = m[A10] * Det3_234_234 - m[A12] * Det3_234_034 +
                            m[A13] * Det3_234_024 - m[A14] * Det3_234_023;
  G4double Det4_1234_1234 = m[A11] * Det3_234_234 - m[A12] * Det3_234_134 +
                            m[A13] * Det3_234_124 - m[A14] * Det3_234_123;

  // Find the 5x5 det:

  G4double det = m[A00] * Det4_1234_1234 - m[A01] * Det4_1234_0234 +
                 m[A02] * Det4_1234_0134 - m[A03] * Det4_1234_0124 +
                 m[A04] * Det4_1234_0123;

  if(det == 0)
  {
    ifail = 1;
    return;
  }

  G4double oneOverDet = 1.0 / det;
  G4double mn1OverDet = -oneOverDet;

  m[A00] = Det4_1234_1234 * oneOverDet;
  m[A01] = Det4_1234_0234 * mn1OverDet;
  m[A02] = Det4_1234_0134 * oneOverDet;
  m[A03] = Det4_1234_0124 * mn1OverDet;
  m[A04] = Det4_1234_0123 * oneOverDet;

  m[A11] = Det4_0234_0234 * oneOverDet;
  m[A12] = Det4_0234_0134 * mn1OverDet;
  m[A13] = Det4_0234_0124 * oneOverDet;
  m[A14] = Det4_0234_0123 * mn1OverDet;

  m[A22] = Det4_0134_0134 * oneOverDet;
  m[A23] = Det4_0134_0124 * mn1OverDet;
  m[A24] = Det4_0134_0123 * oneOverDet;

  m[A33] = Det4_0124_0124 * oneOverDet;
  m[A34] = Det4_0124_0123 * mn1OverDet;

  m[A44] = Det4_0123_0123 * oneOverDet;

  return;
}

void G4ErrorSymMatrix::invertHaywood6(G4int& ifail)
{
  ifail = 0;

  // Find all NECESSARY 2x2 dets:  (39 of them)

  G4double Det2_34_01 = m[A30] * m[A41] - m[A31] * m[A40];
  G4double Det2_34_02 = m[A30] * m[A42] - m[A32] * m[A40];
  G4double Det2_34_03 = m[A30] * m[A43] - m[A33] * m[A40];
  G4double Det2_34_04 = m[A30] * m[A44] - m[A34] * m[A40];
  G4double Det2_34_12 = m[A31] * m[A42] - m[A32] * m[A41];
  G4double Det2_34_13 = m[A31] * m[A43] - m[A33] * m[A41];
  G4double Det2_34_14 = m[A31] * m[A44] - m[A34] * m[A41];
  G4double Det2_34_23 = m[A32] * m[A43] - m[A33] * m[A42];
  G4double Det2_34_24 = m[A32] * m[A44] - m[A34] * m[A42];
  G4double Det2_34_34 = m[A33] * m[A44] - m[A34] * m[A43];
  G4double Det2_35_01 = m[A30] * m[A51] - m[A31] * m[A50];
  G4double Det2_35_02 = m[A30] * m[A52] - m[A32] * m[A50];
  G4double Det2_35_03 = m[A30] * m[A53] - m[A33] * m[A50];
  G4double Det2_35_04 = m[A30] * m[A54] - m[A34] * m[A50];
  G4double Det2_35_05 = m[A30] * m[A55] - m[A35] * m[A50];
  G4double Det2_35_12 = m[A31] * m[A52] - m[A32] * m[A51];
  G4double Det2_35_13 = m[A31] * m[A53] - m[A33] * m[A51];
  G4double Det2_35_14 = m[A31] * m[A54] - m[A34] * m[A51];
  G4double Det2_35_15 = m[A31] * m[A55] - m[A35] * m[A51];
  G4double Det2_35_23 = m[A32] * m[A53] - m[A33] * m[A52];
  G4double Det2_35_24 = m[A32] * m[A54] - m[A34] * m[A52];
  G4double Det2_35_25 = m[A32] * m[A55] - m[A35] * m[A52];
  G4double Det2_35_34 = m[A33] * m[A54] - m[A34] * m[A53];
  G4double Det2_35_35 = m[A33] * m[A55] - m[A35] * m[A53];
  G4double Det2_45_01 = m[A40] * m[A51] - m[A41] * m[A50];
  G4double Det2_45_02 = m[A40] * m[A52] - m[A42] * m[A50];
  G4double Det2_45_03 = m[A40] * m[A53] - m[A43] * m[A50];
  G4double Det2_45_04 = m[A40] * m[A54] - m[A44] * m[A50];
  G4double Det2_45_05 = m[A40] * m[A55] - m[A45] * m[A50];
  G4double Det2_45_12 = m[A41] * m[A52] - m[A42] * m[A51];
  G4double Det2_45_13 = m[A41] * m[A53] - m[A43] * m[A51];
  G4double Det2_45_14 = m[A41] * m[A54] - m[A44] * m[A51];
  G4double Det2_45_15 = m[A41] * m[A55] - m[A45] * m[A51];
  G4double Det2_45_23 = m[A42] * m[A53] - m[A43] * m[A52];
  G4double Det2_45_24 = m[A42] * m[A54] - m[A44] * m[A52];
  G4double Det2_45_25 = m[A42] * m[A55] - m[A45] * m[A52];
  G4double Det2_45_34 = m[A43] * m[A54] - m[A44] * m[A53];
  G4double Det2_45_35 = m[A43] * m[A55] - m[A45] * m[A53];
  G4double Det2_45_45 = m[A44] * m[A55] - m[A45] * m[A54];

  // Find all NECESSARY 3x3 dets:  (65 of them)

  G4double Det3_234_012 =
    m[A20] * Det2_34_12 - m[A21] * Det2_34_02 + m[A22] * Det2_34_01;
  G4double Det3_234_013 =
    m[A20] * Det2_34_13 - m[A21] * Det2_34_03 + m[A23] * Det2_34_01;
  G4double Det3_234_014 =
    m[A20] * Det2_34_14 - m[A21] * Det2_34_04 + m[A24] * Det2_34_01;
  G4double Det3_234_023 =
    m[A20] * Det2_34_23 - m[A22] * Det2_34_03 + m[A23] * Det2_34_02;
  G4double Det3_234_024 =
    m[A20] * Det2_34_24 - m[A22] * Det2_34_04 + m[A24] * Det2_34_02;
  G4double Det3_234_034 =
    m[A20] * Det2_34_34 - m[A23] * Det2_34_04 + m[A24] * Det2_34_03;
  G4double Det3_234_123 =
    m[A21] * Det2_34_23 - m[A22] * Det2_34_13 + m[A23] * Det2_34_12;
  G4double Det3_234_124 =
    m[A21] * Det2_34_24 - m[A22] * Det2_34_14 + m[A24] * Det2_34_12;
  G4double Det3_234_134 =
    m[A21] * Det2_34_34 - m[A23] * Det2_34_14 + m[A24] * Det2_34_13;
  G4double Det3_234_234 =
    m[A22] * Det2_34_34 - m[A23] * Det2_34_24 + m[A24] * Det2_34_23;
  G4double Det3_235_012 =
    m[A20] * Det2_35_12 - m[A21] * Det2_35_02 + m[A22] * Det2_35_01;
  G4double Det3_235_013 =
    m[A20] * Det2_35_13 - m[A21] * Det2_35_03 + m[A23] * Det2_35_01;
  G4double Det3_235_014 =
    m[A20] * Det2_35_14 - m[A21] * Det2_35_04 + m[A24] * Det2_35_01;
  G4double Det3_235_015 =
    m[A20] * Det2_35_15 - m[A21] * Det2_35_05 + m[A25] * Det2_35_01;
  G4double Det3_235_023 =
    m[A20] * Det2_35_23 - m[A22] * Det2_35_03 + m[A23] * Det2_35_02;
  G4double Det3_235_024 =
    m[A20] * Det2_35_24 - m[A22] * Det2_35_04 + m[A24] * Det2_35_02;
  G4double Det3_235_025 =
    m[A20] * Det2_35_25 - m[A22] * Det2_35_05 + m[A25] * Det2_35_02;
  G4double Det3_235_034 =
    m[A20] * Det2_35_34 - m[A23] * Det2_35_04 + m[A24] * Det2_35_03;
  G4double Det3_235_035 =
    m[A20] * Det2_35_35 - m[A23] * Det2_35_05 + m[A25] * Det2_35_03;
  G4double Det3_235_123 =
    m[A21] * Det2_35_23 - m[A22] * Det2_35_13 + m[A23] * Det2_35_12;
  G4double Det3_235_124 =
    m[A21] * Det2_35_24 - m[A22] * Det2_35_14 + m[A24] * Det2_35_12;
  G4double Det3_235_125 =
    m[A21] * Det2_35_25 - m[A22] * Det2_35_15 + m[A25] * Det2_35_12;
  G4double Det3_235_134 =
    m[A21] * Det2_35_34 - m[A23] * Det2_35_14 + m[A24] * Det2_35_13;
  G4double Det3_235_135 =
    m[A21] * Det2_35_35 - m[A23] * Det2_35_15 + m[A25] * Det2_35_13;
  G4double Det3_235_234 =
    m[A22] * Det2_35_34 - m[A23] * Det2_35_24 + m[A24] * Det2_35_23;
  G4double Det3_235_235 =
    m[A22] * Det2_35_35 - m[A23] * Det2_35_25 + m[A25] * Det2_35_23;
  G4double Det3_245_012 =
    m[A20] * Det2_45_12 - m[A21] * Det2_45_02 + m[A22] * Det2_45_01;
  G4double Det3_245_013 =
    m[A20] * Det2_45_13 - m[A21] * Det2_45_03 + m[A23] * Det2_45_01;
  G4double Det3_245_014 =
    m[A20] * Det2_45_14 - m[A21] * Det2_45_04 + m[A24] * Det2_45_01;
  G4double Det3_245_015 =
    m[A20] * Det2_45_15 - m[A21] * Det2_45_05 + m[A25] * Det2_45_01;
  G4double Det3_245_023 =
    m[A20] * Det2_45_23 - m[A22] * Det2_45_03 + m[A23] * Det2_45_02;
  G4double Det3_245_024 =
    m[A20] * Det2_45_24 - m[A22] * Det2_45_04 + m[A24] * Det2_45_02;
  G4double Det3_245_025 =
    m[A20] * Det2_45_25 - m[A22] * Det2_45_05 + m[A25] * Det2_45_02;
  G4double Det3_245_034 =
    m[A20] * Det2_45_34 - m[A23] * Det2_45_04 + m[A24] * Det2_45_03;
  G4double Det3_245_035 =
    m[A20] * Det2_45_35 - m[A23] * Det2_45_05 + m[A25] * Det2_45_03;
  G4double Det3_245_045 =
    m[A20] * Det2_45_45 - m[A24] * Det2_45_05 + m[A25] * Det2_45_04;
  G4double Det3_245_123 =
    m[A21] * Det2_45_23 - m[A22] * Det2_45_13 + m[A23] * Det2_45_12;
  G4double Det3_245_124 =
    m[A21] * Det2_45_24 - m[A22] * Det2_45_14 + m[A24] * Det2_45_12;
  G4double Det3_245_125 =
    m[A21] * Det2_45_25 - m[A22] * Det2_45_15 + m[A25] * Det2_45_12;
  G4double Det3_245_134 =
    m[A21] * Det2_45_34 - m[A23] * Det2_45_14 + m[A24] * Det2_45_13;
  G4double Det3_245_135 =
    m[A21] * Det2_45_35 - m[A23] * Det2_45_15 + m[A25] * Det2_45_13;
  G4double Det3_245_145 =
    m[A21] * Det2_45_45 - m[A24] * Det2_45_15 + m[A25] * Det2_45_14;
  G4double Det3_245_234 =
    m[A22] * Det2_45_34 - m[A23] * Det2_45_24 + m[A24] * Det2_45_23;
  G4double Det3_245_235 =
    m[A22] * Det2_45_35 - m[A23] * Det2_45_25 + m[A25] * Det2_45_23;
  G4double Det3_245_245 =
    m[A22] * Det2_45_45 - m[A24] * Det2_45_25 + m[A25] * Det2_45_24;
  G4double Det3_345_012 =
    m[A30] * Det2_45_12 - m[A31] * Det2_45_02 + m[A32] * Det2_45_01;
  G4double Det3_345_013 =
    m[A30] * Det2_45_13 - m[A31] * Det2_45_03 + m[A33] * Det2_45_01;
  G4double Det3_345_014 =
    m[A30] * Det2_45_14 - m[A31] * Det2_45_04 + m[A34] * Det2_45_01;
  G4double Det3_345_015 =
    m[A30] * Det2_45_15 - m[A31] * Det2_45_05 + m[A35] * Det2_45_01;
  G4double Det3_345_023 =
    m[A30] * Det2_45_23 - m[A32] * Det2_45_03 + m[A33] * Det2_45_02;
  G4double Det3_345_024 =
    m[A30] * Det2_45_24 - m[A32] * Det2_45_04 + m[A34] * Det2_45_02;
  G4double Det3_345_025 =
    m[A30] * Det2_45_25 - m[A32] * Det2_45_05 + m[A35] * Det2_45_02;
  G4double Det3_345_034 =
    m[A30] * Det2_45_34 - m[A33] * Det2_45_04 + m[A34] * Det2_45_03;
  G4double Det3_345_035 =
    m[A30] * Det2_45_35 - m[A33] * Det2_45_05 + m[A35] * Det2_45_03;
  G4double Det3_345_045 =
    m[A30] * Det2_45_45 - m[A34] * Det2_45_05 + m[A35] * Det2_45_04;
  G4double Det3_345_123 =
    m[A31] * Det2_45_23 - m[A32] * Det2_45_13 + m[A33] * Det2_45_12;
  G4double Det3_345_124 =
    m[A31] * Det2_45_24 - m[A32] * Det2_45_14 + m[A34] * Det2_45_12;
  G4double Det3_345_125 =
    m[A31] * Det2_45_25 - m[A32] * Det2_45_15 + m[A35] * Det2_45_12;
  G4double Det3_345_134 =
    m[A31] * Det2_45_34 - m[A33] * Det2_45_14 + m[A34] * Det2_45_13;
  G4double Det3_345_135 =
    m[A31] * Det2_45_35 - m[A33] * Det2_45_15 + m[A35] * Det2_45_13;
  G4double Det3_345_145 =
    m[A31] * Det2_45_45 - m[A34] * Det2_45_15 + m[A35] * Det2_45_14;
  G4double Det3_345_234 =
    m[A32] * Det2_45_34 - m[A33] * Det2_45_24 + m[A34] * Det2_45_23;
  G4double Det3_345_235 =
    m[A32] * Det2_45_35 - m[A33] * Det2_45_25 + m[A35] * Det2_45_23;
  G4double Det3_345_245 =
    m[A32] * Det2_45_45 - m[A34] * Det2_45_25 + m[A35] * Det2_45_24;
  G4double Det3_345_345 =
    m[A33] * Det2_45_45 - m[A34] * Det2_45_35 + m[A35] * Det2_45_34;

  // Find all NECESSARY 4x4 dets:  (55 of them)

  G4double Det4_1234_0123 = m[A10] * Det3_234_123 - m[A11] * Det3_234_023 +
                            m[A12] * Det3_234_013 - m[A13] * Det3_234_012;
  G4double Det4_1234_0124 = m[A10] * Det3_234_124 - m[A11] * Det3_234_024 +
                            m[A12] * Det3_234_014 - m[A14] * Det3_234_012;
  G4double Det4_1234_0134 = m[A10] * Det3_234_134 - m[A11] * Det3_234_034 +
                            m[A13] * Det3_234_014 - m[A14] * Det3_234_013;
  G4double Det4_1234_0234 = m[A10] * Det3_234_234 - m[A12] * Det3_234_034 +
                            m[A13] * Det3_234_024 - m[A14] * Det3_234_023;
  G4double Det4_1234_1234 = m[A11] * Det3_234_234 - m[A12] * Det3_234_134 +
                            m[A13] * Det3_234_124 - m[A14] * Det3_234_123;
  G4double Det4_1235_0123 = m[A10] * Det3_235_123 - m[A11] * Det3_235_023 +
                            m[A12] * Det3_235_013 - m[A13] * Det3_235_012;
  G4double Det4_1235_0124 = m[A10] * Det3_235_124 - m[A11] * Det3_235_024 +
                            m[A12] * Det3_235_014 - m[A14] * Det3_235_012;
  G4double Det4_1235_0125 = m[A10] * Det3_235_125 - m[A11] * Det3_235_025 +
                            m[A12] * Det3_235_015 - m[A15] * Det3_235_012;
  G4double Det4_1235_0134 = m[A10] * Det3_235_134 - m[A11] * Det3_235_034 +
                            m[A13] * Det3_235_014 - m[A14] * Det3_235_013;
  G4double Det4_1235_0135 = m[A10] * Det3_235_135 - m[A11] * Det3_235_035 +
                            m[A13] * Det3_235_015 - m[A15] * Det3_235_013;
  G4double Det4_1235_0234 = m[A10] * Det3_235_234 - m[A12] * Det3_235_034 +
                            m[A13] * Det3_235_024 - m[A14] * Det3_235_023;
  G4double Det4_1235_0235 = m[A10] * Det3_235_235 - m[A12] * Det3_235_035 +
                            m[A13] * Det3_235_025 - m[A15] * Det3_235_023;
  G4double Det4_1235_1234 = m[A11] * Det3_235_234 - m[A12] * Det3_235_134 +
                            m[A13] * Det3_235_124 - m[A14] * Det3_235_123;
  G4double Det4_1235_1235 = m[A11] * Det3_235_235 - m[A12] * Det3_235_135 +
                            m[A13] * Det3_235_125 - m[A15] * Det3_235_123;
  G4double Det4_1245_0123 = m[A10] * Det3_245_123 - m[A11] * Det3_245_023 +
                            m[A12] * Det3_245_013 - m[A13] * Det3_245_012;
  G4double Det4_1245_0124 = m[A10] * Det3_245_124 - m[A11] * Det3_245_024 +
                            m[A12] * Det3_245_014 - m[A14] * Det3_245_012;
  G4double Det4_1245_0125 = m[A10] * Det3_245_125 - m[A11] * Det3_245_025 +
                            m[A12] * Det3_245_015 - m[A15] * Det3_245_012;
  G4double Det4_1245_0134 = m[A10] * Det3_245_134 - m[A11] * Det3_245_034 +
                            m[A13] * Det3_245_014 - m[A14] * Det3_245_013;
  G4double Det4_1245_0135 = m[A10] * Det3_245_135 - m[A11] * Det3_245_035 +
                            m[A13] * Det3_245_015 - m[A15] * Det3_245_013;
  G4double Det4_1245_0145 = m[A10] * Det3_245_145 - m[A11] * Det3_245_045 +
                            m[A14] * Det3_245_015 - m[A15] * Det3_245_014;
  G4double Det4_1245_0234 = m[A10] * Det3_245_234 - m[A12] * Det3_245_034 +
                            m[A13] * Det3_245_024 - m[A14] * Det3_245_023;
  G4double Det4_1245_0235 = m[A10] * Det3_245_235 - m[A12] * Det3_245_035 +
                            m[A13] * Det3_245_025 - m[A15] * Det3_245_023;
  G4double Det4_1245_0245 = m[A10] * Det3_245_245 - m[A12] * Det3_245_045 +
                            m[A14] * Det3_245_025 - m[A15] * Det3_245_024;
  G4double Det4_1245_1234 = m[A11] * Det3_245_234 - m[A12] * Det3_245_134 +
                            m[A13] * Det3_245_124 - m[A14] * Det3_245_123;
  G4double Det4_1245_1235 = m[A11] * Det3_245_235 - m[A12] * Det3_245_135 +
                            m[A13] * Det3_245_125 - m[A15] * Det3_245_123;
  G4double Det4_1245_1245 = m[A11] * Det3_245_245 - m[A12] * Det3_245_145 +
                            m[A14] * Det3_245_125 - m[A15] * Det3_245_124;
  G4double Det4_1345_0123 = m[A10] * Det3_345_123 - m[A11] * Det3_345_023 +
                            m[A12] * Det3_345_013 - m[A13] * Det3_345_012;
  G4double Det4_1345_0124 = m[A10] * Det3_345_124 - m[A11] * Det3_345_024 +
                            m[A12] * Det3_345_014 - m[A14] * Det3_345_012;
  G4double Det4_1345_0125 = m[A10] * Det3_345_125 - m[A11] * Det3_345_025 +
                            m[A12] * Det3_345_015 - m[A15] * Det3_345_012;
  G4double Det4_1345_0134 = m[A10] * Det3_345_134 - m[A11] * Det3_345_034 +
                            m[A13] * Det3_345_014 - m[A14] * Det3_345_013;
  G4double Det4_1345_0135 = m[A10] * Det3_345_135 - m[A11] * Det3_345_035 +
                            m[A13] * Det3_345_015 - m[A15] * Det3_345_013;
  G4double Det4_1345_0145 = m[A10] * Det3_345_145 - m[A11] * Det3_345_045 +
                            m[A14] * Det3_345_015 - m[A15] * Det3_345_014;
  G4double Det4_1345_0234 = m[A10] * Det3_345_234 - m[A12] * Det3_345_034 +
                            m[A13] * Det3_345_024 - m[A14] * Det3_345_023;
  G4double Det4_1345_0235 = m[A10] * Det3_345_235 - m[A12] * Det3_345_035 +
                            m[A13] * Det3_345_025 - m[A15] * Det3_345_023;
  G4double Det4_1345_0245 = m[A10] * Det3_345_245 - m[A12] * Det3_345_045 +
                            m[A14] * Det3_345_025 - m[A15] * Det3_345_024;
  G4double Det4_1345_0345 = m[A10] * Det3_345_345 - m[A13] * Det3_345_045 +
                            m[A14] * Det3_345_035 - m[A15] * Det3_345_034;
  G4double Det4_1345_1234 = m[A11] * Det3_345_234 - m[A12] * Det3_345_134 +
                            m[A13] * Det3_345_124 - m[A14] * Det3_345_123;
  G4double Det4_1345_1235 = m[A11] * Det3_345_235 - m[A12] * Det3_345_135 +
                            m[A13] * Det3_345_125 - m[A15] * Det3_345_123;
  G4double Det4_1345_1245 = m[A11] * Det3_345_245 - m[A12] * Det3_345_145 +
                            m[A14] * Det3_345_125 - m[A15] * Det3_345_124;
  G4double Det4_1345_1345 = m[A11] * Det3_345_345 - m[A13] * Det3_345_145 +
                            m[A14] * Det3_345_135 - m[A15] * Det3_345_134;
  G4double Det4_2345_0123 = m[A20] * Det3_345_123 - m[A21] * Det3_345_023 +
                            m[A22] * Det3_345_013 - m[A23] * Det3_345_012;
  G4double Det4_2345_0124 = m[A20] * Det3_345_124 - m[A21] * Det3_345_024 +
                            m[A22] * Det3_345_014 - m[A24] * Det3_345_012;
  G4double Det4_2345_0125 = m[A20] * Det3_345_125 - m[A21] * Det3_345_025 +
                            m[A22] * Det3_345_015 - m[A25] * Det3_345_012;
  G4double Det4_2345_0134 = m[A20] * Det3_345_134 - m[A21] * Det3_345_034 +
                            m[A23] * Det3_345_014 - m[A24] * Det3_345_013;
  G4double Det4_2345_0135 = m[A20] * Det3_345_135 - m[A21] * Det3_345_035 +
                            m[A23] * Det3_345_015 - m[A25] * Det3_345_013;
  G4double Det4_2345_0145 = m[A20] * Det3_345_145 - m[A21] * Det3_345_045 +
                            m[A24] * Det3_345_015 - m[A25] * Det3_345_014;
  G4double Det4_2345_0234 = m[A20] * Det3_345_234 - m[A22] * Det3_345_034 +
                            m[A23] * Det3_345_024 - m[A24] * Det3_345_023;
  G4double Det4_2345_0235 = m[A20] * Det3_345_235 - m[A22] * Det3_345_035 +
                            m[A23] * Det3_345_025 - m[A25] * Det3_345_023;
  G4double Det4_2345_0245 = m[A20] * Det3_345_245 - m[A22] * Det3_345_045 +
                            m[A24] * Det3_345_025 - m[A25] * Det3_345_024;
  G4double Det4_2345_0345 = m[A20] * Det3_345_345 - m[A23] * Det3_345_045 +
                            m[A24] * Det3_345_035 - m[A25] * Det3_345_034;
  G4double Det4_2345_1234 = m[A21] * Det3_345_234 - m[A22] * Det3_345_134 +
                            m[A23] * Det3_345_124 - m[A24] * Det3_345_123;
  G4double Det4_2345_1235 = m[A21] * Det3_345_235 - m[A22] * Det3_345_135 +
                            m[A23] * Det3_345_125 - m[A25] * Det3_345_123;
  G4double Det4_2345_1245 = m[A21] * Det3_345_245 - m[A22] * Det3_345_145 +
                            m[A24] * Det3_345_125 - m[A25] * Det3_345_124;
  G4double Det4_2345_1345 = m[A21] * Det3_345_345 - m[A23] * Det3_345_145 +
                            m[A24] * Det3_345_135 - m[A25] * Det3_345_134;
  G4double Det4_2345_2345 = m[A22] * Det3_345_345 - m[A23] * Det3_345_245 +
                            m[A24] * Det3_345_235 - m[A25] * Det3_345_234;

  // Find all NECESSARY 5x5 dets:  (19 of them)

  G4double Det5_01234_01234 =
    m[A00] * Det4_1234_1234 - m[A01] * Det4_1234_0234 +
    m[A02] * Det4_1234_0134 - m[A03] * Det4_1234_0124 + m[A04] * Det4_1234_0123;
  G4double Det5_01235_01234 =
    m[A00] * Det4_1235_1234 - m[A01] * Det4_1235_0234 +
    m[A02] * Det4_1235_0134 - m[A03] * Det4_1235_0124 + m[A04] * Det4_1235_0123;
  G4double Det5_01235_01235 =
    m[A00] * Det4_1235_1235 - m[A01] * Det4_1235_0235 +
    m[A02] * Det4_1235_0135 - m[A03] * Det4_1235_0125 + m[A05] * Det4_1235_0123;
  G4double Det5_01245_01234 =
    m[A00] * Det4_1245_1234 - m[A01] * Det4_1245_0234 +
    m[A02] * Det4_1245_0134 - m[A03] * Det4_1245_0124 + m[A04] * Det4_1245_0123;
  G4double Det5_01245_01235 =
    m[A00] * Det4_1245_1235 - m[A01] * Det4_1245_0235 +
    m[A02] * Det4_1245_0135 - m[A03] * Det4_1245_0125 + m[A05] * Det4_1245_0123;
  G4double Det5_01245_01245 =
    m[A00] * Det4_1245_1245 - m[A01] * Det4_1245_0245 +
    m[A02] * Det4_1245_0145 - m[A04] * Det4_1245_0125 + m[A05] * Det4_1245_0124;
  G4double Det5_01345_01234 =
    m[A00] * Det4_1345_1234 - m[A01] * Det4_1345_0234 +
    m[A02] * Det4_1345_0134 - m[A03] * Det4_1345_0124 + m[A04] * Det4_1345_0123;
  G4double Det5_01345_01235 =
    m[A00] * Det4_1345_1235 - m[A01] * Det4_1345_0235 +
    m[A02] * Det4_1345_0135 - m[A03] * Det4_1345_0125 + m[A05] * Det4_1345_0123;
  G4double Det5_01345_01245 =
    m[A00] * Det4_1345_1245 - m[A01] * Det4_1345_0245 +
    m[A02] * Det4_1345_0145 - m[A04] * Det4_1345_0125 + m[A05] * Det4_1345_0124;
  G4double Det5_01345_01345 =
    m[A00] * Det4_1345_1345 - m[A01] * Det4_1345_0345 +
    m[A03] * Det4_1345_0145 - m[A04] * Det4_1345_0135 + m[A05] * Det4_1345_0134;
  G4double Det5_02345_01234 =
    m[A00] * Det4_2345_1234 - m[A01] * Det4_2345_0234 +
    m[A02] * Det4_2345_0134 - m[A03] * Det4_2345_0124 + m[A04] * Det4_2345_0123;
  G4double Det5_02345_01235 =
    m[A00] * Det4_2345_1235 - m[A01] * Det4_2345_0235 +
    m[A02] * Det4_2345_0135 - m[A03] * Det4_2345_0125 + m[A05] * Det4_2345_0123;
  G4double Det5_02345_01245 =
    m[A00] * Det4_2345_1245 - m[A01] * Det4_2345_0245 +
    m[A02] * Det4_2345_0145 - m[A04] * Det4_2345_0125 + m[A05] * Det4_2345_0124;
  G4double Det5_02345_01345 =
    m[A00] * Det4_2345_1345 - m[A01] * Det4_2345_0345 +
    m[A03] * Det4_2345_0145 - m[A04] * Det4_2345_0135 + m[A05] * Det4_2345_0134;
  G4double Det5_02345_02345 =
    m[A00] * Det4_2345_2345 - m[A02] * Det4_2345_0345 +
    m[A03] * Det4_2345_0245 - m[A04] * Det4_2345_0235 + m[A05] * Det4_2345_0234;
  G4double Det5_12345_01234 =
    m[A10] * Det4_2345_1234 - m[A11] * Det4_2345_0234 +
    m[A12] * Det4_2345_0134 - m[A13] * Det4_2345_0124 + m[A14] * Det4_2345_0123;
  G4double Det5_12345_01235 =
    m[A10] * Det4_2345_1235 - m[A11] * Det4_2345_0235 +
    m[A12] * Det4_2345_0135 - m[A13] * Det4_2345_0125 + m[A15] * Det4_2345_0123;
  G4double Det5_12345_01245 =
    m[A10] * Det4_2345_1245 - m[A11] * Det4_2345_0245 +
    m[A12] * Det4_2345_0145 - m[A14] * Det4_2345_0125 + m[A15] * Det4_2345_0124;
  G4double Det5_12345_01345 =
    m[A10] * Det4_2345_1345 - m[A11] * Det4_2345_0345 +
    m[A13] * Det4_2345_0145 - m[A14] * Det4_2345_0135 + m[A15] * Det4_2345_0134;
  G4double Det5_12345_02345 =
    m[A10] * Det4_2345_2345 - m[A12] * Det4_2345_0345 +
    m[A13] * Det4_2345_0245 - m[A14] * Det4_2345_0235 + m[A15] * Det4_2345_0234;
  G4double Det5_12345_12345 =
    m[A11] * Det4_2345_2345 - m[A12] * Det4_2345_1345 +
    m[A13] * Det4_2345_1245 - m[A14] * Det4_2345_1235 + m[A15] * Det4_2345_1234;

  // Find the determinant

  G4double det = m[A00] * Det5_12345_12345 - m[A01] * Det5_12345_02345 +
                 m[A02] * Det5_12345_01345 - m[A03] * Det5_12345_01245 +
                 m[A04] * Det5_12345_01235 - m[A05] * Det5_12345_01234;

  if(det == 0)
  {
    ifail = 1;
    return;
  }

  G4double oneOverDet = 1.0 / det;
  G4double mn1OverDet = -oneOverDet;

  m[A00] = Det5_12345_12345 * oneOverDet;
  m[A01] = Det5_12345_02345 * mn1OverDet;
  m[A02] = Det5_12345_01345 * oneOverDet;
  m[A03] = Det5_12345_01245 * mn1OverDet;
  m[A04] = Det5_12345_01235 * oneOverDet;
  m[A05] = Det5_12345_01234 * mn1OverDet;

  m[A11] = Det5_02345_02345 * oneOverDet;
  m[A12] = Det5_02345_01345 * mn1OverDet;
  m[A13] = Det5_02345_01245 * oneOverDet;
  m[A14] = Det5_02345_01235 * mn1OverDet;
  m[A15] = Det5_02345_01234 * oneOverDet;

  m[A22] = Det5_01345_01345 * oneOverDet;
  m[A23] = Det5_01345_01245 * mn1OverDet;
  m[A24] = Det5_01345_01235 * oneOverDet;
  m[A25] = Det5_01345_01234 * mn1OverDet;

  m[A33] = Det5_01245_01245 * oneOverDet;
  m[A34] = Det5_01245_01235 * mn1OverDet;
  m[A35] = Det5_01245_01234 * oneOverDet;

  m[A44] = Det5_01235_01235 * oneOverDet;
  m[A45] = Det5_01235_01234 * mn1OverDet;

  m[A55] = Det5_01234_01234 * oneOverDet;

  return;
}

void G4ErrorSymMatrix::invertCholesky5(G4int& ifail)
{
  // Invert by
  //
  // a) decomposing M = G*G^T with G lower triangular
  //    (if M is not positive definite this will fail, leaving this unchanged)
  // b) inverting G to form H
  // c) multiplying H^T * H to get M^-1.
  //
  // If the matrix is pos. def. it is inverted and 1 is returned.
  // If the matrix is not pos. def. it remains unaltered and 0 is returned.

  G4double h10;  // below-diagonal elements of H
  G4double h20, h21;
  G4double h30, h31, h32;
  G4double h40, h41, h42, h43;

  G4double h00, h11, h22, h33, h44;  // 1/diagonal elements of G =
                                     // diagonal elements of H

  G4double g10;  // below-diagonal elements of G
  G4double g20, g21;
  G4double g30, g31, g32;
  G4double g40, g41, g42, g43;

  ifail = 1;  // We start by assuing we won't succeed...

  // Form G -- compute diagonal members of H directly rather than of G
  //-------

  // Scale first column by 1/sqrt(A00)

  h00 = m[A00];
  if(h00 <= 0)
  {
    return;
  }
  h00 = 1.0 / std::sqrt(h00);

  g10 = m[A10] * h00;
  g20 = m[A20] * h00;
  g30 = m[A30] * h00;
  g40 = m[A40] * h00;

  // Form G11 (actually, just h11)

  h11 = m[A11] - (g10 * g10);
  if(h11 <= 0)
  {
    return;
  }
  h11 = 1.0 / std::sqrt(h11);

  // Subtract inter-column column dot products from rest of column 1 and
  // scale to get column 1 of G

  g21 = (m[A21] - (g10 * g20)) * h11;
  g31 = (m[A31] - (g10 * g30)) * h11;
  g41 = (m[A41] - (g10 * g40)) * h11;

  // Form G22 (actually, just h22)

  h22 = m[A22] - (g20 * g20) - (g21 * g21);
  if(h22 <= 0)
  {
    return;
  }
  h22 = 1.0 / std::sqrt(h22);

  // Subtract inter-column column dot products from rest of column 2 and
  // scale to get column 2 of G

  g32 = (m[A32] - (g20 * g30) - (g21 * g31)) * h22;
  g42 = (m[A42] - (g20 * g40) - (g21 * g41)) * h22;

  // Form G33 (actually, just h33)

  h33 = m[A33] - (g30 * g30) - (g31 * g31) - (g32 * g32);
  if(h33 <= 0)
  {
    return;
  }
  h33 = 1.0 / std::sqrt(h33);

  // Subtract inter-column column dot product from A43 and scale to get G43

  g43 = (m[A43] - (g30 * g40) - (g31 * g41) - (g32 * g42)) * h33;

  // Finally form h44 - if this is possible inversion succeeds

  h44 = m[A44] - (g40 * g40) - (g41 * g41) - (g42 * g42) - (g43 * g43);
  if(h44 <= 0)
  {
    return;
  }
  h44 = 1.0 / std::sqrt(h44);

  // Form H = 1/G -- diagonal members of H are already correct
  //-------------

  // The order here is dictated by speed considerations

  h43 = -h33 * g43 * h44;
  h32 = -h22 * g32 * h33;
  h42 = -h22 * (g32 * h43 + g42 * h44);
  h21 = -h11 * g21 * h22;
  h31 = -h11 * (g21 * h32 + g31 * h33);
  h41 = -h11 * (g21 * h42 + g31 * h43 + g41 * h44);
  h10 = -h00 * g10 * h11;
  h20 = -h00 * (g10 * h21 + g20 * h22);
  h30 = -h00 * (g10 * h31 + g20 * h32 + g30 * h33);
  h40 = -h00 * (g10 * h41 + g20 * h42 + g30 * h43 + g40 * h44);

  // Change this to its inverse = H^T*H
  //------------------------------------

  m[A00] = h00 * h00 + h10 * h10 + h20 * h20 + h30 * h30 + h40 * h40;
  m[A01] = h10 * h11 + h20 * h21 + h30 * h31 + h40 * h41;
  m[A11] = h11 * h11 + h21 * h21 + h31 * h31 + h41 * h41;
  m[A02] = h20 * h22 + h30 * h32 + h40 * h42;
  m[A12] = h21 * h22 + h31 * h32 + h41 * h42;
  m[A22] = h22 * h22 + h32 * h32 + h42 * h42;
  m[A03] = h30 * h33 + h40 * h43;
  m[A13] = h31 * h33 + h41 * h43;
  m[A23] = h32 * h33 + h42 * h43;
  m[A33] = h33 * h33 + h43 * h43;
  m[A04] = h40 * h44;
  m[A14] = h41 * h44;
  m[A24] = h42 * h44;
  m[A34] = h43 * h44;
  m[A44] = h44 * h44;

  ifail = 0;
  return;
}

void G4ErrorSymMatrix::invertCholesky6(G4int& ifail)
{
  // Invert by
  //
  // a) decomposing M = G*G^T with G lower triangular
  //    (if M is not positive definite this will fail, leaving this unchanged)
  // b) inverting G to form H
  // c) multiplying H^T * H to get M^-1.
  //
  // If the matrix is pos. def. it is inverted and 1 is returned.
  // If the matrix is not pos. def. it remains unaltered and 0 is returned.

  G4double h10;  // below-diagonal elements of H
  G4double h20, h21;
  G4double h30, h31, h32;
  G4double h40, h41, h42, h43;
  G4double h50, h51, h52, h53, h54;

  G4double h00, h11, h22, h33, h44, h55;  // 1/diagonal elements of G =
                                          // diagonal elements of H

  G4double g10;  // below-diagonal elements of G
  G4double g20, g21;
  G4double g30, g31, g32;
  G4double g40, g41, g42, g43;
  G4double g50, g51, g52, g53, g54;

  ifail = 1;  // We start by assuing we won't succeed...

  // Form G -- compute diagonal members of H directly rather than of G
  //-------

  // Scale first column by 1/sqrt(A00)

  h00 = m[A00];
  if(h00 <= 0)
  {
    return;
  }
  h00 = 1.0 / std::sqrt(h00);

  g10 = m[A10] * h00;
  g20 = m[A20] * h00;
  g30 = m[A30] * h00;
  g40 = m[A40] * h00;
  g50 = m[A50] * h00;

  // Form G11 (actually, just h11)

  h11 = m[A11] - (g10 * g10);
  if(h11 <= 0)
  {
    return;
  }
  h11 = 1.0 / std::sqrt(h11);

  // Subtract inter-column column dot products from rest of column 1 and
  // scale to get column 1 of G

  g21 = (m[A21] - (g10 * g20)) * h11;
  g31 = (m[A31] - (g10 * g30)) * h11;
  g41 = (m[A41] - (g10 * g40)) * h11;
  g51 = (m[A51] - (g10 * g50)) * h11;

  // Form G22 (actually, just h22)

  h22 = m[A22] - (g20 * g20) - (g21 * g21);
  if(h22 <= 0)
  {
    return;
  }
  h22 = 1.0 / std::sqrt(h22);

  // Subtract inter-column column dot products from rest of column 2 and
  // scale to get column 2 of G

  g32 = (m[A32] - (g20 * g30) - (g21 * g31)) * h22;
  g42 = (m[A42] - (g20 * g40) - (g21 * g41)) * h22;
  g52 = (m[A52] - (g20 * g50) - (g21 * g51)) * h22;

  // Form G33 (actually, just h33)

  h33 = m[A33] - (g30 * g30) - (g31 * g31) - (g32 * g32);
  if(h33 <= 0)
  {
    return;
  }
  h33 = 1.0 / std::sqrt(h33);

  // Subtract inter-column column dot products from rest of column 3 and
  // scale to get column 3 of G

  g43 = (m[A43] - (g30 * g40) - (g31 * g41) - (g32 * g42)) * h33;
  g53 = (m[A53] - (g30 * g50) - (g31 * g51) - (g32 * g52)) * h33;

  // Form G44 (actually, just h44)

  h44 = m[A44] - (g40 * g40) - (g41 * g41) - (g42 * g42) - (g43 * g43);
  if(h44 <= 0)
  {
    return;
  }
  h44 = 1.0 / std::sqrt(h44);

  // Subtract inter-column column dot product from M54 and scale to get G54

  g54 = (m[A54] - (g40 * g50) - (g41 * g51) - (g42 * g52) - (g43 * g53)) * h44;

  // Finally form h55 - if this is possible inversion succeeds

  h55 = m[A55] - (g50 * g50) - (g51 * g51) - (g52 * g52) - (g53 * g53) -
        (g54 * g54);
  if(h55 <= 0)
  {
    return;
  }
  h55 = 1.0 / std::sqrt(h55);

  // Form H = 1/G -- diagonal members of H are already correct
  //-------------

  // The order here is dictated by speed considerations

  h54 = -h44 * g54 * h55;
  h43 = -h33 * g43 * h44;
  h53 = -h33 * (g43 * h54 + g53 * h55);
  h32 = -h22 * g32 * h33;
  h42 = -h22 * (g32 * h43 + g42 * h44);
  h52 = -h22 * (g32 * h53 + g42 * h54 + g52 * h55);
  h21 = -h11 * g21 * h22;
  h31 = -h11 * (g21 * h32 + g31 * h33);
  h41 = -h11 * (g21 * h42 + g31 * h43 + g41 * h44);
  h51 = -h11 * (g21 * h52 + g31 * h53 + g41 * h54 + g51 * h55);
  h10 = -h00 * g10 * h11;
  h20 = -h00 * (g10 * h21 + g20 * h22);
  h30 = -h00 * (g10 * h31 + g20 * h32 + g30 * h33);
  h40 = -h00 * (g10 * h41 + g20 * h42 + g30 * h43 + g40 * h44);
  h50 = -h00 * (g10 * h51 + g20 * h52 + g30 * h53 + g40 * h54 + g50 * h55);

  // Change this to its inverse = H^T*H
  //------------------------------------

  m[A00] =
    h00 * h00 + h10 * h10 + h20 * h20 + h30 * h30 + h40 * h40 + h50 * h50;
  m[A01] = h10 * h11 + h20 * h21 + h30 * h31 + h40 * h41 + h50 * h51;
  m[A11] = h11 * h11 + h21 * h21 + h31 * h31 + h41 * h41 + h51 * h51;
  m[A02] = h20 * h22 + h30 * h32 + h40 * h42 + h50 * h52;
  m[A12] = h21 * h22 + h31 * h32 + h41 * h42 + h51 * h52;
  m[A22] = h22 * h22 + h32 * h32 + h42 * h42 + h52 * h52;
  m[A03] = h30 * h33 + h40 * h43 + h50 * h53;
  m[A13] = h31 * h33 + h41 * h43 + h51 * h53;
  m[A23] = h32 * h33 + h42 * h43 + h52 * h53;
  m[A33] = h33 * h33 + h43 * h43 + h53 * h53;
  m[A04] = h40 * h44 + h50 * h54;
  m[A14] = h41 * h44 + h51 * h54;
  m[A24] = h42 * h44 + h52 * h54;
  m[A34] = h43 * h44 + h53 * h54;
  m[A44] = h44 * h44 + h54 * h54;
  m[A05] = h50 * h55;
  m[A15] = h51 * h55;
  m[A25] = h52 * h55;
  m[A35] = h53 * h55;
  m[A45] = h54 * h55;
  m[A55] = h55 * h55;

  ifail = 0;
  return;
}

void G4ErrorSymMatrix::invert4(G4int& ifail)
{
  ifail = 0;

  // Find all NECESSARY 2x2 dets:  (14 of them)

  G4double Det2_12_01 = m[A10] * m[A21] - m[A11] * m[A20];
  G4double Det2_12_02 = m[A10] * m[A22] - m[A12] * m[A20];
  G4double Det2_12_12 = m[A11] * m[A22] - m[A12] * m[A21];
  G4double Det2_13_01 = m[A10] * m[A31] - m[A11] * m[A30];
  G4double Det2_13_02 = m[A10] * m[A32] - m[A12] * m[A30];
  G4double Det2_13_03 = m[A10] * m[A33] - m[A13] * m[A30];
  G4double Det2_13_12 = m[A11] * m[A32] - m[A12] * m[A31];
  G4double Det2_13_13 = m[A11] * m[A33] - m[A13] * m[A31];
  G4double Det2_23_01 = m[A20] * m[A31] - m[A21] * m[A30];
  G4double Det2_23_02 = m[A20] * m[A32] - m[A22] * m[A30];
  G4double Det2_23_03 = m[A20] * m[A33] - m[A23] * m[A30];
  G4double Det2_23_12 = m[A21] * m[A32] - m[A22] * m[A31];
  G4double Det2_23_13 = m[A21] * m[A33] - m[A23] * m[A31];
  G4double Det2_23_23 = m[A22] * m[A33] - m[A23] * m[A32];

  // Find all NECESSARY 3x3 dets:   (10 of them)

  G4double Det3_012_012 =
    m[A00] * Det2_12_12 - m[A01] * Det2_12_02 + m[A02] * Det2_12_01;
  G4double Det3_013_012 =
    m[A00] * Det2_13_12 - m[A01] * Det2_13_02 + m[A02] * Det2_13_01;
  G4double Det3_013_013 =
    m[A00] * Det2_13_13 - m[A01] * Det2_13_03 + m[A03] * Det2_13_01;
  G4double Det3_023_012 =
    m[A00] * Det2_23_12 - m[A01] * Det2_23_02 + m[A02] * Det2_23_01;
  G4double Det3_023_013 =
    m[A00] * Det2_23_13 - m[A01] * Det2_23_03 + m[A03] * Det2_23_01;
  G4double Det3_023_023 =
    m[A00] * Det2_23_23 - m[A02] * Det2_23_03 + m[A03] * Det2_23_02;
  G4double Det3_123_012 =
    m[A10] * Det2_23_12 - m[A11] * Det2_23_02 + m[A12] * Det2_23_01;
  G4double Det3_123_013 =
    m[A10] * Det2_23_13 - m[A11] * Det2_23_03 + m[A13] * Det2_23_01;
  G4double Det3_123_023 =
    m[A10] * Det2_23_23 - m[A12] * Det2_23_03 + m[A13] * Det2_23_02;
  G4double Det3_123_123 =
    m[A11] * Det2_23_23 - m[A12] * Det2_23_13 + m[A13] * Det2_23_12;

  // Find the 4x4 det:

  G4double det = m[A00] * Det3_123_123 - m[A01] * Det3_123_023 +
                 m[A02] * Det3_123_013 - m[A03] * Det3_123_012;

  if(det == 0)
  {
    ifail = 1;
    return;
  }

  G4double oneOverDet = 1.0 / det;
  G4double mn1OverDet = -oneOverDet;

  m[A00] = Det3_123_123 * oneOverDet;
  m[A01] = Det3_123_023 * mn1OverDet;
  m[A02] = Det3_123_013 * oneOverDet;
  m[A03] = Det3_123_012 * mn1OverDet;

  m[A11] = Det3_023_023 * oneOverDet;
  m[A12] = Det3_023_013 * mn1OverDet;
  m[A13] = Det3_023_012 * oneOverDet;

  m[A22] = Det3_013_013 * oneOverDet;
  m[A23] = Det3_013_012 * mn1OverDet;

  m[A33] = Det3_012_012 * oneOverDet;

  return;
}

void G4ErrorSymMatrix::invertHaywood4(G4int& ifail)
{
  invert4(ifail);  // For the 4x4 case, the method we use for invert is already
                   // the Haywood method.
}
