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
//      GEANT 4 class  header file 
// ------------------------------------------------------------
//
// Class Description:
//
// Simplified version of CLHEP HepSymMatrix class
//
// History:
// - Created:   P. Arce    May 2007
// --------------------------------------------------------------------

#ifndef G4ErrorSymMatrix_hh
#define G4ErrorSymMatrix_hh

#ifdef GNUPRAGMA
#pragma interface
#endif

#define DISABLE_ALLOC

#include <vector>

//#include "G4ErrorMatrix.hh"

class G4ErrorMatrix;

/**
 * @author
 * @ingroup matrix
 */
class G4ErrorSymMatrix  {
public:
   inline G4ErrorSymMatrix();
   // Default constructor. Gives 0x0 symmetric matrix.
   // Another G4ErrorSymMatrix can be assigned to it.

   explicit G4ErrorSymMatrix(int p);
   G4ErrorSymMatrix(int p, int);
   // Constructor. Gives p x p symmetric matrix.
   // With a second argument, the matrix is initialized. 0 means a zero
   // matrix, 1 means the identity matrix.

   G4ErrorSymMatrix(const G4ErrorSymMatrix &m1);
   // Copy constructor.

   // Constructor from DiagMatrix

   virtual ~G4ErrorSymMatrix();
   // Destructor.

   inline int num_row() const;
   inline int num_col() const;
   // Returns number of rows/columns.

   const double & operator()(int row, int col) const; 
   double & operator()(int row, int col);
   // Read and write a G4ErrorSymMatrix element.
   // ** Note that indexing starts from (1,1). **

   const double & fast(int row, int col) const;
   double & fast(int row, int col);
   // fast element access.
   // Must be row>=col;
   // ** Note that indexing starts from (1,1). **

   void assign(const G4ErrorMatrix &m2);
   // Assigns m2 to s, assuming m2 is a symmetric matrix.

   void assign(const G4ErrorSymMatrix &m2);
   // Another form of assignment. For consistency.

   G4ErrorSymMatrix & operator*=(double t);
   // Multiply a G4ErrorSymMatrix by a floating number.

   G4ErrorSymMatrix & operator/=(double t); 
   // Divide a G4ErrorSymMatrix by a floating number.

   G4ErrorSymMatrix & operator+=( const G4ErrorSymMatrix &m2);
   G4ErrorSymMatrix & operator-=( const G4ErrorSymMatrix &m2);
   // Add or subtract a G4ErrorSymMatrix.

   G4ErrorSymMatrix & operator=( const G4ErrorSymMatrix &m2);
   // Assignment operators. Notice that there is no G4ErrorSymMatrix = Matrix.

   G4ErrorSymMatrix operator- () const;
   // unary minus, ie. flip the sign of each element.

   G4ErrorSymMatrix T() const;
   // Returns the transpose of a G4ErrorSymMatrix (which is itself).

   G4ErrorSymMatrix apply(double (*f)(double, int, int)) const;
   // Apply a function to all elements of the matrix.

   G4ErrorSymMatrix similarity(const G4ErrorMatrix &m1) const;
   G4ErrorSymMatrix similarity(const G4ErrorSymMatrix &m1) const;
   // Returns m1*s*m1.T().

   G4ErrorSymMatrix similarityT(const G4ErrorMatrix &m1) const;
   // temporary. test of new similarity.
   // Returns m1.T()*s*m1.

   G4ErrorSymMatrix sub(int min_row, int max_row) const;
   // Returns a sub matrix of a G4ErrorSymMatrix.
   void sub(int row, const G4ErrorSymMatrix &m1);
   // Sub matrix of this G4ErrorSymMatrix is replaced with m1.
   G4ErrorSymMatrix sub(int min_row, int max_row);
   // SGI CC bug. I have to have both with/without const. I should not need
   // one without const.

   inline G4ErrorSymMatrix inverse(int &ifail) const;
   // Invert a Matrix. The matrix is not changed
   // Returns 0 when successful, otherwise non-zero.

   void invert(int &ifail);
   // Invert a Matrix.
   // N.B. the contents of the matrix are replaced by the inverse.
   // Returns ierr = 0 when successful, otherwise non-zero. 
   // This method has less overhead then inverse().

   double determinant() const;
   // calculate the determinant of the matrix.

   double trace() const;
   // calculate the trace of the matrix (sum of diagonal elements).

   class G4ErrorSymMatrix_row {
   public:
      inline G4ErrorSymMatrix_row(G4ErrorSymMatrix&,int);
      inline double & operator[](int);
   private:
      G4ErrorSymMatrix& _a;
      int _r;
   };
   class G4ErrorSymMatrix_row_const {
   public:
      inline G4ErrorSymMatrix_row_const(const G4ErrorSymMatrix&,int);
      inline const double & operator[](int) const;
   private:
      const G4ErrorSymMatrix& _a;
      int _r;
   };
   // helper class to implement m[i][j]

   inline G4ErrorSymMatrix_row operator[] (int);
   inline G4ErrorSymMatrix_row_const operator[] (int) const;
   // Read or write a matrix element.
   // While it may not look like it, you simply do m[i][j] to get an
   // element. 
   // ** Note that the indexing starts from [0][0]. **

   // Special-case inversions for 5x5 and 6x6 symmetric positive definite:
   // These set ifail=0 and invert if the matrix was positive definite;
   // otherwise ifail=1 and the matrix is left unaltered.
   void invertCholesky5 (int &ifail);  
   void invertCholesky6 (int &ifail);

   // Inversions for 5x5 and 6x6 forcing use of specific methods:  The
   // behavior (though not the speed) will be identical to invert(ifail).
   void invertHaywood4 (int & ifail);  
   void invertHaywood5 (int &ifail);  
   void invertHaywood6 (int &ifail);
   void invertBunchKaufman (int &ifail);  

protected:
   inline int num_size() const;
  
private:
   friend class G4ErrorSymMatrix_row;
   friend class G4ErrorSymMatrix_row_const;
   friend class G4ErrorMatrix;

   friend void tridiagonal(G4ErrorSymMatrix *a,G4ErrorMatrix *hsm);
   friend double condition(const G4ErrorSymMatrix &m);
   friend void diag_step(G4ErrorSymMatrix *t,int begin,int end);
   friend void diag_step(G4ErrorSymMatrix *t,G4ErrorMatrix *u,int begin,int end);
   friend G4ErrorMatrix diagonalize(G4ErrorSymMatrix *s);
   friend void house_with_update2(G4ErrorSymMatrix *a,G4ErrorMatrix *v,int row,int col);

   friend G4ErrorSymMatrix operator+(const G4ErrorSymMatrix &m1, 
				  const G4ErrorSymMatrix &m2);
   friend G4ErrorSymMatrix operator-(const G4ErrorSymMatrix &m1, 
				  const G4ErrorSymMatrix &m2);
   friend G4ErrorMatrix operator*(const G4ErrorSymMatrix &m1, const G4ErrorSymMatrix &m2);
   friend G4ErrorMatrix operator*(const G4ErrorSymMatrix &m1, const G4ErrorMatrix &m2);
   friend G4ErrorMatrix operator*(const G4ErrorMatrix &m1, const G4ErrorSymMatrix &m2);
   // Multiply a Matrix by a Matrix or Vector.
   
   // Returns v * v.T();

#ifdef DISABLE_ALLOC
   std::vector<double > m;
#else
   std::vector<double,Alloc<double,25> > m;
#endif
   int nrow;
   int size;				     // total number of elements

   static double posDefFraction5x5;
   static double adjustment5x5;
   static const  double CHOLESKY_THRESHOLD_5x5;
   static const  double CHOLESKY_CREEP_5x5;

   static double posDefFraction6x6;
   static double adjustment6x6;
   static const double CHOLESKY_THRESHOLD_6x6;
   static const double CHOLESKY_CREEP_6x6;

   void invert4  (int & ifail);
   void invert5  (int & ifail);
   void invert6  (int & ifail);

};

//
// Operations other than member functions for Matrix, G4ErrorSymMatrix, DiagMatrix
// and Vectors implemented in Matrix.cc and Matrix.icc (inline).
//

std::ostream& operator<<(std::ostream &s, const G4ErrorSymMatrix &q);
// Write out Matrix, G4ErrorSymMatrix, DiagMatrix and Vector into ostream.

G4ErrorMatrix operator*(const G4ErrorMatrix &m1, const G4ErrorSymMatrix &m2);
G4ErrorMatrix operator*(const G4ErrorSymMatrix &m1, const G4ErrorMatrix &m2);
G4ErrorMatrix operator*(const G4ErrorSymMatrix &m1, const G4ErrorSymMatrix &m2);
G4ErrorSymMatrix operator*(double t, const G4ErrorSymMatrix &s1);
G4ErrorSymMatrix operator*(const G4ErrorSymMatrix &s1, double t);
// Multiplication operators.
// Note that m *= m1 is always faster than m = m * m1

G4ErrorSymMatrix operator/(const G4ErrorSymMatrix &m1, double t);
// s = s1 / t. (s /= t is faster if you can use it.)

G4ErrorMatrix operator+(const G4ErrorMatrix &m1, const G4ErrorSymMatrix &s2);
G4ErrorMatrix operator+(const G4ErrorSymMatrix &s1, const G4ErrorMatrix &m2);
G4ErrorSymMatrix operator+(const G4ErrorSymMatrix &s1, const G4ErrorSymMatrix &s2);
// Addition operators

G4ErrorMatrix operator-(const G4ErrorMatrix &m1, const G4ErrorSymMatrix &s2);
G4ErrorMatrix operator-(const G4ErrorSymMatrix &m1, const G4ErrorMatrix &m2);
G4ErrorSymMatrix operator-(const G4ErrorSymMatrix &s1, const G4ErrorSymMatrix &s2);
// subtraction operators

G4ErrorSymMatrix dsum(const G4ErrorSymMatrix &s1, const G4ErrorSymMatrix &s2);
// Direct sum of two symmetric matrices;

double condition(const G4ErrorSymMatrix &m);
// Find the conditon number of a symmetric matrix.

void diag_step(G4ErrorSymMatrix *t, int begin, int end);
void diag_step(G4ErrorSymMatrix *t, G4ErrorMatrix *u, int begin, int end);
// Implicit symmetric QR step with Wilkinson Shift

G4ErrorMatrix diagonalize(G4ErrorSymMatrix *s);
// Diagonalize a symmetric matrix.
// It returns the matrix U so that s_old = U * s_diag * U.T()

void house_with_update2(G4ErrorSymMatrix *a, G4ErrorMatrix *v, int row=1, int col=1);
// Finds and does Householder reflection on matrix.

void tridiagonal(G4ErrorSymMatrix *a, G4ErrorMatrix *hsm);
G4ErrorMatrix tridiagonal(G4ErrorSymMatrix *a);
// Does a Householder tridiagonalization of a symmetric matrix.


#include "G4ErrorSymMatrix.icc"

#endif
