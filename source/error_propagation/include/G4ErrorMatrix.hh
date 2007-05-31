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
// Simplified version of CLHEP HepMatrix class
//
// History:
// - Created:   P. Arce    May 2007
// --------------------------------------------------------------------

#ifndef G4ErrorMatrix_hh
#define G4ErrorMatrix_hh

#ifdef GNUPRAGMA
#pragma interface
#endif
#define DISABLE_ALLOC

#include <vector>


class G4ErrorSymMatrix;

typedef std::vector<double >::iterator mIter;
typedef std::vector<double >::const_iterator mcIter;

/**
 * @author
 * @ingroup matrix
 */
class G4ErrorMatrix {
public:
   inline G4ErrorMatrix();
   // Default constructor. Gives 0 x 0 matrix. Another G4ErrorMatrix can be
   // assigned to it.

   G4ErrorMatrix(int p, int q);
   // Constructor. Gives an unitialized p x q matrix.
   G4ErrorMatrix(int p, int q, int i);
   // Constructor. Gives an initialized p x q matrix. 
   // If i=0, it is initialized to all 0. If i=1, the diagonal elements
   // are set to 1.0.

   G4ErrorMatrix(const G4ErrorMatrix &m1);
   // Copy constructor.

   G4ErrorMatrix(const G4ErrorSymMatrix &m1);
   // Constructors from G4ErrorSymG4ErrorMatrix, DiagG4ErrorMatrix and Vector.

   virtual ~G4ErrorMatrix();
   // Destructor.

   inline virtual int num_row() const;
   // Returns the number of rows.

   inline virtual int num_col() const;
   // Returns the number of columns.

   inline virtual const double & operator()(int row, int col) const;
   inline virtual double & operator()(int row, int col);
   // Read or write a matrix element. 
   // ** Note that the indexing starts from (1,1). **

   G4ErrorMatrix & operator *= (double t);
   // Multiply a G4ErrorMatrix by a floating number.

   G4ErrorMatrix & operator /= (double t); 
   // Divide a G4ErrorMatrix by a floating number.

   G4ErrorMatrix & operator += ( const G4ErrorMatrix &m2);
   G4ErrorMatrix & operator += ( const G4ErrorSymMatrix &m2);
   G4ErrorMatrix & operator -= ( const G4ErrorMatrix &m2);
   G4ErrorMatrix & operator -= ( const G4ErrorSymMatrix &m2);
   // Add or subtract a G4ErrorMatrix. 
   // When adding/subtracting Vector, G4ErrorMatrix must have num_col of one.

   G4ErrorMatrix & operator = ( const G4ErrorMatrix &m2);
   G4ErrorMatrix & operator = ( const G4ErrorSymMatrix &m2);
   // Assignment operators.

   G4ErrorMatrix operator- () const;
   // unary minus, ie. flip the sign of each element.

   G4ErrorMatrix apply(double (*f)(double, int, int)) const;
   // Apply a function to all elements of the matrix.

   G4ErrorMatrix T() const;
   // Returns the transpose of a G4ErrorMatrix.

   G4ErrorMatrix sub(int min_row, int max_row, int min_col, int max_col) const;
   // Returns a sub matrix of a G4ErrorMatrix.
   // WARNING: rows and columns are numbered from 1
   void sub(int row, int col, const G4ErrorMatrix &m1);
   // Sub matrix of this G4ErrorMatrix is replaced with m1.
   // WARNING: rows and columns are numbered from 1

  //   friend inline void swap(G4ErrorMatrix &m1, G4ErrorMatrix &m2);
   // Swap m1 with m2.

   inline G4ErrorMatrix inverse(int& ierr) const;
   // Invert a G4ErrorMatrix. G4ErrorMatrix must be square and is not changed.
   // Returns ierr = 0 (zero) when successful, otherwise non-zero.

   virtual void invert(int& ierr);
   // Invert a G4ErrorMatrix. G4ErrorMatrix must be square.
   // N.B. the contents of the matrix are replaced by the inverse.
   // Returns ierr = 0 (zero) when successful, otherwise non-zero. 
   // This method has less overhead then inverse().

   double determinant() const;
   // calculate the determinant of the matrix.

   double trace() const;
   // calculate the trace of the matrix (sum of diagonal elements).

   class G4ErrorMatrix_row {
   typedef std::vector<double >::const_iterator mcIter;
   public:
      inline G4ErrorMatrix_row(G4ErrorMatrix&,int);
      double & operator[](int);
   private:
      G4ErrorMatrix& _a;
      int _r;
   };
   class G4ErrorMatrix_row_const {
   public:
      inline G4ErrorMatrix_row_const (const G4ErrorMatrix&,int);
      const double & operator[](int) const;
   private:
      const G4ErrorMatrix& _a;
      int _r;
   };
   // helper classes for implementing m[i][j]

   inline G4ErrorMatrix_row operator[] (int);
   inline const G4ErrorMatrix_row_const operator[] (int) const;
   // Read or write a matrix element.
   // While it may not look like it, you simply do m[i][j] to get an
   // element. 
   // ** Note that the indexing starts from [0][0]. **

protected:
   virtual inline int num_size() const;
   virtual void invertHaywood4(int& ierr);
   virtual void invertHaywood5(int& ierr);
   virtual void invertHaywood6(int& ierr);

public:
   static void error(const char *s);

private:
   friend class G4ErrorMatrix_row;
   friend class G4ErrorMatrix_row_const;
   friend class G4ErrorSymMatrix;
   // Friend classes.

   friend G4ErrorMatrix operator+(const G4ErrorMatrix &m1, const G4ErrorMatrix &m2);
   friend G4ErrorMatrix operator-(const G4ErrorMatrix &m1, const G4ErrorMatrix &m2);
   friend G4ErrorMatrix operator*(const G4ErrorMatrix &m1, const G4ErrorMatrix &m2);
   friend G4ErrorMatrix operator*(const G4ErrorMatrix &m1, const G4ErrorSymMatrix &m2);
   friend G4ErrorMatrix operator*(const G4ErrorSymMatrix &m1, const G4ErrorMatrix &m2);
   friend G4ErrorMatrix operator*(const G4ErrorSymMatrix &m1, const G4ErrorSymMatrix &m2);
   // Multiply a G4ErrorMatrix by a G4ErrorMatrix or Vector.

   // solve the system of linear eq
   friend G4ErrorMatrix qr_solve(G4ErrorMatrix *, const G4ErrorMatrix &b);
   friend void tridiagonal(G4ErrorSymMatrix *a,G4ErrorMatrix *hsm);
   friend void row_house(G4ErrorMatrix *,const G4ErrorMatrix &, double,
			 int, int, int, int);
   friend void back_solve(const G4ErrorMatrix &R, G4ErrorMatrix *b);
   friend void col_givens(G4ErrorMatrix *A, double c,
			  double s, int k1, int k2, 
			  int rowmin, int rowmax);
   //    Does a column Givens update.
   friend void row_givens(G4ErrorMatrix *A, double c,
			  double s, int k1, int k2, 
			  int colmin, int colmax);
   friend void col_house(G4ErrorMatrix *,const G4ErrorMatrix &, double,
			 int, int, int, int);
   friend void house_with_update(G4ErrorMatrix *a,int row,int col);
   friend void house_with_update(G4ErrorMatrix *a,G4ErrorMatrix *v,int row,int col);
   friend void house_with_update2(G4ErrorSymMatrix *a,G4ErrorMatrix *v,
				  int row,int col); 

   int dfact_matrix(double &det, int *ir);
   // factorize the matrix. If successful, the return code is 0. On
   // return, det is the determinant and ir[] is row-interchange
   // matrix. See CERNLIB's DFACT routine.

   int dfinv_matrix(int *ir);
   // invert the matrix. See CERNLIB DFINV.

#ifdef DISABLE_ALLOC
   std::vector<double > m;
#else
   std::vector<double,Alloc<double,25> > m;
#endif
   int nrow, ncol;
   int size;
};

// Operations other than member functions for G4ErrorMatrix
// implemented in G4ErrorMatrix.cc and G4ErrorMatrix.icc (inline).

G4ErrorMatrix operator*(const G4ErrorMatrix &m1, const G4ErrorMatrix &m2);
G4ErrorMatrix operator*(double t, const G4ErrorMatrix &m1);
G4ErrorMatrix operator*(const G4ErrorMatrix &m1, double t);
// Multiplication operators
// Note that m *= m1 is always faster than m = m * m1.

G4ErrorMatrix operator/(const G4ErrorMatrix &m1, double t);
// m = m1 / t. (m /= t is faster if you can use it.)

G4ErrorMatrix operator+(const G4ErrorMatrix &m1, const G4ErrorMatrix &m2);
// m = m1 + m2;
// Note that m += m1 is always faster than m = m + m1.

G4ErrorMatrix operator-(const G4ErrorMatrix &m1, const G4ErrorMatrix &m2);
// m = m1 - m2;
// Note that m -= m1 is always faster than m = m - m1.

G4ErrorMatrix dsum(const G4ErrorMatrix&, const G4ErrorMatrix&);
// Direct sum of two matrices. The direct sum of A and B is the matrix 
//        A 0
//        0 B

std::ostream& operator<<(std::ostream &s, const G4ErrorMatrix &q);
// Read in, write out G4ErrorMatrix into a stream.
 
//
// Specialized linear algebra functions
//

G4ErrorMatrix qr_solve(const G4ErrorMatrix &A, const G4ErrorMatrix &b);
G4ErrorMatrix qr_solve(G4ErrorMatrix *A, const G4ErrorMatrix &b);
// Works like backsolve, except matrix does not need to be upper
// triangular. For nonsquare matrix, it solves in the least square sense.

G4ErrorMatrix qr_inverse(const G4ErrorMatrix &A);
G4ErrorMatrix qr_inverse(G4ErrorMatrix *A);
// Finds the inverse of a matrix using QR decomposition.  Note, often what
// you really want is solve or backsolve, they can be much quicker than
// inverse in many calculations.


void qr_decomp(G4ErrorMatrix *A, G4ErrorMatrix *hsm);
G4ErrorMatrix qr_decomp(G4ErrorMatrix *A);
// Does a QR decomposition of a matrix.

void back_solve(const G4ErrorMatrix &R, G4ErrorMatrix *b);
// Solves R*x = b where R is upper triangular.  Also has a variation that
// solves a number of equations of this form in one step, where b is a matrix
// with each column a different vector. See also solve.

void col_house(G4ErrorMatrix *a, const G4ErrorMatrix &v, double vnormsq,
	       int row, int col, int row_start, int col_start);
void col_house(G4ErrorMatrix *a, const G4ErrorMatrix &v, int row, int col,
	       int row_start, int col_start);
// Does a column Householder update.

void col_givens(G4ErrorMatrix *A, double c, double s,
		int k1, int k2, int row_min=1, int row_max=0);
// do a column Givens update

void row_givens(G4ErrorMatrix *A, double c, double s,
		int k1, int k2, int col_min=1, int col_max=0);
// do a row Givens update

void givens(double a, double b, double *c, double *s);
// algorithm 5.1.5 in Golub and Van Loan

// Returns a Householder vector to zero elements.

void house_with_update(G4ErrorMatrix *a, int row=1, int col=1);
void house_with_update(G4ErrorMatrix *a, G4ErrorMatrix *v, int row=1, int col=1);
// Finds and does Householder reflection on matrix.

void row_house(G4ErrorMatrix *a, const G4ErrorMatrix &v, double vnormsq,
	       int row, int col, int row_start, int col_start);
void row_house(G4ErrorMatrix *a, const G4ErrorMatrix &v, int row, int col,
	       int row_start, int col_start);
// Does a row Householder update.


#include "G4ErrorMatrix.icc"

#endif
