// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Integrator.hh,v 1.3 1999-11-16 17:30:58 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Class description:
//
// Class collecting template integrator methods for generic funtions.

// History:
//
// 04.09.99 V.Grichine, first implementation based on G4SimpleIntegration class
//                      and H.-P. Wellisch, G.Cosmo, and E.Cherniaev advises
// 08.09.99 V.Grichine, methods involving orthogonal polynomials
// 


#ifndef G4INTEGRATOR_HH
#define G4INTEGRATOR_HH 1

#include "globals.hh"
#include <math.h>

class G4Integrator
{
  public:

   G4Integrator(){;}

  ~G4Integrator(){;}

  // Simpson integration method

 template <class T, class F> G4double 
 Simpson( T& typeT, F f, G4double a, G4double b, G4int n ) ;

 template <class T, class F> G4double 
 Simpson( T* ptrT, F f, G4double a, G4double b, G4int n ) ;

 G4double Simpson( G4double (*f)(G4double), 
                   G4double a, G4double b, G4int n ) ;

  // Adaptive Gauss method

 template <class T, class F> G4double 
 AdaptiveGauss( T& typeT, F f, G4double a, G4double b, G4double e ) ;

 template <class T, class F> G4double 
 AdaptiveGauss( T* ptrT, F f, G4double a, G4double b, G4double e ) ;

 G4double AdaptiveGauss( G4double (*f)(G4double), 
                         G4double a, G4double b, G4double e ) ;

  // Integration methods involving orthogohol polynomials

  // Methods involving Legendre polynomials  


 template <class T, class F>    
 G4double Legendre( T& typeT, F f, G4double a, G4double b, G4int n) ;

 template <class T, class F>    
 G4double Legendre( T* ptrT, F f, G4double a, G4double b, G4int n) ;

 G4double Legendre( G4double (*f)(G4double), G4double a, G4double b, G4int n) ;

  // Legendre10 is very fast and accurate anough

 template <class T, class F>  
 G4double Legendre10( T& typeT, F f,G4double a, G4double b) ;

 template <class T, class F>  
 G4double Legendre10( T* ptrT, F f,G4double a, G4double b) ;

 G4double Legendre10( G4double (*f)(G4double), G4double a, G4double b) ;

  // Legendre96 is very accurate and fast enough

 template <class T, class F>  
 G4double Legendre96( T& typeT, F f,G4double a, G4double b) ;

 template <class T, class F>  
 G4double Legendre96( T* ptrT, F f,G4double a, G4double b) ;

 G4double Legendre96( G4double (*f)(G4double), G4double a, G4double b) ;

  // Methods involving Chebyshev  polynomials  
                          
 template <class T, class F>    
 G4double Chebyshev( T& typeT, F f, G4double a, G4double b, G4int n) ;

 template <class T, class F>    
 G4double Chebyshev( T* ptrT, F f, G4double a, G4double b, G4int n) ;

 G4double Chebyshev( G4double (*f)(G4double), G4double a, G4double b, G4int n) ;

  // Method involving Laguerre polynomials

 template <class T, class F>    
 G4double Laguerre( T& typeT, F f, G4double alpha, G4int n) ;

 template <class T, class F>    
 G4double Laguerre( T* ptrT, F f, G4double alpha, G4int n) ;

 G4double Laguerre( G4double (*f)(G4double), G4double alpha,  G4int n) ;

  // Method involving Hermite polynomials

 template <class T, class F>    
 G4double Hermite( T& typeT, F f, G4int n) ;

 template <class T, class F>    
 G4double Hermite( T* ptrT, F f, G4int n) ;

 G4double Hermite( G4double (*f)(G4double), G4int n) ;

  // Method involving Jacobi polynomials

 template <class T, class F>    
 G4double Jacobi( T& typeT, F f, G4double alpha, G4double beta, G4int n) ;

 template <class T, class F>    
 G4double Jacobi( T* ptrT, F f, G4double alpha, G4double beta, G4int n) ;

 G4double Jacobi( G4double (*f)(G4double), G4double alpha, 
                                           G4double beta, G4int n) ;


 protected:

  // Auxiliary function for adaptive Gauss method

 template <class T, class F> G4double 
 Gauss( T& typeT, F f, G4double a, G4double b ) ;

 template <class T, class F> G4double 
 Gauss( T* ptrT, F f, G4double a, G4double b ) ;

 G4double Gauss( G4double (*f)(G4double), G4double a, G4double b) ;

  //

 template <class T, class F> void
 AdaptGauss( T& typeT, F f, G4double a, G4double b, 
                            G4double e, G4double& sum, G4int& n) ;

 template <class T, class F> void 
 AdaptGauss( T* typeT, F f, G4double a, G4double b, 
                            G4double e, G4double& sum, G4int& n ) ;

 void AdaptGauss( G4double (*f)(G4double), G4double a, G4double b, 
                  G4double e, G4double& sum, G4int& n ) ;

 G4double GammaLogarithm(G4double xx) ;

 
} ;

#include "G4Integrator.icc"

#endif
