// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Integrator.hh,v 1.4 2000-06-15 17:31:23 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Class description:
//
// Template class collecting integrator methods for generic funtions.

// History:
//
// 04.09.99 V.Grichine, first implementation based on G4SimpleIntegration class
//                      H.P.Wellisch, G.Cosmo, and E.Cherniaev advises
// 08.09.99 V.Grichine, methods involving orthogonal polynomials
// 


#ifndef G4INTEGRATOR_HH
#define G4INTEGRATOR_HH 1

#include "globals.hh"
#include <math.h>

template <class T, class F>
class G4Integrator
{
  public:  // with description

   G4Integrator(){;}
  ~G4Integrator(){;}

 G4double Simpson( T& typeT, F f, G4double a, G4double b, G4int n ) ;
 G4double Simpson( T* ptrT, F f, G4double a, G4double b, G4int n ) ;
 G4double Simpson( G4double (*f)(G4double), 
                   G4double a, G4double b, G4int n ) ;
  // Simpson integration method

 G4double AdaptiveGauss( T& typeT, F f, G4double a, G4double b, G4double e ) ;
 G4double AdaptiveGauss( T* ptrT, F f, G4double a, G4double b, G4double e ) ;
 G4double AdaptiveGauss( G4double (*f)(G4double), 
                         G4double a, G4double b, G4double e ) ;
  // Adaptive Gauss method


  // Integration methods involving orthogohol polynomials

 G4double Legendre( T& typeT, F f, G4double a, G4double b, G4int n) ;
 G4double Legendre( T* ptrT, F f, G4double a, G4double b, G4int n) ;
 G4double Legendre( G4double (*f)(G4double), G4double a, G4double b, G4int n) ;
  //
  // Methods involving Legendre polynomials  

 G4double Legendre10( T& typeT, F f,G4double a, G4double b) ;
 G4double Legendre10( T* ptrT, F f,G4double a, G4double b) ;
 G4double Legendre10( G4double (*f)(G4double), G4double a, G4double b) ;
  //
  // Legendre10 is very fast and accurate enough

 G4double Legendre96( T& typeT, F f,G4double a, G4double b) ;
 G4double Legendre96( T* ptrT, F f,G4double a, G4double b) ;
 G4double Legendre96( G4double (*f)(G4double), G4double a, G4double b) ;
  //
  // Legendre96 is very accurate and fast enough

 G4double Chebyshev( T& typeT, F f, G4double a, G4double b, G4int n) ;
 G4double Chebyshev( T* ptrT, F f, G4double a, G4double b, G4int n) ;
 G4double Chebyshev( G4double (*f)(G4double), G4double a, G4double b, G4int n) ;
  //
  // Methods involving Chebyshev  polynomials  
                          
 G4double Laguerre( T& typeT, F f, G4double alpha, G4int n) ;
 G4double Laguerre( T* ptrT, F f, G4double alpha, G4int n) ;
 G4double Laguerre( G4double (*f)(G4double), G4double alpha,  G4int n) ;
  //
  // Method involving Laguerre polynomials

 G4double Hermite( T& typeT, F f, G4int n) ;
 G4double Hermite( T* ptrT, F f, G4int n) ;
 G4double Hermite( G4double (*f)(G4double), G4int n) ;
  //
  // Method involving Hermite polynomials

 G4double Jacobi( T& typeT, F f, G4double alpha, G4double beta, G4int n) ;
 G4double Jacobi( T* ptrT, F f, G4double alpha, G4double beta, G4int n) ;
 G4double Jacobi( G4double (*f)(G4double), G4double alpha, 
                                           G4double beta, G4int n) ;
  // Method involving Jacobi polynomials


 protected:

  // Auxiliary function for adaptive Gauss method

 G4double Gauss( T& typeT, F f, G4double a, G4double b ) ;
 G4double Gauss( T* ptrT, F f, G4double a, G4double b ) ;
 G4double Gauss( G4double (*f)(G4double), G4double a, G4double b) ;

 void AdaptGauss( T& typeT, F f, G4double a, G4double b, 
                                 G4double e, G4double& sum, G4int& n) ;
 void AdaptGauss( T* typeT, F f, G4double a, G4double b, 
                                 G4double e, G4double& sum, G4int& n ) ;
 void AdaptGauss( G4double (*f)(G4double), G4double a, G4double b, 
                  G4double e, G4double& sum, G4int& n ) ;

 G4double GammaLogarithm(G4double xx) ;

 
} ;

#include "G4Integrator.icc"

#endif
