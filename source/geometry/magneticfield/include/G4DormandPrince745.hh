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
// $Id: G4DormandPrince745.hh 97387 2016-06-02 10:03:42Z gcosmo $
//
//  Class desription: 
//    An implementation of the 5th order embedded RK method from the paper
//    J. R. Dormand and P. J. Prince, “A family of embedded Runge-Kutta formulae,”
//	    Journal of computational and applied …, vol. 6, no. 1, pp. 19–26, 1980.
//
//  DormandPrince7 - 5(4) embedded RK method
//
//  Design & Implementation by Somnath Banerjee
//  Supervision & code review: John Apostolakis
//
// Work supported by the Google Summer of Code 2015.
//
//  History
// ------------------------------------------
//  Created   : 25 May 2015.             - Somnath
//   Revisions : 
//   * 29 June 2015:  Added interpolate() method(s) - Somnath
//   *     May 2016:  Cleanup and first comming in G4 - John Apostolakis

#ifndef Dormand_Prince_745
#define Dormand_Prince_745

#include "G4MagIntegratorStepper.hh"

class G4DormandPrince745 : public G4MagIntegratorStepper
{
  public:
    G4DormandPrince745(G4EquationOfMotion *EqRhs,
					 G4int numberOfVariables = 6,
					 G4bool primary =  true);
    ~G4DormandPrince745();
   
    void Stepper( const G4double y[],
                  const G4double dydx[],
                        G4double h,
                        G4double yout[],
                        G4double yerr[] ) ;

    //For Preparing the Interpolant and calculating the extra stages
    void SetupInterpolation_low( /* const G4double yInput[],
                                    const G4double dydx[],
                                    const G4double Step */ );
    
    //For calculating the output at the tau fraction of Step
    void Interpolate_low( /* const G4double yInput[],
                             const G4double dydx[],
                             const G4double Step, */ 
                          G4double yOut[],
                          G4double tau );
    
    inline void SetupInterpolation()
                            /* ( const G4double yInput[],
                                 const G4double dydx[],
                                 const G4double Step ) */
    { 
       SetupInterpolation_low( /* yInput, dydx, Step */ );
       // SetupInterpolation_high( /* yInput, dydx, Step */ );       
    }
    
    //For calculating the output at the tau fraction of Step
    inline void Interpolate(
                         /* const G4double yInput[],
                            const G4double dydx[],
                            const G4double Step,  */
                                  G4double tau,    
                                  G4double yOut[]
       )
    {
        Interpolate_low(  /* yInput, dydx, Step, */  yOut, tau);
        // Interpolate_high(  /* yInput, dydx, Step, */  yOut, tau);        
    }
    
    void SetupInterpolation_high( /* const G4double yInput[],
                               const G4double dydx[],
                               const G4double Step */ );
    
    //For calculating the output at the tau fraction of Step
    void Interpolate_high( /* const G4double yInput[],
                              const G4double dydx[],
                              const G4double Step, */ 
                                 G4double yOut[],
                                 G4double tau );

    G4double  DistChord() const;
    G4double DistChord2() const;
    G4double DistChord3() const;
   
    //  Enabling method, with common code between implementations (and steppers)
    G4double DistLine( G4double yStart[], G4double yMid[], G4double yEnd[] ) const;   
    G4int IntegratorOrder() const {return 4; }
    
    //New copy constructor
    //  G4DormandPrince745(const G4DormandPrince745 &);
    
private :
    
    G4DormandPrince745& operator=(const G4DormandPrince745&);
    
    G4double *ak2, *ak3, *ak4, *ak5, *ak6, *ak7,
      *ak8, *ak9, 	//For additional stages in the interpolant
      *yTemp, *yIn;
    
    G4double fLastStepLength;
    G4double *fLastInitialVector, *fLastFinalVector,
             *fInitialDyDx, *fMidVector, *fMidError;
    // for DistChord calculations
    
    G4DormandPrince745* fAuxStepper;
};
#endif /* defined(__Geant4__G4DormandPrince745__) */
