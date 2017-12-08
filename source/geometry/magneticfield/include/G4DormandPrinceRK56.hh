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
//  Dormand-Prince RK 6(5) non-FSAL implementation by Somnath Banerjee
//  Supervision / code review: John Apostolakis
//
// Sponsored by Google in Google Summer of Code 2015.
// 
// First version: 26 June 2015
///////////////////////////////////////////////////////////////////////////////

#ifndef DORMAND_PRINCE_RK56_H
#define DORMAND_PRINCE_RK56_H

#include "G4MagIntegratorStepper.hh"

class G4DormandPrinceRK56 : public G4MagIntegratorStepper
{
    
public:
    //constructor
    G4DormandPrinceRK56( G4EquationOfMotion *EqRhs,
               G4int numberOfVariables = 6,
               G4bool primary= true ) ;
    
    //destructor
    ~G4DormandPrinceRK56() ;
    
    //Stepper
    void Stepper( const G4double y[],
                  const G4double dydx[],
                 		G4double h,
                 		G4double yout[],
                 		G4double yerr[] ) ;
    
    G4double  DistChord()   const;
    G4int IntegratorOrder() const { return 5; }
    
    G4DormandPrinceRK56(const G4DormandPrinceRK56&);
    G4DormandPrinceRK56& operator=(const G4DormandPrinceRK56&);
    
    //For Preparing the Interpolant and calculating the extra stages
    void SetupInterpolate_low( const G4double yInput[],
                              const G4double dydx[],
                              const G4double Step );
    
    //For calculating the output at the tau fraction of Step
    void Interpolate_low( const G4double yInput[],
                          const G4double dydx[],
                          const G4double Step,
                                G4double yOut[],
                                G4double tau );

    void SetupInterpolation()
      { SetupInterpolate( fLastInitialVector, fLastDyDx, fLastStepLength); }
   
    inline void SetupInterpolate( const G4double yInput[],
                                  const G4double dydx[],
                                  const G4double Step ){
        SetupInterpolate_low( yInput, dydx, Step);
    }
    
    //For calculating the output at the tau fraction of Step
    inline void Interpolate( const G4double yInput[],
                             const G4double dydx[],
                             const G4double Step,
                                   G4double yOut[],
                                   G4double tau ){
        Interpolate_low( yInput, dydx, Step, yOut, tau);
    }

    void Interpolate( G4double tau, G4double yOut[])
    { Interpolate( fLastInitialVector, fLastDyDx, fLastStepLength, yOut, tau ); }

    void SetupInterpolate_high( const G4double yInput[],
                                const G4double dydx[],
                                const G4double Step );
    
    //For calculating the output at the tau fraction of Step
    void Interpolate_high( const G4double yInput[],
                           const G4double dydx[],
                           const G4double Step,
                                 G4double yOut[],
                                 G4double tau );
    
private:
    G4double *ak2, *ak3, *ak4, *ak5, *ak6, *ak7, *ak8, *ak9; // for storing intermediate 'k' values in stepper
    G4double *ak10_low, *ak10, *ak11, * ak12;    //For the additional stages of Interpolant
    G4double *yTemp, *yIn;
    
    G4double fLastStepLength;
    G4double *fLastInitialVector, *fLastFinalVector,
             *fLastDyDx, *fMidVector, *fMidError;
      // for DistChord calculations
    
    G4DormandPrinceRK56* fAuxStepper;
};

#endif /* G4DormandPrinceRK56 */

