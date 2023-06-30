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
// G4DormandPrinceRK56
//
// Class description:
//
// Dormand-Prince RK 6(5) non-FSAL method

// Created: Somnath Banerjee, Google Summer of Code 2015, 26 June 2015
// Supervision: John Apostolakis, CERN
// --------------------------------------------------------------------
#ifndef G4DORMAND_PRINCE_RK56_HH
#define G4DORMAND_PRINCE_RK56_HH

#include "G4MagIntegratorStepper.hh"

class G4DormandPrinceRK56 : public G4MagIntegratorStepper
{
  public:

    G4DormandPrinceRK56( G4EquationOfMotion* EqRhs,
                         G4int numberOfVariables = 6,
                         G4bool primary = true ) ;
    
    ~G4DormandPrinceRK56() override ;
    
    G4DormandPrinceRK56(const G4DormandPrinceRK56&) = delete;
    G4DormandPrinceRK56& operator=(const G4DormandPrinceRK56&) = delete;
    
    void Stepper( const G4double y[],
                  const G4double dydx[],
                        G4double h,
                        G4double yout[],
                        G4double yerr[] ) override ;
    
    G4double  DistChord()   const override;
    G4int IntegratorOrder() const override { return 5; }
    
    void SetupInterpolate_low( const G4double yInput[],
                               const G4double dydx[],
                               const G4double Step );
      // For preparing the Interpolant and calculating the extra stages
    
    void Interpolate_low( const G4double yInput[],
                          const G4double dydx[],
                          const G4double Step,
                                G4double yOut[],
                                G4double tau );
      // For calculating the output at the tau fraction of Step

    inline void SetupInterpolation()
    {
      SetupInterpolate( fLastInitialVector, fLastDyDx, fLastStepLength);
    }
   
    inline void SetupInterpolate( const G4double yInput[],
                                  const G4double dydx[],
                                  const G4double Step )
    {
      SetupInterpolate_low( yInput, dydx, Step);
    }
    
    inline void Interpolate( const G4double yInput[],
                             const G4double dydx[],
                             const G4double Step,
                                   G4double yOut[],
                                   G4double tau )
    {
      Interpolate_low( yInput, dydx, Step, yOut, tau);
    }
      // For calculating the output at the tau fraction of Step

    inline void Interpolate( G4double tau, G4double yOut[])
    {
      Interpolate( fLastInitialVector, fLastDyDx, fLastStepLength, yOut, tau );
    }

    void SetupInterpolate_high( const G4double yInput[],
                                const G4double dydx[],
                                const G4double Step );
    
    void Interpolate_high( const G4double yInput[],
                           const G4double dydx[],
                           const G4double Step,
                                 G4double yOut[],
                                 G4double tau );
      // For calculating the output at the tau fraction of Step
    
  private:

    G4double *ak2, *ak3, *ak4, *ak5, *ak6, *ak7, *ak8, *ak9;
      // For storing intermediate 'k' values in stepper
    G4double *ak10_low, *ak10, *ak11, * ak12;
      // For the additional stages of Interpolant
    G4double *yTemp, *yIn;
    
    G4double fLastStepLength = -1.0;
    G4double *fLastInitialVector, *fLastFinalVector,
             *fLastDyDx, *fMidVector, *fMidError;
      // For DistChord calculations
    
    G4DormandPrinceRK56* fAuxStepper = nullptr;
};

#endif /* G4DormandPrinceRK56 */
