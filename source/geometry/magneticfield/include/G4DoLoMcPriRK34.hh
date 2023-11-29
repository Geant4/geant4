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
// G4DoLoMcPriRK34
//
// Class description:
//
//  Dormand-Lockyer-McGorrigan-Prince-6-3-4 non-FSAL method
//  ( 6 stage, 3rd & 4th order embedded RK method )

// Created: Somnath Banerjee, Google Summer of Code 2015, 7 July 2015
// Supervision: John Apostolakis, CERN
// --------------------------------------------------------------------
#ifndef DOLO_MCPRI_RK34_HH
#define DOLO_MCPRI_RK34_HH

#include "G4MagIntegratorStepper.hh"

class G4DoLoMcPriRK34 : public G4MagIntegratorStepper
{
  public:

    G4DoLoMcPriRK34( G4EquationOfMotion* EqRhs,
                     G4int numberOfVariables = 6,
                     G4bool primary = true );
      // Constructor using Equation

    ~G4DoLoMcPriRK34() override;

    G4DoLoMcPriRK34(const G4DoLoMcPriRK34&) = delete;
    G4DoLoMcPriRK34& operator=(const G4DoLoMcPriRK34&) = delete; 
      // Copy constructor and assignment operator not allowed

    void Stepper( const G4double y[],
                  const G4double dydx[],
                        G4double h,
                        G4double yout[],
                        G4double yerr[] ) override ;
    
    void SetupInterpolation();
    void SetupInterpolate( const G4double yInput[],
                           const G4double dydx[],
                           const G4double Step );
      // For Preparing the interpolation and calculating the extra stages
    
    void Interpolate( const G4double yInput[],
                      const G4double dydx[],
                      const G4double Step,
                            G4double yOut[],
                            G4double tau );
      // For calculating the output at the tau fraction of Step

    void Interpolate( G4double tau,
                      G4double yOut[]);
    
    void interpolate(const G4double yInput[],
                     const G4double dydx[],
                           G4double yOut[],
                           G4double Step,
                           G4double tau ) ;

    G4double DistChord() const override;
    G4int IntegratorOrder() const override { return 3; }
    
  private :
    
    G4double *ak2, *ak3, *ak4, *ak5, *ak6, *yTemp, *yIn;
    
    G4double fLastStepLength = -1.0;
    G4double *fLastInitialVector, *fLastFinalVector,
             *fLastDyDx, *fMidVector, *fMidError;
      // for DistChord calculations
    
    G4DoLoMcPriRK34* fAuxStepper = nullptr;
};

#endif
