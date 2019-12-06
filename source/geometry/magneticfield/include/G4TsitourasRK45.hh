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
// G4TsitourasRK45
//
// Class description:
//
// Tsitouras - 5(4) RK stepper (non-FSAL version)
// Implements RK tableau from 'Table 1' of:
//    C. Tsitouras, "Runge-Kutta pairs of order 5(4) satisfying only
//    the first column simplifying assumption"
//    Computers & Mathematics with Applications,
//    vol. 62, no. 2, pp. 770-775, 2011.

// Author: Somnath Banerjee, Google Summer of Code 2015, 11.06.2015
// Supervision: John Apostolakis, CERN
// --------------------------------------------------------------------
#ifndef TSITOURAS_RK45_HH
#define TSITOURAS_RK45_HH

#include "G4MagIntegratorStepper.hh"

class G4TsitourasRK45 : public G4MagIntegratorStepper
{
  public:

    G4TsitourasRK45(G4EquationOfMotion* EqRhs,
                    G4int numberOfVariables = 6,
                    G4bool primary =  true);
   ~G4TsitourasRK45();

    G4TsitourasRK45(const G4TsitourasRK45&) = delete;
    G4TsitourasRK45& operator=(const G4TsitourasRK45&) = delete;
    
    void Stepper( const G4double y[],
                  const G4double dydx[],
                        G4double h,
                        G4double yout[],
                        G4double yerr[] );
    
    void SetupInterpolation( /* const G4double yInput[],
                              const G4double dydx[],
                              const G4double Step */  );
    
    void Interpolate( const G4double yInput[],
                      const G4double dydx[],
                      const G4double Step,
                            G4double yOut[],
                            G4double tau );
      // For calculating the output at the tau fraction of Step

    void interpolate( const G4double yInput[],
                      const G4double dydx[],
                            G4double yOut[],
                            G4double Step,
                            G4double tau);

    G4double DistChord() const;
    inline G4int IntegratorOrder() const { return 4; }
    
  private :
    
    G4double *ak2, *ak3, *ak4, *ak5, *ak6, *ak7, *ak8, *yTemp, *yIn;
    
    G4double fLastStepLength = 0.0;
    G4double *fLastInitialVector, *fLastFinalVector,
             *fLastDyDx, *fMidVector, *fMidError;
      // For DistChord calculations
    
    G4TsitourasRK45* fAuxStepper = nullptr;
};

#endif
