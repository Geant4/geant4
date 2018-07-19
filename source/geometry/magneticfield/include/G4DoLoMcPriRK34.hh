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
//  DoLoMcPri4(3) RK method - header
//  Implements the 6-3-4 non-FSAL method   ( 6 stage, 3rd & 4th order embedded RK method )
//
//  RK tableau / method develoed by 
//    Dormand Lockyer McGorrigan Prince
//     RK4(3)6FD - forced Non-FSAL
//
//  Header, design, implementation by Somnath Banerjee
//     Supported by Google in Google Summer of Code 2015.
//  Supervision / code review: John Apostolakis
//
// First version: 7 July 2015
//
// This code is made available subject to the Geant4 license, a copy of
// which is available at
//   http://geant4.org/license
//  G4DoLoMcPriRK34.hh
//  Geant4
//
//  History
// -----------------------------
//  Created by Somnath on 7 July 2015


#ifndef DoLo_McPri_34
#define DoLo_McPri_34

#include "G4MagIntegratorStepper.hh"

class G4DoLoMcPriRK34 : public G4MagIntegratorStepper
{
public:
    //Constructor using Equation
	G4DoLoMcPriRK34(G4EquationOfMotion *EqRhs,
					 G4int numberOfVariables = 6,
					 G4bool primary =  true);
	~G4DoLoMcPriRK34();

    void Stepper( const G4double y[],
                  const G4double dydx[],
                        G4double h,
                        G4double yout[],
                        G4double yerr[] ) ;
    
    //For Preparing the Interpolant and calculating the extra stages
    void SetupInterpolation();
    void SetupInterpolate( const G4double yInput[],
                              const G4double dydx[],
                              const G4double Step );
    
    //For calculating the output at the tau fraction of Step
    void Interpolate( const G4double yInput[],
                      const G4double dydx[],
                      const G4double Step,
                            G4double yOut[],
                            G4double tau );

    void Interpolate( G4double tau,
                      G4double yOut[]);
    
    void interpolate(const G4double yInput[],
                     const G4double dydx[],
                           G4double yOut[],
                           G4double Step,
                           G4double tau ) ;

    G4double  DistChord()   const;
    G4int IntegratorOrder() const {return 3; }
    
private :
    
    G4DoLoMcPriRK34(const G4DoLoMcPriRK34&);
        //Copy constructor kept private
    G4DoLoMcPriRK34& operator=(const G4DoLoMcPriRK34&); 
        //assignment operator overloaded
    
    G4double *ak2, *ak3, *ak4, *ak5, *ak6, *yTemp, *yIn;
    
    G4double fLastStepLength;
    G4double *fLastInitialVector, *fLastFinalVector,
    *fLastDyDx, *fMidVector, *fMidError;
    // for DistChord calculations
    
    G4DoLoMcPriRK34* fAuxStepper;
};

#endif /* defined(__Geant4__G4DoLoMcPriRK34__) */
