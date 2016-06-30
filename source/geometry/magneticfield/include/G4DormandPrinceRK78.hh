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
//  G4DormandPrinceRK78.hh
//
//  Created by Somnath on 30/06/15.

#ifndef G4Dormand_Prince_RK78_hh
#define G4Dormand_Prince_RK78_hh

#include "G4MagIntegratorStepper.hh"

class G4DormandPrinceRK78 : public G4MagIntegratorStepper
{
  public:
    G4DormandPrinceRK78(G4EquationOfMotion *EqRhs,
                     G4int numberOfVariables = 6,
                     G4bool primary =  true);
    ~G4DormandPrinceRK78();
    
    void Stepper( const G4double y[],
                 const G4double dydx[],
                 G4double h,
                 G4double yout[],
                 G4double yerr[]) ;

    G4double  DistChord()   const;
    G4int IntegratorOrder() const {return 7; }
    
 private :  
    G4DormandPrinceRK78(const G4DormandPrinceRK78&);
    G4DormandPrinceRK78& operator=(const G4DormandPrinceRK78&);
    
    G4double *ak2,   *ak3,  *ak4,  *ak5,  *ak6,  *ak7, *ak8,
             *ak9,   *ak10, *ak11, *ak12, *ak13,
             *yTemp, *yIn;

    G4double fLastStepLength;
    G4double *fLastInitialVector, *fLastFinalVector,
             *fLastDyDx,     *fMidVector,   *fMidError;
    // for DistChord calculations

    G4DormandPrinceRK78* fAuxStepper;
};

#endif /* defined(__G4_DormandPrinceRK78_hh) */
