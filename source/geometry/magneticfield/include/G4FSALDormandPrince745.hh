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
//  DormandPrince7 - 5(4) FSAL implementation by Somnath Banerjee
//  Supervision / code review: John Apostolakis
//
// Sponsored by Google in Google Summer of Code 2015.
//
// First version: 25 May 2015
//
//  $ID: FDormandPrince745.hh $
//
//  History
// ----------------------------------------------------
//  Created : 25 May 2015.              - Somnath
//
//  Added interpolate() method:         - Somnath
//      29 June 2015


#ifndef FFDormand_Prince_745
#define FFDormand_Prince_745

#include "G4VFSALIntegrationStepper.hh"

class G4FSALDormandPrince745 : public G4VFSALIntegrationStepper
{
  public:
	G4FSALDormandPrince745(G4EquationOfMotion *EqRhs,
					 G4int numberOfVariables = 6,
					 G4bool primary =  true);
	~G4FSALDormandPrince745();

	void Stepper( const G4double y[],
                  const G4double dydx[],
                        G4double h,
                        G4double yout[],
                        G4double yerr[],
                        G4double nextDydx[]) ;
    void interpolate( const G4double yInput[],
                     const G4double dydx[],
                     G4double yOut[],
                     G4double Step,
                     G4double tau ) ;
    //For higher order Interpolant
    void SetupInterpolate( const G4double yInput[],
                           const G4double dydx[],
                           const G4double Step );
    
    //For calculating the output at the tau fraction of Step
    void Interpolate( const G4double yInput[],
                      const G4double dydx[],
                      const G4double Step,
                      G4double yOut[],
                      G4double tau );
    

    G4double  DistChord()   const;
    G4int IntegratorOrder() const {return 4; }
     G4bool isFSAL() const{ return true; }
//    G4double *getLastDydx();
    
private :
    
    G4FSALDormandPrince745(const G4FSALDormandPrince745&);
    G4FSALDormandPrince745& operator=(const G4FSALDormandPrince745&);
    
    G4double *ak2, *ak3, *ak4, *ak5, *ak6, *ak7,
       *ak8, *ak9, 	//For additional stages in the interpolant
       *yTemp, *yIn;
    
    //Only for use with DistChord :-
    G4double *pseudoDydx_for_DistChord;
    
    G4double fLastStepLength;
    G4double *fLastInitialVector, *fLastFinalVector,
       *fInitialDyDx, *fLastDyDx, *fMidVector, *fMidError;
    // for DistChord calculations
    
    G4FSALDormandPrince745* fAuxStepper;
};

#endif /* defined(__Geant4__FSALDormandPrince745__) */
