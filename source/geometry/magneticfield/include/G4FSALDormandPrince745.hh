//  DormandPrince7 - 5(4) FSAL implementation by Somnath Banerjee
//  Supervision / code review: John Apostolakis
//
// Sponsored by Google in Google Summer of Code 2015.
//
// First version: 25 May 2015
//
// This code is made available subject to the Geant4 license, a copy of
// which is available at
//   http://geant4.org/license
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
    *fLastDyDx, *fMidVector, *fMidError;
    // for DistChord calculations
    
    G4FSALDormandPrince745* fAuxStepper;
};

#endif /* defined(__Geant4__FSALDormandPrince745__) */
