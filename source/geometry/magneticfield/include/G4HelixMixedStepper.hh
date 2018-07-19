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
// class G4HelixMixedStepper
//
// Class description:
//
// G4HelixMixedStepper split the Method used for Integration in two:
//
// If Stepping Angle ( h / R_curve) < pi/3 : use Stepper for small step
// 
// Else use  HelixExplicitEuler Stepper
//
// Stepper for the small step is G4ClassicalRK4 by default, but
//  it possible to choose other stepper,like G4CashKarpRK45 or G4RKG3_Stepper,
//  by setting StepperNumber : new HelixMixedStepper(EqRhs,N)
//
//  N=2  G4SimpleRunge;            N=3  G4SimpleHeum;
//  N=4  G4ClassicalRK4;      
//  N=6  G4HelixImplicitEuler;     N=7  G4HelixSimpleRunge;
//  N=8  G4CashKarpRK45;           N=9  G4ExactHelixStepper;
//  N=10 G4RKG3_Stepper;           N=13 G4NystromRK4
//  N=23 BogackiShampine23         N=145 TsitourasRK45 
//  N=45 BogackiShampine45         N=745 DormandPrince745 (ie DoPri5)
//
//  For completeness also available are:
//  N=11 G4ExplicitEuler           N=12 G4ImplicitEuler;   -- Likely poor
//  N=5  G4HelixExplicitEuler (testing only)
//  For recommendations see comments in 'SetupStepper' method.
//
//  Note: Like other helix steppers, only applicable in pure magnetic field
//
// History: 
// Derived from ExactHelicalStepper 18/05/07 
//
// -------------------------------------------------------------------

#ifndef G4HELIXMIXEDSTEPPER_HH
#define G4HELIXMIXEDSTEPPER_HH

#include "G4MagHelicalStepper.hh"


class G4HelixMixedStepper: public G4MagHelicalStepper
{

  public:  

  G4HelixMixedStepper(G4Mag_EqRhs *EqRhs,G4int StepperNumber= -1, G4double Angle_threshold= -1.0);
  ~G4HelixMixedStepper();

   void Stepper( const G4double y[],
                  const G4double dydx[],
                        G4double h,
                        G4double yout[],
                        G4double yerr[]  );
      // Step 'integration' for step size 'h'
      // If SteppingAngle=h/R_curve<pi/3 uses RK4Stepper
      // Else Helix Fast Method 
      
  
    void DumbStepper( const G4double y[],
                            G4ThreeVector  Bfld,
                            G4double       h,
                            G4double       yout[]);
   G4double DistChord() const;
      // Estimate maximum distance of curved solution and chord ... 
    

   public:  // with description

    inline void SetVerbose (G4int newvalue){fVerbose=newvalue;}
  
  public:  // without description
    void PrintCalls();
    G4MagIntegratorStepper* SetupStepper(G4Mag_EqRhs* EqRhs, G4int StepperName);

    void     SetAngleThreshold( G4double val ){ fAngle_threshold= val;}
    G4double GetAngleThreshold(){ return fAngle_threshold; }
  
    G4int IntegratorOrder() const { return 4; }
  private:

      // Mixed Integration RK4 for 'small' steps
      G4MagIntegratorStepper* fRK4Stepper;
      G4int  fStepperNumber; //  Int ID of RK stepper 
   
      // Threshold angle (in radians ) - above it Helical stepper is used
        G4double                fAngle_threshold;
   private:
    // Used for statistic = how many calls to different steppers
       G4int fVerbose;
       G4int fNumCallsRK4;
       G4int fNumCallsHelix;
    
};

#endif /* G4HELIXMIXEDSTEPPER_HH */
