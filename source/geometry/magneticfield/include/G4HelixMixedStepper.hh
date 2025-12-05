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
// G4HelixMixedStepper
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
// it possible to choose other stepper,like G4CashKarpRK45 or G4RKG3_Stepper,
// by setting StepperNumber : new HelixMixedStepper(EqRhs,N)
//
//  N=2  G4SimpleRunge;            N=3  G4SimpleHeum;
//  N=4  G4ClassicalRK4;      
//  N=6  G4HelixImplicitEuler;     N=7  G4HelixSimpleRunge;
//  N=8  G4CashKarpRK45;           N=9  G4ExactHelixStepper;
//  N=10 G4RKG3_Stepper;           N=13 G4NystromRK4
//  N=23 BogackiShampine23         N=145 TsitourasRK45 
//  N=45 BogackiShampine45         N=745 DormandPrince745 (ie DoPri5)
//
// For completeness also available are:
//  N=11 G4ExplicitEuler           N=12 G4ImplicitEuler;   -- Likely poor
//  N=5  G4HelixExplicitEuler (testing only)
//  For recommendations see comments in 'SetupStepper' method.
//
// Note: Like other helix steppers, only applicable in pure magnetic field.

// Author: Tatiana Nikitina (CERN), 18.05.2007
// -------------------------------------------------------------------
#ifndef G4HELIXMIXEDSTEPPER_HH
#define G4HELIXMIXEDSTEPPER_HH

#include "G4MagHelicalStepper.hh"

/**
 * @brief G4HelixMixedStepper is a concrete class for particle motion in
 * magnetic field which splits the method used for Integration in two:
 * if the stepping angle ( h / R_curve) is less than pi/3, use a RK stepper
 * for small step, else use G4HelixExplicitEuler stepper.
 * Like other helix steppers, it is only applicable in pure magnetic field.
 */

class G4HelixMixedStepper : public G4MagHelicalStepper
{
  public:  

    /**
     * Constructor for G4ExactHelixStepper.
     *  @param[in] EqRhs Pointer to the standard equation of motion.
     *  @param[in] StepperNumber Identified for selecting the stepper type;
     *             default (-1) is DormandPrince745.
     *  @param[in] Angle_threshold The stepping angle threshold; default (-1)
     *             is (1/3)*pi.
     */
    G4HelixMixedStepper(G4Mag_EqRhs* EqRhs,
                        G4int StepperNumber = -1,
                        G4double Angle_threshold = -1.0);

    /**
     * Default Destructor.
     */
    ~G4HelixMixedStepper() override;

    /**
     * The integration stepper. The stepsize is fixed, with the step size
     * given by 'hstep'. Integrates ODE starting values yInput[0 to 6].
     * Outputs yout[] and its estimated error yerr[].
     * If SteppingAngle = h/R_curve < pi/3, uses default RK stepper else
     * use Helix fast method.
     *  @param[in] y Starting values array of integration variables.
     *  @param[in] dydx Derivatives array.
     *  @param[in] h The given step size.
     *  @param[out] yout Integration output.
     *  @param[out] yerr The estimated error.
     */
    void Stepper( const G4double y[],
                  const G4double dydx[],
                        G4double h,
                        G4double yout[],
                        G4double yerr[] ) override;

    /**
     * Same as Stepper() function above, but should perform a 'dump' step
     * without error calculation. Assuming a constant field, the solution is
     * a helix.
     *  @param[in] y Starting values array of integration variables.
     *  @param[in] Bfld The field vector.
     *  @param[in] h The given step size.
     *  @param[out] yout Integration output.
     */
    void DumbStepper( const G4double y[],
                            G4ThreeVector Bfld,
                            G4double h,
                            G4double yout[] ) override;

    /**
     * Estimates the maximum distance of curved solution and chord.
     */
    G4double DistChord() const override;
    
    /**
     * Sets the verbosity level.
     */
    inline void SetVerbose (G4int newvalue) { fVerbose = newvalue; }
  
    /**
     * Setter and getter for the stepping angle threshold.
     */
    inline void SetAngleThreshold( G4double val ) { fAngle_threshold = val; }
    inline G4double GetAngleThreshold() { return fAngle_threshold; }

    /**
     * Returns the order, 4, of integration.
     */
    inline G4int IntegratorOrder() const override { return 4; }

    /**
     * Returns the stepper type-ID, "kHelixMixedStepper".
     */
    inline G4StepperType StepperType() const override { return kHelixMixedStepper; }

    /**
     * Logger function for the number of calls.
     */
    void PrintCalls();

    /**
     * Sets the chosen stepper and equation of motion.
     */
    G4MagIntegratorStepper* SetupStepper(G4Mag_EqRhs* EqRhs, G4int StepperName);

  private:

    /** Mixed Integration RK4 for 'small' steps. */
    G4MagIntegratorStepper* fRK4Stepper = nullptr;

    /** Int ID of Runge-Kutta stepper. */ 
    G4int fStepperNumber = -1;

    /** Threshold angle (in radians ); above it, the Helical stepper is used. */
    G4double fAngle_threshold = -1.0;

    /** Verbosity level. */ 
    G4int fVerbose = 0;

    /** Used for statistic, i.e. how many calls to different steppers. */
    G4int fNumCallsRK4 = 0;
    G4int fNumCallsHelix = 0;
};

#endif
