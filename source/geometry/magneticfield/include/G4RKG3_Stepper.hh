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
// G4RKG3_Stepper
//
// Class description:
//
// Runga-Kutta integrator stepper from Geant-3.

// Authors: John Apostolakis & Vladimir Grichine (CERN), 30.01.1997
// -------------------------------------------------------------------
#ifndef G4RKG3_STEPPER_HH
#define G4RKG3_STEPPER_HH

#include "G4Types.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ThreeVector.hh"

class G4Mag_EqRhs;

/**
 * @brief G4RKG3_Stepper implements a Runga-Kutta integrator stepper
 * used in Geant-3.
 */

class G4RKG3_Stepper : public G4MagIntegratorStepper
{
  public:

    /**
     * Constructor for G4RKG3_Stepper.
     *  @param[in] EqRhs Pointer to the provided equation of motion.
     */
    G4RKG3_Stepper(G4Mag_EqRhs* EqRhs);

    /**
     * Default Destructor.
     */
    ~G4RKG3_Stepper() override = default;

    /**
     * Copy constructor and assignment operator not allowed.
     */
    G4RKG3_Stepper(const G4RKG3_Stepper&) = delete;
    G4RKG3_Stepper& operator= (const G4RKG3_Stepper&) = delete;

    /**
     * The stepper for the Runge Kutta integration.
     * Method provided, even if less efficient.
     * The stepsize is fixed, with the step size given by 'h'.
     * Integrates ODE starting values yInput[0 to 6].
     * Outputs yOut[] and its estimated error yErr[].
     *  @param[in] yIn Starting values array of integration variables.
     *  @param[in] dydx Derivatives array.
     *  @param[in] h The given step size.
     *  @param[out] yOut Integration output.
     *  @param[out] yErr The estimated error.
     */
    void Stepper( const G4double yIn[],
                  const G4double dydx[],
                        G4double h,
                        G4double yOut[],
                        G4double yErr[] ) override;

    /**
     * Returns the distance from chord line.
     */
    G4double DistChord() const override ;
 
    /**
     * Integrator of Runge-Kutta Stepper from Geant-3 with only two field
     * evaluation per Step. It is used in propagating the initial Step
     * by small substeps after solution error and delta geometry
     * considerations. B[3] is magnetic field which is passed from substep
     * to substep.
     */
    void StepNoErr( const G4double tIn[8],
                    const G4double dydx[6],
                          G4double Step,
                          G4double tOut[8],
                          G4double B[3] );
    /**
     * Integrator for Runge-Kutta from Geant-3 with evaluation of error in
     * solution and delta geometry based on naive similarity with the case
     * of uniform magnetic field.
     * B1[3] is in input and is the first magnetic field values
     * B2[3] is the output and is the final magnetic field values.
     */
    void StepWithEst( const G4double tIn[8],
                      const G4double dydx[6],
                            G4double Step,
                            G4double tOut[8],
                            G4double& alpha2,
                            G4double& beta2,
                      const G4double B1[3],
                            G4double B2[3] );

    /**
     * Returns the order, 4, of integration.
     */
    inline G4int IntegratorOrder() const override { return 4; }

    /**
     * Returns the stepper type-ID, "kRKG3Stepper".
     */
    inline G4StepperType StepperType() const override { return kRKG3Stepper; }

  private:

    G4ThreeVector fyInitial,
                  fyMidPoint,
                  fyFinal;
    G4ThreeVector fpInitial;
    G4ThreeVector BfldIn;
    G4double      hStep = 0.0;
};

#endif
