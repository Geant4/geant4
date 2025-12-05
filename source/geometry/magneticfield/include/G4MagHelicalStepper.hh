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
// G4MagHelicalStepper
//
// Class description:
//
// Abstract base class for integrator of particle's equation of motion,
// used in tracking in space dependent magnetic field
//
// It is used for a set of steppers which use the helix as a sort of
// 'first order' solution.
//   - Most obtain an error by breaking up the step in two
//   - G4ExactHelicalStepper does not provide an error estimate

// Author: John Apostolakis (CERN), 05.11.1998
// --------------------------------------------------------------------
#ifndef G4MAGHELICALSTEPPER_HH
#define G4MAGHELICALSTEPPER_HH

#include <CLHEP/Units/PhysicalConstants.h>

#include "G4Types.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4Mag_EqRhs.hh"
#include "G4ThreeVector.hh"

/**
 * @brief G4MagHelicalStepper is an abstract base class for integrator of
 * particle's equation of motion, used in tracking in space dependent magnetic
 * field, and for a set of steppers which use the helix as 'first order'
 * solution.
 */

class G4MagHelicalStepper : public G4MagIntegratorStepper
{
  public:

    /**
     * Constructor for G4MagHelicalStepper.
     *  @param[in] EqRhs Pointer to the provided equation of motion.
     */
    G4MagHelicalStepper(G4Mag_EqRhs *EqRhs);

    /**
     * Default Destructor.
     */
    ~G4MagHelicalStepper() override = default;
  
    /**
     * Copy constructor and assignment operator not allowed.
     */
    G4MagHelicalStepper(const G4MagHelicalStepper&) = delete;
    G4MagHelicalStepper& operator=(const G4MagHelicalStepper&) = delete;
 
    /**
     * The stepper for the Runge Kutta integration.
     * The stepsize is fixed, with the step size given by 'h'.
     * Integrates ODE starting values y[0 to 6].
     * Outputs yout[] and its estimated error yerr[].
     *  @param[in] y Starting values array of integration variables.
     *  @param[in] dydx Derivatives array.
     *  @param[in] h The given step size.
     *  @param[out] yout Integration output.
     *  @param[out] yerr The estimated error.
     */
    void Stepper( const G4double y[], // VIRTUAL for ExactHelix
                  const G4double dydx[],
                        G4double h,
                        G4double yout[],
                        G4double yerr[] ) override;
  
    /**
     * Same as Stepper() function above, but should perform a 'dump' step
     * without error calculation. To be implemented in concrete derived classes.
     */
    virtual void DumbStepper( const G4double y[],
                                    G4ThreeVector Bfld,
                                    G4double h,
                                    G4double yout[] ) = 0;
  
    /**
     * Estimates the maximum distance of curved solution and chord.
     */
    G4double DistChord()const override ;

  protected:

    /**
     * Performs a linear Step in regions without magnetic field.
     */
    inline void LinearStep( const G4double yIn[],
                                  G4double h,
                                  G4double yHelix[]) const;

    /**
     * A first order Step along a helix inside the field.
     */
    void AdvanceHelix( const G4double yIn[],
                       const G4ThreeVector& Bfld,
                             G4double h,
                             G4double yHelix[], G4double yHelix2[] = nullptr);

    /**
     * Evaluates the field at a certain point.
     */
    inline void MagFieldEvaluate( const G4double y[], G4ThreeVector& Bfield );
  
    /**
     * Evaluates inverse of Curvature of Track.
     */
    inline G4double GetInverseCurve( const G4double Momentum,
                                     const G4double Bmag );

    // Store and use the parameters of track : 
    // radius of curve, Stepping angle, Radius of projected helix

    /**
     * Modifiers and accessors for storing and using the parameters of a track: 
     * radius of curve, Stepping angle, Radius of projected helix.
     */
    inline void SetAngCurve(const G4double Ang);
    inline G4double GetAngCurve()const;
    inline void SetCurve(const G4double Curve);
    inline G4double GetCurve()const;
    inline void SetRadHelix(const G4double Rad);
    inline G4double GetRadHelix()const;

  private:

    /** As in G4Mag_EqRhs.hh/cc where it is not used. */
    static const G4double fUnitConstant;

    G4Mag_EqRhs* fPtrMagEqOfMot = nullptr;
 
    /** Data stored in order to find the chord. */
    G4double fAngCurve = 0.0;
    G4double frCurve = 0.0;
    G4double frHelix = 0.0;
    G4ThreeVector yInitial, yMidPoint, yFinal;
};

#include  "G4MagHelicalStepper.icc"

#endif
