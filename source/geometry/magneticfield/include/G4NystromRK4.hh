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
// G4NystromRK4
//
// Class description:
//
// Integrates the equations of the motion of a particle in a magnetic field
// using 4th Runge-Kutta-Nystrom method with errors estimation 
// (ATL-SOFT-PUB-2009-01)
// Current form can be used only for 'pure' magnetic field.
// Notes: 1) field must be time-independent.
//        2) time is not integrated

// Author: Igor Gavrilenko (CERN), 15.05.2009 (as G4AtlasRK4)
// Adaptations: John Apostolakis (CERN), 05.11.2009
// -------------------------------------------------------------------
#ifndef G4NYSTROMRK4_HH
#define G4NYSTROMRK4_HH

#include "G4MagIntegratorStepper.hh"
#include "G4Mag_EqRhs.hh"
#include "G4CachedMagneticField.hh"
#include "G4ThreeVector.hh"

#include <memory>

/**
 * @brief G4NystromRK4 integrates the equations of the motion of a particle
 * in a magnetic field using 4th Runge-Kutta-Nystrom method with errors
 * estimation. The current form can be used only for 'pure' magnetic field.
 */

class G4NystromRK4 : public G4MagIntegratorStepper
{
  public: 

    /**
     * Constructor for G4NystromRK4. Can be used only for Magnetic Fields
     * and for 6 variables (x,p).
     *  @param[in] EquationMotion Pointer to the provided equation of motion.
     *  @param[in] distanceConstField Distance value for constant field.
     */
    G4NystromRK4(G4Mag_EqRhs* EquationMotion, 
                 G4double distanceConstField = 0.0); 

    /**
     * Default Destructor.
     */
   ~G4NystromRK4() override = default;
   
    /**
     * The stepper for the Runge Kutta integration.
     * The stepsize is fixed, with the step size given by 'hstep'.
     * Integrates ODE starting values y[0 to 6].
     * Outputs yOut[] and its estimated error yError[].
     * Provides error via analytical method.
     *  @param[in] y Starting values array of integration variables.
     *  @param[in] dydx Derivatives array.
     *  @param[in] hstep The given step size.
     *  @param[out] yOut Integration output.
     *  @param[out] yError The estimated error.
     */
    void Stepper(const G4double y[],
                 const G4double dydx[],
                       G4double hstep,
                       G4double yOut[],
                       G4double yError[]) override;

    /**
     * Setter and getter for the distance value for constant field.
     */
    void SetDistanceForConstantField(G4double length); 
    G4double GetDistanceForConstantField() const; 
   
    /**
     * Returns the order, 4, of integration.
     */
    inline G4int IntegratorOrder() const override;

    /**
     * Returns the distance from chord line.
     */
    G4double DistChord() const override; 

    /**
     * Returns the stepper type-ID, "kNystromRK4".
     */
    inline G4StepperType StepperType() const override;
  
  private:

    /**
     * Private accessors for field data.
     */
    inline void GetFieldValue(const G4double point[4], G4double field[3]);
    inline G4double GetFCof();
    G4CachedMagneticField* GetField();
    const G4CachedMagneticField* GetField() const;

  private:

    G4double fMomentum = 0.0;
    G4double fMomentum2 = 0.0;
    G4double fInverseMomentum = 0.0;
    G4double fCoefficient = 0.0;
    G4ThreeVector fInitialPoint;
    G4ThreeVector fMidPoint;
    G4ThreeVector fEndPoint;

    std::unique_ptr<G4CachedMagneticField> fCachedField;
};

#include "G4NystromRK4.icc"

#endif
