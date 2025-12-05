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
// G4BorisScheme
//
// Class description:
//
// Implementation of the Boris algorithm for advancing 
// charged particles in an electromagnetic field.

// Author: Divyansh Tiwari (CERN, Google Summer of Code 2022), 05.11.2022
// Supervision: John Apostolakis (CERN), Renee Fatemi, Soon Yung Jun (FNAL)
// --------------------------------------------------------------------
#ifndef G4BORIS_SCHEME_HH
#define G4BORIS_SCHEME_HH

#include "G4Types.hh"

#include <CLHEP/Units/PhysicalConstants.h>

class G4EquationOfMotion;

/**
 * @brief The G4BorisScheme class implements of the Boris algorithm for
 * advancing charged particles in an electromagnetic field.
 */

class G4BorisScheme
{
  public:

    /**
     * Default Constructor.
     */
    G4BorisScheme() = default;

    /**
     * Constructor for the equation of motion.
     *  @param[in] equation Pointer to the equation of motion algorithm.
     *  @param[in] nvar The number of integration variables.
     */
    G4BorisScheme( G4EquationOfMotion* equation, G4int nvar = 6 );

    /**
     * Default Destructor.
     */
    ~G4BorisScheme() = default;

    /**
     * Does one step, updating velocity and position.
     *  @param[in] restMass Particle mass.
     *  @param[in] charge Particle charge.
     *  @param[in] yIn Initial position.
     *  @param[out] yOut Updated position.
     *  @param[in] hstep Proposed step.
     */
    void DoStep(G4double restMass, G4double charge, const G4double yIn[], 
                G4double yOut[], G4double hstep) const;

    /**
     * Adopts the Boris Scheme Stepping to estimate the integration error.
     * Uses two half-steps (comparing to a full step) to obtain output and
     * error estimate.
     *  @param[in] yIn Initial position.
     *  @param[in] restMass Particle mass.
     *  @param[in] charge Particle charge.
     *  @param[in] hstep Proposed step.
     *  @param[out] yOut Updated position.
     *  @param[out] yErr The estimated error.
     */
    void StepWithErrorEstimate(const G4double yIn[], G4double restMass,
                               G4double charge, G4double hstep,
                               G4double yOut[], G4double yErr[]) const;

    /**
     * Adopts the Boris Scheme Stepping to estimate the integration error.
     * Uses two half-steps (comparing to a full step) to obtain output and
     * error estimate. Same as above, but also returns the mid-point evaluation.
     *  @param[in] yIn Initial position.
     *  @param[in] restMass Particle mass.
     *  @param[in] charge Particle charge.
     *  @param[in] hstep Proposed step.
     *  @param[out] yMid tThe mid-point evaluation.
     *  @param[out] yOut Updated position.
     *  @param[out] yErr The estimated error.
     */
    void StepWithMidAndErrorEstimate(const G4double yIn[], G4double restMass,
                                     G4double charge, G4double hstep,
                    G4double yMid[], G4double yOut[], G4double yErr[]) const;

    /**
     * Auxiliary methods returning a pointer to the equation of motion
     * and the number of integration variables.
     */
    inline G4EquationOfMotion* GetEquationOfMotion() const;
    inline G4int GetNumberOfVariables() const;
   
  private:

    /**
     * Internal methods for updating position and velocity, used in DoStep().
     */
    void UpdatePosition(const G4double restMass, const G4double charge,
                        const G4double yIn[], G4double yOut[], G4double hstep) const;
    void UpdateVelocity(const G4double restMass, const G4double charge,
                        const G4double yIn[], G4double yOut[], G4double hstep) const;

    /**
     * Utility to mem-copy 'src' array data to 'dst'.
     */
    void copy(G4double dst[], const G4double src[]) const;

  private:

    G4EquationOfMotion* fEquation = nullptr;
    G4int fnvar = 8;
    static constexpr G4double c_l = CLHEP::c_light/CLHEP::m*CLHEP::second;
};

#include "G4BorisScheme.icc"

#endif
