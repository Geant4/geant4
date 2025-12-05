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
// G4ModifiedMidpoint
//
// Class description:
//
// Modified midpoint method implementation, based on Boost odeint.

// Author: Dmitry Sorokin (CERN, Google Summer of Code 2016), 07.10.2016
// Supervision: John Apostolakis (CERN)
// --------------------------------------------------------------------
#ifndef G4MODIFIED_MIDPOINT_HH
#define G4MODIFIED_MIDPOINT_HH

#include "G4Types.hh"
#include "G4EquationOfMotion.hh"
#include "G4FieldTrack.hh"

/**
 * @brief G4ModifiedMidpoint implements a midpoint method adapted from
 * Boost odeint.
 */

class G4ModifiedMidpoint
{
  public:

    /**
     * Constructor for G4ModifiedMidpoint.
     *  @param[in] equation Pointer to the provided equation of motion.
     *  @param[in] nvar The number of integration variables.
     *  @param[in] steps The minimum number of steps.
     */
    G4ModifiedMidpoint( G4EquationOfMotion* equation,
                        G4int nvar = 6, G4int steps = 2 );

    /**
     * Default Destructor.
     */
    ~G4ModifiedMidpoint() = default;

    /**
     * Computes one step.
     *  @param[in] yIn Starting values array of integration variables.
     *  @param[in] dydxIn Derivatives array in input.
     *  @param[out] yOut Integration output.
     *  @param[in] hstep The given step size.
     */
    void DoStep( const G4double yIn[], const G4double dydxIn[],
                 G4double yOut[], G4double hstep) const;

    /**
     * Computes one step, as above but using also intermediate values.
     *  @param[in] yIn Starting values array of integration variables.
     *  @param[in] dydxIn Derivatives array in input.
     *  @param[out] yOut Integration output.
     *  @param[in] hstep The given step size.
     *  @param[in] yMid Mid point integration variables.
     *  @param[in] derivs Intermediate derivatives.
     */
    void DoStep( const G4double yIn[], const G4double dydxIn[],
                 G4double yOut[], G4double hstep, G4double yMid[],
                 G4double derivs[][G4FieldTrack::ncompSVEC]) const;

    /**
     * Setter and getter for steps.
     */
    inline void SetSteps(G4int steps);
    inline G4int GetSteps() const;

    /**
     * Setter and getter for the equation of motion.
     */
    inline void SetEquationOfMotion(G4EquationOfMotion* equation);
    inline G4EquationOfMotion* GetEquationOfMotion() const;

    /**
     * Returns the number of integration variables.
     */
    inline G4int GetNumberOfVariables() const;

  private:

    /**
     * Utility for copying array content from 'src' to 'dst'.
     */
    void copy(G4double dst[], const G4double src[]) const;

  private:

    G4EquationOfMotion* fEquation = nullptr;
    G4int fnvar = 0;
    G4int fsteps = 0;
};

#include "G4ModifiedMidpoint.icc"

#endif
