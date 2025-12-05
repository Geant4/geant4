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
// G4EquationOfMotion
//
// Class description:
//
// Abstract Base Class for the right hand size of the equation of
// motion of a particle in a field.

// Author: John Apostolakis (CERN), 1998
// -------------------------------------------------------------------
#ifndef G4EQUATIONOFMOTION_HH
#define G4EQUATIONOFMOTION_HH

#include "G4Types.hh"
#include "G4Field.hh"   // required in inline method implementations
#include "G4FieldParameters.hh"

#include "G4ChargeState.hh"

/**
 * @brief G4EquationOfMotion is the abstract base class for the right
 * hand size of the equation of motion of a particle in a field.
 */

class G4EquationOfMotion 
{
  public:

    /**
     * Constructor for G4EquationOfMotion.
     *  @param[in] Field Pointer to the field.
     */
    G4EquationOfMotion( G4Field* Field );

    /**
     * Default virtual Destructor.
     */
    virtual ~G4EquationOfMotion() = default;

    /**
     * Calculates the value of the derivative, given the value of the field.
     *  @param[in] y Coefficients array.
     *  @param[in] Field Field value.
     *  @param[out] dydx Derivatives array.
     */
    virtual void EvaluateRhsGivenB( const G4double y[],
                                    const G4double B[3],
                                          G4double dydx[] ) const = 0;

    /**
     * Sets the charge, momentum and mass of the current particle.
     * Used to set the equation's coefficients.
     *  @param[in] particleCharge Magnetic charge and moments in e+ units.
     *  @param[in] MomentumXc Particle momentum.
     *  @param[in] mass Particle mass.
     */
    virtual void SetChargeMomentumMass(G4ChargeState particleCharge,
                                       G4double MomentumXc,
                                       G4double MassXc2) = 0;

    /**
     * Returns the equation type-ID, "kUserEquation".
     */
    virtual G4EquationType GetEquationType() const { return kUserEquation; }

    /**
     * Calculates the value of the derivative 'dydx' at 'y'.
     * Calls the virtual function above.
     *  @param[in] y Coefficients array.
     *  @param[out] dydx Derivatives array.
     */
    inline void RightHandSide( const G4double y[],
                                     G4double dydx[] ) const;

    /**
     * Calculates the value of the derivative 'dydx' at 'y' as above,
     * but also returns the value of B.
     *  @param[in] y Coefficients array.
     *  @param[out] dydx Derivatives array.
     *  @param[out] Field Field value.
     */
    inline void EvaluateRhsReturnB( const G4double y[],
                                          G4double dydx[],
                                          G4double Field[] ) const;

    /**
     * Returns the 'Field' value at the given time 'Point'.
     *  @param[in] Point The time point (x,y,z,t).
     *  @param[out] Field The returned field value.
     */
    inline void GetFieldValue( const G4double Point[4],
                                     G4double Field[] ) const;

    /**
     * Accessors and modifier for the field.
     */
    inline const G4Field* GetFieldObj() const;
    inline G4Field* GetFieldObj();
    inline void SetFieldObj(G4Field* pField);

  private:

    G4Field* itsField = nullptr;
};

#include "G4EquationOfMotion.icc"

#endif
