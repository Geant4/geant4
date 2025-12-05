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
// G4MonopoleEq
//
// Class description:
//
// This is the right-hand side of equation of motion for monopole
// in a combined electric and magnetic field:
// d(p_c)/ds=g{c-energyB_ - p_c x E}/pc.

// Author: Vladimir Grichine (CERN), 17.11.2009
// -------------------------------------------------------------------
#ifndef G4EQMAGELECTRICFIELD_HH
#define G4EQMAGELECTRICFIELD_HH

#include "G4ChargeState.hh"
#include "G4EquationOfMotion.hh"
#include "G4ElectroMagneticField.hh"

/**
 * @brief G4MonopoleEq defines the right-hand side of equation of motion
 * for monopole in a combined electric and magnetic field:
 * d(p_c)/ds=g{c-energyB_ - p_c x E}/pc.
 */

class G4MonopoleEq : public G4EquationOfMotion
{
  public:

    /**
     * Constructor for G4MonopoleEq.
     *  @param[in] emField Pointer to the field.
     */
    G4MonopoleEq(G4ElectroMagneticField* emField);

    /**
     * Default Destructor.
     */
    ~G4MonopoleEq() override = default;

    /**
     * Sets the charge, momentum and mass of the current particle.
     * Used to set the equation's coefficients.
     *  @param[in] particleCharge Magnetic charge and moments in e+ units.
     *  @param[in] MomentumXc Particle momentum.
     *  @param[in] mass Particle mass.
     */
    void  SetChargeMomentumMass(G4ChargeState particleCharge,
                                G4double MomentumXc,
                                G4double mass) override;

    /**
     * Calculates the value of the derivative, given the value of the field.
     *  @param[in] y Coefficients array.
     *  @param[in] Field Field value.
     *  @param[out] dydx Derivatives array.
     */
    void EvaluateRhsGivenB(const G4double y[],
                           const G4double Field[],
                                 G4double dydx[] ) const override;

    /**
     * Returns the equation type-ID, "kEqMonopole".
     */
    inline G4EquationType GetEquationType() const override { return kEqMonopole; }

  private:

    G4double fElectroMagCof = 0.0;
    G4double fMassCof = 0.0;
};

#endif
