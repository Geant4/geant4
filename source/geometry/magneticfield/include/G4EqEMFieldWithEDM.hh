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
// G4EqEMFieldWithEDM
//
// Class description:
//
// This is the right-hand side of equation of motion in a combined
// electric and magnetic field, with spin tracking for both MDM and
// EDM terms.

// Author: Kevin Lynch (Boston Univ.), 19.02.2009 - Based on G4EqEMFieldWithSpin
// -------------------------------------------------------------------
#ifndef G4EQEMFIELDWITHEDM_HH
#define G4EQEMFIELDWITHEDM_HH

#include "G4ChargeState.hh"
#include "G4EquationOfMotion.hh"

class G4ElectroMagneticField;

/**
 * @brief G4EqEMFieldWithEDM implements the right-hand side of equation of
 * motion in a combined electric and magnetic field, with spin tracking for
 * both MDM and EDM terms.
 */

class G4EqEMFieldWithEDM : public G4EquationOfMotion
{
  public:

    /**
     * Constructor for G4EqEMFieldWithEDM.
     *  @param[in] emField Pointer to the electromagnetic field.
     */
    G4EqEMFieldWithEDM(G4ElectroMagneticField* emField);

    /**
     * Default Destructor.
     */
    ~G4EqEMFieldWithEDM() override = default;
  
    /**
     * Sets the charge, momentum and mass of the current particle.
     * Used to set the equation's coefficients.
     *  @param[in] particleCharge Magnetic charge and moments in e+ units.
     *  @param[in] MomentumXc Particle momentum.
     *  @param[in] mass Particle mass.
     */
    void SetChargeMomentumMass(G4ChargeState particleCharge, // in e+ units
                               G4double MomentumXc,
                               G4double mass) override;

    /**
     * Calculates the value of the derivative, given the value of the
     * electromagnetic field.
     *  @param[in] y Coefficients array.
     *  @param[in] Field Field value.
     *  @param[out] dydx Derivatives array.
     */
    void EvaluateRhsGivenB(const G4double y[],
                           const G4double Field[],
                                 G4double dydx[] ) const override;

    /**
     * Setter and getter for magnetic anomaly.
     */
    inline void SetAnomaly(G4double a) { anomaly = a; }
    inline G4double GetAnomaly() const { return anomaly; }

    /**
     * Setter and getter for EDM eta parameter.
     */
    inline void SetEta(G4double n) { eta = n; }
    inline G4double GetEta() const { return eta; }

    /**
     * Returns the equation type-ID, "kEqEMfieldWithEDM".
     */
    inline G4EquationType GetEquationType() const override { return kEqEMfieldWithEDM; }

  private:

    G4double charge{0.}, mass{0.}, magMoment{0.}, spin{0.};

    G4double fElectroMagCof{0.} ;
    G4double fMassCof{0.};

    G4double omegac{0.}, anomaly{0.0011659208}, eta{0.};
    G4double beta{0.}, gamma{0.};
};

#endif
