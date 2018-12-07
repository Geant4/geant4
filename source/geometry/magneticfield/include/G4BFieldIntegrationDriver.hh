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
//
//
//
// class G4BFieldIntegrationDriver
//
// Class description:
//
// Specialized integration driver for pure magnetic field

// History:
// - Created. D.Sorokin
// --------------------------------------------------------------------

#ifndef G4BFieldIntegrationDriver_HH
#define G4BFieldIntegrationDriver_HH

#include "G4IntegrationDriver.hh"
#include "G4HelixExplicitEuler.hh"

template <class T>
class G4BFieldIntegrationDriver : public G4IntegrationDriver<T> {
public:
    G4BFieldIntegrationDriver(  G4double hminimum,
                                T*       stepper,
                                G4int    numberOfComponents = 6,
                                G4int    statisticsVerbosity = 1);

    G4BFieldIntegrationDriver(const G4BFieldIntegrationDriver &) = delete;
    const G4BFieldIntegrationDriver& operator =(const G4BFieldIntegrationDriver &) = delete;

    virtual G4bool QuickAdvance(
        G4FieldTrack& fieldTrack,
        const G4double dydx[],
        G4double hstep,
        G4double inverseCurvatureRadius,
        G4double& dchord_step,
        G4double& dyerr) override;

    virtual void SetEquationOfMotion(G4EquationOfMotion* equation) override;

    virtual G4double GetInverseCurvatureRadius(const G4FieldTrack& track,
                                               G4double field[]) const override;

private:
    G4double fallbackThreshold;
    G4Mag_EqRhs* fequation;
    G4HelixExplicitEuler fallbackStepper;
};

#include "G4BFieldIntegrationDriver.icc"

#endif
