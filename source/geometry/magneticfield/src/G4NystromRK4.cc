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
// History:
// - Created:      I.Gavrilenko    15.05.2009   (as G4AtlasRK4)
// - Adaptations:  J.Apostolakis  May-Nov 2009
// -------------------------------------------------------------------

#include "G4NystromRK4.hh"

#include "G4Exception.hh"
#include "G4SystemOfUnits.hh"
#include "G4FieldUtils.hh"
#include "G4LineSection.hh"

using namespace field_utils;

namespace {

G4bool notEquals(G4double p1, G4double p2)
{
    return std::fabs(p1 - p2) > perMillion * p2;
  }
   constexpr G4int INTEGRATED_COMPONENTS = 6;
} // namespace


G4NystromRK4::G4NystromRK4(G4Mag_EqRhs* equation, G4double distanceConstField)
    : G4MagIntegratorStepper(equation, INTEGRATED_COMPONENTS),
      fMomentum(0),
      fMomentum2(0),
      fInverseMomentum(0),
      fCoefficient(0)
{
    if (distanceConstField > 0)
    {
        SetDistanceForConstantField(distanceConstField);
    }
}

void G4NystromRK4::Stepper(const G4double y[],
                           const G4double dydx[],
                           G4double hstep, 
                           G4double yOut[], 
                           G4double yError[])
{
    fInitialPoint = { y[0], y[1], y[2] };

    G4double field[3];

    constexpr G4double one_sixth= 1./6.;
    const G4double S5 = 0.5 * hstep;
    const G4double S4 = .25 * hstep;
    const G4double S6 = hstep * one_sixth;
   
    const G4double momentum2 = getValue2(y, Value3D::Momentum);
    if (notEquals(momentum2, fMomentum2))
    {
        fMomentum = std::sqrt(momentum2);
        fMomentum2 = momentum2;
        fInverseMomentum  = 1. / fMomentum;
        fCoefficient = GetFCof() * fInverseMomentum;
    }
    
    // Point 1
    const G4double K1[3] = { 
        fInverseMomentum * dydx[3], 
        fInverseMomentum * dydx[4], 
        fInverseMomentum * dydx[5]
    };
    
    // Point2
    G4double p[4] = {
        y[0] + S5 * (dydx[0] + S4 * K1[0]),
        y[1] + S5 * (dydx[1] + S4 * K1[1]),
        y[2] + S5 * (dydx[2] + S4 * K1[2]),
        y[7]
    };

    GetFieldValue(p, field);
  
    const G4double A2[3] = {
        dydx[0] + S5 * K1[0],
        dydx[1] + S5 * K1[1],
        dydx[2] + S5 * K1[2]
    };

    const G4double K2[3] = {
        (A2[1] * field[2] - A2[2] * field[1]) * fCoefficient,
        (A2[2] * field[0] - A2[0] * field[2]) * fCoefficient,
        (A2[0] * field[1] - A2[1] * field[0]) * fCoefficient
    };

    fMidPoint = { p[0], p[1], p[2] };

    // Point 3 with the same magnetic field
    const G4double A3[3] = {
        dydx[0] + S5 * K2[0],
        dydx[1] + S5 * K2[1],
        dydx[2] + S5 * K2[2]
    };

    const G4double K3[3] = {
        (A3[1] * field[2] - A3[2] * field[1]) * fCoefficient,
        (A3[2] * field[0] - A3[0] * field[2]) * fCoefficient,
        (A3[0] * field[1] - A3[1] * field[0]) * fCoefficient
    };

    // Point 4
    p[0] = y[0] + hstep * (dydx[0] + S5 * K3[0]);
    p[1] = y[1] + hstep * (dydx[1] + S5 * K3[1]);
    p[2] = y[2] + hstep * (dydx[2] + S5 * K3[2]);             

    GetFieldValue(p, field);
  
    const G4double A4[3] = {
        dydx[0] + hstep * K3[0],
        dydx[1] + hstep * K3[1],
        dydx[2] + hstep * K3[2]
    };

    const G4double K4[3] = {
        (A4[1] * field[2] - A4[2] * field[1]) * fCoefficient,
        (A4[2] * field[0] - A4[0] * field[2]) * fCoefficient,
        (A4[0] * field[1] - A4[1] * field[0]) * fCoefficient
    };
  
    // New position
    yOut[0] = y[0] + hstep * (dydx[0] + S6 * (K1[0] + K2[0] + K3[0]));
    yOut[1] = y[1] + hstep * (dydx[1] + S6 * (K1[1] + K2[1] + K3[1]));
    yOut[2] = y[2] + hstep * (dydx[2] + S6 * (K1[2] + K2[2] + K3[2]));
    // New direction
    yOut[3] = dydx[0] + S6 * (K1[0] + K4[0] + 2. * (K2[0] + K3[0]));
    yOut[4] = dydx[1] + S6 * (K1[1] + K4[1] + 2. * (K2[1] + K3[1]));
    yOut[5] = dydx[2] + S6 * (K1[2] + K4[2] + 2. * (K2[2] + K3[2]));
    // Pass Energy, time unchanged -- time is not integrated !!
    yOut[6] = y[6]; 
    yOut[7] = y[7];

    fEndPoint = { yOut[0], yOut[1], yOut[2] };

    // Errors
    yError[3] = hstep * std::fabs(K1[0] - K2[0] - K3[0] + K4[0]);
    yError[4] = hstep * std::fabs(K1[1] - K2[1] - K3[1] + K4[1]);
    yError[5] = hstep * std::fabs(K1[2] - K2[2] - K3[2] + K4[2]);
    yError[0] = hstep * yError[3];
    yError[1] = hstep * yError[4];
    yError[2] = hstep * yError[5];
    yError[3] *= fMomentum;
    yError[4] *= fMomentum;
    yError[5] *= fMomentum;

    // Normalize momentum
    const G4double normF = fMomentum / getValue(yOut, Value3D::Momentum);
    yOut[3] *= normF; 
    yOut[4] *= normF; 
    yOut[5] *= normF; 

    // My trial code:
    // G4double endMom2 = yOut[3]*yOut[3]+yOut[4]*yOut[4]+yOut[5]*yOut[5];
    // G4double normF = std::sqrt( startMom2 / endMom2 );  
}
        
G4double G4NystromRK4::DistChord() const 
{
    return G4LineSection::Distline(fMidPoint, fInitialPoint, fEndPoint);
}

void G4NystromRK4::SetDistanceForConstantField(G4double length)
{
    if (!GetField())
    {
        G4Exception("G4NystromRK4::SetDistanceForConstantField","Nystrom 001",
                    JustWarning, "Provided field is not G4CachedMagneticField. Changing field type.");

        fCachedField = std::unique_ptr<G4CachedMagneticField>(
            new G4CachedMagneticField(
                dynamic_cast<G4MagneticField*>(GetEquationOfMotion()->GetFieldObj()), 
                length));

        GetEquationOfMotion()->SetFieldObj(fCachedField.get());
    }

    GetField()->SetConstDistance(length);
}

G4double G4NystromRK4::GetDistanceForConstantField() const
{
    if (!GetField())
    {
        return 0;
    }

    return GetField()->GetConstDistance(); 
}

G4CachedMagneticField* G4NystromRK4::GetField()
{
    return dynamic_cast<G4CachedMagneticField*>(GetEquationOfMotion()->GetFieldObj());
}

const G4CachedMagneticField* G4NystromRK4::GetField() const
{
    return const_cast<G4NystromRK4*>(this)->GetField();
}
