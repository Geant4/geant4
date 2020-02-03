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
// G4BFieldIntegrationDriver implementation
//
// Specialized integration driver for pure magnetic field
//
// Author: D.Sorokin
// --------------------------------------------------------------------

#include "G4BFieldIntegrationDriver.hh"

#include "G4FieldTrack.hh"
#include "G4FieldUtils.hh"
#include "G4Exception.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "templates.hh"


namespace {

G4Mag_EqRhs* toMagneticEquation(G4EquationOfMotion* equation)
{
    auto e = dynamic_cast<G4Mag_EqRhs*>(equation);

    if (!e) 
    {
        G4Exception("G4BFieldIntegrationDriver::G4BFieldIntegrationDriver",
                    "GeomField0003", FatalErrorInArgument,
                    "Works only with G4Mag_EqRhs");
    }

    return e;
}

} // namespace


G4BFieldIntegrationDriver::G4BFieldIntegrationDriver(
    std::unique_ptr<G4VIntegrationDriver> smallStepDriver, 
    std::unique_ptr<G4VIntegrationDriver> largeStepDriver)
    : fSmallStepDriver(std::move(smallStepDriver)),
      fLargeStepDriver(std::move(largeStepDriver)),
      fCurrDriver(fSmallStepDriver.get()),
      fEquation(toMagneticEquation(fCurrDriver->GetEquationOfMotion()))
{
    if (fSmallStepDriver->GetEquationOfMotion()
     != fLargeStepDriver->GetEquationOfMotion())
    {
        G4Exception("G4BFieldIntegrationDriver Constructor:",
                    "GeomField1001", FatalException, "different EoM");  
    }
}

G4double G4BFieldIntegrationDriver::AdvanceChordLimited(G4FieldTrack& yCurrent, 
                                                        G4double stepMax, 
                                                        G4double epsStep, 
                                                        G4double chordDistance)
{
    const G4double radius = CurvatureRadius(yCurrent);

    G4VIntegrationDriver* driver = nullptr;
    if (chordDistance < 2 * radius)
    {
        stepMax = std::min(stepMax, twopi * radius);
        driver = fSmallStepDriver.get();
        ++fSmallDriverSteps;
    } else
    {
        driver = fLargeStepDriver.get();
        ++fLargeDriverSteps;
    }

    if (driver != fCurrDriver)
    {
        driver->OnComputeStep();
    }

    fCurrDriver = driver;

    return fCurrDriver->AdvanceChordLimited(yCurrent, stepMax,
                                            epsStep, chordDistance);
}

void
G4BFieldIntegrationDriver::SetEquationOfMotion(G4EquationOfMotion* equation)
{
    fEquation = toMagneticEquation(equation);
    fSmallStepDriver->SetEquationOfMotion(equation);
    fLargeStepDriver->SetEquationOfMotion(equation);
}

G4double
G4BFieldIntegrationDriver::CurvatureRadius(const G4FieldTrack& track) const
{
    G4double field[G4Field::MAX_NUMBER_OF_COMPONENTS];
    
    GetFieldValue(track, field);
    
    const G4double Bmag2 =   field[0] * field[0]
                           + field[1] * field[1]
                           + field[2] * field[2] ;
    if (Bmag2 == 0.0 )
    {
       return DBL_MAX;
    }

    const G4double momentum2 = track.GetMomentum().mag2();
    const G4double fCof_inv = eplus / std::abs(fEquation->FCof());

    return std::sqrt(momentum2 / Bmag2) * fCof_inv;
}

void
G4BFieldIntegrationDriver::GetFieldValue(const G4FieldTrack& track,
                                         G4double Field[] ) const
{
    G4ThreeVector pos= track.GetPosition();
    G4double positionTime[4]= { pos.x(), pos.y(), pos.z(),
                                track.GetLabTimeOfFlight() } ;
    
    fEquation->GetFieldValue(positionTime, Field);
}

void G4BFieldIntegrationDriver::PrintStatistics() const
{
    const auto totSteps = fSmallDriverSteps + fLargeDriverSteps;
    const auto toFraction = [&](double value) { return value / totSteps * 100; };

    G4cout << "============= G4BFieldIntegrationDriver statistics ===========\n"
           << "total steps " << totSteps << " "
           << "smallDriverSteps " << toFraction(fSmallDriverSteps) << " "
           << "largeDriverSteps " << toFraction(fLargeDriverSteps) << "\n"
           << "======================================\n";
}
