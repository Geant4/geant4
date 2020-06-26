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
// G4BFieldIntegrationDriver
//
// Class description:
//
// Specialized integration driver for pure magnetic field

// Author: D.Sorokin
// --------------------------------------------------------------------
#ifndef G4BFIELD_INTEGRATION_DRIVER_HH
#define G4BFIELD_INTEGRATION_DRIVER_HH

#include "G4VIntegrationDriver.hh"
#include "G4Mag_EqRhs.hh"

#include <memory>

class G4BFieldIntegrationDriver : public G4VIntegrationDriver
{
  public:

    G4BFieldIntegrationDriver(
        std::unique_ptr<G4VIntegrationDriver> smallStepDriver, 
        std::unique_ptr<G4VIntegrationDriver> largeStepDriver);

    G4BFieldIntegrationDriver(const G4BFieldIntegrationDriver &) = delete;
    const G4BFieldIntegrationDriver& operator =(const G4BFieldIntegrationDriver &) = delete;

    virtual G4double AdvanceChordLimited(G4FieldTrack& track,
                                         G4double hstep,
                                         G4double eps,
                                         G4double chordDistance) override;

    virtual G4bool AccurateAdvance(G4FieldTrack& track,
                                   G4double hstep,
                                   G4double eps,
                                   G4double hinitial = 0) override
    {
        return fCurrDriver->AccurateAdvance(track, hstep, eps, hinitial);
    }

    virtual G4bool DoesReIntegrate() const override
    {
       return fCurrDriver->DoesReIntegrate();
    }
   
    //[[deprecated("will be removed")]]
    virtual void GetDerivatives(const G4FieldTrack& track,
                                G4double dydx[]) const override
    {
        fCurrDriver->GetDerivatives(track, dydx);
    }

    //[[deprecated("will be removed")]]
    virtual void GetDerivatives(const G4FieldTrack& track,
                                G4double dydx[],
                                G4double field[]) const override
    {
        fCurrDriver->GetDerivatives(track, dydx, field);
    }

    virtual void SetEquationOfMotion(G4EquationOfMotion* equation) override;

    virtual G4EquationOfMotion* GetEquationOfMotion() override
    {
        return fCurrDriver->GetEquationOfMotion();
    }

    //[[deprecated("use GetEquationOfMotion() instead of GetStepper()->GetEquationOfMotion()")]]
    virtual const G4MagIntegratorStepper* GetStepper() const override
    {
        return fCurrDriver->GetStepper();
    }

    virtual G4MagIntegratorStepper* GetStepper() override
    {
        return fCurrDriver->GetStepper();
    }

    virtual G4double ComputeNewStepSize(G4double errMaxNorm,
                                        G4double hstepCurrent) override
    {
        return fCurrDriver->ComputeNewStepSize(errMaxNorm, hstepCurrent);
    }

    virtual void SetVerboseLevel(G4int level) override
    {
        fSmallStepDriver->SetVerboseLevel(level);
        fLargeStepDriver->SetVerboseLevel(level);
    }

    virtual G4int GetVerboseLevel() const override
    {
        return fCurrDriver->GetVerboseLevel();
    }

    virtual void OnComputeStep() override
    {
        fSmallStepDriver->OnComputeStep();
        fLargeStepDriver->OnComputeStep();
    }

    virtual void OnStartTracking() override
    {
        fSmallStepDriver->OnStartTracking();
        fLargeStepDriver->OnStartTracking();
    }

    virtual void  StreamInfo( std::ostream& os ) const override
    {
       os << "Small Step Driver Info: " << std::endl;
       fSmallStepDriver->StreamInfo(os);
       os << "Large Step Driver Info: " << std::endl;        
       fLargeStepDriver->StreamInfo(os);
    }
    // Write out the parameters / state of the driver
   
    void PrintStatistics() const;

  private:

    G4double CurvatureRadius(const G4FieldTrack& track) const;

    void GetFieldValue(const G4FieldTrack& track, 
                             G4double      Field[] ) const;
   
    std::unique_ptr<G4VIntegrationDriver> fSmallStepDriver;
    std::unique_ptr<G4VIntegrationDriver> fLargeStepDriver;
    G4VIntegrationDriver* fCurrDriver = nullptr;
    G4Mag_EqRhs* fEquation = nullptr;

    G4int fSmallDriverSteps = 0;
    G4int fLargeDriverSteps = 0;
};

#endif
