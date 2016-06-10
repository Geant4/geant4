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
//---------------------------------------------------------------------------
//
// ClassName:   G4OpticalPhysics
//
// Author:      P.Gumplinger 30.09.2009
//
// Modified:    P.Gumplinger 29.09.2011
//              (based on code from I. Hrivnacova)
//
//---------------------------------------------------------------------------
//
// This class provides construction of default optical physics
//

#ifndef G4OpticalPhysics_h
#define G4OpticalPhysics_h 1

#include "globals.hh"

#include "G4VPhysicsConstructor.hh"

#include "G4OpticalProcessIndex.hh"
#include "G4OpticalPhysicsMessenger.hh"

#include "G4OpticalSurface.hh"

#include <vector>

class G4VProcess;
class G4EmSaturation;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4OpticalPhysics : public G4VPhysicsConstructor
{
  public:

    G4OpticalPhysics(G4int verbose = 0, const G4String& name = "Optical");
    virtual ~G4OpticalPhysics();

  protected:

    // construct particle and physics
    virtual void ConstructParticle();
    virtual void ConstructProcess();

  private:

    /// Not implemented
    G4OpticalPhysics(const G4OpticalPhysics& right);
    /// Not implemented
    G4OpticalPhysics& operator=(const G4OpticalPhysics& right);

  public:

    // configure G4OpticalPhysics builder
    void Configure(G4OpticalProcessIndex, G4bool );

    void SetMaxNumPhotonsPerStep(G4int );
    void SetMaxBetaChangePerStep(G4double );
    void SetScintillationYieldFactor(G4double );
    void SetScintillationExcitationRatio(G4double );

    void SetWLSTimeProfile(G4String );
    void SetScintillationByParticleType(G4bool );
    //void AddScintillationSaturation(G4EmSaturation* );

    void SetTrackSecondariesFirst(G4OpticalProcessIndex, G4bool );
    void SetFiniteRiseTime(G4bool );

  private:

    // methods
    void PrintStatistics() const;

    // messenger
    G4OpticalPhysicsMessenger* fMessenger;

    // The vector of process configuration
    std::vector<G4bool>         fProcessUse;

    // The vector of track secondaries options;
    // the option to track secondaries before finishing their parent track
    std::vector<G4bool>         fProcessTrackSecondariesFirst;

    /// max number of Cerenkov photons per step
    G4int                       fMaxNumPhotons;

    /// max change of beta per step
    G4double                    fMaxBetaChange;

    /// scintillation yield factor
    G4double                    fYieldFactor;

    /// scintillation excitation ratio
    G4double                    fExcitationRatio;

    /// the WLS process time profile
    G4String                    fProfile;

    /// option to set a finite rise-time; Note: the G4Scintillation
    /// process expects the user to have set the constant material
    /// property FAST/SLOWSCINTILLATIONRISETIME
    G4bool                      fFiniteRiseTime;

    /// option to  allow for the light yield to be a function of
    /// particle type and deposited energy in case of non-linear
    /// light emission in scintillators
    G4bool                      fScintillationByParticleType;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif // G4OpticalPhysics_h
