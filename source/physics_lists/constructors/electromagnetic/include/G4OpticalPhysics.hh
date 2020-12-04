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

#include "G4VPhysicsConstructor.hh"
#include "G4OpticalParameters.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4OpticalPhysics : public G4VPhysicsConstructor
{
  public:

    G4OpticalPhysics(G4int verbose = 0, const G4String& name = "Optical");
    virtual ~G4OpticalPhysics();
    virtual void PrintStatistics() const;

  protected:

    // construct particle and physics
    virtual void ConstructParticle();
    virtual void ConstructProcess();

  private:

    G4OpticalPhysics(const G4OpticalPhysics& right) = delete;
    G4OpticalPhysics& operator=(const G4OpticalPhysics& right) = delete;

  public:

    // DEPRECATED
    // the methods below are kept for backwards compatibility, and are to be
    // removed in future. Please use the methods in 
    // G4OpticalParameters instead.

    void Configure(G4OpticalProcessIndex, G4bool);
    void SetTrackSecondariesFirst(G4OpticalProcessIndex, G4bool);

    // Cerenkov
    void SetMaxNumPhotonsPerStep(G4int);
    void SetMaxBetaChangePerStep(G4double);
    void SetCerenkovStackPhotons(G4bool);
    void SetCerenkovTrackSecondariesFirst(G4bool);
    void SetCerenkovVerbosity(G4int);

    // Scintillation
    void SetScintillationYieldFactor(G4double);
    void SetScintillationExcitationRatio(G4double);
    void SetScintillationByParticleType(G4bool);
    void SetScintillationTrackInfo(G4bool);
    void SetScintillationTrackSecondariesFirst(G4bool);
    void SetFiniteRiseTime(G4bool);
    void SetScintillationStackPhotons(G4bool);
    void SetScintillationVerbosity(G4int);
    void SetScintillationEnhancedTimeConstants(G4bool);
    //void AddScintillationSaturation(G4EmSaturation* );

    // WLS
    void SetWLSTimeProfile(G4String);
    void SetWLSVerbosity(G4int);

    // WLS2
    void SetWLS2TimeProfile(G4String);
    void SetWLS2Verbosity(G4int);

    //boundary
    void SetBoundaryVerbosity(G4int);
    void SetInvokeSD(G4bool);

    void SetAbsorptionVerbosity(G4int);
    void SetRayleighVerbosity(G4int);
    void SetMieVerbosity(G4int);

  private:

    void PrintWarning(G4ExceptionDescription&) const;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif // G4OpticalPhysics_h
