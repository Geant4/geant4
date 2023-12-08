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
//------------------ G4XrayReflection physics process -----------------------
//
// History:
// 19-09-23 H. Burkhardt, initial implementation
//
// --------------------------------------------------------------------------

#ifndef G4XrayReflection_h
#define G4XrayReflection_h 1

#include "G4Gamma.hh"
#include "G4VEmProcess.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4ParticleDefinition;
class G4VEmModel;
class G4PropagatorInField;

class G4XrayReflection : public G4VDiscreteProcess
{
  public:  // with description
    explicit G4XrayReflection(const G4String& processName = "XrayReflection",
                              G4ProcessType type = fElectromagnetic);
    ~G4XrayReflection() override = default;

    G4double GetMeanFreePath(const G4Track& track, G4double previousStepSize,
                             G4ForceCondition* condition) override;

    G4VParticleChange* PostStepDoIt(const G4Track& track, const G4Step& Step) override;

    G4bool IsApplicable(const G4ParticleDefinition&) override;
    void BuildPhysicsTable(const G4ParticleDefinition&) override;

    void ProcessDescription(std::ostream&) const override;
    void DumpInfo() const override { ProcessDescription(G4cout); }

    void SetSurfaceRoughness(const G4double value);

    G4double Reflectivity(const G4double GamEner, const G4double SinIncidentAngle,
                          const G4Material* theMat) const;

    G4int ReadHenkeXrayData(std::string ElName, std::vector<G4double>& Ephot,
                            std::vector<G4double>& f1, std::vector<G4double>& f2);
    void SaveHenkeDataAsMaterialProperty();  // Save Henke data as Material propoerty

  private:
    G4VPhysicalVolume* fLastVolume = nullptr;
    G4ThreeVector fSurfaceNormal;
    static G4double fSurfaceRoughness;  // same value in all threads
};

#endif
