//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//
//      ---------- G4NuclearStopping physics process -----
//                by Vladimir Ivanchenko, 27 April 2004
//
// Modifications:
//

// ------------------------------------------------------------

// Class Description:
// Nuclear process of charged hadrons and ions is the implementation
// of non-ionizing energy loss
// The physics model is described in CERN-OPEN-99-121 and CERN-OPEN-99-300.
// The user may select parametrisation tables for nuclear stopping powers
// The list of available tables:
//     "ICRU_49" (default), "Ziegler1977", "Ziegler1985"
// Further documentation available from http://www.ge.infn.it/geant4/lowE
// and in the Physics Reference Manual

// ------------------------------------------------------------

#ifndef G4NuclearStopping_h
#define G4NuclearStopping_h 1

#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4VContinuousProcess.hh"
#include "G4ParticleChangeForLoss.hh"
#include "G4hNuclearStoppingModel.hh"
#include <map>

class G4Track;
class G4Step;
class G4Material;

class G4NuclearStopping : public G4VContinuousProcess
{
public: // With description

  G4NuclearStopping(const G4String& processName = "NucEloss",
                          G4ProcessType aType = fElectromagnetic);
  // The ionisation process for hadrons/ions to be include in the
  // UserPhysicsList

  ~G4NuclearStopping();
  // Destructor

  G4bool IsApplicable(const G4ParticleDefinition&);
  // True for all charged hadrons/ions

  virtual G4VParticleChange* AlongStepDoIt(const G4Track&,
                                           const G4Step&);

  virtual void BuildPhysicsTable(const G4ParticleDefinition&);

  virtual void PrintInfoDefinition();
  // Print out of the class parameters

  void SetHighEnergy(G4double energy) {highEnergy = energy;} ;
  // For higher energies nuclear stopping power neglected

  void SetLowEnergy(G4double energy) {lowEnergy = energy;} ;
  // For lower energies nuclear stopping power is constant

  void SetNuclearStoppingPowerModel(const G4String& dedxTable);
  // This method defines the electron ionisation parametrisation method
  // via the name of the table. Default is "ICRU_49".

  G4double StoppingPower(const G4Track&);

  void AddSaturationFactor(const G4Material* material, G4double val);

  void SetFluctuations(G4bool val) {fluctuations = val;};

protected:

  virtual G4double GetContinuousStepLimit(
                       const G4Track& aTrack,
                             G4double  previousStepSize,
                             G4double  currentMinimumStep,
                             G4double& currentSafety);

  G4ParticleChangeForLoss    fParticleChange;

private:

  // hide assignment operator
  G4NuclearStopping & operator=(const G4NuclearStopping &right);
  G4NuclearStopping(const G4NuclearStopping&);

  // name of parametrisation table of electron stopping power
  G4String                    theTable;

  // interval of parametrisation of nuclear stopping power
  G4double                    lowEnergy;
  G4double                    highEnergy;

  G4hNuclearStoppingModel*    model;
  G4bool                      fluctuations;
  G4bool                      initialised;

  std::map<const G4Material*,double,std::less<const G4Material*> > factors;

  //cash
  G4double                    preStepDEDX;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4bool G4NuclearStopping::IsApplicable(
                                const G4ParticleDefinition& particle)
{
  return(particle.GetPDGCharge() != 0.0
      && particle.GetPDGMass() > proton_mass_c2*0.1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4NuclearStopping::GetContinuousStepLimit(
          const G4Track& aTrack, G4double, G4double, G4double&)
{
  preStepDEDX = StoppingPower(aTrack);
  G4double step = DBL_MAX;
  if(preStepDEDX > 0.0) step = aTrack.GetKineticEnergy()/preStepDEDX;
  return step;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4NuclearStopping::StoppingPower(const G4Track& aTrack)
{
  G4double dedx = 0.0;
  G4double e = aTrack.GetKineticEnergy();
  if(e < lowEnergy) e = lowEnergy;
  if(e < highEnergy) dedx = model->TheValue(aTrack.GetDefinition(),
                                            aTrack.GetMaterial(), e);
  return dedx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif

