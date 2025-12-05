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
// GEM de-excitation model
// by V. Ivanchenko (July 2019)
//

#ifndef G4GEMChannelVI_h
#define G4GEMChannelVI_h 1

#include "G4VEvaporationChannel.hh"
#include "G4VSIntegration.hh"

class G4PairingCorrection;
class G4VCoulombBarrier;
class G4LevelManager;
class G4NuclearLevelData;
class G4InterfaceToXS;
class G4ParticleDefinition;
class G4Pow;

class G4GEMChannelVI : public G4VEvaporationChannel, public G4VSIntegration
{
public:

  explicit G4GEMChannelVI(G4int theA, G4int theZ);

  ~G4GEMChannelVI() override;

  void Initialise() override;

  G4double ProbabilityDensityFunction(G4double ekin) override;

  G4double GetEmissionProbability(G4Fragment* theNucleus) override;

  G4Fragment* EmittedFragment(G4Fragment* theNucleus) override;

  const G4String& ModelName() const override;

  G4double GetCurrentXS() { return recentXS; };

  void Dump() const override;

  G4GEMChannelVI(const G4GEMChannelVI & right) = delete;
  const G4GEMChannelVI & operator=(const G4GEMChannelVI & right) = delete;
  G4bool operator==(const G4GEMChannelVI & right) const = delete;
  G4bool operator!=(const G4GEMChannelVI & right) const = delete;

private: 

  G4double CrossSection(G4double ekin);

  G4double CorrectExcitation(G4double energy, const G4LevelManager*);

  G4NuclearLevelData* nData;
  const G4VCoulombBarrier* cBarrier;
  const G4PairingCorrection* pairingCorrection;
  const G4LevelManager* lManagerEvap{nullptr};
  const G4LevelManager* lManagerRes{nullptr};
  G4InterfaceToXS* fXSection{nullptr};
  G4Pow* g4pow;
  const G4ParticleDefinition* fProton;
  const G4ParticleDefinition* fNeutron;
  
  G4double fEvapMass;     // ground state mass of the evaporated fragment 
  G4double fEvapMass2;    // ground state mass of the evaporated fragment square
  G4double fMass{0.0};    // mass of the initial fragment
  G4double fResMass{0.0}; // ground state mass of the residual fragment
  G4double fResA13{0.0};  //
  G4double fFragExc{0.0}; // excitation energy of the evaporated fragment
  G4double fEvapExc{0.0}; // excitation energy of the evaporated fragment
  G4double fResExc{0.0};  // excitation energy of the residual fragment
  G4double bCoulomb{0.0};
  G4double fDeltaEvap{0.0};
  G4double fE0{0.0};
  G4double fE1{0.0};
  G4double a0{0.0};
  G4double a1{0.0};
  G4double delta0{0.0};
  G4double delta1{0.0};
  G4double recentXS{0.0};
  G4double fEnergyLimitXS{0.0};
  G4double xsfactor{1.0};
  G4double fTolerance;
  G4double fCoeff;

  G4int evapA;
  G4int evapZ;
  G4int resA{0};
  G4int resZ{0};
  G4int fragA{0};
  G4int fragZ{0};
  G4int fVerbose{1};
  G4int nProbEvap{1};
  G4int nProbRes{1};
  G4int indexC{7};
  G4int secID;

  G4String fModelName;
};

#endif
