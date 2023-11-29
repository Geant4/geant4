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
// -------------------------------------------------------------------
//
//      GEANT4 class file
//
//      CERN, Geneva, Switzerland
//
//      File name:     G4PhotonEvaporation
//
//      Author:        Vladimir Ivantchenko
//
//      Creation date: 22 October 2015 
//
// -------------------------------------------------------------------
//
// This is gamma deexcitation model based on the nuclear levels data
//

#ifndef G4PHOTONEVAPORATION_HH
#define G4PHOTONEVAPORATION_HH 1

#include "globals.hh"
#include "G4VEvaporationChannel.hh"
#include "G4NuclearLevelData.hh"
#include "G4LevelManager.hh"
#include "G4Fragment.hh"

const G4int MAXDEPOINT = 10;
const G4int MAXGRDATA = 300;

class G4GammaTransition;

class G4PhotonEvaporation : public G4VEvaporationChannel {

public:

  explicit G4PhotonEvaporation(G4GammaTransition* ptr=nullptr);

  ~G4PhotonEvaporation() override;

  void Initialise() override;

  // one photon or e- emission
  G4Fragment* EmittedFragment(G4Fragment* theNucleus) override;

  // returns "false", emitted gamma and e- are added to the results
  G4bool 
  BreakUpChain(G4FragmentVector* theResult, G4Fragment* theNucleus) override;

  // emitted gamma, e-, and residual fragment are added to the results
  G4FragmentVector* BreakItUp(const G4Fragment& theNucleus);

  // compute emission probability for both continum and discrete cases
  // must be called before any method above
  G4double GetEmissionProbability(G4Fragment* theNucleus) override;

  // methods for unit tests
  G4double ComputeInverseXSection(G4Fragment* theNucleus, 
                                  G4double kinEnergy) override;
  G4double ComputeProbability(G4Fragment* theNucleus, 
			      G4double kinEnergy) override;

  G4double GetFinalLevelEnergy(G4int Z, G4int A, G4double energy);

  G4double GetUpperLevelEnergy(G4int Z, G4int A);

  void SetGammaTransition(G4GammaTransition*);

  void SetICM(G4bool) override;

  void RDMForced (G4bool) override;
  
  inline void SetVerboseLevel(G4int verbose);

  inline G4int GetVacantShellNumber() const;

  G4PhotonEvaporation(const G4PhotonEvaporation & right) = delete;
  const G4PhotonEvaporation & operator = 
    (const G4PhotonEvaporation & right) = delete;
 
private:

  void InitialiseGRData();

  G4Fragment* GenerateGamma(G4Fragment* nucleus);

  inline void InitialiseLevelManager(G4int Z, G4int A);

  G4NuclearLevelData*   fNuclearLevelData;
  const G4LevelManager* fLevelManager;
  G4GammaTransition*    fTransition;

  // fPolarization stores polarization tensor for consecutive
  // decays of a nucleus 
  G4NuclearPolarization* fPolarization;

  G4int    fVerbose;
  G4int    theZ;
  G4int    theA;
  G4int    fPoints;
  G4int    fCode;
  G4int    vShellNumber;
  size_t   fIndex;

  G4int fSecID;  // Creator model ID for the secondaries created by this model

  static G4float GREnergy[MAXGRDATA];
  static G4float GRWidth[MAXGRDATA];

  G4double fCummProbability[MAXDEPOINT]; 

  G4double fLevelEnergyMax;
  G4double fExcEnergy;
  G4double fProbability;
  G4double fStep;
  G4double fMaxLifeTime;

  G4double fTolerance;

  G4bool   fICM;
  G4bool   fRDM;
  G4bool   fSampleTime;
  G4bool   fCorrelatedGamma;
  G4bool   isInitialised;
};

inline void G4PhotonEvaporation::SetVerboseLevel(G4int verbose)
{
  fVerbose = verbose;
}

inline void 
G4PhotonEvaporation::InitialiseLevelManager(G4int Z, G4int A)
{
  if(Z != theZ || A != theA) {
    theZ = Z;
    theA = A;
    fIndex = 0;
    fLevelManager = fNuclearLevelData->GetLevelManager(theZ, theA);
    fLevelEnergyMax = fLevelManager ? fLevelManager->MaxLevelEnergy() : 0.0;
  }
}

inline G4int G4PhotonEvaporation::GetVacantShellNumber() const 
{ 
  return vShellNumber;
}

#endif
