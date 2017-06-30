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
// $Id: G4PhotonEvaporation.hh 103899 2017-05-03 08:26:53Z gcosmo $
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
//Modifications:
//
// 
// -------------------------------------------------------------------
//
// This is a new class which has different design and uses different data 
// structure than the old one
//

#ifndef G4PHOTONEVAPORATION_HH
#define G4PHOTONEVAPORATION_HH 1

#include "globals.hh"
#include "G4VEvaporationChannel.hh"
#include "G4NuclearLevelData.hh"
#include "G4LevelManager.hh"
#include "G4Fragment.hh"

const G4int MAXDEPOINT = 10;
const G4int MAXGRDATA  = 300;

class G4GammaTransition;

class G4PhotonEvaporation : public G4VEvaporationChannel {

public:

  explicit G4PhotonEvaporation(G4GammaTransition* ptr=nullptr);

  virtual ~G4PhotonEvaporation();

  virtual void Initialise() final;

  // one photon or e- emission
  virtual G4Fragment* EmittedFragment(G4Fragment* theNucleus) final;

  // returns "false", emitted gamma and e- are added to the results
  virtual G4bool 
  BreakUpChain(G4FragmentVector* theResult, G4Fragment* theNucleus) final;

  // emitted gamma, e-, and residual fragment are added to the results
  G4FragmentVector* BreakItUp(const G4Fragment& theNucleus);

  // compute emission probability for both continum and discrete cases
  // must be called before any method above
  virtual G4double GetEmissionProbability(G4Fragment* theNucleus) final;

  virtual G4double GetFinalLevelEnergy(G4int Z, G4int A, G4double energy) final;

  virtual G4double GetUpperLevelEnergy(G4int Z, G4int A) final;

  void SetGammaTransition(G4GammaTransition*);

  void SetMaxHalfLife(G4double);

  virtual void SetICM(G4bool);

  virtual void RDMForced (G4bool);
  
  inline void SetVerboseLevel(G4int verbose);

  inline G4int GetVacantShellNumber() const;
 
private:

  G4Fragment* GenerateGamma(G4Fragment* nucleus);

  inline void InitialiseLevelManager(G4int Z, G4int A);

  G4PhotonEvaporation(const G4PhotonEvaporation & right) = delete;
  const G4PhotonEvaporation & operator = (const G4PhotonEvaporation & right) = delete;

  G4NuclearLevelData*   fNuclearLevelData;
  const G4LevelManager* fLevelManager;
  G4GammaTransition*    fTransition;

  G4int    fVerbose;
  G4int    theZ;
  G4int    theA;
  G4int    fPoints;
  G4int    fCode;
  G4int    vShellNumber;
  size_t   fIndex;

  static G4float GREnergy[MAXGRDATA];
  static G4float GRWidth[MAXGRDATA];

  G4double fCummProbability[MAXDEPOINT]; 

  G4double fLevelEnergyMax;
  G4double fExcEnergy;
  G4double fProbability;
  G4double fStep;
  G4double fMaxLifeTime;

  G4double LevelDensity;
  G4double Tolerance;

  G4bool   fICM;
  G4bool   fRDM;
  G4bool   fSampleTime;
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
