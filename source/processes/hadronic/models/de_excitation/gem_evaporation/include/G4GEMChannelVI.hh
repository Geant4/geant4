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

class G4PairingCorrection;
class G4VCoulombBarrier;
class G4LevelManager;
class G4GEMProbabilityVI;

class G4GEMChannelVI : public G4VEvaporationChannel
{
public:

  explicit G4GEMChannelVI(G4int theA, G4int theZ);

  ~G4GEMChannelVI() override;

  void Initialise() override;

  G4double GetEmissionProbability(G4Fragment* theNucleus) override;

  G4Fragment* EmittedFragment(G4Fragment* theNucleus) override;

  void Dump() const override;

  G4GEMChannelVI(const G4GEMChannelVI & right) = delete;
  const G4GEMChannelVI & operator=(const G4GEMChannelVI & right) = delete;
  G4bool operator==(const G4GEMChannelVI & right) const = delete;
  G4bool operator!=(const G4GEMChannelVI & right) const = delete;

private: 

  const G4VCoulombBarrier* cBarrier;
  const G4PairingCorrection* pairingCorrection;
  G4GEMProbabilityVI* fProbability;
  
  G4double fEvapMass;
  G4double fEvapMass2;
  G4double fMass{0.0};
  G4double fResMass{0.0};
  G4double fExc{0.0};
  G4double bCoulomb{0.0};
  G4double fCoeff;

  G4int A;
  G4int Z;
  G4int resA{0};
  G4int resZ{0};
  G4int fragA{0};
  G4int fragZ{0};
  G4int fVerbose{1};
  G4int nProb{1};
  G4int secID;
  G4int indexC;
  
  // evaporation fragment data
  struct evapData {
    G4double exc{0.0};   // excitation
    G4double ekin1{0.0}; // min kinetic energy
    G4double ekin2{0.0}; // max kinetic energy
    G4double prob{0.0};  // probability
  };
  evapData fEData[10]; 
};

#endif
