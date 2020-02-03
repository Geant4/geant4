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

  ~G4GEMChannelVI() final;
    
  G4double GetEmissionProbability(G4Fragment* theNucleus) final;

  G4Fragment* EmittedFragment(G4Fragment* theNucleus) final;

  void Dump() const final;

private: 

  G4GEMChannelVI(const G4GEMChannelVI & right);  
  const G4GEMChannelVI & operator=(const G4GEMChannelVI & right);
  G4bool operator==(const G4GEMChannelVI & right) const;
  G4bool operator!=(const G4GEMChannelVI & right) const;

  const G4VCoulombBarrier* cBarrier;
  const G4PairingCorrection* pairingCorrection;
  G4GEMProbabilityVI* fProbability;
  
  G4int A;
  G4int Z;
  G4int resA;
  G4int resZ;
  G4int fragA;
  G4int fragZ;

  G4double mass;
  G4double resMass;
  G4double evapMass;
  G4double evapMass2;
};

#endif
