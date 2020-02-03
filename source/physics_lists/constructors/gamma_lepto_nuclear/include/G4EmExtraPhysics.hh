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
//---------------------------------------------------------------------------
//
// ClassName:   G4EmExtraPhysics
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 10.11.2005 V.Ivanchenko edit to provide a standard
// 19.06.2006 V.Ivanchenko add mu-nuclear process
// 16.10.2012 A.Ribon: renamed G4EmExtraBertiniPhysics as G4EmExtraPhysics
// 31.01.2018 V. Grichine: add neutrino-electron process and xsc
//
//----------------------------------------------------------------------------
//

#ifndef G4EmExtraPhysics_h
#define G4EmExtraPhysics_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

#include "G4EmMessenger.hh"

class G4CascadeInterface;
class G4PhotoNuclearProcess;

class G4EmExtraPhysics : public G4VPhysicsConstructor
{
public:

  G4EmExtraPhysics(G4int ver = 1);

  // obsolete
  G4EmExtraPhysics(const G4String& name);

  virtual ~G4EmExtraPhysics();

  void ConstructParticle();
  void ConstructProcess();

  void Synch(G4bool val);
  void SynchAll(G4bool val);
  void GammaNuclear(G4bool val);
  void LENDGammaNuclear(G4bool val);
  void ElectroNuclear(G4bool val);
  void MuonNuclear(G4bool val);
  void GammaToMuMu(G4bool val);
  void PositronToMuMu(G4bool val);
  void PositronToHadrons(G4bool val);
  void GammaToMuMuFactor(G4double val);
  void PositronToMuMuFactor(G4double val);
  void PositronToHadronsFactor(G4double val);

  void NeutrinoActivated(G4bool val);
  void NuETotXscActivated(G4bool val);
  void SetNuEleCcBias(G4double bf);
  void SetNuEleNcBias(G4double bf);
  void SetNuNucleusBias(G4double bf);
  void SetNuDetectorName(const G4String& dn);

private:

  void ConstructGammaElectroNuclear();

  void ConstructLENDGammaNuclear(G4CascadeInterface* cascade,
                                 G4PhotoNuclearProcess* gnuc);

  G4bool gnActivated;
  G4bool eActivated;
  G4bool gLENDActivated;
  G4bool munActivated;
  G4bool synActivated;
  G4bool synActivatedForAll;
  G4bool gmumuActivated;
  G4bool pmumuActivated;
  G4bool phadActivated;
  G4bool fNuActivated;
  G4bool fNuETotXscActivated;

  G4double gmumuFactor;
  G4double pmumuFactor;
  G4double phadFactor;
  G4double fNuEleCcBias;
  G4double fNuEleNcBias;
  G4double fNuNucleusBias;

  G4String fNuDetectorName;

  G4EmMessenger* theMessenger;
  G4int verbose;
};

#endif





