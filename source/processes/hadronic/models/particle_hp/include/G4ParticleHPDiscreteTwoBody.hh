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
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#ifndef G4ParticleHPDiscreteTwoBody_h
#define G4ParticleHPDiscreteTwoBody_h 1

// 101110 Bug fix in MF=6, LAW=2 case; contribution from E. Mendoza, D. Cano-Ott (CIEMAT)

#include "G4InterpolationManager.hh"
#include "G4ParticleHPInterpolator.hh"
#include "G4ParticleHPLegendreTable.hh"
#include "G4VParticleHPEnergyAngular.hh"
#include "G4ios.hh"
#include "globals.hh"

#include <fstream>

class G4ParticleHPDiscreteTwoBody : public G4VParticleHPEnergyAngular
{
public:
  G4ParticleHPDiscreteTwoBody();
  ~G4ParticleHPDiscreteTwoBody() override;

  void Init(std::istream& aDataFile) override;

  G4ReactionProduct* Sample(G4double anEnergy, G4double massCode, G4double mass) override;
  G4double MeanEnergyOfThisInteraction() override { return -1.0; }

private:
  G4int nEnergy{0};
  G4ParticleHPLegendreTable* theCoeff{nullptr};
  G4bool bCheckDiffCoeffRepr{true};
  // for example ENDF-VII0_proton/Inelastic/F01/4_9_Beryllium has 0
  // for energy 7.5E+07 and 12 for energy 1.e+08

  G4InterpolationManager theManager;  // knows the interpolation between stores
  G4ParticleHPInterpolator theInt;
};

#endif
