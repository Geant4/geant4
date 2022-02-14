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
// G4ForwardXrayTR
//
// Class for description
//
// Class for forward X-ray transition radiation generated
// by relativistic charged particle crossed interface between material 1
// and material 2 (1 -> 2)

// History:
// 22.09.97, V. Grichine (Vladimir.Grichine@cern.ch)
// 26.01.00, V.Grichine, new constructor and protected DM for fast sim. models
// 10.03.03, V.Ivanchenko migrade to "cut per region"
// 03.06.03, V.Ivanchenko fix compilation warnings

#ifndef G4FORWARDXRAYTR_H
#define G4FORWARDXRAYTR_H

#include "globals.hh"
#include "G4Track.hh"
#include "G4TransitionRadiation.hh"
#include "G4VParticleChange.hh"

class G4ParticleDefinition;
class G4PhysicsTable;
class G4PhysicsLogVector;

class G4ForwardXrayTR : public G4TransitionRadiation
{
 public:
  explicit G4ForwardXrayTR(const G4String& matName1, const G4String& matName2,
                           const G4String& processName = "XrayTR");

  explicit G4ForwardXrayTR(const G4String& processName = "XrayTR");

  ~G4ForwardXrayTR();

  G4ForwardXrayTR(const G4ForwardXrayTR& right) = delete;
  G4ForwardXrayTR& operator=(const G4ForwardXrayTR& right) = delete;

  ///////////////////////    Methods    /////////////////////////////////

  void ProcessDescription(std::ostream&) const override;
  void DumpInfo() const override { ProcessDescription(G4cout); };

  void BuildXrayTRtables();

  G4double GetMeanFreePath(const G4Track&, G4double,
                           G4ForceCondition* condition) override;

  G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                  const G4Step& aStep) override;

  G4double GetEnergyTR(G4int iMat, G4int jMat, G4int iTkin) const;

  G4double GetThetaTR(G4int iMat, G4int jMat, G4int iTkin) const;

  ///////////////////// Angle distribution  /////////////////////////////

  G4double SpectralAngleTRdensity(G4double energy,
                                  G4double varAngle) const override;

  G4double AngleDensity(G4double energy, G4double varAngle) const;

  G4double EnergyInterval(G4double energy1, G4double energy2,
                          G4double varAngle) const;

  G4double AngleSum(G4double varAngle1, G4double varAngle2) const;

  /////////////////////////  Energy distribution ///////////////////////////////

  G4double SpectralDensity(G4double energy, G4double x) const;

  G4double AngleInterval(G4double energy, G4double varAngle1,
                         G4double varAngle2) const;

  G4double EnergySum(G4double energy1, G4double energy2) const;

  ///////////////////////////   Access functions  ////////////////////////////

  G4PhysicsTable* GetAngleDistrTable();
  G4PhysicsTable* GetEnergyDistrTable();

  static G4int GetSympsonNumber();
  static G4int GetBinTR();

  static G4double GetMinProtonTkin();
  static G4double GetMaxProtonTkin();
  static G4int GetTotBin();

 protected:  // for access from X-ray TR fast simulation models
  static constexpr G4double fTheMinEnergyTR =
    1. * CLHEP::keV;  //  static min TR energy
  static constexpr G4double fTheMaxEnergyTR =
    100. * CLHEP::keV;                              //  static max TR energy
  static constexpr G4double fTheMaxAngle = 1.0e-3;  //  max theta of TR quanta
  static constexpr G4double fTheMinAngle = 5.0e-6;  //  min theta of TR quanta
  static constexpr G4double fMinProtonTkin =
    100. * CLHEP::GeV;  // min Tkin of proton in tables
  static constexpr G4double fMaxProtonTkin =
    100. * CLHEP::TeV;  // max Tkin of proton in tables
  static constexpr G4double fPlasmaCof =
    4.0 * CLHEP::pi * CLHEP::fine_structure_const * CLHEP::hbarc *
    CLHEP::hbarc * CLHEP::hbarc /
    CLHEP::electron_mass_c2;  // physical consts for plasma energy
  static constexpr G4double fCofTR = CLHEP::fine_structure_const / CLHEP::pi;

  static constexpr G4int fSympsonNumber =
    100;                                // Accuracy of Sympson integration
  static constexpr G4int fBinTR  = 50;  //  number of bins in TR vectors
  static constexpr G4int fTotBin = 50;  // number of bins in log scale

  const std::vector<G4double>* fGammaCutInKineticEnergy;
  // TR photon cut in energy array

  G4ParticleDefinition* fPtrGamma;  // pointer to TR photon

  G4PhysicsTable* fAngleDistrTable;
  G4PhysicsTable* fEnergyDistrTable;

  G4PhysicsLogVector* fProtonEnergyVector;

  G4double fMinEnergyTR;   //  min TR energy in material
  G4double fMaxEnergyTR;   //  max TR energy in material
  G4double fMaxThetaTR;    //  max theta of TR quanta
  G4double fGamma;         // current Lorentz factor
  G4double fGammaTkinCut;  // Tkin cut of TR photon in current mat.
  G4double fSigma1;        // plasma energy Sq of matter1
  G4double fSigma2;        // plasma energy Sq of matter2

  G4int secID = -1;  // creator modelID
};

#endif  // G4FORWARDXRAYTR_H
