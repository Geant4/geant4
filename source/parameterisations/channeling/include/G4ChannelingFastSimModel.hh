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
// Author:      Alexei Sytov
// Co-author:   Gianfranco Paternò (modifications & testing)
// On the base of the CRYSTALRAD realization of channeling model:
// A. I. Sytov, V. V. Tikhomirov, and L. Bandiera PRAB 22, 064601 (2019)

#ifndef G4ChannelingFastSimModel_h
#define G4ChannelingFastSimModel_h 1

#include "G4VFastSimulationModel.hh"
#include "G4Step.hh"
#include "G4TouchableHandle.hh"
#include <vector>
#include <CLHEP/Units/SystemOfUnits.h>
#include <CLHEP/Units/PhysicalConstants.h>

#include "G4ChannelingFastSimCrystalData.hh"
#include <unordered_map>
#include "G4BaierKatkov.hh"
#include "G4LogicalVolume.hh"
#include "G4ParticleTable.hh"

/** \file G4ChannelingFastSimModel.hh
* \brief Definition of the G4ChannelingFastSimModel class
* FastSimulation Channeling model: calculates charge particle trajectories
* in oriented crystals in the field of crystal planes/axes either straight or bent.
* It is also possible to simulate radiation using Baier-Katkov method.
*/

class G4ChannelingFastSimModel : public G4VFastSimulationModel
{
public:
  // Constructor, destructor
  G4ChannelingFastSimModel (const G4String&, G4Region*);
  G4ChannelingFastSimModel (const G4String&);
  ~G4ChannelingFastSimModel ();

  /// -- IsApplicable
  G4bool IsApplicable(const G4ParticleDefinition&) override;
  /// -- ModelTrigger
  G4bool ModelTrigger(const G4FastTrack &) override;
  /// -- User method DoIt
  void DoIt(const G4FastTrack&, G4FastStep&) override;

  ///special functions
  void Input(const G4Material* crystal,
             const G4String &lattice)
            {Input(crystal,lattice,"");}

  void Input(const G4Material* crystal,
             const G4String &lattice,
             const G4String &filePath);

  void RadiationModelActivate();

  G4ChannelingFastSimCrystalData* GetCrystalData() {return fCrystalData;}

  G4BaierKatkov* GetRadiationModel() {return fBaierKatkov;}

  G4bool GetIfRadiationModelActive(){return fRad;}

  ///set cuts
  void SetLowKineticEnergyLimit(G4double ekinetic, const G4String& particleName)
   {fLowEnergyLimit[particleTable->FindParticle(particleName)->
              GetParticleDefinitionID()] = ekinetic;}
  void SetLindhardAngleNumberHighLimit(G4double angleNumber, const G4String& particleName)
   {fLindhardAngleNumberHighLimit[particleTable->FindParticle(particleName)->
              GetParticleDefinitionID()]=angleNumber;}
  void SetHighAngleLimit(G4double anglemax, const G4String& particleName)
   {fHighAngleLimit[particleTable->FindParticle(particleName)->
                       GetParticleDefinitionID()] = anglemax;}

  void SetDefaultLowKineticEnergyLimit(G4double ekinetic)
           {fDefaultLowEnergyLimit=ekinetic;}
  void SetDefaultLindhardAngleNumberHighLimit(G4double angleNumber)
           {fDefaultLindhardAngleNumberHighLimit=angleNumber;}
  void SetDefaultHighAngleLimit(G4double anglemax)
           {fDefaultHighAngleLimit=anglemax;}


  /// get the maximal number of photons that can be produced per fastStep
  /// Caution: is redundant, if the radiation model is not activated
  void SetMaxPhotonsProducedPerStep(G4double nPhotons)
           {fMaxPhotonsProducedPerStep=nPhotons;}

  ///get cuts
  G4double GetLowKineticEnergyLimit(G4int particleDefinitionID)
               {return (fLowEnergyLimit.count(particleDefinitionID) == 1)
                        ? fLowEnergyLimit[particleDefinitionID]
                        : fDefaultLowEnergyLimit;}
  G4double GetLindhardAngleNumberHighLimit(G4int particleDefinitionID)
               {return (fLindhardAngleNumberHighLimit.count(particleDefinitionID) == 1)
                        ? fLindhardAngleNumberHighLimit[particleDefinitionID]
                        : fDefaultLindhardAngleNumberHighLimit;}
  G4double GetHighAngleLimit(G4int particleDefinitionID)
               {return (fHighAngleLimit.count(particleDefinitionID) == 1)
                              ? fHighAngleLimit[particleDefinitionID]
                              : fDefaultHighAngleLimit;}

  /// get the maximal number of photons that can be produced per fastStep
  G4int GetMaxPhotonsProducedPerStep(){return fMaxPhotonsProducedPerStep;}

private:

  G4ChannelingFastSimCrystalData* fCrystalData{nullptr};
  G4BaierKatkov* fBaierKatkov{nullptr};

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  ///flag of radiation model
  G4bool fRad = false;

  /// maps of cuts (angular cuts are chosen as std::max of
  /// fHighAngleLimit and calculated Lindhard angle)
  std::unordered_map<G4int, G4double> fLowEnergyLimit;
  std::unordered_map<G4int, G4double> fLindhardAngleNumberHighLimit;
  std::unordered_map<G4int, G4double> fHighAngleLimit;

  G4double fDefaultLowEnergyLimit = 200*CLHEP::MeV;
  G4double fDefaultLindhardAngleNumberHighLimit = 100.;
  G4double fDefaultHighAngleLimit = 0.;

  /// the maximal number of photons that can be produced per fastStep
  G4int fMaxPhotonsProducedPerStep=1000.;

};
#endif




