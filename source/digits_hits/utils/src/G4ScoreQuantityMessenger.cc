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
//
// ---------------------------------------------------------------------
// Modifications
// 08-Oct-2010 T.Aso remove unit of G4PSPassageCellCurrent.
//  24-Mar-2011  T.Aso  Add StepChecker for debugging.
//  24-Mar-2011  T.Aso  Size and segmentation for replicated cylinder.
//  01-Jun-2012  T.Aso  Support weighted/dividedByArea options
//                      in flatCurrent and flatFulx commands.
//  27-Mar-2013  T.Aso  Unit option in the kineticEnergy filter was
//                     supported.
//
// ---------------------------------------------------------------------

#include "G4ScoreQuantityMessenger.hh"
#include "G4ScoringManager.hh"
#include "G4VScoringMesh.hh"
#include "G4VPrimitiveScorer.hh"

#include "G4PSCellCharge3D.hh"
#include "G4PSCellFlux3D.hh"
#include "G4PSCellFluxForCylinder3D.hh"
#include "G4PSPassageCellFlux3D.hh"
#include "G4PSPassageCellFluxForCylinder3D.hh"
#include "G4PSEnergyDeposit3D.hh"
#include "G4PSDoseDeposit3D.hh"
#include "G4PSDoseDepositForCylinder3D.hh"
#include "G4PSNofStep3D.hh"
#include "G4PSNofSecondary3D.hh"
//
#include "G4PSTrackLength3D.hh"
#include "G4PSPassageCellCurrent3D.hh"
#include "G4PSPassageTrackLength3D.hh"
#include "G4PSFlatSurfaceCurrent3D.hh"
#include "G4PSFlatSurfaceFlux3D.hh"
#include "G4PSSphereSurfaceCurrent3D.hh"
#include "G4PSSphereSurfaceFlux3D.hh"
#include "G4PSCylinderSurfaceCurrent3D.hh"
#include "G4PSCylinderSurfaceFlux3D.hh"
#include "G4PSVolumeFlux3D.hh"
#include "G4PSNofCollision3D.hh"
#include "G4PSPopulation3D.hh"
#include "G4PSTrackCounter3D.hh"
#include "G4PSTermination3D.hh"
#include "G4PSMinKinEAtGeneration3D.hh"

#include "G4PSCellCharge.hh"
#include "G4PSCellFlux.hh"
#include "G4PSPassageCellFlux.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSDoseDeposit.hh"
#include "G4PSNofStep.hh"
#include "G4PSNofSecondary.hh"
//
#include "G4PSTrackLength.hh"
#include "G4PSPassageCellCurrent.hh"
#include "G4PSPassageTrackLength.hh"
#include "G4PSFlatSurfaceCurrent.hh"
#include "G4PSFlatSurfaceFlux.hh"
#include "G4PSSphereSurfaceCurrent.hh"
#include "G4PSSphereSurfaceFlux.hh"
#include "G4PSCylinderSurfaceCurrent.hh"
#include "G4PSCylinderSurfaceFlux.hh"
#include "G4PSNofCollision.hh"
#include "G4PSPopulation.hh"
#include "G4PSTrackCounter.hh"
#include "G4PSTermination.hh"
#include "G4PSMinKinEAtGeneration.hh"

//
// For debug purpose
#include "G4PSStepChecker3D.hh"

#include "G4SDChargedFilter.hh"
#include "G4SDNeutralFilter.hh"
#include "G4SDKineticEnergyFilter.hh"
#include "G4SDParticleFilter.hh"
#include "G4SDParticleWithEnergyFilter.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4Tokenizer.hh"
#include "G4UnitsTable.hh"

G4ScoreQuantityMessenger::G4ScoreQuantityMessenger(G4ScoringManager* SManager)
  : fSMan(SManager)
{
  QuantityCommands();
  FilterCommands();
}

void G4ScoreQuantityMessenger::QuantityCommands()
{
  G4UIparameter* param;

  //
  // Quantity commands
  quantityDir = new G4UIdirectory("/score/quantity/");
  quantityDir->SetGuidance("Scoring quantity of the mesh.");
  //
  qTouchCmd = new G4UIcmdWithAString("/score/quantity/touch", this);
  qTouchCmd->SetGuidance(
    "Assign previously defined quantity to the current quantity.");
  qTouchCmd->SetParameterName("qname", false);
  //
  qGetUnitCmd = new G4UIcmdWithoutParameter("/score/quantity/get/unit", this);
  qGetUnitCmd->SetGuidance("Print output unit of the current quantity.");
  //
  qSetUnitCmd = new G4UIcmdWithAString("/score/quantity/set/unit", this);
  qSetUnitCmd->SetGuidance("Set output unit of the current quantity.");
  qSetUnitCmd->SetParameterName("unit", false);

  // Primitive Scorers
  qeDepCmd = new G4UIcommand("/score/quantity/energyDeposit", this);
  qeDepCmd->SetGuidance("Energy deposit scorer.");
  qeDepCmd->SetGuidance("[usage] /score/quantity/energyDeposit qname unit");
  qeDepCmd->SetGuidance("  qname  :(String) scorer name");
  qeDepCmd->SetGuidance("  unit   :(String) unit");
  param = new G4UIparameter("qname", 's', false);
  qeDepCmd->SetParameter(param);
  param = new G4UIparameter("unit", 's', true);
  param->SetDefaultUnit("MeV");
  qeDepCmd->SetParameter(param);
  //
  qCellChgCmd = new G4UIcommand("/score/quantity/cellCharge", this);
  qCellChgCmd->SetGuidance("Cell charge scorer.");
  qCellChgCmd->SetGuidance("[usage] /score/quantity/cellCharge qname unit");
  qCellChgCmd->SetGuidance("  qname  :(String) scorer name");
  qCellChgCmd->SetGuidance("  unit   :(String) unit");
  param = new G4UIparameter("qname", 's', false);
  qCellChgCmd->SetParameter(param);
  param = new G4UIparameter("unit", 's', true);
  param->SetDefaultUnit("e+");
  qCellChgCmd->SetParameter(param);
  //
  qCellFluxCmd = new G4UIcommand("/score/quantity/cellFlux", this);
  qCellFluxCmd->SetGuidance("Cell flux scorer.");
  qCellFluxCmd->SetGuidance("[usage] /score/quantity/cellFlux qname unit");
  qCellFluxCmd->SetGuidance("  qname  :(String) scorer name");
  qCellFluxCmd->SetGuidance("  unit   :(String) unit");
  param = new G4UIparameter("qname", 's', false);
  qCellFluxCmd->SetParameter(param);
  param = new G4UIparameter("unit", 's', true);
  param->SetDefaultValue("percm2");
  qCellFluxCmd->SetParameter(param);
  //
  qPassCellFluxCmd = new G4UIcommand("/score/quantity/passageCellFlux", this);
  qPassCellFluxCmd->SetGuidance("Passage cell flux scorer");
  qPassCellFluxCmd->SetGuidance(
    "[usage] /score/quantity/passageCellFlux qname unit");
  qPassCellFluxCmd->SetGuidance("  qname  :(String) scorer name");
  qPassCellFluxCmd->SetGuidance("  unit   :(String) unit");
  param = new G4UIparameter("qname", 's', false);
  qPassCellFluxCmd->SetParameter(param);
  param = new G4UIparameter("unit", 's', true);
  param->SetDefaultValue("percm2");
  qPassCellFluxCmd->SetParameter(param);
  //
  qdoseDepCmd = new G4UIcommand("/score/quantity/doseDeposit", this);
  qdoseDepCmd->SetGuidance("Dose deposit scorer.");
  qdoseDepCmd->SetGuidance("[usage] /score/quantity/doseDeposit qname unit");
  qdoseDepCmd->SetGuidance("  qname  :(String) scorer name");
  qdoseDepCmd->SetGuidance("  unit   :(String) unit");
  param = new G4UIparameter("qname", 's', false);
  qdoseDepCmd->SetParameter(param);
  param = new G4UIparameter("unit", 's', true);
  param->SetDefaultUnit("Gy");
  qdoseDepCmd->SetParameter(param);
  //
  qnOfStepCmd = new G4UIcommand("/score/quantity/nOfStep", this);
  qnOfStepCmd->SetGuidance("Number of step scorer.");
  qnOfStepCmd->SetGuidance("[usage] /score/quantity/nOfStep qname");
  qnOfStepCmd->SetGuidance("[usage] /score/quantity/nOfStep qname  bflag");
  qnOfStepCmd->SetGuidance("  qname  :(String) scorer name");
  qnOfStepCmd->SetGuidance("  bflag  :(Bool) Skip zero step ");
  qnOfStepCmd->SetGuidance("          at geometry boundary if true");
  param = new G4UIparameter("qname", 's', false);
  qnOfStepCmd->SetParameter(param);
  param = new G4UIparameter("bflag", 'b', true);
  param->SetDefaultValue("false");
  qnOfStepCmd->SetParameter(param);
  //
  qnOfSecondaryCmd = new G4UIcommand("/score/quantity/nOfSecondary", this);
  qnOfSecondaryCmd->SetGuidance("Number of secondary scorer.");
  qnOfSecondaryCmd->SetGuidance("[usage] /score/quantity/nOfSecondary qname");
  qnOfSecondaryCmd->SetGuidance("  qname  :(String) scorer name");
  param = new G4UIparameter("qname", 's', false);
  qnOfSecondaryCmd->SetParameter(param);
  //
  qTrackLengthCmd = new G4UIcommand("/score/quantity/trackLength", this);
  qTrackLengthCmd->SetGuidance("Track length scorer.");
  qTrackLengthCmd->SetGuidance(
    "[usage] /score/quantity/trackLength qname wflag kflag vflag unit");
  qTrackLengthCmd->SetGuidance("  qname  :(String) scorer name");
  qTrackLengthCmd->SetGuidance("  wflag  :(Bool) weighted");
  qTrackLengthCmd->SetGuidance("  kflag  :(Bool) multiply kinetic energy");
  qTrackLengthCmd->SetGuidance("  vflag  :(Bool) divide by velocity");
  qTrackLengthCmd->SetGuidance("  unit   :(String) unit");
  param = new G4UIparameter("qname", 's', false);
  qTrackLengthCmd->SetParameter(param);
  param = new G4UIparameter("wflag", 'b', true);
  param->SetDefaultValue("false");
  qTrackLengthCmd->SetParameter(param);
  param = new G4UIparameter("kflag", 'b', true);
  param->SetDefaultValue("false");
  qTrackLengthCmd->SetParameter(param);
  param = new G4UIparameter("vflag", 'b', true);
  param->SetDefaultValue("false");
  qTrackLengthCmd->SetParameter(param);
  param = new G4UIparameter("unit", 's', true);
  param->SetDefaultValue("mm");
  qTrackLengthCmd->SetParameter(param);
  //
  qPassCellCurrCmd =
    new G4UIcommand("/score/quantity/passageCellCurrent", this);
  qPassCellCurrCmd->SetGuidance("Passage cell current scorer.");
  qPassCellCurrCmd->SetGuidance(
    "[usage] /score/quantity/passageCellCurrent qname wflag");
  qPassCellCurrCmd->SetGuidance("  qname  :(String) scorer name");
  qPassCellCurrCmd->SetGuidance("  wflag  :(Bool) weighted");
  param = new G4UIparameter("qname", 's', false);
  qPassCellCurrCmd->SetParameter(param);
  param = new G4UIparameter("wflag", 'b', true);
  param->SetDefaultValue("true");
  qPassCellCurrCmd->SetParameter(param);
  //
  qPassTrackLengthCmd =
    new G4UIcommand("/score/quantity/passageTrackLength", this);
  qPassTrackLengthCmd->SetGuidance("Passage track length scorer.");
  qPassTrackLengthCmd->SetGuidance(
    "[usage] /score/quantity/passageTrackLength qname wflag unit");
  qPassTrackLengthCmd->SetGuidance("  qname  :(String) scorer name");
  qPassTrackLengthCmd->SetGuidance("  wflag  :(Bool) weighted");
  qPassTrackLengthCmd->SetGuidance("  unit   :(Bool) unit");
  param = new G4UIparameter("qname", 's', false);
  qPassTrackLengthCmd->SetParameter(param);
  param = new G4UIparameter("wflag", 'b', true);
  param->SetDefaultValue("true");
  qPassTrackLengthCmd->SetParameter(param);
  param = new G4UIparameter("unit", 's', true);
  param->SetDefaultUnit("mm");
  qPassTrackLengthCmd->SetParameter(param);
  //
  qFlatSurfCurrCmd =
    new G4UIcommand("/score/quantity/flatSurfaceCurrent", this);
  qFlatSurfCurrCmd->SetGuidance("Flat surface current Scorer.");
  qFlatSurfCurrCmd->SetGuidance(
    "[usage] /score/quantity/flatSurfaceCurrent qname dflag wflag aflag unit");
  qFlatSurfCurrCmd->SetGuidance("  qname  :(String) scorer name");
  qFlatSurfCurrCmd->SetGuidance("  dflag  :(Int) direction flag");
  qFlatSurfCurrCmd->SetGuidance("         : 0 = Both In and Out");
  qFlatSurfCurrCmd->SetGuidance("         : 1 = In only");
  qFlatSurfCurrCmd->SetGuidance("         : 2 = Out only");
  qFlatSurfCurrCmd->SetGuidance("  wflag  :(Bool) weighted");
  qFlatSurfCurrCmd->SetGuidance("  aflag  :(Bool) divide by area");
  qFlatSurfCurrCmd->SetGuidance("  unit   :(String) unit");
  param = new G4UIparameter("qname", 's', false);
  qFlatSurfCurrCmd->SetParameter(param);
  param = new G4UIparameter("dflag", 'i', true);
  param->SetDefaultValue("0");
  qFlatSurfCurrCmd->SetParameter(param);
  param = new G4UIparameter("wflag", 'b', true);
  param->SetDefaultValue("true");
  qFlatSurfCurrCmd->SetParameter(param);
  param = new G4UIparameter("aflag", 'b', true);
  param->SetDefaultValue("true");
  qFlatSurfCurrCmd->SetParameter(param);
  param = new G4UIparameter("unit", 's', true);
  param->SetDefaultValue("percm2");
  qFlatSurfCurrCmd->SetParameter(param);
  //
  qFlatSurfFluxCmd = new G4UIcommand("/score/quantity/flatSurfaceFlux", this);
  qFlatSurfFluxCmd->SetGuidance("Flat surface flux scorer.");
  qFlatSurfFluxCmd->SetGuidance(
    "[usage] /score/quantity/flatSurfaceFlux qname dflag unit");
  qFlatSurfFluxCmd->SetGuidance("  qname  :(String) scorer name");
  qFlatSurfFluxCmd->SetGuidance("  dflag  :(Int) direction flag");
  qFlatSurfFluxCmd->SetGuidance("         : 0 = Both In and Out");
  qFlatSurfFluxCmd->SetGuidance("         : 1 = In only");
  qFlatSurfFluxCmd->SetGuidance("         : 2 = Out only");
  qFlatSurfFluxCmd->SetGuidance("  wflag  :(Bool) weighted");
  qFlatSurfFluxCmd->SetGuidance("  aflag  :(Bool) divide by area");
  qFlatSurfFluxCmd->SetGuidance("  unit   :(String) unit");
  param = new G4UIparameter("qname", 's', false);
  qFlatSurfFluxCmd->SetParameter(param);
  param = new G4UIparameter("dflag", 'i', true);
  param->SetDefaultValue("0");
  qFlatSurfFluxCmd->SetParameter(param);
  param = new G4UIparameter("wflag", 'b', true);
  param->SetDefaultValue("true");
  qFlatSurfFluxCmd->SetParameter(param);
  param = new G4UIparameter("aflag", 'b', true);
  param->SetDefaultValue("true");
  qFlatSurfFluxCmd->SetParameter(param);
  param = new G4UIparameter("unit", 's', true);
  param->SetDefaultValue("percm2");
  qFlatSurfFluxCmd->SetParameter(param);
  //

  qVolFluxCmd = new G4UIcommand("/score/quantity/volumeFlux", this);
  qVolFluxCmd->SetGuidance("Volume flux scorer.");
  qVolFluxCmd->SetGuidance(
    "This scorer scores the number of particles getting into the volume "
    "without normalized by the surface area.");
  qVolFluxCmd->SetGuidance(
    "[usage] /score/quantity/volumeFlux qname divcos dflag");
  qVolFluxCmd->SetGuidance("  qname  :(String) scorer name");
  qVolFluxCmd->SetGuidance("  divcos :(Bool) divide by cos(theta), where theta "
                           "is the incident angle (default : false)");
  qVolFluxCmd->SetGuidance(
    "  dflag  :(Int) direction, 1 : inward (default), 2 : outward");
  param = new G4UIparameter("qname", 's', false);
  qVolFluxCmd->SetParameter(param);
  param = new G4UIparameter("divcos", 'b', true);
  param->SetDefaultValue(0);
  qVolFluxCmd->SetParameter(param);
  param = new G4UIparameter("dflag", 'i', true);
  param->SetParameterRange("dflag>=1 && dflag<=2");
  param->SetDefaultValue(1);
  qVolFluxCmd->SetParameter(param);

  qNofCollisionCmd = new G4UIcommand("/score/quantity/nOfCollision", this);
  qNofCollisionCmd->SetGuidance("Number of collision scorer.");
  qNofCollisionCmd->SetGuidance(
    "[usage] /score/quantity/nOfCollision qname wflag");
  qNofCollisionCmd->SetGuidance("  qname  :(String) scorer name");
  param = new G4UIparameter("qname", 's', false);
  qNofCollisionCmd->SetParameter(param);
  param = new G4UIparameter("wflag", 'b', true);
  param->SetDefaultValue("false");
  qNofCollisionCmd->SetParameter(param);
  //
  qPopulationCmd = new G4UIcommand("/score/quantity/population", this);
  qPopulationCmd->SetGuidance("Population scorer.");
  qPopulationCmd->SetGuidance("[usage] /score/quantity/population qname wflag");
  qPopulationCmd->SetGuidance("  qname  :(String) scorer name");
  qPopulationCmd->SetGuidance("  wflag  :(Bool) weighted");
  param = new G4UIparameter("qname", 's', false);
  qPopulationCmd->SetParameter(param);
  param = new G4UIparameter("wflag", 'b', true);
  param->SetDefaultValue("false");
  qPopulationCmd->SetParameter(param);

  //
  qTrackCountCmd = new G4UIcommand("/score/quantity/nOfTrack", this);
  qTrackCountCmd->SetGuidance("Number of track scorer.");
  qTrackCountCmd->SetGuidance(
    "[usage] /score/quantity/nOfTrack qname dflag wflag");
  qTrackCountCmd->SetGuidance("  qname  :(String) scorer name");
  qTrackCountCmd->SetGuidance("  dflag  :(Int) direction");
  qTrackCountCmd->SetGuidance("         : 0 = Both In and Out");
  qTrackCountCmd->SetGuidance("         : 1 = In only");
  qTrackCountCmd->SetGuidance("         : 2 = Out only");
  qTrackCountCmd->SetGuidance("  wflag  :(Bool) weighted");
  param = new G4UIparameter("qname", 's', false);
  qTrackCountCmd->SetParameter(param);
  param = new G4UIparameter("dflag", 'i', true);
  param->SetDefaultValue("0");
  qTrackCountCmd->SetParameter(param);
  param = new G4UIparameter("wflag", 'b', true);
  param->SetDefaultValue("false");
  qTrackCountCmd->SetParameter(param);

  //
  qTerminationCmd = new G4UIcommand("/score/quantity/nOfTerminatedTrack", this);
  qTerminationCmd->SetGuidance("Number of terminated tracks scorer.");
  qTerminationCmd->SetGuidance(
    "[usage] /score/quantity/nOfTerminatedTrack qname wflag");
  qTerminationCmd->SetGuidance("  qname  :(String) scorer name");
  qTerminationCmd->SetGuidance("  wflag  :(Bool) weighted");
  param = new G4UIparameter("qname", 's', false);
  qTerminationCmd->SetParameter(param);
  param = new G4UIparameter("wflag", 'b', true);
  param->SetDefaultValue("false");
  qTerminationCmd->SetParameter(param);

  //
  qMinKinEAtGeneCmd =
    new G4UIcommand("/score/quantity/minKinEAtGeneration", this);
  qMinKinEAtGeneCmd->SetGuidance("Min Kinetic Energy at Generation");
  qMinKinEAtGeneCmd->SetGuidance(
    "[usage] /score/quantity/minKinEAtGeneration qname unit");
  qMinKinEAtGeneCmd->SetGuidance("  qname  :(String) scorer name");
  qMinKinEAtGeneCmd->SetGuidance("  unit   :(String) unit name");
  param = new G4UIparameter("qname", 's', false);
  qMinKinEAtGeneCmd->SetParameter(param);
  param = new G4UIparameter("unit", 's', true);
  param->SetDefaultUnit("MeV");
  qMinKinEAtGeneCmd->SetParameter(param);
  //
  qStepCheckerCmd = new G4UIcommand("/score/quantity/stepChecker", this);
  qStepCheckerCmd->SetGuidance("Display a comment when this PS is invoked");
  qStepCheckerCmd->SetGuidance("[usage] /score/quantity/stepChecker qname");
  qStepCheckerCmd->SetGuidance("  qname  :(String) scorer name");
  param = new G4UIparameter("qname", 's', false);
  qStepCheckerCmd->SetParameter(param);
}

void G4ScoreQuantityMessenger::FilterCommands()
{
  G4UIparameter* param;

  //
  // Filter commands
  filterDir = new G4UIdirectory("/score/filter/");
  filterDir->SetGuidance("  Scoring filter commands.");
  //
  fchargedCmd = new G4UIcmdWithAString("/score/filter/charged", this);
  fchargedCmd->SetGuidance("Charged particle filter.");
  fchargedCmd->SetParameterName("fname", false);
  //
  fneutralCmd = new G4UIcmdWithAString("/score/filter/neutral", this);
  fneutralCmd->SetGuidance("Neutral particle filter.");
  fneutralCmd->SetParameterName("fname", false);
  //
  fkinECmd = new G4UIcommand("/score/filter/kineticEnergy", this);
  fkinECmd->SetGuidance("Kinetic energy filter.");
  fkinECmd->SetGuidance(
    "[usage] /score/filter/kineticEnergy fname Elow Ehigh unit");
  fkinECmd->SetGuidance("  fname     :(String) Filter Name ");
  fkinECmd->SetGuidance("  Elow      :(Double) Lower edge of kinetic energy");
  fkinECmd->SetGuidance("  Ehigh     :(Double) Higher edge of kinetic energy");
  fkinECmd->SetGuidance("  unit      :(String) unit of given kinetic energy");
  param = new G4UIparameter("fname", 's', false);
  fkinECmd->SetParameter(param);
  param = new G4UIparameter("elow", 'd', true);
  param->SetDefaultValue("0.0");
  fkinECmd->SetParameter(param);
  param = new G4UIparameter("ehigh", 'd', true);
  fkinECmd->SetParameter(param);
  G4String smax = DtoS(DBL_MAX);
  param->SetDefaultValue(smax);
  param = new G4UIparameter("unit", 's', true);
  param->SetDefaultUnit("keV");
  fkinECmd->SetParameter(param);
  //
  fparticleCmd = new G4UIcommand("/score/filter/particle", this);
  fparticleCmd->SetGuidance("Particle filter.");
  fparticleCmd->SetGuidance("[usage] /score/filter/particle fname p0 .. pn");
  fparticleCmd->SetGuidance("  fname     :(String) Filter Name ");
  fparticleCmd->SetGuidance("  p0 .. pn  :(String) particle names");
  param = new G4UIparameter("fname", 's', false);
  fparticleCmd->SetParameter(param);
  param = new G4UIparameter("particlelist", 's', false);
  param->SetDefaultValue("");
  fparticleCmd->SetParameter(param);
  //
  //
  //
  fparticleKinECmd =
    new G4UIcommand("/score/filter/particleWithKineticEnergy", this);
  fparticleKinECmd->SetGuidance("Particle with kinetic energy filter.");
  fparticleKinECmd->SetGuidance(
    "[usage] /score/filter/particleWithKineticEnergy fname Elow Ehigh unit p0 "
    ".. pn");
  fparticleKinECmd->SetGuidance("  fname     :(String) Filter Name ");
  fparticleKinECmd->SetGuidance(
    "  Elow      :(Double) Lower edge of kinetic energy");
  fparticleKinECmd->SetGuidance(
    "  Ehigh     :(Double) Higher edge of kinetic energy");
  fparticleKinECmd->SetGuidance(
    "  unit      :(String) unit of given kinetic energy");
  fparticleKinECmd->SetGuidance("  p0 .. pn  :(String) particle names");
  param = new G4UIparameter("fname", 's', false);
  fparticleKinECmd->SetParameter(param);
  param = new G4UIparameter("elow", 'd', false);
  param->SetDefaultValue("0.0");
  fparticleKinECmd->SetParameter(param);
  param = new G4UIparameter("ehigh", 'd', true);
  param->SetDefaultValue(smax);
  fparticleKinECmd->SetParameter(param);
  param = new G4UIparameter("unit", 's', true);
  param->SetDefaultUnit("keV");
  fparticleKinECmd->SetParameter(param);
  param = new G4UIparameter("particlelist", 's', false);
  param->SetDefaultValue("");
  fparticleKinECmd->SetParameter(param);
  //
  //
}

G4ScoreQuantityMessenger::~G4ScoreQuantityMessenger()
{
  delete quantityDir;
  delete qTouchCmd;
  delete qGetUnitCmd;
  delete qSetUnitCmd;

  //
  delete qCellChgCmd;
  delete qCellFluxCmd;
  delete qPassCellFluxCmd;
  delete qeDepCmd;
  delete qdoseDepCmd;
  delete qnOfStepCmd;
  delete qnOfSecondaryCmd;
  //
  delete qTrackLengthCmd;
  delete qPassCellCurrCmd;
  delete qPassTrackLengthCmd;
  delete qFlatSurfCurrCmd;
  delete qFlatSurfFluxCmd;
  delete qVolFluxCmd;

  delete qNofCollisionCmd;
  delete qPopulationCmd;
  delete qTrackCountCmd;
  delete qTerminationCmd;
  delete qMinKinEAtGeneCmd;
  //
  delete qStepCheckerCmd;
  //
  delete filterDir;
  delete fchargedCmd;
  delete fneutralCmd;
  delete fkinECmd;
  delete fparticleCmd;
  delete fparticleKinECmd;
}

void G4ScoreQuantityMessenger::SetNewValue(G4UIcommand* command,
                                           G4String newVal)
{
  using MeshShape = G4VScoringMesh::MeshShape;

  G4ExceptionDescription ed;

  //
  // Get Current Mesh
  //
  G4VScoringMesh* mesh = fSMan->GetCurrentMesh();
  if(mesh == nullptr)
  {
    ed << "ERROR : No mesh is currently open. Open/create a mesh first. "
          "Command ignored.";
    command->CommandFailed(ed);
    return;
  }
  // Mesh type
  auto shape = mesh->GetShape();
  // Tokens
  G4TokenVec token;
  FillTokenVec(newVal, token);
  //
  // Commands for Current Mesh
  if(command == qTouchCmd)
  {
    mesh->SetCurrentPrimitiveScorer(newVal);
  }
  else if(command == qGetUnitCmd)
  {
    G4cout << "Unit:  " << mesh->GetCurrentPSUnit() << G4endl;
  }
  else if(command == qSetUnitCmd)
  {
    mesh->SetCurrentPSUnit(newVal);
  }
  else if(command == qCellChgCmd)
  {
    if(CheckMeshPS(mesh, token[0], command))
    {
      G4PSCellCharge* ps = nullptr;
      if(shape == MeshShape::realWorldLogVol || shape == MeshShape::probe)
      {
        ps = new G4PSCellCharge(token[0], mesh->GetCopyNumberLevel());
      }
      else
      {
        ps = new G4PSCellCharge3D(token[0]);
      }
      ps->SetUnit(token[1]);
      mesh->SetPrimitiveScorer(ps);
    }
  }
  else if(command == qCellFluxCmd)
  {
    if(CheckMeshPS(mesh, token[0], command))
    {
      G4PSCellFlux* ps = nullptr;
      if(shape == MeshShape::box)
      {
        ps = new G4PSCellFlux3D(token[0]);
      }
      else if(shape == MeshShape::cylinder)
      {
        auto  pps =
          new G4PSCellFluxForCylinder3D(token[0]);
        pps->SetCylinderSize(mesh->GetSize(),mesh->GetStartAngle(),mesh->GetAngleSpan());   
        G4int nSeg[3];
        mesh->GetNumberOfSegments(nSeg);
        pps->SetNumberOfSegments(nSeg);
        ps = pps;
      }
      else if(shape == MeshShape::realWorldLogVol)
      {
        ed << "Cell flux for real world volume is not yet supported. Command "
              "ignored.";
        command->CommandFailed(ed);
        return;
      }
      else if(shape == MeshShape::probe)
      {
        ps = new G4PSCellFlux(token[0]);
      }
      ps->SetUnit(token[1]);
      mesh->SetPrimitiveScorer(ps);
    }
  }
  else if(command == qPassCellFluxCmd)
  {
    if(CheckMeshPS(mesh, token[0], command))
    {
      G4PSPassageCellFlux* ps = nullptr;
      if(shape == MeshShape::box)
      {
        ps = new G4PSPassageCellFlux3D(token[0]);
      }
      else if(shape == MeshShape::cylinder)
      {
        auto  pps =
          new G4PSPassageCellFluxForCylinder3D(token[0]);
        pps->SetCylinderSize(mesh->GetSize(),mesh->GetStartAngle(),mesh->GetAngleSpan());   
        G4int nSeg[3];
        mesh->GetNumberOfSegments(nSeg);
        pps->SetNumberOfSegments(nSeg);
        ps = pps;
      }
      else if(shape == MeshShape::realWorldLogVol)
      {
        ed << "Passing cell flux for real world volume is not yet supported. "
              "Command ignored.";
        command->CommandFailed(ed);
        return;
      }
      else if(shape == MeshShape::probe)
      {
        ps = new G4PSPassageCellFlux(token[0]);
      }
      ps->SetUnit(token[1]);
      mesh->SetPrimitiveScorer(ps);
    }
  }
  else if(command == qeDepCmd)
  {
    if(CheckMeshPS(mesh, token[0], command))
    {
      G4PSEnergyDeposit* ps = nullptr;
      if(shape == MeshShape::realWorldLogVol || shape == MeshShape::probe)
      {
        ps = new G4PSEnergyDeposit(token[0], mesh->GetCopyNumberLevel());
      }
      else
      {
        ps = new G4PSEnergyDeposit3D(token[0]);
      }
      ps->SetUnit(token[1]);
      mesh->SetPrimitiveScorer(ps);
    }
  }
  else if(command == qdoseDepCmd)
  {
    if(CheckMeshPS(mesh, token[0], command))
    {
      G4PSDoseDeposit* ps = nullptr;
      if(shape == MeshShape::box)
      {
        ps = new G4PSDoseDeposit3D(token[0]);
      }
      else if(shape == MeshShape::cylinder)
      {
        auto  pps =
          new G4PSDoseDepositForCylinder3D(token[0]);
        pps->SetUnit(token[1]);
        pps->SetCylinderSize(mesh->GetSize(),mesh->GetStartAngle(),mesh->GetAngleSpan());   
        G4int nSeg[3];
        mesh->GetNumberOfSegments(nSeg);
        pps->SetNumberOfSegments(nSeg);
        ps = pps;
      }
      else if(shape == MeshShape::realWorldLogVol || shape == MeshShape::probe)
      {
        ps = new G4PSDoseDeposit(token[0], mesh->GetCopyNumberLevel());
      }
      ps->SetUnit(token[1]);
      mesh->SetPrimitiveScorer(ps);
    }
  }
  else if(command == qnOfStepCmd)
  {
    if(CheckMeshPS(mesh, token[0], command))
    {
      G4PSNofStep* ps = nullptr;
      if(shape == MeshShape::realWorldLogVol || shape == MeshShape::probe)
      {
        ps = new G4PSNofStep(token[0], mesh->GetCopyNumberLevel());
      }
      else
      {
        ps = new G4PSNofStep3D(token[0]);
      }
      ps->SetBoundaryFlag(StoB(token[1]));
      mesh->SetPrimitiveScorer(ps);
    }
  }
  else if(command == qnOfSecondaryCmd)
  {
    if(CheckMeshPS(mesh, token[0], command))
    {
      G4PSNofSecondary* ps = nullptr;
      if(shape == MeshShape::realWorldLogVol || shape == MeshShape::probe)
      {
        ps = new G4PSNofSecondary(token[0], mesh->GetCopyNumberLevel());
      }
      else
      {
        ps = new G4PSNofSecondary3D(token[0]);
      }
      mesh->SetPrimitiveScorer(ps);
    }
  }
  else if(command == qTrackLengthCmd)
  {
    if(CheckMeshPS(mesh, token[0], command))
    {
      G4PSTrackLength* ps = nullptr;
      if(shape == MeshShape::realWorldLogVol || shape == MeshShape::probe)
      {
        ps = new G4PSTrackLength(token[0], mesh->GetCopyNumberLevel());
      }
      else
      {
        ps = new G4PSTrackLength3D(token[0]);
      }
      ps->Weighted(StoB(token[1]));
      ps->MultiplyKineticEnergy(StoB(token[2]));
      ps->DivideByVelocity(StoB(token[3]));
      ps->SetUnit(token[4]);
      mesh->SetPrimitiveScorer(ps);
    }
  }
  else if(command == qPassCellCurrCmd)
  {
    if(CheckMeshPS(mesh, token[0], command))
    {
      G4PSPassageCellCurrent* ps = nullptr;
      if(shape == MeshShape::realWorldLogVol || shape == MeshShape::probe)
      {
        ps = new G4PSPassageCellCurrent(token[0], mesh->GetCopyNumberLevel());
      }
      else
      {
        ps = new G4PSPassageCellCurrent3D(token[0]);
      }
      ps->Weighted(StoB(token[1]));
      mesh->SetPrimitiveScorer(ps);
    }
  }
  else if(command == qPassTrackLengthCmd)
  {
    if(CheckMeshPS(mesh, token[0], command))
    {
      G4PSPassageTrackLength* ps = nullptr;
      if(shape == MeshShape::realWorldLogVol || shape == MeshShape::probe)
      {
        ps = new G4PSPassageTrackLength(token[0], mesh->GetCopyNumberLevel());
      }
      else
      {
        ps = new G4PSPassageTrackLength3D(token[0]);
      }
      ps->Weighted(StoB(token[1]));
      ps->SetUnit(token[2]);
      mesh->SetPrimitiveScorer(ps);
    }
  }
  else if(command == qFlatSurfCurrCmd)
  {
    if(CheckMeshPS(mesh, token[0], command))
    {
      G4PSFlatSurfaceCurrent* ps = nullptr;
      if(shape == MeshShape::realWorldLogVol || shape == MeshShape::probe)
      {
        ps = new G4PSFlatSurfaceCurrent(token[0], StoI(token[1]),
                                        mesh->GetCopyNumberLevel());
      }
      else
      {
        ps = new G4PSFlatSurfaceCurrent3D(token[0], StoI(token[1]));
      }
      ps->Weighted(StoB(token[2]));
      ps->DivideByArea(StoB(token[3]));
      if(StoB(token[3]))
      {
        ps->SetUnit(token[4]);
      }
      else
      {
        ps->SetUnit("");
      }
      mesh->SetPrimitiveScorer(ps);
    }
  }
  else if(command == qFlatSurfFluxCmd)
  {
    if(CheckMeshPS(mesh, token[0], command))
    {
      G4PSFlatSurfaceFlux* ps = nullptr;
      if(shape == MeshShape::realWorldLogVol || shape == MeshShape::probe)
      {
        ps = new G4PSFlatSurfaceFlux(token[0], StoI(token[1]),
                                     mesh->GetCopyNumberLevel());
      }
      else
      {
        ps = new G4PSFlatSurfaceFlux3D(token[0], StoI(token[1]));
      }
      ps->Weighted(StoB(token[2]));
      ps->DivideByArea(StoB(token[3]));
      if(StoB(token[3]))
      {
        ps->SetUnit(token[4]);
      }
      else
      {
        ps->SetUnit("");
      }
      mesh->SetPrimitiveScorer(ps);
    }
  }
  else if(command == qVolFluxCmd)
  {
    if(CheckMeshPS(mesh, token[0], command))
    {
      G4PSVolumeFlux* ps = nullptr;
      if(shape == MeshShape::realWorldLogVol || shape == MeshShape::probe)
      {
        ps = new G4PSVolumeFlux(token[0], StoI(token[2]),
                                mesh->GetCopyNumberLevel());
      }
      else
      {
        ps = new G4PSVolumeFlux3D(token[0], StoI(token[2]));
      }
      ps->SetDivCos(StoI(token[1]) != 0);
      mesh->SetPrimitiveScorer(ps);
    }
  }
  else if(command == qNofCollisionCmd)
  {
    if(CheckMeshPS(mesh, token[0], command))
    {
      G4PSNofCollision* ps = nullptr;
      if(shape == MeshShape::realWorldLogVol || shape == MeshShape::probe)
      {
        ps = new G4PSNofCollision3D(token[0], mesh->GetCopyNumberLevel());
      }
      else
      {
        ps = new G4PSNofCollision3D(token[0]);
      }
      ps->Weighted(StoB(token[1]));
      mesh->SetPrimitiveScorer(ps);
    }
  }
  else if(command == qPopulationCmd)
  {
    if(CheckMeshPS(mesh, token[0], command))
    {
      G4PSPopulation* ps = nullptr;
      if(shape == MeshShape::realWorldLogVol || shape == MeshShape::probe)
      {
        ps = new G4PSPopulation(token[0], mesh->GetCopyNumberLevel());
      }
      else
      {
        ps = new G4PSPopulation3D(token[0]);
      }
      ps->Weighted(StoB(token[1]));
      mesh->SetPrimitiveScorer(ps);
    }
  }
  else if(command == qTrackCountCmd)
  {
    if(CheckMeshPS(mesh, token[0], command))
    {
      G4PSTrackCounter* ps = nullptr;
      if(shape == MeshShape::realWorldLogVol || shape == MeshShape::probe)
      {
        ps = new G4PSTrackCounter(token[0], StoI(token[1]),
                                  mesh->GetCopyNumberLevel());
      }
      else
      {
        ps = new G4PSTrackCounter3D(token[0], StoI(token[1]));
      }
      ps->Weighted(StoB(token[2]));
      mesh->SetPrimitiveScorer(ps);
    }
  }
  else if(command == qTerminationCmd)
  {
    if(CheckMeshPS(mesh, token[0], command))
    {
      G4PSTermination* ps = nullptr;
      if(shape == MeshShape::realWorldLogVol || shape == MeshShape::probe)
      {
        ps = new G4PSTermination(token[0], mesh->GetCopyNumberLevel());
      }
      else
      {
        ps = new G4PSTermination3D(token[0]);
      }
      ps->Weighted(StoB(token[1]));
      mesh->SetPrimitiveScorer(ps);
    }
  }
  else if(command == qMinKinEAtGeneCmd)
  {
    if(CheckMeshPS(mesh, token[0], command))
    {
      G4PSMinKinEAtGeneration* ps = nullptr;
      if(shape == MeshShape::realWorldLogVol || shape == MeshShape::probe)
      {
        ps = new G4PSMinKinEAtGeneration(token[0], mesh->GetCopyNumberLevel());
      }
      else
      {
        ps = new G4PSMinKinEAtGeneration3D(token[0]);
      }
      ps->SetUnit(token[1]);
      mesh->SetPrimitiveScorer(ps);
    }
  }
  else if(command == qStepCheckerCmd)
  {
    if(CheckMeshPS(mesh, token[0], command))
    {
      G4PSStepChecker* ps = nullptr;
      if(shape == MeshShape::realWorldLogVol || shape == MeshShape::probe)
      {
        ps = new G4PSStepChecker(token[0], mesh->GetCopyNumberLevel());
      }
      else
      {
        ps = new G4PSStepChecker3D(token[0]);
      }
      mesh->SetPrimitiveScorer(ps);
    }

    //
    // Filters
    //
  }
  else if(command == fchargedCmd)
  {
    if(!mesh->IsCurrentPrimitiveScorerNull())
    {
      mesh->SetFilter(new G4SDChargedFilter(token[0]));
    }
    else
    {
      ed << "WARNING[" << fchargedCmd->GetCommandPath()
         << "] : Current quantity is not set. Set or touch a quantity first.";
      command->CommandFailed(ed);
    }
  }
  else if(command == fneutralCmd)
  {
    if(!mesh->IsCurrentPrimitiveScorerNull())
    {
      mesh->SetFilter(new G4SDNeutralFilter(token[0]));
    }
    else
    {
      ed << "WARNING[" << fneutralCmd->GetCommandPath()
         << "] : Current quantity is not set. Set or touch a quantity first.";
      command->CommandFailed(ed);
    }
  }
  else if(command == fkinECmd)
  {
    if(!mesh->IsCurrentPrimitiveScorerNull())
    {
      G4String& name   = token[0];
      G4double elow    = StoD(token[1]);
      G4double ehigh   = StoD(token[2]);
      G4double unitVal = G4UnitDefinition::GetValueOf(token[3]);
      mesh->SetFilter(
        new G4SDKineticEnergyFilter(name, elow * unitVal, ehigh * unitVal));
    }
    else
    {
      ed << "WARNING[" << fkinECmd->GetCommandPath()
         << "] : Current quantity is not set. Set or touch a quantity first.";
      command->CommandFailed(ed);
    }
  }
  else if(command == fparticleKinECmd)
  {
    if(!mesh->IsCurrentPrimitiveScorerNull())
    {
      FParticleWithEnergyCommand(mesh, token);
    }
    else
    {
      ed << "WARNING[" << fparticleKinECmd->GetCommandPath()
         << "] : Current quantity is not set. Set or touch a quantity first.";
      command->CommandFailed(ed);
    }
  }
  else if(command == fparticleCmd)
  {
    if(!mesh->IsCurrentPrimitiveScorerNull())
    {
      FParticleCommand(mesh, token);
    }
    else
    {
      ed << "WARNING[" << fparticleCmd->GetCommandPath()
         << "] : Current quantity is not set. Set or touch a quantity first.";
      command->CommandFailed(ed);
    }
  }
}

G4String G4ScoreQuantityMessenger::GetCurrentValue(G4UIcommand* /*command*/)
{
  G4String val;

  return val;
}

void G4ScoreQuantityMessenger::FillTokenVec(G4String newValues,
                                            G4TokenVec& token)
{
  G4Tokenizer next(newValues);
  G4String val;
  while(!(val = next()).empty())
  {  // Loop checking 12.18.2015 M.Asai
    token.push_back(val);
  }
}

void G4ScoreQuantityMessenger::FParticleCommand(G4VScoringMesh* mesh,
                                                G4TokenVec& token)
{
  //
  // Filter name
  G4String name = token[0];
  //
  // particle list
  std::vector<G4String> pnames;
  for(G4int i = 1; i < (G4int) token.size(); i++)
  {
    pnames.push_back(token[i]);
  }
  //
  // Attach Filter
  mesh->SetFilter(new G4SDParticleFilter(name, pnames));
}

void G4ScoreQuantityMessenger::FParticleWithEnergyCommand(G4VScoringMesh* mesh,
                                                          G4TokenVec& token)
{
  G4String& name   = token[0];
  G4double elow    = StoD(token[1]);
  G4double ehigh   = StoD(token[2]);
  G4double unitVal = G4UnitDefinition::GetValueOf(token[3]);
  auto  filter =
    new G4SDParticleWithEnergyFilter(name, elow * unitVal, ehigh * unitVal);
  for(G4int i = 4; i < (G4int) token.size(); i++)
  {
    filter->add(token[i]);
  }
  mesh->SetFilter(filter);
}

G4bool G4ScoreQuantityMessenger::CheckMeshPS(G4VScoringMesh* mesh,
                                             G4String& psname,
                                             G4UIcommand* command)
{
  if(!mesh->FindPrimitiveScorer(psname))
  {
    return true;
  }
  
  G4ExceptionDescription ed;
  ed << "WARNING[" << qTouchCmd->GetCommandPath() << "] : Quantity name, \""
     << psname << "\", is already existing.";
  command->CommandFailed(ed);
  mesh->SetNullToCurrentPrimitiveScorer();
  return false;
}
