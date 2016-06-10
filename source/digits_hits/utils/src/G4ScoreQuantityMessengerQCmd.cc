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
// $Id: G4ScoreQuantityMessengerQCmd.cc 73819 2013-09-12 15:52:52Z gcosmo $
//
// ---------------------------------------------------------------------
// Modifications
// 08-Oct-2010 T.Aso remove unit of G4PSPassageCellCurrent.
// 01-Jun-2012 T.Aso Support weighted/dividedByArea options 
//                      in flatCurrent and flatFulx commands.
// 29-Mar-2013 T.Aso Support weighted options in the nOfTrack command.
//                   Support a boundary flag option in the nOfStep command
//                  for skipping stepLength=0 steps.
// 
// ---------------------------------------------------------------------

#include "G4ScoreQuantityMessenger.hh"
#include "G4ScoringManager.hh"
#include "G4VScoringMesh.hh"

#include "G4PSCellCharge3D.hh"
#include "G4PSCellFlux3D.hh"
#include "G4PSPassageCellFlux3D.hh"
#include "G4PSEnergyDeposit3D.hh"
#include "G4PSDoseDeposit3D.hh"
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
#include "G4PSNofCollision3D.hh"
#include "G4PSPopulation3D.hh"
#include "G4PSTrackCounter3D.hh"
#include "G4PSTermination3D.hh"

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
#include "G4Tokenizer.hh"
#include "G4UnitsTable.hh"

void G4ScoreQuantityMessenger::QuantityCommands()
{
  G4UIparameter* param;

  //
  // Quantity commands
  quantityDir = new G4UIdirectory("/score/quantity/");
  quantityDir->SetGuidance("Scoring quantity of the mesh.");
  //
  qTouchCmd= new G4UIcmdWithAString("/score/quantity/touch",this);
  qTouchCmd->SetGuidance("Assign previously defined quantity to the current quantity.");
  qTouchCmd->SetParameterName("qname",false);
  //
  qGetUnitCmd = new G4UIcmdWithoutParameter("/score/quantity/get/unit",this);
  qGetUnitCmd->SetGuidance("Print output unit of the current quantity.");
  //
  qSetUnitCmd = new G4UIcmdWithAString("/score/quantity/set/unit",this);
  qSetUnitCmd->SetGuidance("Set output unit of the current quantity.");
  qSetUnitCmd->SetParameterName("unit",false);

  // Primitive Scorers
  qeDepCmd = new G4UIcommand("/score/quantity/energyDeposit",this);
  qeDepCmd->SetGuidance("Energy deposit scorer.");
  qeDepCmd->
  SetGuidance("[usage] /score/quantiy/energyDeposit qname unit");
  qeDepCmd->SetGuidance("  qname  :(String) scorer name");
  qeDepCmd->SetGuidance("  unit   :(String) unit");
  param = new G4UIparameter("qname",'s',false);
  qeDepCmd->SetParameter(param);
  param = new G4UIparameter("unit",'s',true);
  param->SetDefaultValue("MeV");
  qeDepCmd->SetParameter(param);
  //
  qCellChgCmd  = new G4UIcommand("/score/quantity/cellCharge",this);
  qCellChgCmd->SetGuidance("Cell charge scorer.");
  qCellChgCmd->
  SetGuidance("[usage] /score/quantiy/cellCharge qname unit");
  qCellChgCmd->SetGuidance("  qname  :(String) scorer name");
  qCellChgCmd->SetGuidance("  unit   :(String) unit");
  param = new G4UIparameter("qname",'s',false);
  qCellChgCmd->SetParameter(param);
  param = new G4UIparameter("unit",'s',true);
  param->SetDefaultValue("e+");
  qCellChgCmd->SetParameter(param);
  //
  qCellFluxCmd = new G4UIcommand("/score/quantity/cellFlux",this);
  qCellFluxCmd->SetGuidance("Cell flux scorer.");
  qCellFluxCmd->
  SetGuidance("[usage] /score/quantiy/cellFlux qname unit");
  qCellFluxCmd->SetGuidance("  qname  :(String) scorer name");
  qCellFluxCmd->SetGuidance("  unit   :(String) unit");
  param = new G4UIparameter("qname",'s',false);
  qCellFluxCmd->SetParameter(param);
  param = new G4UIparameter("unit",'s',true);
  param->SetDefaultValue("percm2");
  qCellFluxCmd->SetParameter(param);
  //
  qPassCellFluxCmd = new G4UIcommand("/score/quantity/passageCellFlux",this);
  qPassCellFluxCmd->SetGuidance("Passage cell flux scorer");
  qPassCellFluxCmd->
  SetGuidance("[usage] /score/quantiy/passageCellFlux qname unit");
  qPassCellFluxCmd->SetGuidance("  qname  :(String) scorer name");
  qPassCellFluxCmd->SetGuidance("  unit   :(String) unit");
  param = new G4UIparameter("qname",'s',false);
  qPassCellFluxCmd->SetParameter(param);
  param = new G4UIparameter("unit",'s',true);
  param->SetDefaultValue("percm2");
  qPassCellFluxCmd->SetParameter(param);
  //
  qdoseDepCmd = new G4UIcommand("/score/quantity/doseDeposit",this);
  qdoseDepCmd->SetGuidance("Dose deposit scorer.");
  qdoseDepCmd->
  SetGuidance("[usage] /score/quantiy/doseDeposit qname unit");
  qdoseDepCmd->SetGuidance("  qname  :(String) scorer name");
  qdoseDepCmd->SetGuidance("  unit   :(String) unit");
  param = new G4UIparameter("qname",'s',false);
  qdoseDepCmd->SetParameter(param);
  param = new G4UIparameter("unit",'s',true);
  param->SetDefaultValue("Gy");
  qdoseDepCmd->SetParameter(param);
  //
  qnOfStepCmd = new G4UIcommand("/score/quantity/nOfStep",this);
  qnOfStepCmd->SetGuidance("Number of step scorer.");
  qnOfStepCmd->
  SetGuidance("[usage] /score/quantiy/nOfStep qname");
  qnOfStepCmd->
  SetGuidance("[usage] /score/quantiy/nOfStep qname  bflag"); 
  qnOfStepCmd->SetGuidance("  qname  :(String) scorer name");
  qnOfStepCmd->SetGuidance("  bflag  :(Bool) Skip zero step ");
  qnOfStepCmd->SetGuidance("          at geometry boundary if true");
  param = new G4UIparameter("qname",'s',false);
  qnOfStepCmd->SetParameter(param);
  param = new G4UIparameter("bflag",'b',true);  
  param->SetDefaultValue("false");
  qnOfStepCmd->SetParameter(param);
  //
  qnOfSecondaryCmd = new G4UIcommand("/score/quantity/nOfSecondary",this);
  qnOfSecondaryCmd->SetGuidance("Number of secondary scorer.");
  qnOfSecondaryCmd->
  SetGuidance("[usage] /score/quantiy/nOfSecondary qname"); 
  qnOfSecondaryCmd->SetGuidance("  qname  :(String) scorer name");
  param = new G4UIparameter("qname",'s',false);
  qnOfSecondaryCmd->SetParameter(param);
  //
  qTrackLengthCmd = new G4UIcommand("/score/quantity/trackLength",this);
  qTrackLengthCmd->SetGuidance("Track length scorer.");
  qTrackLengthCmd->
  SetGuidance("[usage] /score/quantiy/trackLength qname wflag kflag vflag unit");
  qTrackLengthCmd->SetGuidance("  qname  :(String) scorer name");
  qTrackLengthCmd->SetGuidance("  wflag  :(Bool) weighted");
  qTrackLengthCmd->SetGuidance("  kflag  :(Bool) multiply kinetic energy");
  qTrackLengthCmd->SetGuidance("  vflag  :(Bool) divide by velocity");
  qTrackLengthCmd->SetGuidance("  unit   :(String) unit");
  param = new G4UIparameter("qname",'s',false);
  qTrackLengthCmd->SetParameter(param);
  param = new G4UIparameter("wflag",'b',true);
  param->SetDefaultValue("false");
  qTrackLengthCmd->SetParameter(param);
  param = new G4UIparameter("kflag",'b',true);
  param->SetDefaultValue("false");
  qTrackLengthCmd->SetParameter(param);
  param = new G4UIparameter("vflag",'b',true);
  param->SetDefaultValue("false");
  qTrackLengthCmd->SetParameter(param);
  param = new G4UIparameter("unit",'s',true);
  param->SetDefaultValue("mm");
  qTrackLengthCmd->SetParameter(param);
  //
  qPassCellCurrCmd = new G4UIcommand("/score/quantity/passageCellCurrent",this);
  qPassCellCurrCmd->SetGuidance("Passage cell current scorer.");
  qPassCellCurrCmd->
  SetGuidance("[usage] /score/quantiy/passageCellCurrent qname wflag");
  qPassCellCurrCmd->SetGuidance("  qname  :(String) scorer name");
  qPassCellCurrCmd->SetGuidance("  wflag  :(Bool) weighted");
  param = new G4UIparameter("qname",'s',false);
  qPassCellCurrCmd->SetParameter(param);
  param = new G4UIparameter("wflag",'b',true);
  param->SetDefaultValue("true");
  qPassCellCurrCmd->SetParameter(param);
  //
  qPassTrackLengthCmd = new G4UIcommand("/score/quantity/passageTrackLength",this);
  qPassTrackLengthCmd->SetGuidance("Passage track length scorer.");
  qPassTrackLengthCmd->
  SetGuidance("[usage] /score/quantiy/passageTrackLength qname wflag unit");
  qPassTrackLengthCmd->SetGuidance("  qname  :(String) scorer name");
  qPassTrackLengthCmd->SetGuidance("  wflag  :(Bool) weighted");
  qPassTrackLengthCmd->SetGuidance("  unit   :(Bool) unit");
  param = new G4UIparameter("qname",'s',false);
  qPassTrackLengthCmd->SetParameter(param);
  param = new G4UIparameter("wflag",'b',true);
  param->SetDefaultValue("true");
  qPassTrackLengthCmd->SetParameter(param);
  param = new G4UIparameter("unit",'s',true);
  param->SetDefaultValue("mm");
  qPassTrackLengthCmd->SetParameter(param);
  //
  qFlatSurfCurrCmd = new G4UIcommand("/score/quantity/flatSurfaceCurrent",this);
  qFlatSurfCurrCmd->SetGuidance("Flat surface current Scorer.");
  qFlatSurfCurrCmd->
  SetGuidance("[usage] /score/quantiy/flatSurfaceCurrent qname dflag wflag aflag unit");
  qFlatSurfCurrCmd->SetGuidance("  qname  :(String) scorer name");
  qFlatSurfCurrCmd->SetGuidance("  dflag  :(Int) direction flag");
  qFlatSurfCurrCmd->SetGuidance("         : 0 = Both In and Out");
  qFlatSurfCurrCmd->SetGuidance("         : 1 = In only");
  qFlatSurfCurrCmd->SetGuidance("         : 2 = Out only");
  qFlatSurfCurrCmd->SetGuidance("  wflag  :(Bool) weighted");
  qFlatSurfCurrCmd->SetGuidance("  aflag  :(Bool) divide by area");
  qFlatSurfCurrCmd->SetGuidance("  unit   :(Bool) unit");
  param = new G4UIparameter("qname",'s',false);
  qFlatSurfCurrCmd->SetParameter(param);
  param = new G4UIparameter("dflag",'i',true);
  param->SetDefaultValue("0");
  qFlatSurfCurrCmd->SetParameter(param);
  param = new G4UIparameter("wflag",'b',true);
  param->SetDefaultValue("true");
  qFlatSurfCurrCmd->SetParameter(param);
  param = new G4UIparameter("aflag",'b',true);
  param->SetDefaultValue("true");
  qFlatSurfCurrCmd->SetParameter(param);
  param = new G4UIparameter("unit",'s',true);
  param->SetDefaultValue("percm2");
  qFlatSurfCurrCmd->SetParameter(param);
  //
  qFlatSurfFluxCmd = new G4UIcommand("/score/quantity/flatSurfaceFlux",this);
  qFlatSurfFluxCmd->SetGuidance("Flat surface flux scorer.");
  qFlatSurfFluxCmd->
  SetGuidance("[usage] /score/quantiy/flatSurfaceFlux qname dflag unit");
  qFlatSurfFluxCmd->SetGuidance("  qname  :(String) scorer name");
  qFlatSurfFluxCmd->SetGuidance("  dflag  :(Int) direction flag");
  qFlatSurfFluxCmd->SetGuidance("         : 0 = Both In and Out");
  qFlatSurfFluxCmd->SetGuidance("         : 1 = In only");
  qFlatSurfFluxCmd->SetGuidance("         : 2 = Out only");
  qFlatSurfFluxCmd->SetGuidance("  wflag  :(Bool) weighted");
  qFlatSurfFluxCmd->SetGuidance("  aflag  :(Bool) divide by area");
  qFlatSurfFluxCmd->SetGuidance("  unit   :(String) unit");
  param = new G4UIparameter("qname",'s',false);
  qFlatSurfFluxCmd->SetParameter(param);
  param = new G4UIparameter("dflag",'i',true);
  param->SetDefaultValue("0");
  qFlatSurfFluxCmd->SetParameter(param);
  param = new G4UIparameter("wflag",'b',true);
  param->SetDefaultValue("true");
  qFlatSurfFluxCmd->SetParameter(param);
  param = new G4UIparameter("aflag",'b',true);
  param->SetDefaultValue("true");
  qFlatSurfFluxCmd->SetParameter(param);
  param = new G4UIparameter("unit",'s',true);
  param->SetDefaultValue("percm2");
  qFlatSurfFluxCmd->SetParameter(param);
  //
//  qSphereSurfCurrCmd = new G4UIcommand("/score/quantity/sphereSurfaceCurrent",this);
//  qSphereSurfCurrCmd->SetGuidance("Sphere surface current Scorer.");
//  qSphereSurfCurrCmd->
//      SetGuidance("[usage] /score/quantiy/sphereSurfaceCurrent qname dflag wflag aflag unit");
//  qSphereSurfCurrCmd->SetGuidance("  qname  :(String) scorer name");
//  qSphereSurfCurrCmd->SetGuidance("  dflag  :(Int) direction flag");
//  qSphereSurfCurrCmd->SetGuidance("         : 0 = Both In and Out");
//  qSphereSurfCurrCmd->SetGuidance("         : 1 = In only");
//  qSphereSurfCurrCmd->SetGuidance("         : 2 = Out only");
//  qSphereSurfCurrCmd->SetGuidance("  wflag  :(Bool) Weighted");
//  qSphereSurfCurrCmd->SetGuidance("  aflag  :(Bool) DivideByArea");
//  qSphereSurfCurrCmd->SetGuidance("  unit   :(String) unit");
//  param = new G4UIparameter("qname",'s',false);
//  qSphereSurfCurrCmd->SetParameter(param);
//  param = new G4UIparameter("dflag",'i',true);
//  param->SetDefaultValue("0");
//  qSphereSurfCurrCmd->SetParameter(param);
//  param = new G4UIparameter("wflag",'b',true);
//  param->SetDefaultValue("true");
//  qSphereSurfCurrCmd->SetParameter(param);
//  param = new G4UIparameter("aflag",'b',true);
//  param->SetDefaultValue("true");
//  qSphereSurfCurrCmd->SetParameter(param);
//  param = new G4UIparameter("unit",'s',true);
//  param->SetDefaultValue("percm2");
//  qSphereSurfCurrCmd->SetParameter(param);

  //
//  qSphereSurfFluxCmd = new G4UIcommand("/score/quantity/sphereSurfaceFlux",this);
//  qSphereSurfFluxCmd->SetGuidance("Sphere surface Flux Scorer.");
//  qSphereSurfFluxCmd->
//  SetGuidance("[usage] /score/quantiy/sphereSurfaceFlux qname dflag unit");
//  qSphereSurfFluxCmd->SetGuidance("  qname  :(String) scorer name");
//  qSphereSurfFluxCmd->SetGuidance("  dflag  :(Int) direction flag");
//  qSphereSurfFluxCmd->SetGuidance("         : 0 = Both In and Out");
//  qSphereSurfFluxCmd->SetGuidance("         : 1 = In only");
//  qSphereSurfFluxCmd->SetGuidance("         : 2 = Out only");
//  qSphereSurfFluxCmd->SetGuidance("  wflag  :(Bool) weighted");
//  qSphereSurfFluxCmd->SetGuidance("  aflag  :(Bool) divide by area");
//  qSphereSurfFluxCmd->SetGuidance("  unit   :(String) unit");
//  param = new G4UIparameter("qname",'s',false);
//  qSphereSurfFluxCmd->SetParameter(param);
//  param = new G4UIparameter("dflag",'i',true);
//  param->SetDefaultValue("0");
//  qSphereSurfFluxCmd->SetParameter(param);
//  param = new G4UIparameter("wflag",'b',true);
//  param->SetDefaultValue("true");
//  qSphereSurfFluxCmd->SetParameter(param);
//  param = new G4UIparameter("aflag",'b',true);
//  param->SetDefaultValue("true");
//  qSphereSurfFluxCmd->SetParameter(param);
//  param = new G4UIparameter("unit",'s',true);
//  param->SetDefaultValue("percm2");
//  qSphereSurfFluxCmd->SetParameter(param);

  //
//  qCylSurfCurrCmd = new G4UIcommand("/score/quantity/cylinderSurfaceCurrent",this);
//  qCylSurfCurrCmd->SetGuidance("Cylinder surface current Scorer.");
//  qCylSurfCurrCmd->
//  SetGuidance("[usage] /score/quantiy/cylinderSurfaceCurrent qname dflag wflag aflag unit");
//  qCylSurfCurrCmd->SetGuidance("  qname  :(String) scorer name");
//  qCylSurfCurrCmd->SetGuidance("  dflag  :(Int) direction flag");
//  qCylSurfCurrCmd->SetGuidance("         : 0 = Both In and Out");
//  qCylSurfCurrCmd->SetGuidance("         : 1 = In only");
//  qCylSurfCurrCmd->SetGuidance("         : 2 = Out only");
//  qCylSurfCurrCmd->SetGuidance("  wflag  :(Bool) Weighted");
//  qCylSurfCurrCmd->SetGuidance("  aflag  :(Bool) DivideByArea");
//  qCylSurfCurrCmd->SetGuidance("  unit   :(String) unit");
//  param = new G4UIparameter("qname",'s',false);
//  qCylSurfCurrCmd->SetParameter(param);
//  param = new G4UIparameter("dflag",'i',true);
//  param->SetDefaultValue("0");
//  qCylSurfCurrCmd->SetParameter(param);
//  param = new G4UIparameter("wflag",'b',true);
//  param->SetDefaultValue("true");
//  qCylSurfCurrCmd->SetParameter(param);
//  param = new G4UIparameter("aflag",'b',true);
//  param->SetDefaultValue("true");
//  qCylSurfCurrCmd->SetParameter(param);
//  param = new G4UIparameter("unit",'s',true);
//  param->SetDefaultValue("percm2");
//  qCylSurfCurrCmd->SetParameter(param);
//
//  qCylSurfFluxCmd = new G4UIcommand("/score/quantity/cylinderSurfaceFlux",this);
//  qCylSurfFluxCmd->SetGuidance("Cylinder surface Flux Scorer.");
//  qCylSurfFluxCmd->
//  SetGuidance("[usage] /score/quantiy/cylinderSurfaceFlux qname dflag unit");
//  qCylSurfFluxCmd->SetGuidance("  qname  :(String) scorer name");
//  qCylSurfFluxCmd->SetGuidance("  dflag  :(Int) direction flag");
//  qCylSurfFluxCmd->SetGuidance("         : 0 = Both In and Out");
//  qCylSurfFluxCmd->SetGuidance("         : 1 = In only");
//  qCylSurfFluxCmd->SetGuidance("         : 2 = Out only");
//  qCylSurfFluxCmd->SetGuidance("  wflag  :(Bool) weighted");
//  qCylSurfFluxCmd->SetGuidance("  aflag  :(Bool) divide by area");
//  qCylSurfFluxCmd->SetGuidance("  unit   :(String) unit");
//  param = new G4UIparameter("qname",'s',false);
//  qCylSurfFluxCmd->SetParameter(param);
//  param = new G4UIparameter("dflag",'i',true);
//  param->SetDefaultValue("0");
//  qCylSurfFluxCmd->SetParameter(param);
//  param = new G4UIparameter("wflag",'b',true);
//  param->SetDefaultValue("true");
//  qCylSurfFluxCmd->SetParameter(param);
//  param = new G4UIparameter("aflag",'b',true);
//  param->SetDefaultValue("true");
//  qCylSurfFluxCmd->SetParameter(param);
//  param = new G4UIparameter("unit",'s',true);
//  param->SetDefaultValue("percm2");
//  qCylSurfFluxCmd->SetParameter(param);
//
  //
  qNofCollisionCmd = new G4UIcommand("/score/quantity/nOfCollision",this);
  qNofCollisionCmd->SetGuidance("Number of collision scorer.");
  qNofCollisionCmd->
  SetGuidance("[usage] /score/quantiy/nOfCollision qname wflag");
  qNofCollisionCmd->SetGuidance("  qname  :(String) scorer name");
  param = new G4UIparameter("qname",'s',false);
  qNofCollisionCmd->SetParameter(param);
  param = new G4UIparameter("wflag",'b',true);
  param->SetDefaultValue("false");
  qNofCollisionCmd->SetParameter(param);
  //
  qPopulationCmd = new G4UIcommand("/score/quantity/population",this);
  qPopulationCmd->SetGuidance("Population scorer.");
  qPopulationCmd->
   SetGuidance("[usage] /score/quantiy/population qname wflag");
  qPopulationCmd->SetGuidance("  qname  :(String) scorer name");
  qPopulationCmd->SetGuidance("  wflag  :(Bool) weighted");
  param = new G4UIparameter("qname",'s',false);
  qPopulationCmd->SetParameter(param);
  param = new G4UIparameter("wflag",'b',true);
  param->SetDefaultValue("false");
  qPopulationCmd->SetParameter(param);

  //
  qTrackCountCmd = new G4UIcommand("/score/quantity/nOfTrack",this);
  qTrackCountCmd->SetGuidance("Number of track scorer.");
  qTrackCountCmd->
  SetGuidance("[usage] /score/quantiy/nOfTrack qname dflag wflag");
  qTrackCountCmd->SetGuidance("  qname  :(String) scorer name");
  qTrackCountCmd->SetGuidance("  dflag  :(Int) direction");
  qTrackCountCmd->SetGuidance("         : 0 = Both In and Out");
  qTrackCountCmd->SetGuidance("         : 1 = In only");
  qTrackCountCmd->SetGuidance("         : 2 = Out only");
  qTrackCountCmd->SetGuidance("  wflag  :(Bool) weighted");
  param = new G4UIparameter("qname",'s',false);
  qTrackCountCmd->SetParameter(param);
  param = new G4UIparameter("dflag",'i',true);
  param->SetDefaultValue("0");
  qTrackCountCmd->SetParameter(param);
  param = new G4UIparameter("wflag",'b',true);
  param->SetDefaultValue("false");
  qTrackCountCmd->SetParameter(param);

  //
  qTerminationCmd = new G4UIcommand("/score/quantity/nOfTerminatedTrack",this);
  qTerminationCmd->SetGuidance("Number of terminated tracks scorer.");
  qTerminationCmd->
      SetGuidance("[usage] /score/quantiy/nOfTerminatedTrack qname wflag");
  qTerminationCmd->SetGuidance("  qname  :(String) scorer name");
  qTerminationCmd->SetGuidance("  wflag  :(Bool) weighted");
  param = new G4UIparameter("qname",'s',false);
  qTerminationCmd->SetParameter(param);
  param = new G4UIparameter("wflag",'b',true);
  param->SetDefaultValue("false");
  qTerminationCmd->SetParameter(param);

  //
  qMinKinEAtGeneCmd = new G4UIcommand("/score/quantity/minKinEAtGeneration",this);
  qMinKinEAtGeneCmd->SetGuidance("Min Kinetic Energy at Generation");
  qMinKinEAtGeneCmd->
  SetGuidance("[usage] /score/quantiy/minKinEAtGeneration qname unit");
  qMinKinEAtGeneCmd->SetGuidance("  qname  :(String) scorer name");
  qMinKinEAtGeneCmd->SetGuidance("  unit   :(String) unit name");
  param = new G4UIparameter("qname",'s',false);
  qMinKinEAtGeneCmd->SetParameter(param);
  param = new G4UIparameter("unit",'s',true);
  param->SetDefaultValue("MeV");
  qMinKinEAtGeneCmd->SetParameter(param);
  //
  qStepCheckerCmd = new G4UIcommand("/score/quantity/stepChecker",this);
  qStepCheckerCmd->SetGuidance("Display a comment when this PS is invoked");
  qStepCheckerCmd->
  SetGuidance("[usage] /score/quantiy/stepChecker qname");
  qStepCheckerCmd->SetGuidance("  qname  :(String) scorer name");
  param = new G4UIparameter("qname",'s',false);
  qStepCheckerCmd->SetParameter(param);

}

