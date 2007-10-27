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
// $Id: G4ScoringMessenger.cc,v 1.23 2007-10-27 00:28:44 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ---------------------------------------------------------------------

#include "G4ScoringMessenger.hh"
#include "G4ScoringManager.hh"
#include "G4VScoringMesh.hh"
#include "G4ScoringBox.hh"

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

G4ScoringMessenger::G4ScoringMessenger(G4ScoringManager* SManager)
:fSMan(SManager)
{
  G4UIparameter* param;

  scoreDir = new G4UIdirectory("/score/");
  scoreDir->SetGuidance("Interactive scoring commands.");

  listCmd = new G4UIcmdWithoutParameter("/score/list",this);
  listCmd->SetGuidance("List scoring worlds.");

  dumpCmd = new G4UIcmdWithoutParameter("/score/dump",this);
  dumpCmd->SetGuidance("Dump scorer results ");

  verboseCmd = new G4UIcmdWithAnInteger("/score/verbose",this);
  verboseCmd->SetGuidance("Verbosity");

  meshCreateDir = new G4UIdirectory("/score/create/");
  meshCreateDir->SetGuidance("Interactive scoring commands.");
  //
  // Mesh commands
  meshBoxCreateCmd = new G4UIcmdWithAString("/score/create/boxMesh",this);
  meshBoxCreateCmd->SetGuidance("Create scoring mesh.");
  meshBoxCreateCmd->SetParameterName("MeshName",false);
  //
  meshTubsCreateCmd = new G4UIcmdWithAString("/score/create/tubsMesh",this);
  meshTubsCreateCmd->SetGuidance("Create scoring mesh.");
  meshTubsCreateCmd->SetParameterName("MeshName",false);
  //
  meshSphereCreateCmd = new G4UIcmdWithAString("/score/create/sphereMesh",this);
  meshSphereCreateCmd->SetGuidance("Create scoring mesh.");
  meshSphereCreateCmd->SetParameterName("MeshName",false);
  //
  meshOpnCmd = new G4UIcmdWithAString("/score/open",this);
  meshOpnCmd->SetGuidance("Open scoring mesh.");
  meshOpnCmd->SetParameterName("MeshName",false);
  //
  meshClsCmd = new G4UIcmdWithoutParameter("/score/close",this);
  meshClsCmd->SetGuidance("Close scoring mesh.");
  //
  meshActCmd = new G4UIcmdWithABool("/score/mesh/activate",this);
  meshActCmd->SetGuidance("Activate scoring mesh.");
  meshActCmd->SetParameterName("MeshName",false);
  //
  mBoxSizeCmd = new G4UIcmdWith3VectorAndUnit("/score/mesh/boxsize",this);
  mBoxSizeCmd->SetGuidance("Define Size of scoring mesh.");
  mBoxSizeCmd->SetParameterName("Di","Dj","Dk",false,false);
  mBoxSizeCmd->SetDefaultUnit("mm");
  //
  //   Division command
  mBinCmd = new G4UIcommand("/score/mesh/nbin",this);
  mBinCmd->SetGuidance("Define segmentation of scoring mesh.");
  mBinCmd->SetGuidance("[usage] /score/mesh/nbin");
  mBinCmd->SetGuidance("  Ni  :(int) Number of bins i ");
  mBinCmd->SetGuidance("  Nj  :(int) Number of bins j ");
  mBinCmd->SetGuidance("  Nk  :(int) Number of bins k ");
  mBinCmd->SetGuidance("  Axis:(int) Axis of division ");
  mBinCmd->SetGuidance("  P1..Pn-1  :(double) \"paramter from P1 to Pn-1 for division.\"");
  param = new G4UIparameter("Ni",'i',false);
  param->SetDefaultValue("1");
  mBinCmd->SetParameter(param);
  param = new G4UIparameter("Nj",'i',false);
  param->SetDefaultValue("1");
  mBinCmd->SetParameter(param);
  param = new G4UIparameter("Nk",'i',false);
  param->SetDefaultValue("1");
  mBinCmd->SetParameter(param);
  param = new G4UIparameter("Axis",'i',true);
  param->SetDefaultValue("3");
  mBinCmd->SetParameter(param);
  //
  //   Placement command
  mTransDir = new G4UIdirectory("/score/mesh/translate/");
  mTransDir->SetGuidance("Placement of scoring mesh");
  //
  mTResetCmd = new G4UIcmdWithoutParameter("/score/mesh/translate/reset",this);
  mTResetCmd->SetGuidance("Reset translation of scoring mesh placement.");
  //
  mTXyzCmd = new G4UIcmdWith3VectorAndUnit("/score/mesh/translate/xyz",this);
  mTXyzCmd->SetGuidance("Translation the current scoring mesh to the position.");
  mTXyzCmd->SetParameterName("X","Y","Z",false,false);
  mTXyzCmd->SetDefaultUnit("mm");
  //
  mRotDir = new G4UIdirectory("/score/mesh/rotate/");
  mRotDir->SetGuidance("Placement of scoring mesh");
  //
  mRResetCmd = new G4UIcmdWithoutParameter("/score/mesh/rotate/reset",this);
  mRResetCmd->SetGuidance("Reset rotation of scoring mesh placement.");
  //
  mRotXCmd = new G4UIcmdWithADoubleAndUnit("/score/mesh/rotate/rotX",this);
  mRotXCmd->SetGuidance("Add rotation to the current scoring mesh in X.");
  mRotXCmd->SetParameterName("Rx",false);
  mRotXCmd->SetDefaultUnit("deg");
  //
  mRotYCmd = new G4UIcmdWithADoubleAndUnit("/score/mesh/rotate/rotY",this);
  mRotYCmd->SetGuidance("Add rotation to the current scoring mesh in Y.");
  mRotYCmd->SetParameterName("Ry",false);
  mRotYCmd->SetDefaultUnit("deg");
  //
  mRotZCmd = new G4UIcmdWithADoubleAndUnit("/score/mesh/rotate/rotZ",this);
  mRotZCmd->SetGuidance("Add rotation to the current scoring mesh in Z.");
  mRotZCmd->SetParameterName("Rz",false);
  mRotZCmd->SetDefaultUnit("deg");
  //

  // Draw Scoring result
  drawCmd = new G4UIcommand("/score/draw",this);
  drawCmd->SetGuidance("Draw scorer results ");
  param = new G4UIparameter("meshName",'s',false);
  drawCmd->SetParameter(param);
  param = new G4UIparameter("psName",'s',false);
  drawCmd->SetParameter(param);
  param = new G4UIparameter("proj",'i',true);
  param->SetDefaultValue(111);
  drawCmd->SetParameter(param);

  // Dump scoring result
  dumpToFileCmd = new G4UIcommand("/score/dumpToFile", this);
  dumpToFileCmd->SetGuidance("Dump scorer results to file ");
  param = new G4UIparameter("meshName", 's', false);
  dumpToFileCmd->SetParameter(param);
  param = new G4UIparameter("psName", 's', false);
  dumpToFileCmd->SetParameter(param);
  param = new G4UIparameter("fileName", 's', false);
  dumpToFileCmd->SetParameter(param);
  param = new G4UIparameter("option", 's', true);
  param->SetDefaultValue("csv");
  dumpToFileCmd->SetParameter(param);

  //
  // Quantity commands
  quantityDir = new G4UIdirectory("/score/quantity/");
  quantityDir->SetGuidance("Scoring quantity of the mesh");
  //
  qTouchCmd= new G4UIcmdWithAString("/score/quantity/touch",this);
  qTouchCmd->SetGuidance("Assign previously defined quantity to current quantity.");
  qTouchCmd->SetParameterName("qname",false);
  //
  qeDepCmd = new G4UIcmdWithAString("/score/quantity/eDep",this);
  qeDepCmd->SetGuidance("Energy Deposit Scorer");
  qeDepCmd->SetParameterName("qname",false);
  //
  qCellChgCmd  = new G4UIcmdWithAString("/score/quantity/cellCharge",this);
  qCellChgCmd->SetGuidance("Cell Charge Scorer");
  qCellChgCmd->SetParameterName("qname",false);
  //
  qCellFluxCmd = new G4UIcmdWithAString("/score/quantity/cellFlux",this);
  qCellFluxCmd->SetGuidance("Cell Flux Scorer");
  qCellFluxCmd->SetParameterName("qname",false);
  //
  qPassCellFluxCmd = new G4UIcmdWithAString("/score/quantity/passageCellFlux",this);
  qPassCellFluxCmd->SetGuidance("Passage Cell Flux Scorer");
  qPassCellFluxCmd->SetParameterName("qname",false);
  //
  qdoseDepCmd = new G4UIcmdWithAString("/score/quantity/doseDeposit",this);
  qdoseDepCmd->SetGuidance("Dose Deposit Scorer");
  qdoseDepCmd->SetParameterName("qname",false);
  //
  qnOfStepCmd = new G4UIcmdWithAString("/score/quantity/nOfStep",this);
  qnOfStepCmd->SetGuidance("Number of Step Scorer ");
  qnOfStepCmd->SetParameterName("qname",false);
  //
  qnOfSecondaryCmd = new G4UIcmdWithAString("/score/quantity/nOfSecondary",this);
  qnOfSecondaryCmd->SetGuidance("Number of Secondary Scorer ");
  qnOfSecondaryCmd->SetParameterName("qname",false);
  //
  qTrackLengthCmd = new G4UIcommand("/score/quantity/trackLength",this);
  qTrackLengthCmd->SetGuidance("TrackLength Scorer");
  qTrackLengthCmd->
      SetGuidance("[usage] /score/quantiy/trackLength qname wflag kflag vflag ");
  qTrackLengthCmd->SetGuidance("  qname  :(String) scorer name");
  qTrackLengthCmd->SetGuidance("  wflag  :(Bool) Weighted");
  qTrackLengthCmd->SetGuidance("  kflag  :(Bool) MultiplyKineticEnergy");
  qTrackLengthCmd->SetGuidance("  vflag  :(Bool) DivideByVelocity");
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
  //
  qPassCellCurrCmd = new G4UIcommand("/score/quantity/passageCellCurrent",this);
  qPassCellCurrCmd->SetGuidance("PassageCellCurrent Scorer");
  qPassCellCurrCmd->
      SetGuidance("[usage] /score/quantiy/passageCellCurrent qname wflag");
  qPassCellCurrCmd->SetGuidance("  qname  :(String) scorer name");
  qPassCellCurrCmd->SetGuidance("  wflag  :(Bool) Weighted");
  param = new G4UIparameter("qname",'s',false);
  qPassCellCurrCmd->SetParameter(param);
  param = new G4UIparameter("wflag",'b',true);
  param->SetDefaultValue("true");
  qPassCellCurrCmd->SetParameter(param);
  //
  qPassTrackLengthCmd = new G4UIcommand("/score/quantity/passageTrackLength",this);
  qPassTrackLengthCmd->SetGuidance("PassageTrackLength Scorer");
  qPassTrackLengthCmd->
      SetGuidance("[usage] /score/quantiy/passageTrackLength qname wflag");
  qPassTrackLengthCmd->SetGuidance("  qname  :(String) scorer name");
  qPassTrackLengthCmd->SetGuidance("  wflag  :(Bool) Weighted");
  param = new G4UIparameter("qname",'s',false);
  qPassTrackLengthCmd->SetParameter(param);
  param = new G4UIparameter("wflag",'b',true);
  param->SetDefaultValue("true");
  qPassTrackLengthCmd->SetParameter(param);
  //
  qFlatSurfCurrCmd = new G4UIcommand("/score/quantity/flatSurfCurrent",this);
  qFlatSurfCurrCmd->SetGuidance("Flat surface current Scorer");
  qFlatSurfCurrCmd->
      SetGuidance("[usage] /score/quantiy/flatSurfCurrent qname dflag wflag aflag");
  qFlatSurfCurrCmd->SetGuidance("  qname  :(String) scorer name");
  qFlatSurfCurrCmd->SetGuidance("  dflag  :(Int) direction flag");
  qFlatSurfCurrCmd->SetGuidance("         : 0 = Both In and Out");
  qFlatSurfCurrCmd->SetGuidance("         : 1 = In only");
  qFlatSurfCurrCmd->SetGuidance("         : 2 = Out only");
  qFlatSurfCurrCmd->SetGuidance("  wflag  :(Bool) Weighted");
  qFlatSurfCurrCmd->SetGuidance("  aflag  :(Bool) DivideByArea");
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
  //
  qFlatSurfFluxCmd = new G4UIcommand("/score/quantity/flatSurfFlux",this);
  qFlatSurfFluxCmd->SetGuidance("Flat surface Flux Scorer");
  qFlatSurfFluxCmd->
      SetGuidance("[usage] /score/quantiy/flatSurfFlux qname dflag");
  qFlatSurfFluxCmd->SetGuidance("  qname  :(String) scorer name");
  qFlatSurfFluxCmd->SetGuidance("  dflag  :(Int) direction flag");
  qFlatSurfFluxCmd->SetGuidance("         : 0 = Both In and Out");
  qFlatSurfFluxCmd->SetGuidance("         : 1 = In only");
  qFlatSurfFluxCmd->SetGuidance("         : 2 = Out only");
  param = new G4UIparameter("qname",'s',false);
  qFlatSurfFluxCmd->SetParameter(param);
  param = new G4UIparameter("dflag",'i',true);
  param->SetDefaultValue("0");
  qFlatSurfFluxCmd->SetParameter(param);
  //
  qSphereSurfCurrCmd = new G4UIcommand("/score/quantity/sphereSurfCurrent",this);
  qSphereSurfCurrCmd->SetGuidance("Sphere surface current Scorer");
  qSphereSurfCurrCmd->
      SetGuidance("[usage] /score/quantiy/sphereSurfCurrent qname dflag wflag aflag");
  qSphereSurfCurrCmd->SetGuidance("  qname  :(String) scorer name");
  qSphereSurfCurrCmd->SetGuidance("  dflag  :(Int) direction flag");
  qSphereSurfCurrCmd->SetGuidance("         : 0 = Both In and Out");
  qSphereSurfCurrCmd->SetGuidance("         : 1 = In only");
  qSphereSurfCurrCmd->SetGuidance("         : 2 = Out only");
  qSphereSurfCurrCmd->SetGuidance("  wflag  :(Bool) Weighted");
  qSphereSurfCurrCmd->SetGuidance("  aflag  :(Bool) DivideByArea");
  param = new G4UIparameter("qname",'s',false);
  qSphereSurfCurrCmd->SetParameter(param);
  param = new G4UIparameter("dflag",'i',true);
  param->SetDefaultValue("0");
  qSphereSurfCurrCmd->SetParameter(param);
  param = new G4UIparameter("wflag",'b',true);
  param->SetDefaultValue("true");
  qSphereSurfCurrCmd->SetParameter(param);
  param = new G4UIparameter("aflag",'b',true);
  param->SetDefaultValue("true");
  qSphereSurfCurrCmd->SetParameter(param);

  //
  qSphereSurfFluxCmd = new G4UIcommand("/score/quantity/sphereSurfFlux",this);
  qSphereSurfFluxCmd->SetGuidance("Sphere surface Flux Scorer");
  qSphereSurfFluxCmd->
      SetGuidance("[usage] /score/quantiy/sphereSurfFlux qname dflag");
  qSphereSurfFluxCmd->SetGuidance("  qname  :(String) scorer name");
  qSphereSurfFluxCmd->SetGuidance("  dflag  :(Int) direction flag");
  qSphereSurfFluxCmd->SetGuidance("         : 0 = Both In and Out");
  qSphereSurfFluxCmd->SetGuidance("         : 1 = In only");
  qSphereSurfFluxCmd->SetGuidance("         : 2 = Out only");
  param = new G4UIparameter("qname",'s',false);
  qSphereSurfFluxCmd->SetParameter(param);
  param = new G4UIparameter("dflag",'i',true);
  param->SetDefaultValue("0");
  qSphereSurfFluxCmd->SetParameter(param);

  //
  qCylSurfCurrCmd = new G4UIcommand("/score/quantity/cylSurfCurrent",this);
  qCylSurfCurrCmd->SetGuidance("Cylinder surface current Scorer");
  qCylSurfCurrCmd->
      SetGuidance("[usage] /score/quantiy/cylSurfCurrent qname dflag wflag aflag");
  qCylSurfCurrCmd->SetGuidance("  qname  :(String) scorer name");
  qCylSurfCurrCmd->SetGuidance("  dflag  :(Int) direction flag");
  qCylSurfCurrCmd->SetGuidance("         : 0 = Both In and Out");
  qCylSurfCurrCmd->SetGuidance("         : 1 = In only");
  qCylSurfCurrCmd->SetGuidance("         : 2 = Out only");
  qCylSurfCurrCmd->SetGuidance("  wflag  :(Bool) Weighted");
  qCylSurfCurrCmd->SetGuidance("  aflag  :(Bool) DivideByArea");
  param = new G4UIparameter("qname",'s',false);
  qCylSurfCurrCmd->SetParameter(param);
  param = new G4UIparameter("dflag",'i',true);
  param->SetDefaultValue("0");
  qCylSurfCurrCmd->SetParameter(param);
  param = new G4UIparameter("wflag",'b',true);
  param->SetDefaultValue("true");
  qCylSurfCurrCmd->SetParameter(param);
  param = new G4UIparameter("aflag",'b',true);
  param->SetDefaultValue("true");
  qCylSurfCurrCmd->SetParameter(param);
  //
  qCylSurfFluxCmd = new G4UIcommand("/score/quantity/cylSurfFlux",this);
  qCylSurfFluxCmd->SetGuidance("Cylinder surface Flux Scorer");
  qCylSurfFluxCmd->
      SetGuidance("[usage] /score/quantiy/cylSurfFlux qname dflag");
  qCylSurfFluxCmd->SetGuidance("  qname  :(String) scorer name");
  qCylSurfFluxCmd->SetGuidance("  dflag  :(Int) direction flag");
  qCylSurfFluxCmd->SetGuidance("         : 0 = Both In and Out");
  qCylSurfFluxCmd->SetGuidance("         : 1 = In only");
  qCylSurfFluxCmd->SetGuidance("         : 2 = Out only");
  param = new G4UIparameter("qname",'s',false);
  qCylSurfFluxCmd->SetParameter(param);
  param = new G4UIparameter("dflag",'i',true);
  param->SetDefaultValue("0");
  qCylSurfFluxCmd->SetParameter(param);
  //
  qNofCollisionCmd = new G4UIcommand("/score/quantity/nOfCollision",this);
  qNofCollisionCmd->SetGuidance("Number of Collision Scorer");
  qNofCollisionCmd->
      SetGuidance("[usage] /score/quantiy/nOfCollision qname wflag");
  qNofCollisionCmd->SetGuidance("  qname  :(String) scorer name");
  qNofCollisionCmd->SetGuidance("  wflag  :(Bool) Weighted");
  param = new G4UIparameter("qname",'s',false);
  qNofCollisionCmd->SetParameter(param);
  param = new G4UIparameter("wflag",'b',true);
  param->SetDefaultValue("false");
  qNofCollisionCmd->SetParameter(param);
  //
  qPopulationCmd = new G4UIcommand("/score/quantity/population",this);
  qPopulationCmd->SetGuidance("Population Scorer");
  qPopulationCmd->
      SetGuidance("[usage] /score/quantiy/population qname wflag");
  qPopulationCmd->SetGuidance("  qname  :(String) scorer name");
  qPopulationCmd->SetGuidance("  wflag  :(Bool) Weighted");
  param = new G4UIparameter("qname",'s',false);
  qPopulationCmd->SetParameter(param);
  param = new G4UIparameter("wflag",'b',true);
  param->SetDefaultValue("false");
  qPopulationCmd->SetParameter(param);

  //
  qTrackCountCmd = new G4UIcommand("/score/quantity/trackCounter",this);
  qTrackCountCmd->SetGuidance("Number of Track Counter Scorer");
  qTrackCountCmd->
      SetGuidance("[usage] /score/quantiy/trackCounter qname dflag wflag");
  qTrackCountCmd->SetGuidance("  qname  :(String) scorer name");
  qTrackCountCmd->SetGuidance("  dflag  :(Int) Direction");
  qTrackCountCmd->SetGuidance("         : 0 = Both In and Out");
  qTrackCountCmd->SetGuidance("         : 1 = In only");
  qTrackCountCmd->SetGuidance("         : 2 = Out only");
  qTrackCountCmd->SetGuidance("  wflag  :(Bool) Weighted");
  param = new G4UIparameter("qname",'s',false);
  qTrackCountCmd->SetParameter(param);
  param = new G4UIparameter("dflag",'i',true);
  param->SetDefaultValue("0");
  qTrackCountCmd->SetParameter(param);
  param = new G4UIparameter("wflag",'b',true);
  param->SetDefaultValue("false");
  qTrackCountCmd->SetParameter(param);

  //
  qTerminationCmd = new G4UIcommand("/score/quantity/termination",this);
  qTerminationCmd->SetGuidance("Number of Terminated tracks Scorer");
  qTerminationCmd->
      SetGuidance("[usage] /score/quantiy/termination qname wflag");
  qTerminationCmd->SetGuidance("  qname  :(String) scorer name");
  qTerminationCmd->SetGuidance("  wflag  :(Bool) Weighted");
  param = new G4UIparameter("qname",'s',false);
  qTerminationCmd->SetParameter(param);
  param = new G4UIparameter("wflag",'b',true);
  param->SetDefaultValue("false");
  qTerminationCmd->SetParameter(param);

  //
  // Filter commands 
  filterDir = new G4UIdirectory("/score/filter/");
  filterDir->SetGuidance("Filter for scoring");
  //
  fchargedCmd = new G4UIcmdWithAString("/score/filter/charged",this);
  fchargedCmd->SetGuidance("Charge filter ( charged )");
  fchargedCmd->SetParameterName("fname",false);
  //
  fneutralCmd = new G4UIcmdWithAString("/score/filter/neutral",this);
  fneutralCmd->SetGuidance("Charge filter ( Neutral )");
  fneutralCmd->SetParameterName("fname",false);
  //
  fkinECmd = new G4UIcommand("/score/filter/kinE",this);
  fkinECmd->SetGuidance("Kinetic Energy Filter");
  fkinECmd->SetGuidance("[usage] /score/filter/kinE fname Elow Ehigh unit");
  fkinECmd->SetGuidance("  fname     :(String) Filter Name ");
  fkinECmd->SetGuidance("  Elow      :(Double) Lower edge of kinetic energy");
  fkinECmd->SetGuidance("  Ehigh     :(Double) Higher edge of kinetic energy");
  fkinECmd->SetGuidance("  unit      :(String) unit of given kinetic energy");
  param = new G4UIparameter("fname",'s',false);
  fkinECmd->SetParameter(param);
  param = new G4UIparameter("elow",'d',true);
  param->SetDefaultValue("0.0");
  fkinECmd->SetParameter(param);
  param = new G4UIparameter("ehigh",'d',false);
  fkinECmd->SetParameter(param);
  G4String smax = DtoS(DBL_MAX);
  param->SetDefaultValue(smax);
  param = new G4UIparameter("unit",'s',false);
  param->SetDefaultValue("keV");
  fkinECmd->SetParameter(param);
  //
  fparticleCmd = new G4UIcommand("/score/filter/particle",this);
  fparticleCmd->SetGuidance("Touch particle filter into current quantity");
  fparticleCmd->SetGuidance("[usage] /score/filter/particle fname p0 .. pn");
  fparticleCmd->SetGuidance("  fname     :(String) Filter Name ");
  fparticleCmd->SetGuidance("  p0 .. pn  :(String) particle names");
  param = new G4UIparameter("fname",'s',false);
  fparticleCmd->SetParameter(param);
  param = new G4UIparameter("particlelist",'s',false);
  param->SetDefaultValue("");
  fparticleCmd->SetParameter(param);
  //
  //
  //
  fparticleKinECmd = new G4UIcommand("/score/filter/particleWithKinE",this);
  fparticleKinECmd->SetGuidance("Particle with kinetic energy filter");
  fparticleKinECmd->SetGuidance("[usage] /score/filter/particleWithKinE fname Elow Ehigh unit p0 .. pn");
  fparticleKinECmd->SetGuidance("  fname     :(String) Filter Name ");
  fparticleKinECmd->SetGuidance("  Elow      :(Double) Lower edge of kinetic energy");
  fparticleKinECmd->SetGuidance("  Ehigh     :(Double) Higher edge of kinetic energy");
  fparticleKinECmd->SetGuidance("  unit      :(String) unit of given kinetic energy");
  fparticleKinECmd->SetGuidance("  p0 .. pn  :(String) particle names");
  param = new G4UIparameter("fname",'s',false);
  fparticleKinECmd->SetParameter(param);
  param = new G4UIparameter("elow",'d',false);
  param->SetDefaultValue("0.0");
  fparticleKinECmd->SetParameter(param);
  param = new G4UIparameter("ehigh",'d',true);
  param->SetDefaultValue(smax);
  fparticleKinECmd->SetParameter(param);
  param = new G4UIparameter("unit",'s',true);
  param->SetDefaultValue("keV");
  fparticleKinECmd->SetParameter(param);
  param = new G4UIparameter("particlelist",'s',false);
  param->SetDefaultValue("");
  fparticleKinECmd->SetParameter(param);
  //
  //

}

G4ScoringMessenger::~G4ScoringMessenger()
{
    delete scoreDir;
    delete listCmd;
    delete verboseCmd;
    //
    delete           meshCreateDir;
    delete           meshBoxCreateCmd;
    delete           meshTubsCreateCmd;
    delete           meshSphereCreateCmd;
    //
    delete          meshDir;
    delete          meshOpnCmd;
    //
    delete    meshClsCmd;
    delete    meshActCmd;
    //
    delete  mBoxSizeCmd;
    delete  mTubsSizeCmd;
    delete  mSphereSizeCmd;
    //
    delete      mBinCmd;
    //
    delete   mTransDir;
    delete   mTResetCmd;
    delete   mTXyzCmd;
    delete   mRotDir;
    delete   mRResetCmd;
    delete   mRotXCmd;
    delete   mRotYCmd;
    delete   mRotZCmd;
    //
    delete     dumpCmd;
    delete     drawCmd;
    delete     dumpToFileCmd;
    //
    delete         quantityDir;
    delete         qTouchCmd;
    //
    delete    qCellChgCmd;
    delete    qCellFluxCmd;
    delete    qPassCellFluxCmd;
    delete    qeDepCmd;
    delete    qdoseDepCmd;
    delete    qnOfStepCmd;
    delete    qnOfSecondaryCmd;
    //
    delete          qTrackLengthCmd;
    delete          qPassCellCurrCmd;
    delete          qPassTrackLengthCmd;
    delete          qFlatSurfCurrCmd;
    delete          qFlatSurfFluxCmd;
    delete          qSphereSurfCurrCmd;
    delete          qSphereSurfFluxCmd;
    delete          qCylSurfCurrCmd;
    delete          qCylSurfFluxCmd;
    delete          qNofCollisionCmd;
    delete          qPopulationCmd;
    delete          qTrackCountCmd;
    delete          qTerminationCmd;
    //
    delete   filterDir;
    delete   fchargedCmd;
    delete   fneutralCmd;
    delete   fkinECmd;
    delete   fparticleCmd;
    delete   fparticleKinECmd;
}

void G4ScoringMessenger::SetNewValue(G4UIcommand * command,G4String newVal)
{
  if(command==listCmd) { 
      fSMan->List(); 
  } else if(command==dumpCmd) { 
      fSMan->Dump(); 
  } else if(command==drawCmd) { 
      G4Tokenizer next(newVal);
      G4String meshName = next();
      G4String psName = next();
      G4int axflg = StoI(next());
      fSMan->DrawMesh(meshName,psName,axflg);
  } else if(command==dumpToFileCmd) { 
      G4Tokenizer next(newVal);
      G4String meshName = next();
      G4String psName = next();
      G4String fileName = next();
      G4String option = next("\n");
      fSMan->DumpToFile(meshName, psName, fileName, option);
  } else if(command==verboseCmd) { 
      fSMan->SetVerboseLevel(verboseCmd->GetNewIntValue(newVal)); 
  } else if(command==meshBoxCreateCmd) {
      G4VScoringMesh*  mesh = fSMan->FindMesh(newVal);
      if ( !mesh ){
	  mesh = new G4ScoringBox(newVal);
	  fSMan->RegisterScoringMesh(mesh);
      }else{
	  /////////////////////G4Exception("G4ScroingMessenger:: Mesh has already existed. Error!");
          G4cerr << "Scoring mesh <" << newVal << "> already exists. Command ignored." << G4endl;
      }
  } else if(command==meshOpnCmd) {
      G4VScoringMesh* currentmesh = fSMan->GetCurrentMesh(); 
      if ( currentmesh ){
 	  /////////////////////G4Exception("G4ScroingMessenger:: Close current mesh first!. Error!");
          G4cerr << "Mesh <" << currentmesh->GetWorldName() << "> is still open. Close it first. Command ignored." << G4endl;
      } else {
	G4VScoringMesh* mesh = fSMan->FindMesh(newVal); 
	if ( !mesh ){
 	  /////////////////////G4Exception("G4ScroingMessenger:: Mesh has not existed. Error!");
          G4cerr << "Scoring mesh <" << newVal << "> does not exist. Command ignored." << G4endl;
	} else {
	  fSMan->SetCurrentMesh(mesh);
	}
      }
  } else if(command==meshClsCmd) {
      fSMan->CloseCurrentMesh();
  } else {
      // Tokens
      G4TokenVec token;
      FillTokenVec(newVal,token);
      //
      // Get Current Mesh
      //
      G4VScoringMesh* mesh = fSMan->GetCurrentMesh();
      //
      // Commands for Current Mesh
      if ( mesh ){
	  // 
	  // Mesh Geometry
	  //
	  if(command==meshActCmd) {
	      mesh->Activate(meshActCmd->GetNewBoolValue(newVal)); 
	  } else if(command==mBoxSizeCmd) {
	      MeshShape shape = mesh->GetShape();
	      if ( shape == boxMesh ){
		  G4ThreeVector size = mBoxSizeCmd->GetNew3VectorValue(newVal);
		  G4double vsize[3];
		  vsize[0] = size.x();
		  vsize[1] = size.y();
		  vsize[2] = size.z();
		  mesh->SetSize(vsize);
	      } else {
		  ////////////////////G4Exception("G4ScroingMessenger:: Mesh is not Box type. Error!");
                 G4cerr << "This mesh is not Box. Command ignored." << G4endl;
	      }
	  } else if(command==mBinCmd) {
	      MeshBinCommand(mesh,token);
	  } else if(command==mTResetCmd) {
	      G4double centerPosition[3] ={ 0., 0., 0.};
	      mesh->SetCenterPosition(centerPosition);
	  } else if(command==mTXyzCmd) {
	      G4ThreeVector xyz = mTXyzCmd->GetNew3VectorValue(newVal);
	      G4double centerPosition[3];
	      centerPosition[0] = xyz.x();
	      centerPosition[1] = xyz.y();
	      centerPosition[2] = xyz.z();
	      mesh->SetCenterPosition(centerPosition);
	  } else if(command==mRResetCmd) {
	  } else if(command==mRotXCmd) {
	      G4double value = mRotXCmd->GetNewDoubleValue(newVal);
	      mesh->RotateX(value);
	  } else if(command==mRotYCmd) {
	      G4double value = mRotYCmd->GetNewDoubleValue(newVal);
	      mesh->RotateY(value);
	  } else if(command==mRotZCmd) {
	      G4double value = mRotZCmd->GetNewDoubleValue(newVal);
	      mesh->RotateZ(value);
	  } else if(command==qTouchCmd) {
	      mesh->SetCurrentPrimitiveScorer(newVal);
	  //
	  // Quantity
          //    
	  } else if(command== qCellChgCmd) {
	    if(!mesh->FindPrimitiveScorer(newVal)) {
	      mesh->SetPrimitiveScorer(new G4PSCellCharge3D(newVal));
	    } else {
	      G4cout << " Quantity name, \"" << newVal << "\", is already existing." << G4endl;
	      mesh->SetNullToCurrentPrimitiveScorer();
	    }
	  } else if(command== qCellFluxCmd) {
	    if(!mesh->FindPrimitiveScorer(newVal)) { 
	      mesh->SetPrimitiveScorer(new G4PSCellFlux3D(newVal));
	    } else {
	      G4cout << " Quantity name, \"" << newVal << "\", is already existing." << G4endl;
	      mesh->SetNullToCurrentPrimitiveScorer();
	    }
	  } else if(command== qPassCellFluxCmd) {
	    if(!mesh->FindPrimitiveScorer(newVal))  {
	      mesh->SetPrimitiveScorer(new G4PSPassageCellFlux3D(newVal));
	    } else {
	      G4cout << " Quantity name, \"" << newVal << "\", is already existing." << G4endl;
	      mesh->SetNullToCurrentPrimitiveScorer();
	    }
	  } else if(command==qeDepCmd) {
	    if(!mesh->FindPrimitiveScorer(newVal))  {
	      mesh->SetPrimitiveScorer(new G4PSEnergyDeposit3D(newVal));
	    } else {
	      G4cout << " Quantity name, \"" << newVal << "\", is already existing." << G4endl;
	      mesh->SetNullToCurrentPrimitiveScorer();
	    }
	  } else if(command== qdoseDepCmd) {
	    if(!mesh->FindPrimitiveScorer(newVal))  {
	      mesh->SetPrimitiveScorer(new G4PSDoseDeposit3D(newVal));
	    } else {
	      G4cout << " Quantity name, \"" << newVal << "\", is already existing." << G4endl;
	      mesh->SetNullToCurrentPrimitiveScorer();
	    }
	  } else if(command== qnOfStepCmd) {
	    if(!mesh->FindPrimitiveScorer(newVal))  {
	      mesh->SetPrimitiveScorer(new G4PSNofStep3D(newVal));
	    } else {
	      G4cout << " Quantity name, \"" << newVal << "\", is already existing." << G4endl;
	      mesh->SetNullToCurrentPrimitiveScorer();
	    }
	  } else if(command== qnOfSecondaryCmd) {
	    if(!mesh->FindPrimitiveScorer(newVal))  {
	      mesh->SetPrimitiveScorer(new G4PSNofSecondary3D(newVal));
	    } else {
	      G4cout << " Quantity name, \"" << newVal << "\", is already existing." << G4endl;
	      mesh->SetNullToCurrentPrimitiveScorer();
	    }
	  } else if(command== qTrackLengthCmd) {
	    if(!mesh->FindPrimitiveScorer(newVal)) {
	      G4PSTrackLength3D* ps = new G4PSTrackLength3D(token[0]);
	      ps->Weighted(StoB(token[1]));
	      ps->MultiplyKineticEnergy(StoB(token[2]));
	      ps->DivideByVelocity(StoB(token[3]));
	      mesh->SetPrimitiveScorer(ps);
	    } else {
	      G4cout << " Quantity name, \"" << newVal << "\", is already existing." << G4endl;
	      mesh->SetNullToCurrentPrimitiveScorer();
	    }
          } else if(command== qPassCellCurrCmd){
	    if(!mesh->FindPrimitiveScorer(newVal)) {
	      G4PSPassageCellCurrent* ps = new G4PSPassageCellCurrent3D(token[0]);
	      ps->Weighted(StoB(token[1]));
	      mesh->SetPrimitiveScorer(ps);
	    } else {
	      G4cout << " Quantity name, \"" << newVal << "\", is already existing." << G4endl;
	      mesh->SetNullToCurrentPrimitiveScorer();
	    }
          } else if(command== qPassTrackLengthCmd){
	    if(!mesh->FindPrimitiveScorer(newVal)) {
	      G4PSPassageTrackLength* ps = new G4PSPassageTrackLength3D(token[0]);
	      ps->Weighted(StoB(token[1]));
	      mesh->SetPrimitiveScorer(ps);
	    } else {
	      G4cout << " Quantity name, \"" << newVal << "\", is already existing." << G4endl;
	      mesh->SetNullToCurrentPrimitiveScorer();
	    }
          } else if(command== qFlatSurfCurrCmd){
	    if(!mesh->FindPrimitiveScorer(newVal)) {
	      G4PSFlatSurfaceCurrent3D* ps = 
		new G4PSFlatSurfaceCurrent3D(token[0],StoI(token[1]));
	      ps->Weighted(StoB(token[2]));
	      ps->DivideByArea(StoB(token[3]));
	      mesh->SetPrimitiveScorer(ps);
	    } else {
	      G4cout << " Quantity name, \"" << newVal << "\", is already existing." << G4endl;
	      mesh->SetNullToCurrentPrimitiveScorer();
	    }
          } else if(command== qFlatSurfFluxCmd){
	    if(!mesh->FindPrimitiveScorer(newVal)) {
	      mesh->SetPrimitiveScorer(
				       new G4PSFlatSurfaceFlux3D(token[0],StoI(token[1])));
	    } else {
	      G4cout << " Quantity name, \"" << newVal << "\", is already existing." << G4endl;
	      mesh->SetNullToCurrentPrimitiveScorer();
	    }
          } else if(command== qSphereSurfCurrCmd){
	    if(!mesh->FindPrimitiveScorer(newVal)) {
	      G4PSSphereSurfaceCurrent3D* ps = 
		new G4PSSphereSurfaceCurrent3D(token[0],StoI(token[1]));
	      ps->Weighted(StoB(token[2]));
	      ps->DivideByArea(StoB(token[3]));
	      mesh->SetPrimitiveScorer(ps);
	    } else {
	      G4cout << " Quantity name, \"" << newVal << "\", is already existing." << G4endl;
	      mesh->SetNullToCurrentPrimitiveScorer();
	    }
	  } else if(command== qSphereSurfFluxCmd){
	    if(!mesh->FindPrimitiveScorer(newVal)) {
	      mesh->SetPrimitiveScorer(
				       new G4PSSphereSurfaceFlux3D(token[0], StoI(token[1])));
	    } else {
	      G4cout << " Quantity name, \"" << newVal << "\", is already existing." << G4endl;
	      mesh->SetNullToCurrentPrimitiveScorer();
	    }
          } else if(command== qCylSurfCurrCmd){
	    if(!mesh->FindPrimitiveScorer(newVal)) {
	      G4PSCylinderSurfaceCurrent3D* ps = 
		new G4PSCylinderSurfaceCurrent3D(token[0],StoI(token[1]));
	      ps->Weighted(StoB(token[2]));
	      ps->DivideByArea(StoB(token[3]));
	      mesh->SetPrimitiveScorer(ps);
	    } else {
	      G4cout << " Quantity name, \"" << newVal << "\", is already existing." << G4endl;
	      mesh->SetNullToCurrentPrimitiveScorer();
	    }
          } else if(command== qCylSurfFluxCmd){
	    if(!mesh->FindPrimitiveScorer(newVal)) {
	      mesh->SetPrimitiveScorer(
				       new G4PSCylinderSurfaceFlux3D(token[0], StoI(token[1])));
	    } else {
	      G4cout << " Quantity name, \"" << newVal << "\", is already existing." << G4endl;
	      mesh->SetNullToCurrentPrimitiveScorer();
	    }
          } else if(command== qNofCollisionCmd){
	    if(!mesh->FindPrimitiveScorer(newVal)) {
	      G4PSNofCollision3D* ps =new G4PSNofCollision3D(token[0]); 
	      ps->Weighted(StoB(token[1]));
	      mesh->SetPrimitiveScorer(ps);
	    } else {
	      G4cout << " Quantity name, \"" << newVal << "\", is already existing." << G4endl;
	      mesh->SetNullToCurrentPrimitiveScorer();
	    }
          } else if(command== qPopulationCmd){
	    if(!mesh->FindPrimitiveScorer(newVal)) {
	      G4PSPopulation3D* ps =new G4PSPopulation3D(token[0]); 
	      ps->Weighted(StoB(token[1]));
	      mesh->SetPrimitiveScorer(ps);
	    } else {
	      G4cout << " Quantity name, \"" << newVal << "\", is already existing." << G4endl;
	      mesh->SetNullToCurrentPrimitiveScorer();
	    }
          } else if(command== qTrackCountCmd){
	    if(!mesh->FindPrimitiveScorer(newVal)) {
	      G4PSTrackCounter3D* ps =new G4PSTrackCounter3D(token[0],StoI(token[1])); 
	      mesh->SetPrimitiveScorer(ps);
	    } else {
	      G4cout << " Quantity name, \"" << newVal << "\", is already existing." << G4endl;
	      mesh->SetNullToCurrentPrimitiveScorer();
	    }
          } else if(command== qTerminationCmd){
	    if(!mesh->FindPrimitiveScorer(newVal)) {
	      G4PSTermination3D* ps =new G4PSTermination3D(token[0]); 
	      ps->Weighted(StoB(token[1]));
	      mesh->SetPrimitiveScorer(ps);
	    } else {
	      G4cout << " Quantity name, \"" << newVal << "\", is already existing." << G4endl;
	      mesh->SetNullToCurrentPrimitiveScorer();
	    }

	    //
	    // Filters 
	    // 
	  }else if(command== fchargedCmd){
	    if(!mesh->IsCurrentPrimitiveScorerNull()) {
	      mesh->SetFilter(new G4SDChargedFilter(token[0])); 
	    } else {
	      G4cout << " Filter, \"" << token[0] << "\", is not registered." << G4endl;
	    }
          }else if(command== fneutralCmd){
	    if(!mesh->IsCurrentPrimitiveScorerNull()) {
	      mesh->SetFilter(new G4SDNeutralFilter(token[0])); 
	    } else {
	      G4cout << " Filter, \"" << token[0] << "\", is not registered." << G4endl;
	    }
          }else if(command== fkinECmd){
	    if(!mesh->IsCurrentPrimitiveScorerNull()) {
	      G4String& name = token[0];
	      G4double elow  = StoD(token[1]);
	      G4double ehigh = StoD(token[2]);
	      mesh->SetFilter(new G4SDKineticEnergyFilter(name,elow,ehigh));
	    } else {
	      G4cout << " Filter, \"" << token[0] << "\", is not registered." << G4endl;
	    }
          }else if(command== fparticleKinECmd){
	    if(!mesh->IsCurrentPrimitiveScorerNull()) {
	      FParticleWithEnergyCommand(mesh,token); 
	    } else {
	      G4cout << " Filter, \"" << token[0] << "\", is not registered." << G4endl;
	    }
	  } else if(command==fparticleCmd) {
	    if(!mesh->IsCurrentPrimitiveScorerNull()) {
	      FParticleCommand(mesh,token);
	    } else {
	      G4cout << " Filter, \"" << token[0] << "\", is not registered." << G4endl;
	    }
	  }
      }else{
///////////////////	  G4Exception("G4ScroingMessenger:: Current Mesh has not opened. Error!");
        G4cerr << "No mesh is currently open. Open/create a mesh first. Command ignored." << G4endl;
      }
  }
}

G4String G4ScoringMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String val;
  if(command==verboseCmd)
  { val = verboseCmd->ConvertToString(fSMan->GetVerboseLevel()); }

  return val;
}

void G4ScoringMessenger::FillTokenVec(G4String newValues, G4TokenVec& token){

    G4Tokenizer next(newValues);
    G4String val;
    while ( !(val = next()).isNull() ) {
	//if ( val.first('"')==0 ) {
	//     val.remove(0,1);
	//}
	//if ( (val.last('"')) == (G4int)(val.length()-1) ){
	//    val.remove(val.length()-1,1);
	//}
	token.push_back(val);

//      G4cout << "@GetToken:"
//	       << val
//	       << "::"
//	       << G4endl;
    }
}


void G4ScoringMessenger::MeshBinCommand(G4VScoringMesh* mesh,G4TokenVec& token){
    G4int Ni = StoI(token[0]);
    G4int Nj = StoI(token[1]);
    G4int Nk = StoI(token[2]);
    G4int nSegment[3];
    nSegment[0] = Ni;
    nSegment[1] = Nj;
    nSegment[2] = Nk;
    //
    mesh->SetNumberOfSegments(nSegment);
    //
    //
    /*
    G4int iAxis = 3;
    G4String Axis = token[3];
    if ( ! Axis.isNull() ){
       iAxis = StoI(Axis);
       //
       //==== Implementation for variable bin size Here
       //
       //  .........
    }
    */
}

void G4ScoringMessenger::FParticleCommand(G4VScoringMesh* mesh, G4TokenVec& token){
    //
    // Filter name
    G4String name = token[0];
    //
    // particle list
    std::vector<G4String> pnames;
    for ( G4int i = 1; i<(G4int)token.size(); i++){
	pnames.push_back(token[i]);
    }
    //
    // Attach Filter
    mesh->SetFilter(new G4SDParticleFilter(name,pnames));
}    

void G4ScoringMessenger::FParticleWithEnergyCommand(G4VScoringMesh* mesh,G4TokenVec& token){
    G4String& name = token[0];
    G4double  elow = StoD(token[1]);
    G4double  ehigh= StoD(token[2]);
    G4double  unitVal = G4UnitDefinition::GetValueOf(token[3]);
    G4SDParticleWithEnergyFilter* filter = 
	new G4SDParticleWithEnergyFilter(name,elow*unitVal,ehigh*unitVal);
    for ( G4int i = 4; i < (G4int)token.size(); i++){
	filter->add(token[i]);
    }
    mesh->SetFilter(filter);
}
 
