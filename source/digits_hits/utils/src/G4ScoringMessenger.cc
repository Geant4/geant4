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
// $Id: G4ScoringMessenger.cc,v 1.11 2007-08-29 05:14:40 taso Exp $
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

#include "G4SDParticleFilter.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcommand.hh"
#include "G4Tokenizer.hh"

G4ScoringMessenger::G4ScoringMessenger(G4ScoringManager* SManager)
:fSMan(SManager)
{
  scoreDir = new G4UIdirectory("/score/");
  scoreDir->SetGuidance("Interactive scoring commands.");

  listCmd = new G4UIcmdWithoutParameter("/score/list",this);
  listCmd->SetGuidance("List scoring worlds.");

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
  mBinCmd->SetGuidance("  P1..Pn-1  :(double) paramter from P1 to Pn-1 for division.");
  G4UIparameter* param;
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

  // Dump Scoring result
  dumpCmd = new G4UIcmdWithAString("/score/dump",this);
  dumpCmd->SetGuidance("Dump scorer results ");
  dumpCmd->SetParameterName("qname",true);

  //
  // Quantity commands
  quantityDir = new G4UIdirectory("/score/quantity/");
  quantityDir->SetGuidance("Scoring quantity of the mesh");
  //
  qTouchCmd= new G4UIcmdWithAString("/score/quantity/touch",this);
  qTouchCmd->SetGuidance("Assign previously defined quantity to current quantity.");
  qTouchCmd->SetParameterName("qName",false);
  //
  qeDepCmd = new G4UIcmdWithAString("/score/quantity/eDep",this);
  qeDepCmd->SetGuidance("Energy Deposit Scorer");
  qeDepCmd->SetParameterName("name",false);
  qCellChgCmd  = new G4UIcmdWithAString("/score/quantity/cellCharge",this);
  qCellChgCmd->SetGuidance("Cell Charge Scorer");
  qCellChgCmd->SetParameterName("name",false);
  qCellFluxCmd = new G4UIcmdWithAString("/score/quantity/cellFlux",this);
  qCellFluxCmd->SetGuidance("Cell Flux Scorer");
  qCellFluxCmd->SetParameterName("name",false);
  qPassCellFluxCmd = new G4UIcmdWithAString("/score/quantity/passageCellFlux",this);
  qPassCellFluxCmd->SetGuidance("Passage Cell Flux Scorer");
  qPassCellFluxCmd->SetParameterName("name",false);
  qdoseDepCmd = new G4UIcmdWithAString("/score/quantity/doseDeposit",this);
  qdoseDepCmd->SetGuidance("Dose Deposit Scorer");
  qdoseDepCmd->SetParameterName("name",false);
  qnOfStepCmd = new G4UIcmdWithAString("/score/quantity/nOfStep",this);
  qnOfStepCmd->SetGuidance("Number of Step Scorer ");
  qnOfStepCmd->SetParameterName("name",false);
  qnOfSecondaryCmd = new G4UIcmdWithAString("/score/quantity/nOfSecondary",this);
  qnOfSecondaryCmd->SetGuidance("Number of Secondary Scorer ");
  qnOfSecondaryCmd->SetParameterName("name",false);
  //
  //
  // Filter commands 
  filterDir = new G4UIdirectory("/score/filter/");
  filterDir->SetGuidance("Filter for scoring");
  //
  fparticleCmd = new G4UIcommand("/score/filter/particle",this);
  fparticleCmd->SetGuidance("Touch particle filter into current quantity");
  fparticleCmd->SetGuidance("[usage] /score/filter/particle fname p0 .. pn");
  fparticleCmd->SetGuidance("  fname     :(String) Filter Name ");
  fparticleCmd->SetGuidance("  p0 .. pn  :(String) particle names ");
  G4UIparameter* particleParam;
  particleParam = new G4UIparameter("name",'s',false);
  fparticleCmd->SetParameter(particleParam);
  particleParam = new G4UIparameter("particlelist",'s',false);
  fparticleCmd->SetParameter(particleParam);
  //
  //
}

G4ScoringMessenger::~G4ScoringMessenger()
{
  delete listCmd;
  delete scoreDir;

  delete meshBoxCreateCmd;
  delete mBoxSizeCmd;

  delete meshTubsCreateCmd;
  delete meshSphereCreateCmd;

  delete meshClsCmd;
  delete meshActCmd;

  delete mBinCmd;

  delete mTResetCmd;
  delete mTXyzCmd;
  delete mTransDir;
  delete mRResetCmd;
  delete mRotXCmd;
  delete mRotYCmd;
  delete mRotZCmd;
  delete mRotDir;
  delete meshDir;
 
  delete   qCellChgCmd;
  delete   qCellFluxCmd;
  delete   qPassCellFluxCmd;
  delete   qeDepCmd;
  delete   qdoseDepCmd;
  delete   qnOfStepCmd;
  delete   qnOfSecondaryCmd;

  delete qTouchCmd;
  delete quantityDir;

  delete fparticleCmd;

  delete filterDir;

}

void G4ScoringMessenger::SetNewValue(G4UIcommand * command,G4String newVal)
{
  if(command==listCmd) { 
      fSMan->List(); 
  } else if(command==dumpCmd) { 
      fSMan->Dump(); 
  } else if(command==verboseCmd) { 
      fSMan->SetVerboseLevel(verboseCmd->GetNewIntValue(newVal)); 
  } else if(command==meshBoxCreateCmd) {
      G4VScoringMesh*  mesh = fSMan->FindMesh(newVal);
      if ( !mesh ){
	  mesh = new G4ScoringBox(newVal);
	  fSMan->RegisterScoringMesh(mesh);
      }else{
	  G4Exception("G4ScroingMessenger:: Mesh has already existed. Error!");
      }
  } else if(command==meshOpnCmd) {
      G4VScoringMesh* currentmesh = fSMan->GetCurrentMesh(); 
      if ( currentmesh ){
 	  G4Exception("G4ScroingMessenger:: Close current mesh first!. Error!");
      }
      G4VScoringMesh* mesh = fSMan->FindMesh(newVal); 
      if ( !mesh ){
 	  G4Exception("G4ScroingMessenger:: Mesh has not existed. Error!");
      }
  } else {
      //
      // Get Current Mesh
      //
      G4VScoringMesh* mesh = fSMan->GetCurrentMesh();
      //
      // Commands for Current Mesh
      if ( mesh ){
	  if(command==meshClsCmd) {
	      fSMan->CloseCurrentMesh();
	  } else if(command==meshActCmd) {
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
		  G4Exception("G4ScroingMessenger:: Mesh is not Box type. Error!");
	      }
	  } else if(command==mBinCmd) {
	      MeshBinCommand(newVal);
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
	  } else if(command== qCellChgCmd) {
	      mesh->SetPrimitiveScorer(new G4PSCellCharge3D(newVal));
	  } else if(command== qCellFluxCmd) {
	      mesh->SetPrimitiveScorer(new G4PSCellFlux3D(newVal));
	  } else if(command== qPassCellFluxCmd) {
	      mesh->SetPrimitiveScorer(new G4PSPassageCellFlux3D(newVal));
	  } else if(command==qeDepCmd) {
	      mesh->SetPrimitiveScorer(new G4PSEnergyDeposit3D(newVal));
	  } else if(command== qdoseDepCmd) {
	      mesh->SetPrimitiveScorer(new G4PSDoseDeposit3D(newVal));
	  } else if(command== qnOfStepCmd) {
	      mesh->SetPrimitiveScorer(new G4PSNofStep3D(newVal));
	  } else if(command== qnOfSecondaryCmd) {
	      mesh->SetPrimitiveScorer(new G4PSNofSecondary3D(newVal));
	  } else if(command==fparticleCmd) {
	      FParticleCommand(newVal);
	  }
      }else{
	  G4Exception("G4ScroingMessenger:: Current Mesh has not opened. Error!");
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

void G4ScoringMessenger::MeshBinCommand(G4String newVal){
    G4Tokenizer next(newVal);
    G4int Ni = StoI(next());
    G4int Nj = StoI(next());
    G4int Nk = StoI(next());
    G4int nSegment[3];
    nSegment[0] = Ni;
    nSegment[1] = Nj;
    nSegment[2] = Nk;
    //
    G4VScoringMesh* mesh = fSMan->GetCurrentMesh();
    if ( !mesh ){
	G4Exception("G4ScroingMessenger:: Current Mesh has not opened. Error!");
    }
    mesh->SetNumberOfSegments(nSegment);
    //
    //
    /*
    G4int iAxis = 3;
    G4String Axis = next();
    if ( Axis.isNull() ) {
    } else {
	iAxis = StoI(Axis);
	//
	//==== Implementation for variable bin size Here
	//
	//  .........
    }
    */
}

void G4ScoringMessenger::FParticleCommand(G4String newVal){
    G4Tokenizer next(newVal);
    //
    G4String name = next();
    //
    G4String p = next();
    std::vector<G4String> pnames;
    pnames.push_back(p);
    //G4cout << " XXX " << name << "   XXX  " << p << G4endl;
    /*
    if ( ! p.isNull() ) {
	pnames.push_back(p);
	do {
	    G4String p = next();
	    pnames.push_back(p);
	} while (!p.isNull());
    }
    */
    G4VScoringMesh* mesh = fSMan->GetCurrentMesh();
    if ( !mesh ){
	G4Exception("G4ScroingMessenger:: Current Mesh has not opened. Error!");
    }
    mesh->SetFilter(new G4SDParticleFilter(name,pnames));

}    

