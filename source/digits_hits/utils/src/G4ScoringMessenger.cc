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
// $Id: G4ScoringMessenger.cc,v 1.7 2007-08-28 10:14:26 taso Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ---------------------------------------------------------------------

#include "G4ScoringMessenger.hh"
#include "G4ScoringManager.hh"
#include "G4VScoringMesh.hh"
#include "G4ScoringBox.hh"

#include "G4PSEnergyDeposit3D.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
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
  meshBoxCreateCmd->SetGuidance("Open scoring mesh.");
  meshBoxCreateCmd->SetParameterName("MeshName",false);
  //
  meshTubsCreateCmd = new G4UIcmdWithAString("/score/create/tubsMesh",this);
  meshTubsCreateCmd->SetGuidance("Open scoring mesh.");
  meshTubsCreateCmd->SetParameterName("MeshName",false);
  //
  meshSphereCreateCmd = new G4UIcmdWithAString("/score/create/sphereMesh",this);
  meshSphereCreateCmd->SetGuidance("Open scoring mesh.");
  meshSphereCreateCmd->SetParameterName("MeshName",false);
  //
  meshOpnCmd = new G4UIcmdWithAString("/score/open",this);
  meshOpnCmd->SetGuidance("Open scoring mesh.");
  meshOpnCmd->SetParameterName("MeshName",false);
  //
  meshClsCmd = new G4UIcmdWithoutParameter("/score/close",this);
  meshClsCmd->SetGuidance("Close scoring mesh.");
  //
  meshDelCmd = new G4UIcmdWithAString("/score/mesh/delete",this);
  meshDelCmd->SetGuidance("Delete scoring mesh.");
  meshDelCmd->SetParameterName("MeshName",false);
  //
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
  param = new G4UIparameter("Nj",'i',false);
  param->SetDefaultValue("1");
  param = new G4UIparameter("Nk",'i',false);
  param->SetDefaultValue("1");
  param = new G4UIparameter("Axis",'i',false);
  param->SetDefaultValue("2");
  //
  //   Placement command
  mTransDir = new G4UIdirectory("/score/mesh/translate/");
  mTransDir->SetGuidance("Placement of scoring mesh");
  //
  mTSetCmd = new G4UIcmdWith3VectorAndUnit("/score/mesh/translate/set",this);
  mTSetCmd->SetGuidance("Set translation of scoring mesh placement.");
  mTSetCmd->SetParameterName("X","Y","Z",false,false);
  mTSetCmd->SetDefaultUnit("mm");
  mTAddCmd = new G4UIcmdWith3VectorAndUnit("/score/mesh/translate/add",this);
  mTAddCmd->SetGuidance("Add translation to the current scoring mesh position.");
  mTAddCmd->SetParameterName("DX","DY","DZ",false,false);
  mTAddCmd->SetDefaultUnit("mm");
  //
  mRotDir = new G4UIdirectory("/score/mesh/rotate/");
  mRotDir->SetGuidance("Placement of scoring mesh");
  //
  mRSetCmd = new G4UIcmdWith3VectorAndUnit("/score/mesh/rotate/set",this);
  mRSetCmd->SetGuidance("Set rotation of scoring mesh placement.");
  mRSetCmd->SetParameterName("Rx","Ry","Rz",false,false);
  mRSetCmd->SetDefaultUnit("deg");
  mRAddCmd = new G4UIcmdWith3VectorAndUnit("/score/mesh/rotate/add",this);
  mRAddCmd->SetGuidance("Add rotation to the current scoring mesh position.");
  mRAddCmd->SetParameterName("DRx","DRy","DRz",false,false);
  mRAddCmd->SetDefaultUnit("deg");

  //
  // Quantity commands
  quantityDir = new G4UIdirectory("/score/quantity/");
  quantityDir->SetGuidance("Scoring quantity of the mesh");
  //
  qeDepCmd = new G4UIcmdWithoutParameter("/score/quantity/eDep",this);
  qeDepCmd->SetGuidance("Energy Deposit Scorer");
  //
  //
  // Filter commands 
  filterDir = new G4UIdirectory("/score/filter/");
  filterDir->SetGuidance("Filter for scoring");
  //
  filterOpnCmd = new G4UIcmdWithAString("/score/filter/open",this);
  filterOpnCmd->SetGuidance("Open filter of scoring mesh.");
  filterOpnCmd->SetParameterName("FilterName",false);
  filterClsCmd = new G4UIcmdWithoutParameter("/score/filter/close",this);
  filterClsCmd->SetGuidance("Close filter of scoring mesh.");
  filterDelCmd = new G4UIcmdWithAString("/score/filter/delete",this);
  filterDelCmd->SetGuidance("Delete filter of scoring mesh.");
  filterDelCmd->SetParameterName("FilterName",false);
  fAttachCmd= new G4UIcmdWithAString("/score/filter/attach",this);
  fAttachCmd->SetGuidance("Attach previously defined filter into a mesh or quantity,");
  fAttachCmd->SetParameterName("FilterName",false);
  //

  fdumpCmd = new G4UIcmdWithAString("/score/dump",this);
  fdumpCmd->SetGuidance("Dump scorer results ");

}

G4ScoringMessenger::~G4ScoringMessenger()
{
  delete listCmd;
  delete scoreDir;

  delete meshBoxCreateCmd;
  delete meshTubsCreateCmd;
  delete meshSphereCreateCmd;
  delete meshClsCmd;
  delete meshDelCmd;
  delete mBoxSizeCmd;
  delete mBinCmd;
  delete mBinDir;
  delete mTSetCmd;
  delete mTAddCmd;
  delete mTransDir;
  delete mRSetCmd;
  delete mRAddCmd;
  delete mRotDir;
  delete meshDir;

  delete qeDepCmd;
  delete quantityDir;

  delete filterOpnCmd;
  delete filterClsCmd;
  delete filterDelCmd;
  delete fAttachCmd;
  delete filterDir;

}

void G4ScoringMessenger::SetNewValue(G4UIcommand * command,G4String newVal)
{
  if(command==listCmd) { 
      fSMan->List(); 
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
      G4VScoringMesh* mesh = fSMan->FindMesh(newVal); 
      if ( !mesh ){
 	  G4Exception("G4ScroingMessenger:: Mesh has not existed. Error!");
      }
  } else if(command==meshClsCmd) {
      G4VScoringMesh* mesh = fSMan->GetCurrentMesh();
      if ( !mesh ){
 	  G4Exception("G4ScroingMessenger:: Mesh has not existed. Error!");
      }else{
	  fSMan->CloseCurrentMesh();
      }
  } else if(command==meshDelCmd) {
  } else if(command==mBoxSizeCmd) {
      G4VScoringMesh* mesh = fSMan->GetCurrentMesh();
      if ( !mesh ){
 	  G4Exception("G4ScroingMessenger:: Current Mesh has not opened. Error!");
      }
      MeshShape shape = mesh->GetShape();
      if ( shape == boxMesh ){
	  G4ThreeVector size = mBoxSizeCmd->GetNew3VectorValue(newVal);
	  G4double vsize[3];
	  vsize[0] = size.x();
	  vsize[1] = size.y();
	  vsize[2] = size.z();
	  mesh->SetSize(vsize);
      } else {
 	  G4Exception("G4ScroingMessenger:: Current Mesh is not Box type. Error!");
      }
  } else if(command==mBinCmd) {
      MeshBinCommand(newVal);
  } else if(command==mTSetCmd) {
  } else if(command==mTAddCmd) {
  } else if(command==mRSetCmd) {
  } else if(command==mRAddCmd) {
  } else if(command==qeDepCmd) {
      G4VScoringMesh* mesh = fSMan->GetCurrentMesh();
      if ( !mesh ){
 	  G4Exception("G4ScroingMessenger:: Current Mesh has not opened. Error!");
      }
      mesh->SetPrimitiveScorer(new G4PSEnergyDeposit3D(newVal));
  } else if(command==filterOpnCmd) {
  } else if(command==filterClsCmd) {
  } else if(command==filterDelCmd) {
  } else if(command==fAttachCmd) {
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
