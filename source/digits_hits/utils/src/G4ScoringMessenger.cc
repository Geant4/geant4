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
// $Id: G4ScoringMessenger.cc,v 1.3 2007-08-14 16:51:10 taso Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ---------------------------------------------------------------------

#include "G4ScoringMessenger.hh"
#include "G4ScoringManager.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcommand.hh"

G4ScoringMessenger::G4ScoringMessenger(G4ScoringManager* SManager):fSMan(SManager)
{
  scoreDir = new G4UIdirectory("/score/");
  scoreDir->SetGuidance("Interactive scoring commands.");

  listCmd = new G4UIcmdWithoutParameter("/score/list",this);
  listCmd->SetGuidance("List scoring worlds.");

  verboseCmd = new G4UIcmdWithAnInteger("/score/verbose",this);
  verboseCmd->SetGuidance("Verbosity");

  //
  // Mesh commands
  meshDir  = new G4UIdirectory("/score/mesh/");
  meshDir->SetGuidance("Scoring mesh commands.");
  //
  meshOpnCmd = new G4UIcmdWithAString("/score/mesh/open",this);
  meshOpnCmd->SetGuidance("Open scoring mesh.");
  meshOpnCmd->SetParameterName("MeshName",false);
  meshClsCmd = new G4UIcmdWithoutParameter("/score/mesh/close",this);
  meshClsCmd->SetGuidance("Close scoring mesh.");
  meshDelCmd = new G4UIcmdWithAString("/score/mesh/delete",this);
  meshDelCmd->SetGuidance("Delete scoring mesh.");
  meshDelCmd->SetParameterName("MeshName",false);
  //
  //   Shape commands
  mShapeDir = new G4UIdirectory("/score/mesh/shape/");
  mShapeDir->SetGuidance("Shape commands for scoring mesh");
  //
  mSBoxCmd = new G4UIcmdWith3VectorAndUnit("/score/mesh/shape/box",this);
  mSBoxCmd->SetGuidance("Define Box type scoring mesh.");
  mSBoxCmd->SetParameterName("dX","dY","dZ",false,false);
  mSBoxCmd->SetDefaultUnit("mm");
  mSTubsCmd= new G4UIcmdWith3VectorAndUnit("/score/mesh/shape/tubs",this);
  mSTubsCmd->SetGuidance("Define Tubs type scoring mesh.");
  mSTubsCmd->SetParameterName("Rin","Rout","dZ",false,false);
  mSTubsCmd->SetDefaultUnit("mm");
  mSSphereCmd = new G4UIcommand("/score/mesh/shape/sphere",this);
  mSSphereCmd->SetGuidance("Define Sphere type scoring mesh.");
  //mSSphereCmd->SetParameterName("Rin","Rout",false,false);
  //mSSphereCmd->SetDefaultUnit("mm");
  //
  //   Division command
  mBinDir = new G4UIdirectory("/score/mesh/bin/");
  mShapeDir->SetGuidance("Binning of scoring mesh");
  //
  mBinCmd = new G4UIcommand("/score/mesh/bin/numberOfBin",this);
  mBinCmd->SetGuidance("Define segmentation of scoring mesh.");
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
  quantityOpnCmd = new G4UIcmdWithAString("/score/quantity/open",this);
  quantityOpnCmd->SetGuidance("Open quantity of scoring mesh.");
  quantityOpnCmd->SetParameterName("QuantityName",false);
  quantityClsCmd = new G4UIcmdWithoutParameter("/score/quantity/close",this);
  quantityClsCmd->SetGuidance("Close quantity of scoring mesh.");
  quantityDelCmd = new G4UIcmdWithAString("/score/quantity/delete",this);
  quantityDelCmd->SetGuidance("Delete quantity of scoring mesh.");
  quantityDelCmd->SetParameterName("QuantityName",false);
  //
  //    Type Command
  qTypeCmd = new G4UIcommand("/score/quantity/type",this);
  qTypeCmd->SetGuidance("Quantity type of the scorer.");
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
  //    Type Command
  fTypeCmd = new G4UIcommand("/score/filter/type",this);
  fTypeCmd->SetGuidance("Filter type of the scorer.");

}

G4ScoringMessenger::~G4ScoringMessenger()
{
  delete listCmd;
  delete scoreDir;

  delete meshOpnCmd;
  delete meshClsCmd;
  delete meshDelCmd;
  delete mSBoxCmd;
  delete mSTubsCmd;
  delete mSSphereCmd;
  delete mShapeDir;
  delete mBinCmd;
  delete mBinDir;
  delete mTSetCmd;
  delete mTAddCmd;
  delete mTransDir;
  delete mRSetCmd;
  delete mRAddCmd;
  delete mRotDir;
  delete meshDir;

  delete quantityOpnCmd;
  delete quantityClsCmd;
  delete quantityDelCmd;
  delete qTypeCmd;
  delete quantityDir;

  delete filterOpnCmd;
  delete filterClsCmd;
  delete filterDelCmd;
  delete fTypeCmd;
  delete fAttachCmd;
  delete filterDir;

}

void G4ScoringMessenger::SetNewValue(G4UIcommand * command,G4String newVal)
{
  if(command==listCmd) { 
      fSMan->List(); 
  } else if(command==verboseCmd) { 
      fSMan->SetVerboseLevel(verboseCmd->GetNewIntValue(newVal)); 
  } else if(command==meshOpnCmd) {
  } else if(command==meshClsCmd) {
  } else if(command==meshDelCmd) {
  } else if(command==mSBoxCmd) {
  } else if(command==mSTubsCmd) {
  } else if(command==mSSphereCmd) {
  } else if(command==mBinCmd) {
  } else if(command==mTSetCmd) {
  } else if(command==mTAddCmd) {
  } else if(command==mRSetCmd) {
  } else if(command==mRAddCmd) {
  } else if(command==quantityOpnCmd) {
  } else if(command==quantityClsCmd) {
  } else if(command==quantityDelCmd) {
  } else if(command==qTypeCmd) {

  } else if(command==filterOpnCmd) {
  } else if(command==filterClsCmd) {
  } else if(command==filterDelCmd) {
  } else if(command==fTypeCmd) {
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

