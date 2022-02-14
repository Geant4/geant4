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
//  Gorad (Geant4 Open-source Radiation Analysis and Design)
//
//  Author : Makoto Asai (SLAC National Accelerator Laboratory)
//
//  Development of Gorad is funded by NASA Johnson Space Center (JSC)
//  under the contract NNJ15HK11B.
//
// ********************************************************************
//
// GRGeomBiasMessenger.cc
//   A messenger class that handles the UI commands for geometry
//   imprtance biasing.
//
// History
//   September 8th, 2020 : first implementation
//
// ********************************************************************

#include "GRGeomBiasMessenger.hh"

#include "GRDetectorConstruction.hh"
#include "GRPhysicsList.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIparameter.hh"
#include "G4Tokenizer.hh"
#include "G4UnitsTable.hh"

GRGeomBiasMessenger::GRGeomBiasMessenger(
  GRDetectorConstruction* det,GRPhysicsList* phys,G4int verboseLvl)
: detector(det),physics(phys),verboseLevel(verboseLvl)
{
  G4UIparameter* param = nullptr;

  biasDir = new G4UIdirectory("/gorad/bias/");
  biasDir->SetGuidance("Gorad biasing commands");

  geoBiasCmd = new G4UIcommand("/gorad/bias/geomImportance",this);
  geoBiasCmd->SetGuidance("Geometry importance biasing");
  geoBiasCmd->SetGuidance("This command defines the number of layers and radius of the outermost sphere for geometry importance biasing.");
  geoBiasCmd->SetGuidance("If radius is set to -1 (default), actual radius is set to 80% of the world volume");
  geoBiasCmd->SetGuidance("Note: There must be at least two layers.");
  param = new G4UIparameter("nLayer",'i',false);
  param->SetParameterRange("nLayer > 1");
  geoBiasCmd->SetParameter(param);
  param = new G4UIparameter("radius",'d',true);
  param->SetParameterRange("radius == -1.0 || radius > 0.");
  param->SetDefaultValue(-1.0);
  geoBiasCmd->SetParameter(param);
  param = new G4UIparameter("unit",'s',true);
  param->SetDefaultUnit("mm");
  geoBiasCmd->SetParameter(param);
  geoBiasCmd->SetToBeBroadcasted(false);
  geoBiasCmd->AvailableForStates(G4State_PreInit);

  geoBiasLocCmd = new G4UIcmdWith3VectorAndUnit("/gorad/bias/geomImpLocate",this);
  geoBiasLocCmd->SetGuidance("Position of the center of the outermost sphere.");
  geoBiasLocCmd->SetGuidance("By default, the sphere is located at the origin of the world volume.");
  geoBiasLocCmd->SetGuidance("This command has to follow /gorad/bias/geomImportance command.");
  geoBiasLocCmd->SetParameterName("x0","y0","z0",false);
  geoBiasLocCmd->SetDefaultUnit("mm");
  geoBiasLocCmd->SetToBeBroadcasted(false);
  geoBiasLocCmd->AvailableForStates(G4State_PreInit);
  
  geoBiasInRadCmd = new G4UIcmdWithADoubleAndUnit("/gorad/bias/geomImpInnerRadius",this);
  geoBiasInRadCmd->SetGuidance("Radius of the innermost sphere.");
  geoBiasInRadCmd->SetGuidance("By default it is defined as 1/n of radius of the outermost sphere.");
  geoBiasInRadCmd->SetParameterName("rT",false);
  geoBiasInRadCmd->SetDefaultUnit("mm");
  geoBiasInRadCmd->SetRange("rT > 0.");
  geoBiasInRadCmd->SetToBeBroadcasted(false);
  geoBiasInRadCmd->AvailableForStates(G4State_PreInit);
  
  geoBiasLocTgtCmd = new G4UIcmdWith3VectorAndUnit("/gorad/bias/geomImpLocTgt",this);
  geoBiasLocTgtCmd->SetGuidance("Position of the center of the innermost sphere.");
  geoBiasLocTgtCmd->SetGuidance("By default, it is located at the center of the outermost sphere.");
  geoBiasLocTgtCmd->SetGuidance("This command has to follow /gorad/bias/geomImportance command.");
  geoBiasLocTgtCmd->SetGuidance("Note: distance between (x0,y,0,z0) and (xT,yT,zT) must be smaller than r0*(nLayer-1)/nLayer.");
  geoBiasLocTgtCmd->SetGuidance("          (smaller than r0-rT if radius of innermost sphere is set)");
  geoBiasLocTgtCmd->SetParameterName("xT","yT","zT",false);
  geoBiasLocTgtCmd->SetDefaultUnit("mm");
  geoBiasLocTgtCmd->SetToBeBroadcasted(false);
  geoBiasLocTgtCmd->AvailableForStates(G4State_PreInit);

  geoBiasFucCmd = new G4UIcmdWithAnInteger("/gorad/bias/geomImpFactor",this);
  geoBiasFucCmd->SetGuidance("Alternate the geometry importance biasing factor.");
  geoBiasFucCmd->SetGuidance("By default the factor is set to 2. We do not recommend the factor to be much larger than 2.");
  geoBiasFucCmd->SetGuidance("This command has to follow /gorad/bias/geomImportance command.");
  geoBiasFucCmd->SetParameterName("factor",false);
  geoBiasFucCmd->SetDefaultValue(2);
  geoBiasFucCmd->SetRange("factor > 0");
  geoBiasFucCmd->SetToBeBroadcasted(false);
  geoBiasFucCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  geoBiasProbCmd = new G4UIcmdWithADouble("/gorad/bias/geomImpProbability",this);
  geoBiasProbCmd->SetGuidance("Reduce the probability of geometry importance biasing to avoid over biasing.");
  geoBiasProbCmd->SetGuidance("By default the probability is set to 1.0 (i.e. 100%).");
  geoBiasProbCmd->SetGuidance("This command has to follow /gorad/bias/geomImportance command.");
  geoBiasProbCmd->SetParameterName("prob",true);
  geoBiasProbCmd->SetDefaultValue(1.);
  geoBiasProbCmd->SetRange("prob > 0. && prob <= 1.0");
  geoBiasProbCmd->SetToBeBroadcasted(false);
  geoBiasProbCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  if(verboseLevel>0) 
  { G4cout << "UI commands /gorad/bias/ instantiated." << G4endl; }
}

GRGeomBiasMessenger::~GRGeomBiasMessenger()
{
  delete geoBiasProbCmd;
  delete geoBiasFucCmd;
  delete geoBiasLocTgtCmd;
  delete geoBiasInRadCmd;
  delete geoBiasLocCmd;
  delete geoBiasCmd;
  delete biasDir;
}

void GRGeomBiasMessenger::SetNewValue(G4UIcommand * command,G4String newVal)
{
  if(command==geoBiasCmd)
  {
    G4Tokenizer next(newVal);
    G4int nL = StoI(next());
    G4String r = next() + " ";
    r += next();
    detector->GeomImp(nL,geoBiasCmd->ConvertToDimensionedDouble(r));
    physics->ApplyGeomImpBias();
  }
  else
  { // following commands have to come after /gorad/bias/geomImportance command.
    G4bool applyGeomImpBias = detector->ApplyGeomImpBias();
    if(!applyGeomImpBias)
    {
      G4ExceptionDescription ed;
      ed << "This command has to follow /gorad/bias/geomImportance command. Command failed.";
      command->CommandFailed(ed);
    }
    else if(command==geoBiasLocCmd)
    {
      detector->GeomImpLocate(geoBiasLocCmd->GetNew3VectorValue(newVal));
    }
    else if(command==geoBiasInRadCmd)
    {
      G4double rt = detector->GeomImpInnerRadius(geoBiasInRadCmd->GetNewDoubleValue(newVal));
      if(rt<0.)
      {
        G4ExceptionDescription ed;
        ed << "Specified radius is too large. It has to be smaller than the outermost sphere "
           << -rt << " (mm)\n" << "command failed.";
        command->CommandFailed(ed);
      }
    }
    else if(command==geoBiasLocTgtCmd)
    {
      G4double dr = detector->GeomImpLocateTgt(geoBiasLocTgtCmd->GetNew3VectorValue(newVal));
      if(dr>0.)
      {
        G4ExceptionDescription ed;
        ed << "Distance between (x0,y,0,z0) and (xT,yT,zT) must be smaller than radius*(nLayer-1)/nLayer, "
           << dr << " (mm)\n" << "command failed."; 
        command->CommandFailed(ed);
      }
    }
    else if(command==geoBiasFucCmd)
    {
      detector->GeomImpFactor(StoI(newVal));
    }
    else if(command==geoBiasProbCmd)
    {
      detector->GeomImpProb(StoD(newVal));
    }
  }
}

G4String GRGeomBiasMessenger::GetCurrentValue(G4UIcommand* /*command*/)
{
  G4String val;
  return val;
}
