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
// GRPrimGenActionMessenger.hh
//   Header file of a messenger that handles primary generator action.
//   Input radiation spectrum file should be in ASCII format and each
//   row should have low-end kinetic energy and differential flux
//   separated by a space.
//
// History
//   September 8th, 2020 : first implementation
//
// ********************************************************************

#include "GRPrimGenActionMessenger.hh"

#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"

GRPrimGenActionMessenger::GRPrimGenActionMessenger(G4int verboseLvl)
: verboseLevel(verboseLvl)
{
  G4UIparameter* param = nullptr;

  srcDir = new G4UIdirectory("/gorad/source/");
  srcDir->SetGuidance("Define primary particle spectrum.");

  defCmd = new G4UIcommand("/gorad/source/define",this);
  defCmd->SetGuidance("Define primary particle spectrum.");
  defCmd->SetGuidance("[usage] /gorad/source/define pName fName srcType radius unit (x0 y0 z0)");
  defCmd->SetGuidance("  pName : Particle name.");
  defCmd->SetGuidance("  fName : File name of the spectrum. Directory may be preceded to the file name.");
  defCmd->SetGuidance("  srcType : defines how the primaries are generated.");
  defCmd->SetGuidance("            Arb : Generate primaries along with the defined spectrum.");
  defCmd->SetGuidance("            LW  : Generate primaries in flat ditribution with track weight representing the spectrum.");
  defCmd->SetGuidance("  radius : Radius of the source sphere. If it is set to -1 (default), radius is set to 98% of the world volume.");
  defCmd->SetGuidance("  unit : Unit of radius value.");
  defCmd->SetGuidance("  x0,y0,z0 : (Optional) location of the center of the sphere (same unit is applied as radius).");
  defCmd->SetGuidance("             By default the sphere is located at the center of the world volume.");
  defCmd->SetGuidance("[note] This command is not mandatory but dwialternatively granular /gps/ commands can be used for defining primary particles.");
  particlePara = new G4UIparameter("pName",'s',false);
  defCmd->SetParameter(particlePara);
  param = new G4UIparameter("fName",'s',false);
  defCmd->SetParameter(param);
  param = new G4UIparameter("srcType",'s',false);
  param->SetParameterCandidates("Arb LW");
  defCmd->SetParameter(param);
  param = new G4UIparameter("radius",'d',true);
  param->SetParameterRange("radius == -1. || radius>0.");
  param->SetDefaultValue(-1.);
  defCmd->SetParameter(param);
  param = new G4UIparameter("unit",'s',true);
  param->SetDefaultUnit("m");
  defCmd->SetParameter(param);
  param = new G4UIparameter("x0",'d',true);
  param->SetDefaultValue(0.);
  defCmd->SetParameter(param);
  param = new G4UIparameter("y0",'d',true);
  param->SetDefaultValue(0.);
  defCmd->SetParameter(param);
  param = new G4UIparameter("z0",'d',true);
  param->SetDefaultValue(0.);
  defCmd->SetParameter(param);
  defCmd->AvailableForStates(G4State_Idle);
  defCmd->SetToBeBroadcasted(false);
}

GRPrimGenActionMessenger::~GRPrimGenActionMessenger()
{
  delete defCmd;
  delete srcDir;
}

#include "G4UImanager.hh"
#include "G4Threading.hh"
#include "G4Tokenizer.hh"
#include "G4UnitsTable.hh"
#include <fstream>
#include "GRDetectorConstruction.hh"

void GRPrimGenActionMessenger::SetNewValue(G4UIcommand * command,G4String newVal)
{
  auto UI = G4UImanager::GetUIpointer();
  if(verboseLevel>0) G4cout << newVal << G4endl;

  G4Tokenizer next(newVal);
  G4String pName = next();
  G4String fName = next();
  std::ifstream specFile;
  specFile.open(fName,std::ios::in);
  if(specFile.fail())
  {
    G4ExceptionDescription ed;
    ed << "ERROR : File <" << fName << "> is not found.";
    command->CommandFailed(ed);
    return;
  }

  G4String srcType = next();
  G4String radius = next();
  G4String unit = next();
  G4double r = StoD(radius);
  if(r<0.)
  {
    r = GRDetectorConstruction::GetWorldSize() * 0.98;
    radius = DtoS(r);
    unit = "mm";
    G4cout<<"################ Radius of the GPS sphere is set to "<<r<<" (mm)"<<G4endl;
  }
  G4String pos = next() + " " + next() + " " + next();

  G4String cmd;
  G4int ec = 0;
  ec = std::max(UI->ApplyCommand("/gps/pos/shape Sphere"),ec);
  cmd = "/gps/pos/centre " + pos + " " + unit;
  ec = std::max(UI->ApplyCommand(cmd),ec);
  cmd = "/gps/pos/radius " + radius + " " + unit;
  ec = std::max(UI->ApplyCommand(cmd),ec);
  ec = std::max(UI->ApplyCommand("/gps/pos/type Surface"),ec);
  ec = std::max(UI->ApplyCommand("/gps/number 1"),ec);
  cmd = "/gps/particle " + pName;
  ec = std::max(UI->ApplyCommand(cmd),ec);
  ec = std::max(UI->ApplyCommand("/gps/ang/type cos"),ec);
  ec = std::max(UI->ApplyCommand("/gps/ang/maxtheta 90.0 deg"),ec);
  ec = std::max(UI->ApplyCommand("/gps/ang/mintheta 0.0 deg"),ec);
  cmd = "/gps/ene/type " + srcType;
  ec = std::max(UI->ApplyCommand(cmd),ec);
  ec = std::max(UI->ApplyCommand("/gps/hist/type arb"),ec);

  enum { bufsize = 128 };
  static char* line = new char[bufsize];
  while(specFile.good())
  {
    specFile.getline(line,bufsize);
    if(line[(size_t)0]=='#') continue;
    G4String valStr = line;
    G4StrUtil::strip(valStr);
    G4StrUtil::rstrip(valStr, 0x0d);
    if(valStr.size()==0) continue;
    cmd = "/gps/hist/point " + valStr;
    ec = std::max(UI->ApplyCommand(cmd),ec);
    if(specFile.eof()) break;
  }
  ec = std::max(UI->ApplyCommand("/gps/hist/inter Lin"),ec);

  if(ec>0)
  {
    G4ExceptionDescription ed;
    ed << "ERROR : Internal error while processing /gorad/source/define command.";
    command->CommandFailed(ec,ed);
  }
}

G4String GRPrimGenActionMessenger::GetCurrentValue(G4UIcommand* /*command*/)
{
  G4String val;
  return val;
}

#include "G4ParticleTable.hh"
void GRPrimGenActionMessenger::UpdateParticleList()
{
  auto particleTable = G4ParticleTable::GetParticleTable();
  G4String candList;
  for(G4int i=0;i<particleTable->entries();i++)
  { candList += particleTable->GetParticleName(i) + " "; }
  candList += "ion";
  particlePara->SetParameterCandidates(candList);
}

