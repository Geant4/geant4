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
// GRDetectorConstructionMessenger.cc
//   A messenger class that handles geometry configuration.
//
// History
//   September 8th, 2020 : first implementation
//
// ********************************************************************

#include "GRDetectorConstructionMessenger.hh"

#include "GRDetectorConstruction.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"

GRDetectorConstructionMessenger::GRDetectorConstructionMessenger(GRDetectorConstruction* dc)
: pDC(dc)
{
  G4UIparameter* para = nullptr;

  geomDir = new G4UIdirectory("/gorad/geometry/");
  geomDir->SetGuidance("GORAD geometry selection");

  selectCmd = new G4UIcmdWithAString("/gorad/geometry/selectGDML",this);
  selectCmd->SetGuidance("Select GDML file");
  selectCmd->SetParameterName("gdml",false);
  selectCmd->AvailableForStates(G4State_PreInit);
  selectCmd->SetToBeBroadcasted(false);

  listSolidCmd = new G4UIcmdWithAnInteger("/gorad/geometry/listSolids",this);
  listSolidCmd->SetGuidance("List all the registered solids");
  listSolidCmd->SetParameterName("level",true);
  listSolidCmd->SetDefaultValue(0);
  listSolidCmd->SetRange("level>=0 && level<=2");
  listSolidCmd->AvailableForStates(G4State_Idle);
  listSolidCmd->SetToBeBroadcasted(false);

  listLogVolCmd = new G4UIcmdWithAnInteger("/gorad/geometry/listLogicalVolumes",this);
  listLogVolCmd->SetGuidance("List all the registered logical volumes");
  listLogVolCmd->SetParameterName("level",true);
  listLogVolCmd->SetDefaultValue(0);
  listLogVolCmd->SetRange("level>=0 && level<=2");
  listLogVolCmd->AvailableForStates(G4State_Idle);
  listLogVolCmd->SetToBeBroadcasted(false);

  listPhysVolCmd = new G4UIcmdWithAnInteger("/gorad/geometry/listPhysicalVolumes",this);
  listPhysVolCmd->SetGuidance("List all the registered physical volumes");
  listPhysVolCmd->SetParameterName("level",true);
  listPhysVolCmd->SetDefaultValue(0);
  listPhysVolCmd->SetRange("level>=0 && level<=2");
  listPhysVolCmd->AvailableForStates(G4State_Idle);
  listPhysVolCmd->SetToBeBroadcasted(false);

  listRegionCmd = new G4UIcmdWithAnInteger("/gorad/geometry/listRegions",this);
  listRegionCmd->SetGuidance("List all the registered regions");
  listRegionCmd->SetParameterName("level",true);
  listRegionCmd->SetDefaultValue(0);
  listRegionCmd->SetRange("level>=0 && level<=2");
  listRegionCmd->AvailableForStates(G4State_Idle);
  listRegionCmd->SetToBeBroadcasted(false);

  createRegionCmd = new G4UIcommand("/gorad/geometry/createRegion",this);
  createRegionCmd->SetGuidance("Create a region and set the root logical volume to it.");
  createRegionCmd->SetGuidance("Region propagates to the daughter volumes. So, only the root logical volume (i.e. top of the hierarchy) should be defined.");
  createRegionCmd->SetGuidance("If two isolated root logical volumes should share the same region, the same region name can be used.");
  createRegionCmd->SetGuidance("Region must not be set to the world volume.");
  para = new G4UIparameter("regionName",'s',false);
  para->SetGuidance("Name of the region to be created");
  createRegionCmd->SetParameter(para);
  para = new G4UIparameter("logVolName",'s',false);
  para->SetGuidance("Name of the root logical volume");
  createRegionCmd->SetParameter(para);
  createRegionCmd->AvailableForStates(G4State_Idle);
  createRegionCmd->SetToBeBroadcasted(false);

  //// This command is fragile for large-scale geometry - temporally disabled
  ////checkOverlapCmd = new G4UIcommand("/gorad/geometry/checkOverlap",this);
  ////checkOverlapCmd->SetGuidance("Check volume overlap with existing volumes");
  ////checkOverlapCmd->SetGuidance(" i.e. with mother volume for protrusion and with other siblings for overlap.");
  ////checkOverlapCmd->SetGuidance(" - This command is valid only for placement and parameterized volumes. If this command is");
  ////checkOverlapCmd->SetGuidance("   used for other physical volume type, e.g. replica, command will be simply ignored.");
  ////checkOverlapCmd->SetGuidance("   If \"**ALL**\" is used as the volume name, all physical volumes are examined (SLOW!!).");
  ////checkOverlapCmd->SetGuidance(" - nSpots specifies number of spots on the surface of the volume to be examined.");
  ////checkOverlapCmd->SetGuidance("   The more spots used, the more chances to detect overlaps, but the more time it takes.");
  ////checkOverlapCmd->SetGuidance(" - maxErr specifies maximum number of errors to be generated (default 1) before quiting.");
  ////para = new G4UIparameter("physVol",'s',true);
  ////para->SetDefaultValue("**ALL**");
  ////checkOverlapCmd->SetParameter(para);
  ////para = new G4UIparameter("nSpots",'i',true);
  ////para->SetDefaultValue(1000);
  ////para->SetGuidance("Number of trial spots on the volume surface");
  ////checkOverlapCmd->SetParameter(para);
  ////para = new G4UIparameter("maxErr",'i',true);
  ////para->SetDefaultValue(1);
  ////para->SetParameterRange("maxErr > 0");
  ////para->SetGuidance("Maxinum number of report to be generated");
  ////checkOverlapCmd->SetParameter(para);
  ////para = new G4UIparameter("tolerance",'d',true);
  ////para->SetDefaultValue(0.);
  ////para->SetParameterRange("tolerance >= 0.");
  ////para->SetGuidance("Tolerance (default 0.)");
  ////checkOverlapCmd->SetParameter(para);
  ////para = new G4UIparameter("unit",'s',true);
  ////para->SetDefaultUnit("mm");
  ////checkOverlapCmd->SetParameter(para);
  ////checkOverlapCmd->AvailableForStates(G4State_Idle);
  ////checkOverlapCmd->SetToBeBroadcasted(false);
  
  materialDir = new G4UIdirectory("/gorad/material/");
  materialDir->SetGuidance("GORAD material commands");

  listMatCmd = new G4UIcmdWithAString("/gorad/material/list",this);
  listMatCmd->SetGuidance("List material property");
  listMatCmd->SetGuidance(" If material name is not specified, this command list all registered materials");
  listMatCmd->SetParameterName("matName",true);
  listMatCmd->SetDefaultValue("**ALL**");
  listMatCmd->AvailableForStates(G4State_Idle);
  listMatCmd->SetToBeBroadcasted(false);

  dumpMatCmd = new G4UIcmdWithoutParameter("/gorad/material/dumpNistMaterials",this);
  dumpMatCmd->SetGuidance("List all pre-defined material names in G4NistManager.");
  dumpMatCmd->SetGuidance(" Note : a material has to be instantiated with /gorad/material/create before setting it to a logical volume");
  dumpMatCmd->AvailableForStates(G4State_Idle);
  dumpMatCmd->SetToBeBroadcasted(false);

  createMatCmd = new G4UIcmdWithAString("/gorad/material/create",this);
  createMatCmd->SetGuidance("Instantiate a material defined in G4NistManager");
  createMatCmd->SetGuidance(" If the material has already existed, this command does nothing.");
  createMatCmd->SetParameterName("matName",false);
  createMatCmd->AvailableForStates(G4State_Idle);
  createMatCmd->SetToBeBroadcasted(false);

  getMatCmd = new G4UIcmdWithAString("/gorad/material/show",this); 
  getMatCmd->SetGuidance("Show the current material of the specified logical volume");
  getMatCmd->SetParameterName("logVol",false);
  getMatCmd->AvailableForStates(G4State_Idle);
  getMatCmd->SetToBeBroadcasted(false);

  setMatCmd = new G4UIcommand("/gorad/material/set",this);
  setMatCmd->SetGuidance("Set the material to the logical volume. The material has to be instantiated in advance.");
  setMatCmd->SetGuidance("  [usage] /gorad/material/set logicalVolumeName materialName");
  para = new G4UIparameter("logVol",'s',false);
  setMatCmd->SetParameter(para);
  para = new G4UIparameter("matName",'s',false);
  setMatCmd->SetParameter(para);
  setMatCmd->AvailableForStates(G4State_Idle);
  setMatCmd->SetToBeBroadcasted(false);
}

GRDetectorConstructionMessenger::~GRDetectorConstructionMessenger()
{
  delete selectCmd;
  delete listSolidCmd;
  delete listLogVolCmd;
  delete listPhysVolCmd;
  delete listRegionCmd;
  delete createRegionCmd;
  ////delete checkOverlapCmd;
  delete geomDir;

  delete listMatCmd;
  delete dumpMatCmd;
  delete createMatCmd;
  delete getMatCmd;
  delete setMatCmd;
  delete materialDir;
}

#include "G4Tokenizer.hh"

void GRDetectorConstructionMessenger::SetNewValue(G4UIcommand* cmd, G4String val)
{
  if(cmd==selectCmd)
  {
    auto valid = pDC->SetGDMLFile(val);
    if(!valid)
    {
      G4ExceptionDescription ed;
      ed << "<" << val << "> is not a valid GDML file.";
      cmd->CommandFailed(ed);
    }
  }
  else if(cmd==listSolidCmd)
  { pDC->ListSolids(listSolidCmd->GetNewIntValue(val)); }
  else if(cmd==listLogVolCmd)
  { pDC->ListLogVols(listLogVolCmd->GetNewIntValue(val)); }
  else if(cmd==listPhysVolCmd)
  { pDC->ListPhysVols(listPhysVolCmd->GetNewIntValue(val)); }
  else if(cmd==listRegionCmd)
  { pDC->ListRegions(listRegionCmd->GetNewIntValue(val)); }
  else if(cmd==createRegionCmd)
  {
    G4Tokenizer next(val);
    G4String regionName = next();
    G4String logVolName = next();
    auto valid = pDC->CreateRegion(regionName,logVolName);
    if(!valid)
    {
      G4ExceptionDescription ed;
      ed << "Logical volume <" << logVolName << "> is not defined. Command ignored.";
      cmd->CommandFailed(ed);
    }
  }
  ////else if(cmd==checkOverlapCmd)
  ////{
    ////G4Tokenizer next(val);
    ////G4String physVolName = next();
    ////G4int nSpots = StoI(next());
    ////G4int maxErr = StoI(next());
    ////G4String tolStr = next();
    ////G4double tol = StoD(tolStr);
    ////if(tol>0.)
    ////{
      ////tolStr += " ";
      ////tolStr += next();
      ////tol = checkOverlapCmd->ConvertToDimensionedDouble(tolStr);
    ////}
    ////auto valid = pDC->CheckOverlap(physVolName,nSpots,maxErr,tol);
    ////if(!valid)
    ////{
      ////G4ExceptionDescription ed;
      ////ed << "Physical volume <" << physVolName << "> is not defined. Command ignored.";
      ////cmd->CommandFailed(ed);
    ////}
  ////}

  else if(cmd==listMatCmd)
  { 
    if(val=="**ALL**") 
    { pDC->ListAllMaterial(); }
    else
    {
      auto valid = pDC->ListMaterial(val);
      if(!valid)
      {
        G4ExceptionDescription ed;
        ed << "<" << val << "> is not defined. If necessary, create it with /gorad/material/create command.";
        cmd->CommandFailed(ed);
      }
    }
  }
  else if(cmd==dumpMatCmd)
  { pDC->DumpNistMaterials(); }
  else if(cmd==createMatCmd)
  { 
    auto valid = pDC->CreateMaterial(val); 
    if(!valid)
    {
      G4ExceptionDescription ed;
      ed << "The material name <" << val << "> is not defined in G4NistManager.";
      cmd->CommandFailed(ed);
    }
  }
  else if(cmd==getMatCmd)
  {
    auto valid = pDC->GetMaterial(val);
    if(!valid)
    {
      G4ExceptionDescription ed;
      ed << "<" << val << "> is not a name of registered logical volume.\n"
         << "Check existing logical volumes with /gorad/geometry/listLogicalVolumes command.";
      cmd->CommandFailed(ed);
    }
  }
  else if(cmd==setMatCmd)
  {
    G4Tokenizer next(val);
    G4String logVolName = next();
    G4String matName = next();
    auto valid = pDC->SetMaterial(logVolName,matName);
    if(valid!=0)
    {
      G4ExceptionDescription ed;
      if(valid==1 || valid==3)
      {
        ed << "<" << logVolName << "> is not a name of registered logical volume.\n"
           << "Check existing logical volumes with /gorad/geometry/listLogicalVolumes command.\n";
      }
      if(valid==2 || valid==3)
      {
        ed << "<" << matName << "> is not defined. If necessary, create it with /gorad/material/create command.";
      }
      cmd->CommandFailed(ed);
    }
  }

}

G4String GRDetectorConstructionMessenger::GetCurrentValue(G4UIcommand* cmd)
{
  G4String val("");
  if(cmd==selectCmd)
  { val = pDC->GetGDMLFile(); }
  return val;
}


