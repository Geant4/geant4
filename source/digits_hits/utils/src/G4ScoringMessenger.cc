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
// ---------------------------------------------------------------------

#include "G4ScoringMessenger.hh"
#include "G4ScoringManager.hh"
#include "G4VScoringMesh.hh"
#include "G4ScoringBox.hh"
#include "G4ScoringCylinder.hh"
#include "G4ScoringRealWorld.hh"
#include "G4ScoringProbe.hh"

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
#include "G4VScoreColorMap.hh"

G4ScoringMessenger::G4ScoringMessenger(G4ScoringManager* SManager)
  : fSMan(SManager)
{
  G4UIparameter* param = nullptr;

  scoreDir = new G4UIdirectory("/score/");
  scoreDir->SetGuidance("Interactive scoring commands.");

  listCmd = new G4UIcmdWithoutParameter("/score/list", this);
  listCmd->SetGuidance("List scoring worlds.");

  dumpCmd = new G4UIcmdWithoutParameter("/score/dump", this);
  dumpCmd->SetGuidance("Dump results of scorers.");

  verboseCmd = new G4UIcmdWithAnInteger("/score/verbose", this);
  verboseCmd->SetGuidance("Verbosity.");
  verboseCmd->SetGuidance("  0) errors or warnings,");
  verboseCmd->SetGuidance("  1) information with 0)");

  meshCreateDir = new G4UIdirectory("/score/create/");
  meshCreateDir->SetGuidance("  Mesh creation commands.");
  //
  // Mesh commands
  meshBoxCreateCmd = new G4UIcmdWithAString("/score/create/boxMesh", this);
  meshBoxCreateCmd->SetGuidance("Create scoring box mesh.");
  meshBoxCreateCmd->SetParameterName("MeshName", false);
  //
  meshCylinderCreateCmd =
    new G4UIcmdWithAString("/score/create/cylinderMesh", this);
  meshCylinderCreateCmd->SetGuidance("Create scoring mesh.");
  meshCylinderCreateCmd->SetParameterName("MeshName", false);

  meshRWLogVolCreateCmd =
    new G4UIcommand("/score/create/realWorldLogVol", this);
  meshRWLogVolCreateCmd->SetGuidance(
    "Define scorers to a logical volume defined in the real world.");
  meshRWLogVolCreateCmd->SetGuidance(
    "  - Name of the specified logical volume is used as the mesh name.");
  meshRWLogVolCreateCmd->SetGuidance(
    "  - /score/mesh commands do not affect for this mesh.");
  meshRWLogVolCreateCmd->SetGuidance(
    "  - If copyNumberLevel is set, the copy number of that-level higher");
  meshRWLogVolCreateCmd->SetGuidance(
    "    in the geometrical hierarchy is used as the index.");
  param = new G4UIparameter("logVol", 's', false);
  meshRWLogVolCreateCmd->SetParameter(param);
  param = new G4UIparameter("copyNumberLevel", 'i', true);
  param->SetParameterRange("copyNumberLevel>=0");
  param->SetDefaultValue(0);
  meshRWLogVolCreateCmd->SetParameter(param);
  //
  probeCreateCmd = new G4UIcommand("/score/create/probe", this);
  probeCreateCmd->SetGuidance("Define scoring probe.");
  probeCreateCmd->SetGuidance(
    "  halfSize defines the half-width of the probing cube.");
  param = new G4UIparameter("pname", 's', false);
  probeCreateCmd->SetParameter(param);
  param = new G4UIparameter("halfSize", 'd', false);
  probeCreateCmd->SetParameter(param);
  param = new G4UIparameter("unit", 's', true);
  param->SetDefaultUnit("mm");
  probeCreateCmd->SetParameter(param);
  param = new G4UIparameter("checkOverlap", 'b', true);
  param->SetDefaultValue(0);
  probeCreateCmd->SetParameter(param);
  //
  meshOpnCmd = new G4UIcmdWithAString("/score/open", this);
  meshOpnCmd->SetGuidance("Open scoring mesh.");
  meshOpnCmd->SetParameterName("MeshName", false);
  //
  meshClsCmd = new G4UIcmdWithoutParameter("/score/close", this);
  meshClsCmd->SetGuidance("Close scoring mesh.");
 
  meshDir = new G4UIdirectory("/score/mesh/");
  meshDir->SetGuidance("    Mesh processing commands.");
  //
  mBoxSizeCmd = new G4UIcmdWith3VectorAndUnit("/score/mesh/boxSize", this);
  mBoxSizeCmd->SetGuidance("Define size of the scoring mesh.");
  mBoxSizeCmd->SetGuidance("Dx  Dy  Dz  unit");
  mBoxSizeCmd->SetParameterName("Di", "Dj", "Dk", false, false);
  mBoxSizeCmd->SetRange("Di>0. && Dj>0. && Dk>0.");
  mBoxSizeCmd->SetDefaultUnit("mm");
  //
  mCylinderSizeCmd = new G4UIcommand("/score/mesh/cylinderSize", this);
  mCylinderSizeCmd->SetGuidance("Define size of the scoring mesh.");
  mCylinderSizeCmd->SetGuidance("R   Dz  unit");
  param = new G4UIparameter("R", 'd', false);
  param->SetParameterRange("R>0");
  mCylinderSizeCmd->SetParameter(param);
  param = new G4UIparameter("Dz", 'd', false);
  param->SetParameterRange("Dz>0");
  mCylinderSizeCmd->SetParameter(param);
  param = new G4UIparameter("unit", 's', true);
  param->SetDefaultUnit("mm");
  mCylinderSizeCmd->SetParameter(param);
  //
  mCylinderRMinCmd =
    new G4UIcmdWithADoubleAndUnit("/score/mesh/cylinderRMin", this);
  mCylinderRMinCmd->SetGuidance("Define the inner radius of the tube mesh.");
  mCylinderRMinCmd->SetGuidance("This command is not needed for cylinder mesh");
  mCylinderRMinCmd->SetParameterName("RMin", false, false);
  mCylinderRMinCmd->SetRange("RMin>=0.");
  mCylinderRMinCmd->SetDefaultUnit("mm");
  //
  mCylinderAngleCmd = new G4UIcommand("/score/mesh/cylinderAngles", this);
  mCylinderAngleCmd->SetGuidance(
    "Define starting angle and span for tube segment mesh.");
  mCylinderAngleCmd->SetGuidance(
    "This command is not needed for cylinder mesh");
  param = new G4UIparameter("startPhi", 'd', false);
  mCylinderAngleCmd->SetParameter(param);
  param = new G4UIparameter("deltaPhi", 'd', false);
  param->SetParameterRange("deltaPhi>0.");
  mCylinderAngleCmd->SetParameter(param);
  param = new G4UIparameter("unit", 's', true);
  param->SetDefaultUnit("deg");
  mCylinderAngleCmd->SetParameter(param);
  //
  //   Division command
  mBinCmd = new G4UIcommand("/score/mesh/nBin", this);
  mBinCmd->SetGuidance("Define segments of the scoring mesh.");
  mBinCmd->SetGuidance("[usage] /score/mesh/nBin");
  mBinCmd->SetGuidance(" In case of boxMesh, parameters are given in");
  mBinCmd->SetGuidance("   Ni  :(int) Number of bins i (in x-axis) ");
  mBinCmd->SetGuidance("   Nj  :(int) Number of bins j (in y-axis) ");
  mBinCmd->SetGuidance("   Nk  :(int) Number of bins k (in z-axis) ");
  mBinCmd->SetGuidance(" In case of cylinderMesh, parameters are given in");
  mBinCmd->SetGuidance("   Nr  :(int) Number of bins in radial axis ");
  mBinCmd->SetGuidance("   Nz  :(int) Number of bins in z axis ");
  mBinCmd->SetGuidance("   Nphi:(int) Number of bins in phi axis ");

  param = new G4UIparameter("Ni", 'i', false);
  param->SetDefaultValue("1");
  param->SetParameterRange("Ni>0");
  mBinCmd->SetParameter(param);
  param = new G4UIparameter("Nj", 'i', false);
  param->SetDefaultValue("1");
  param->SetParameterRange("Nj>0");
  mBinCmd->SetParameter(param);
  param = new G4UIparameter("Nk", 'i', false);
  param->SetDefaultValue("1");
  mBinCmd->SetParameter(param);
  param->SetParameterRange("Nk>0");

  //   Placement command
  mTransDir = new G4UIdirectory("/score/mesh/translate/");
  mTransDir->SetGuidance("Mesh translation commands.");
  //
  mTResetCmd = new G4UIcmdWithoutParameter("/score/mesh/translate/reset", this);
  mTResetCmd->SetGuidance("Reset translated position of the scoring mesh.");
  //
  mTXyzCmd = new G4UIcmdWith3VectorAndUnit("/score/mesh/translate/xyz", this);
  mTXyzCmd->SetGuidance("Translate the scoring mesh.");
  mTXyzCmd->SetParameterName("X", "Y", "Z", false, false);
  mTXyzCmd->SetDefaultUnit("mm");
  //
  mRotDir = new G4UIdirectory("/score/mesh/rotate/");
  mRotDir->SetGuidance("Mesh rotation commands.");

  mRotXCmd = new G4UIcmdWithADoubleAndUnit("/score/mesh/rotate/rotateX", this);
  mRotXCmd->SetGuidance("Rotate the scoring mesh in X axis.");
  mRotXCmd->SetParameterName("Rx", false);
  mRotXCmd->SetDefaultUnit("deg");
  //
  mRotYCmd = new G4UIcmdWithADoubleAndUnit("/score/mesh/rotate/rotateY", this);
  mRotYCmd->SetGuidance("Rotate the scoring mesh in Y axis.");
  mRotYCmd->SetParameterName("Ry", false);
  mRotYCmd->SetDefaultUnit("deg");
  //
  mRotZCmd = new G4UIcmdWithADoubleAndUnit("/score/mesh/rotate/rotateZ", this);
  mRotZCmd->SetGuidance("Rotate the scoring mesh in Z axis.");
  mRotZCmd->SetParameterName("Rz", false);
  mRotZCmd->SetDefaultUnit("deg");
  //
  probeDir = new G4UIdirectory("/score/probe/");
  probeDir->SetGuidance("Probe commands");

  probeMatCmd = new G4UIcmdWithAString("/score/probe/material", this);
  probeMatCmd->SetGuidance("Specify a material to the probe cube.");
  probeMatCmd->SetGuidance("Material name has to be taken from G4NistManager.");
  probeMatCmd->SetGuidance("Once this command is used, the specified material "
                           "overlays the material in the mass geometry");
  probeMatCmd->SetGuidance("with \"Layered Mass Geometry\" mechanism so that "
                           "physics quantities such as energy deposition");
  probeMatCmd->SetGuidance("or dose will be calculated with this material.");
  probeMatCmd->SetGuidance("To switch-off this overlaying, use \"none\".");
  probeMatCmd->SetParameterName("matName", true);
  probeMatCmd->SetDefaultValue("none");
  probeMatCmd->SetToBeBroadcasted(false);

  probeLocateCmd = new G4UIcmdWith3VectorAndUnit("/score/probe/locate", this);
  probeLocateCmd->SetGuidance(
    "Locate a probe in the global coordinate system.");
  probeLocateCmd->SetParameterName("x", "y", "z", false);
  probeLocateCmd->SetDefaultUnit("mm");
  probeLocateCmd->SetToBeBroadcasted(false);

  // Draw Scoring result
  drawCmd = new G4UIcommand("/score/drawProjection", this);
  drawCmd->SetGuidance("Draw projection(s) of scored quantities.");
  drawCmd->SetGuidance(
    "Parameter <proj> specified which projection(s) to be drawn.");
  drawCmd->SetGuidance(
    "  100 : xy-plane, 010 : yz-plane,    001 : zx-plane -- default 111");
  drawCmd->SetGuidance(
    "  100 : N/A,      010 : z_phi-plane, 001 : r_phi-plane -- default 111");
  param = new G4UIparameter("meshName", 's', false);
  drawCmd->SetParameter(param);
  param = new G4UIparameter("psName", 's', false);
  drawCmd->SetParameter(param);
  param = new G4UIparameter("colorMapName", 's', true);
  param->SetDefaultValue("defaultLinearColorMap");
  drawCmd->SetParameter(param);
  param = new G4UIparameter("proj", 'i', true);
  param->SetDefaultValue(111);
  drawCmd->SetParameter(param);
  drawCmd->SetToBeBroadcasted(false);

  // Draw column
  drawColumnCmd = new G4UIcommand("/score/drawColumn", this);
  drawColumnCmd->SetGuidance("Draw a cell column.");
  drawColumnCmd->SetGuidance(" plane = 0 : x-y, 1: y-z, 2: z-x  for box mesh");
  drawColumnCmd->SetGuidance(
    "         0 : z-phi, 1: r-phi, 2: r-z  for cylinder mesh");
  param = new G4UIparameter("meshName", 's', false);
  drawColumnCmd->SetParameter(param);
  param = new G4UIparameter("psName", 's', false);
  drawColumnCmd->SetParameter(param);
  param = new G4UIparameter("plane", 'i', false);
  param->SetParameterRange("plane>=0 && plane<=2");
  drawColumnCmd->SetParameter(param);
  param = new G4UIparameter("column", 'i', false);
  drawColumnCmd->SetParameter(param);
  param = new G4UIparameter("colorMapName", 's', true);
  param->SetDefaultValue("defaultLinearColorMap");
  drawColumnCmd->SetParameter(param);
  drawColumnCmd->SetToBeBroadcasted(false);

  colorMapDir = new G4UIdirectory("/score/colorMap/");
  colorMapDir->SetGuidance("Color map commands.");

  listColorMapCmd =
    new G4UIcmdWithoutParameter("/score/colorMap/listScoreColorMaps", this);
  listColorMapCmd->SetGuidance("List registered score color maps.");
  listColorMapCmd->SetToBeBroadcasted(false);

  floatMinMaxCmd = new G4UIcmdWithAString("/score/colorMap/floatMinMax", this);
  floatMinMaxCmd->SetGuidance(
    "Min/Max of the color map is calculated according to the actual scores.");
  floatMinMaxCmd->SetParameterName("colorMapName", true, false);
  floatMinMaxCmd->SetDefaultValue("defaultLinearColorMap");
  floatMinMaxCmd->SetToBeBroadcasted(false);

  colorMapMinMaxCmd = new G4UIcommand("/score/colorMap/setMinMax", this);
  colorMapMinMaxCmd->SetGuidance("Define min/max value of the color map.");
  param = new G4UIparameter("colorMapMame", 's', true);
  param->SetDefaultValue("defaultLinearColorMap");
  colorMapMinMaxCmd->SetParameter(param);
  param = new G4UIparameter("minValue", 'd', false);
  colorMapMinMaxCmd->SetParameter(param);
  param = new G4UIparameter("maxValue", 'd', false);
  colorMapMinMaxCmd->SetParameter(param);
  colorMapMinMaxCmd->SetToBeBroadcasted(false);

  // Dump a scored quantity
  dumpQtyToFileCmd = new G4UIcommand("/score/dumpQuantityToFile", this);
  dumpQtyToFileCmd->SetGuidance("Dump one scored quantity to file.");
  param = new G4UIparameter("meshName", 's', false);
  dumpQtyToFileCmd->SetParameter(param);
  param = new G4UIparameter("psName", 's', false);
  dumpQtyToFileCmd->SetParameter(param);
  param = new G4UIparameter("fileName", 's', false);
  dumpQtyToFileCmd->SetParameter(param);
  param = new G4UIparameter("option", 's', true);
  dumpQtyToFileCmd->SetParameter(param);
  dumpQtyToFileCmd->SetToBeBroadcasted(false);

  dumpQtyWithFactorCmd = new G4UIcommand("/score/dumpQuantityWithFactor", this);
  dumpQtyWithFactorCmd->SetGuidance("Dump one scored quantity to file.");
  dumpQtyWithFactorCmd->SetGuidance(
    "Each value is multiplied by the specified factor.");
  param = new G4UIparameter("meshName", 's', false);
  dumpQtyWithFactorCmd->SetParameter(param);
  param = new G4UIparameter("psName", 's', false);
  dumpQtyWithFactorCmd->SetParameter(param);
  param = new G4UIparameter("fileName", 's', false);
  dumpQtyWithFactorCmd->SetParameter(param);
  param = new G4UIparameter("factor", 'd', false);
  param->SetParameterRange("factor>0.");
  dumpQtyWithFactorCmd->SetParameter(param);
  param = new G4UIparameter("option", 's', true);
  dumpQtyWithFactorCmd->SetParameter(param);
  dumpQtyWithFactorCmd->SetToBeBroadcasted(false);

  // Dump all scored quantities
  dumpAllQtsToFileCmd = new G4UIcommand("/score/dumpAllQuantitiesToFile", this);
  dumpAllQtsToFileCmd->SetGuidance("Dump all quantities of the mesh to file.");
  param = new G4UIparameter("meshName", 's', false);
  dumpAllQtsToFileCmd->SetParameter(param);
  param = new G4UIparameter("fileName", 's', false);
  dumpAllQtsToFileCmd->SetParameter(param);
  param = new G4UIparameter("option", 's', true);
  dumpAllQtsToFileCmd->SetParameter(param);
  dumpAllQtsToFileCmd->SetToBeBroadcasted(false);

  dumpAllQtsWithFactorCmd =
    new G4UIcommand("/score/dumpAllQuantitiesWithFactor", this);
  dumpAllQtsWithFactorCmd->SetGuidance(
    "Dump all quantities of the mesh to file.");
  dumpAllQtsWithFactorCmd->SetGuidance(
    "Each value is multiplied by the specified factor.");
  param = new G4UIparameter("meshName", 's', false);
  dumpAllQtsWithFactorCmd->SetParameter(param);
  param = new G4UIparameter("fileName", 's', false);
  dumpAllQtsWithFactorCmd->SetParameter(param);
  param = new G4UIparameter("factor", 'd', false);
  param->SetParameterRange("factor>0.");
  dumpAllQtsWithFactorCmd->SetParameter(param);
  param = new G4UIparameter("option", 's', true);
  dumpAllQtsWithFactorCmd->SetParameter(param);
  dumpAllQtsWithFactorCmd->SetToBeBroadcasted(false);

  fill1DCmd = new G4UIcommand("/score/fill1D", this);
  fill1DCmd->SetGuidance("Let a primitive scorer fill 1-D histogram");
  fill1DCmd->SetGuidance("Before using this command, primitive scorer must be "
                         "defined and assigned.");
  fill1DCmd->SetGuidance("Also before using this command, a histogram has to "
                         "be defined by /analysis/h1/create command.");
  fill1DCmd->SetGuidance(
    "This command is available only for real-world volume or probe.");
  fill1DCmd->SetGuidance("Please note that this command has to be applied to "
                         "each copy number of the scoring volume.");
  fill1DCmd->SetGuidance("If same histogram ID is used more than once, more "
                         "than one scorers fill that histogram.");
  param = new G4UIparameter("histID", 'i', false);
  fill1DCmd->SetParameter(param);
  param = new G4UIparameter("meshName", 's', false);
  fill1DCmd->SetParameter(param);
  param = new G4UIparameter("scorerName", 's', false);
  fill1DCmd->SetParameter(param);
  param = new G4UIparameter("copyNo", 'i', true);
  param->SetDefaultValue(0);
  fill1DCmd->SetParameter(param);
}

G4ScoringMessenger::~G4ScoringMessenger()
{
  delete listCmd;
  delete verboseCmd;
  //
  delete meshBoxCreateCmd;
  delete meshCylinderCreateCmd;
  delete meshRWLogVolCreateCmd;
  delete probeCreateCmd;
  delete meshCreateDir;
  //
  delete meshOpnCmd;
  //
  delete meshClsCmd;
  delete meshDir;
  //
  delete mBoxSizeCmd;
  delete mCylinderSizeCmd;
  delete mCylinderRMinCmd;
  delete mCylinderAngleCmd;
  //
  delete mBinCmd;
  //
  delete mTResetCmd;
  delete mTXyzCmd;
  delete mTransDir;
  delete mRotXCmd;
  delete mRotYCmd;
  delete mRotZCmd;
  delete mRotDir;
  //
  delete probeLocateCmd;
  delete probeMatCmd;
  delete probeDir;
  //
  delete dumpCmd;
  delete drawCmd;
  delete drawColumnCmd;
  delete listColorMapCmd;
  delete floatMinMaxCmd;
  delete colorMapMinMaxCmd;
  delete colorMapDir;
  delete dumpQtyToFileCmd;
  delete dumpQtyWithFactorCmd;
  delete dumpAllQtsToFileCmd;
  delete dumpAllQtsWithFactorCmd;
  delete fill1DCmd;
  //
  delete scoreDir;
}

void G4ScoringMessenger::SetNewValue(G4UIcommand* command, G4String newVal)
{
  using MeshShape = G4VScoringMesh::MeshShape;

  if(command == listCmd)
  {
    fSMan->List();
  }
  else if(command == dumpCmd)
  {
    fSMan->Dump();
  }
  else if(command == drawCmd)
  {
    G4Tokenizer next(newVal);
    G4String meshName     = next();
    G4String psName       = next();
    G4String colorMapName = next();
    G4int axflg           = StoI(next());
    fSMan->DrawMesh(meshName, psName, colorMapName, axflg);
  }
  else if(command == drawColumnCmd)
  {
    G4Tokenizer next(newVal);
    G4String meshName     = next();
    G4String psName       = next();
    G4int iPlane          = StoI(next());
    G4int iColumn         = StoI(next());
    G4String colorMapName = next();
    fSMan->DrawMesh(meshName, psName, iPlane, iColumn, colorMapName);
  }
  else if(command == dumpQtyToFileCmd)
  {
    G4Tokenizer next(newVal);
    G4String meshName = next();
    G4String psName   = next();
    G4String fileName = next();
    G4String option   = next("\n");
    auto mesh         = fSMan->FindMesh(meshName);
    if(mesh == nullptr)
    {
      G4ExceptionDescription ed;
      ed << "Mesh name <" << meshName << "> is not found. Command ignored.";
      command->CommandFailed(ed);
      return;
    }
    fSMan->DumpQuantityToFile(meshName, psName, fileName, option);
  }
  else if(command == dumpQtyWithFactorCmd)
  {
    G4Tokenizer next(newVal);
    G4String meshName = next();
    G4String psName   = next();
    G4String fileName = next();
    G4double fac      = StoD(next());
    G4String option   = next("\n");
    auto mesh         = fSMan->FindMesh(meshName);
    if(mesh == nullptr)
    {
      G4ExceptionDescription ed;
      ed << "Mesh name <" << meshName << "> is not found. Command ignored.";
      command->CommandFailed(ed);
      return;
    }
    fSMan->SetFactor(fac);
    fSMan->DumpQuantityToFile(meshName, psName, fileName, option);
    fSMan->SetFactor(1.0);
  }
  else if(command == dumpAllQtsToFileCmd)
  {
    G4Tokenizer next(newVal);
    G4String meshName = next();
    G4String fileName = next();
    G4String option   = next("\n");
    auto mesh         = fSMan->FindMesh(meshName);
    if(mesh == nullptr)
    {
      G4ExceptionDescription ed;
      ed << "Mesh name <" << meshName << "> is not found. Command ignored.";
      command->CommandFailed(ed);
      return;
    }
    fSMan->DumpAllQuantitiesToFile(meshName, fileName, option);
  }
  else if(command == dumpAllQtsWithFactorCmd)
  {
    G4Tokenizer next(newVal);
    G4String meshName = next();
    G4String fileName = next();
    G4double fac      = StoD(next());
    G4String option   = next("\n");
    auto mesh         = fSMan->FindMesh(meshName);
    if(mesh == nullptr)
    {
      G4ExceptionDescription ed;
      ed << "Mesh name <" << meshName << "> is not found. Command ignored.";
      command->CommandFailed(ed);
      return;
    }
    fSMan->SetFactor(fac);
    fSMan->DumpAllQuantitiesToFile(meshName, fileName, option);
    fSMan->SetFactor(1.0);
  }
  else if(command == fill1DCmd)
  {
    Fill1D(command, newVal);
  }
  else if(command == verboseCmd)
  {
    fSMan->SetVerboseLevel(verboseCmd->GetNewIntValue(newVal));
  }
  else if(command == meshBoxCreateCmd)
  {
    G4VScoringMesh* currentmesh = fSMan->GetCurrentMesh();
    if(currentmesh != nullptr)
    {
      G4ExceptionDescription ed;
      ed << "ERROR[" << meshBoxCreateCmd->GetCommandPath() << "] : Mesh <"
         << currentmesh->GetWorldName()
         << "> is still open. Close it first. Command ignored.";
      command->CommandFailed(ed);
    }
    else
    {
      G4VScoringMesh* mesh = fSMan->FindMesh(newVal);
      if(mesh == nullptr)
      {
        mesh = new G4ScoringBox(newVal);
        fSMan->RegisterScoringMesh(mesh);
      }
      else
      {
        G4ExceptionDescription ed;
        ed << "ERROR[" << meshBoxCreateCmd->GetCommandPath()
           << "] : Scoring mesh <" << newVal
           << "> already exists. Command ignored.";
        command->CommandFailed(ed);
      }
    }
  }
  else if(command == meshCylinderCreateCmd)
  {
    G4VScoringMesh* currentmesh = fSMan->GetCurrentMesh();
    if(currentmesh != nullptr)
    {
      G4ExceptionDescription ed;
      ed << "ERROR[" << meshCylinderCreateCmd->GetCommandPath() << "] : Mesh <"
         << currentmesh->GetWorldName()
         << "> is still open. Close it first. Command ignored.";
      command->CommandFailed(ed);
    }
    else
    {
      G4VScoringMesh* mesh = fSMan->FindMesh(newVal);
      if(mesh == nullptr)
      {
        mesh = new G4ScoringCylinder(newVal);
        fSMan->RegisterScoringMesh(mesh);
      }
      else
      {
        G4ExceptionDescription ed;
        ed << "ERROR[" << meshCylinderCreateCmd->GetCommandPath()
           << "] : Scoring mesh <" << newVal
           << "> already exists. Command ignored.";
        command->CommandFailed(ed);
      }
    }
  }
  else if(command == meshRWLogVolCreateCmd)
  {
    auto mesh = fSMan->GetCurrentMesh();
    if(mesh != nullptr)
    {
      G4ExceptionDescription ed;
      ed << "ERROR[" << meshRWLogVolCreateCmd->GetCommandPath() << "] : Mesh <"
         << mesh->GetWorldName()
         << "> is still open. Close it first. Command ignored.";
      command->CommandFailed(ed);
    }
    else
    {
      G4Tokenizer next(newVal);
      G4String meshName = next();
      G4int idx         = StoI(next());
      mesh              = fSMan->FindMesh(meshName);
      if(mesh == nullptr)
      {
        mesh = new G4ScoringRealWorld(meshName);
        mesh->SetCopyNumberLevel(idx);
        fSMan->RegisterScoringMesh(mesh);
      }
      else
      {
        G4ExceptionDescription ed;
        ed << "ERROR[" << meshRWLogVolCreateCmd->GetCommandPath()
           << "] : Scoring mesh <" << meshName
           << "> already exists. Command ignored.";
        command->CommandFailed(ed);
      }
    }
  }
  else if(command == probeCreateCmd)
  {
    auto mesh = fSMan->GetCurrentMesh();
    if(mesh != nullptr)
    {
      G4ExceptionDescription ed;
      ed << "ERROR[" << meshRWLogVolCreateCmd->GetCommandPath() << "] : Mesh <"
         << mesh->GetWorldName()
         << "> is still open. Close it first. Command ignored.";
      command->CommandFailed(ed);
    }
    else
    {
      G4Tokenizer next(newVal);
      G4String qname    = next();
      G4double halfSize = StoD(next());
      halfSize *= G4UIcommand::ValueOf(next());
      G4bool checkOverlap = StoB(next());
      mesh                = fSMan->FindMesh(qname);
      if(mesh == nullptr)
      {
        mesh = new G4ScoringProbe(qname, halfSize, checkOverlap);
        fSMan->RegisterScoringMesh(mesh);
      }
      else
      {
        G4ExceptionDescription ed;
        ed << "ERROR[" << probeCreateCmd->GetCommandPath() << "] : Mesh name <"
           << qname << "> already exists. Use another name.";
        command->CommandFailed(ed);
      }
    }
  }
  else if(command == probeMatCmd || command == probeLocateCmd)
  {
    auto mesh = fSMan->GetCurrentMesh();
    if(mesh == nullptr)
    {
      G4ExceptionDescription ed;
      ed << "ERROR : No mesh is currently open. Open/create a mesh first. "
            "Command ignored.";
      command->CommandFailed(ed);
      return;
    }
    if(mesh->GetShape() != MeshShape::probe)
    {
      G4ExceptionDescription ed;
      ed << "ERROR : Inconsistent mesh type. Close current mesh and open "
            "Scoring Probe.";
      command->CommandFailed(ed);
      return;
    }

    if(command == probeMatCmd)
    {
      G4bool succ = static_cast<G4ScoringProbe*>(mesh)->SetMaterial(newVal);
      if(!succ)
      {
        G4ExceptionDescription ed;
        ed << "Material <" << newVal
           << "> is not defind in G4NistManager. Command is ignored.\n"
           << "Use /material/nist/listMaterials command to see the available "
              "materials.";
        command->CommandFailed(ed);
        return;
      }
    }
    else if(command == probeLocateCmd)
    {
      G4ThreeVector loc = probeLocateCmd->GetNew3VectorValue(newVal);
      static_cast<G4ScoringProbe*>(mesh)->LocateProbe(loc);
    }
  }
  else if(command == listColorMapCmd)
  {
    fSMan->ListScoreColorMaps();
  }
  else if(command == floatMinMaxCmd)
  {
    G4VScoreColorMap* colorMap = fSMan->GetScoreColorMap(newVal);
    if(colorMap != nullptr)
    {
      colorMap->SetFloatingMinMax(true);
    }
    else
    {
      G4ExceptionDescription ed;
      ed << "ERROR[" << floatMinMaxCmd->GetCommandPath() << "] : color map <"
         << newVal << "> is not defined. Command ignored.";
      command->CommandFailed(ed);
    }
  }
  else if(command == colorMapMinMaxCmd)
  {
    G4Tokenizer next(newVal);
    G4String mapName           = next();
    G4double minVal            = StoD(next());
    G4double maxVal            = StoD(next());
    G4VScoreColorMap* colorMap = fSMan->GetScoreColorMap(mapName);
    if(colorMap != nullptr)
    {
      colorMap->SetFloatingMinMax(false);
      colorMap->SetMinMax(minVal, maxVal);
    }
    else
    {
      G4ExceptionDescription ed;
      ed << "ERROR[" << colorMapMinMaxCmd->GetCommandPath() << "] : color map <"
         << newVal << "> is not defined. Command ignored." << G4endl;
      command->CommandFailed(ed);
    }
  }
  else if(command == meshOpnCmd)
  {
    G4VScoringMesh* currentmesh = fSMan->GetCurrentMesh();
    if(currentmesh != nullptr)
    {
      G4ExceptionDescription ed;
      ed << "ERROR[" << meshOpnCmd->GetCommandPath() << "] : Mesh <"
         << currentmesh->GetWorldName()
         << "> is still open. Close it first. Command ignored.";
      command->CommandFailed(ed);
    }
    else
    {
      G4VScoringMesh* mesh = fSMan->FindMesh(newVal);
      if(mesh == nullptr)
      {
        G4ExceptionDescription ed;
        ed << "ERROR[" << meshOpnCmd->GetCommandPath() << "] : Scoring mesh <"
           << newVal << "> does not exist. Command ignored.";
        command->CommandFailed(ed);
      }
      else
      {
        fSMan->SetCurrentMesh(mesh);
      }
    }
  }
  else if(command == meshClsCmd)
  {
    fSMan->CloseCurrentMesh();
  }
  else
  {
    //
    // Get Current Mesh
    //
    G4VScoringMesh* mesh = fSMan->GetCurrentMesh();
    //
    // Commands for Current Mesh
    if(mesh != nullptr)
    {
      MeshShape shape = mesh->GetShape();
      if(shape == MeshShape::realWorldLogVol)
      {
        G4ExceptionDescription ed;
        ed << "ERROR[" << mBinCmd->GetCommandPath()
           << "] : Number of mesh command cannot be set for this type of mesh. "
              "Command ignored.";
        command->CommandFailed(ed);
      }
      else
      {
        // Tokens
        G4TokenVec token;
        FillTokenVec(newVal, token);
        //
        // Mesh Geometry
        if(command == mBoxSizeCmd)
        {
          if(shape == MeshShape::box)
          {
            G4ThreeVector size = mBoxSizeCmd->GetNew3VectorValue(newVal);
            G4double vsize[3];
            vsize[0] = size.x();
            vsize[1] = size.y();
            vsize[2] = size.z();
            mesh->SetSize(vsize);
          }
          else
          {
            G4ExceptionDescription ed;
            ed << "ERROR[" << mBoxSizeCmd->GetCommandPath()
               << "] : This mesh is not Box. Command ignored.";
            command->CommandFailed(ed);
          }
        }
        else if(command == mCylinderSizeCmd || command == mCylinderRMinCmd ||
                command == mCylinderAngleCmd)
        {
          if(shape != MeshShape::cylinder)
          {
            G4ExceptionDescription ed;
            ed << "ERROR[" << command->GetCommandPath()
               << "] : This mesh is not Cylinder. Command ignored.";
            command->CommandFailed(ed);
          }
          else
          {
            if(command == mCylinderSizeCmd)
            {
              G4double vsize[3];
              vsize[0]     = (mesh->GetSize()).x();
              vsize[1]     = StoD(token[0]);
              vsize[2]     = StoD(token[1]);
              G4double unt = mCylinderSizeCmd->ValueOf(token[2]);
              vsize[1] *= unt;
              vsize[2] *= unt;
              mesh->SetSize(vsize);
            }
            else if(command == mCylinderRMinCmd)
            {
              G4double vsize[3];
              vsize[0] = mCylinderRMinCmd->GetNewDoubleValue(newVal);
              vsize[1] = (mesh->GetSize()).y();
              vsize[2] = (mesh->GetSize()).z();
              mesh->SetSize(vsize);
            }
            else if(command == mCylinderAngleCmd)
            {
              G4double stphi = StoD(token[0]);
              G4double spphi = StoD(token[1]);
              G4double unt   = mCylinderAngleCmd->ValueOf(token[2]);
              mesh->SetAngles(stphi * unt, spphi * unt);
            }
          }
        }
        else if(command == mBinCmd)
        {
          MeshBinCommand(mesh, token);
        }
        else if(command == mTResetCmd)
        {
          G4double centerPosition[3] = { 0., 0., 0. };
          mesh->SetCenterPosition(centerPosition);
        }
        else if(command == mTXyzCmd)
        {
          G4ThreeVector xyz = mTXyzCmd->GetNew3VectorValue(newVal);
          G4double centerPosition[3];
          centerPosition[0] = xyz.x();
          centerPosition[1] = xyz.y();
          centerPosition[2] = xyz.z();
          mesh->SetCenterPosition(centerPosition);
        }
        else if(command == mRotXCmd)
        {
          G4double value = mRotXCmd->GetNewDoubleValue(newVal);
          mesh->RotateX(value);
        }
        else if(command == mRotYCmd)
        {
          G4double value = mRotYCmd->GetNewDoubleValue(newVal);
          mesh->RotateY(value);
        }
        else if(command == mRotZCmd)
        {
          G4double value = mRotZCmd->GetNewDoubleValue(newVal);
          mesh->RotateZ(value);
        }
      }
    }
    else
    {
      G4ExceptionDescription ed;
      ed << "ERROR: No mesh is currently open. Open/create a mesh first. "
            "Command ignored.";
      command->CommandFailed(ed);
    }
  }
}

G4String G4ScoringMessenger::GetCurrentValue(G4UIcommand* command)
{
  G4String val;
  if(command == verboseCmd)
  {
    val = verboseCmd->ConvertToString(fSMan->GetVerboseLevel());
  }

  return val;
}

void G4ScoringMessenger::FillTokenVec(G4String newValues, G4TokenVec& token)
{
  G4Tokenizer next(newValues);
  G4String val;
  while(!(val = next()).empty())
  {  // Loop checking 12.18.2015 M.Asai
    token.push_back(val);
  }
}

void G4ScoringMessenger::MeshBinCommand(G4VScoringMesh* mesh, G4TokenVec& token)
{
  G4int Ni = StoI(token[0]);
  G4int Nj = StoI(token[1]);
  G4int Nk = StoI(token[2]);
  G4int nSegment[3];

  if(dynamic_cast<G4ScoringBox*>(mesh) != nullptr)
  {
    G4cout << ".... G4ScoringMessenger::MeshBinCommand - G4ScoringBox"
           << G4endl;
    nSegment[0] = Ni;
    nSegment[1] = Nj;
    nSegment[2] = Nk;
  }
  else if(dynamic_cast<G4ScoringCylinder*>(mesh) != nullptr)
  {
    G4cout << ".... G4ScoringMessenger::MeshBinCommand - G4ScoringCylinder"
           << G4endl;
    nSegment[0] = Nj;
    nSegment[1] = Nk;
    nSegment[2] = Ni;
  }
  else
  {
    G4Exception("G4ScoringMessenger::MeshBinCommand()", "001", FatalException,
                "invalid mesh type");
    return;
  }
  //
  mesh->SetNumberOfSegments(nSegment);
}

#include "G4VPrimitivePlotter.hh"
#include "G4VScoreHistFiller.hh"

void G4ScoringMessenger::Fill1D(G4UIcommand* cmd, G4String newVal)
{
  using MeshShape = G4VScoringMesh::MeshShape;

  G4Tokenizer next(newVal);
  G4int histID      = StoI(next());
  G4String meshName = next();
  G4String primName = next();
  G4int copyNo      = StoI(next());

  auto filler = G4VScoreHistFiller::Instance();
  if(filler == nullptr)
  {
    G4ExceptionDescription ed;
    ed << "G4TScoreHistFiller is not instantiated in this application.";
    cmd->CommandFailed(ed);
    return;
  }
  //
  // To do : check the validity of histID
  //

  auto sm   = G4ScoringManager::GetScoringManagerIfExist();
  auto mesh = sm->FindMesh(meshName);
  if(mesh == nullptr)
  {
    G4ExceptionDescription ed;
    ed << "Mesh name <" << meshName << "> is not found.";
    cmd->CommandFailed(ed);
    return;
  }
  auto shape = mesh->GetShape();
  if(shape != MeshShape::realWorldLogVol && shape != MeshShape::probe)
  {
    G4ExceptionDescription ed;
    ed << "Mesh <" << meshName
       << "> is not real-world logical volume or probe.";
    cmd->CommandFailed(ed);
    return;
  }

  auto prim = mesh->GetPrimitiveScorer(primName);
  if(prim == nullptr)
  {
    G4ExceptionDescription ed;
    ed << "Primitive scorer name <" << primName << "> is not found.";
    cmd->CommandFailed(ed);
    return;
  }
  auto pp = dynamic_cast<G4VPrimitivePlotter*>(prim);
  if(pp == nullptr)
  {
    G4ExceptionDescription ed;
    ed << "Primitive scorer <" << primName
       << "> does not support direct histogram filling.";
    cmd->CommandFailed(ed);
    return;
  }

  pp->Plot(copyNo, histID);
}
