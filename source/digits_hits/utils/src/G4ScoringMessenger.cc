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
// $Id: G4ScoringMessenger.cc 108471 2018-02-15 14:12:06Z gcosmo $
//
// ---------------------------------------------------------------------

#include "G4ScoringMessenger.hh"
#include "G4ScoringManager.hh"
#include "G4VScoringMesh.hh"
#include "G4ScoringBox.hh"
#include "G4ScoringCylinder.hh"

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
#include "G4VScoreColorMap.hh"

G4ScoringMessenger::G4ScoringMessenger(G4ScoringManager* SManager)
:fSMan(SManager)
{
  G4UIparameter* param;

  scoreDir = new G4UIdirectory("/score/");
  scoreDir->SetGuidance("Interactive scoring commands.");

  listCmd = new G4UIcmdWithoutParameter("/score/list",this);
  listCmd->SetGuidance("List scoring worlds.");

  dumpCmd = new G4UIcmdWithoutParameter("/score/dump",this);
  dumpCmd->SetGuidance("Dump results of scorers.");

  verboseCmd = new G4UIcmdWithAnInteger("/score/verbose",this);
  verboseCmd->SetGuidance("Verbosity.");
  verboseCmd->SetGuidance("  0) errors or warnings,");
  verboseCmd->SetGuidance("  1) information with 0)");

  meshDir = new G4UIdirectory("/score/mesh/");
  meshDir->SetGuidance("    Mesh processing commands.");

  meshCreateDir = new G4UIdirectory("/score/create/");
  meshCreateDir->SetGuidance("  Mesh creation commands.");
  //
  // Mesh commands
  meshBoxCreateCmd = new G4UIcmdWithAString("/score/create/boxMesh",this);
  meshBoxCreateCmd->SetGuidance("Create scoring box mesh.");
  meshBoxCreateCmd->SetParameterName("MeshName",false);
  //
  meshCylinderCreateCmd = new G4UIcmdWithAString("/score/create/cylinderMesh",this);
  meshCylinderCreateCmd->SetGuidance("Create scoring mesh.");
  meshCylinderCreateCmd->SetParameterName("MeshName",false);
  //
//  meshSphereCreateCmd = new G4UIcmdWithAString("/score/create/sphereMesh",this);
//  meshSphereCreateCmd->SetGuidance("Create scoring mesh.");
//  meshSphereCreateCmd->SetParameterName("MeshName",false);
  //
  meshOpnCmd = new G4UIcmdWithAString("/score/open",this);
  meshOpnCmd->SetGuidance("Open scoring mesh.");
  meshOpnCmd->SetParameterName("MeshName",false);
  //
  meshClsCmd = new G4UIcmdWithoutParameter("/score/close",this);
  meshClsCmd->SetGuidance("Close scoring mesh.");
  //
//  meshActCmd = new G4UIcmdWithABool("/score/mesh/activate",this);
//  meshActCmd->SetGuidance("Activate scoring mesh.");
//  meshActCmd->SetParameterName("MeshName",false);
  //
  mBoxSizeCmd = new G4UIcmdWith3VectorAndUnit("/score/mesh/boxSize",this);
  mBoxSizeCmd->SetGuidance("Define size of the scoring mesh.");
  mBoxSizeCmd->SetGuidance("Dx  Dy  Dz  unit");
  mBoxSizeCmd->SetParameterName("Di","Dj","Dk",false,false);
  mBoxSizeCmd->SetRange("Di>0. && Dj>0. && Dk>0.");
  mBoxSizeCmd->SetDefaultUnit("mm");
  //
  mCylinderSizeCmd = new G4UIcommand("/score/mesh/cylinderSize",this);
  mCylinderSizeCmd->SetGuidance("Define size of the scoring mesh.");
  mCylinderSizeCmd->SetGuidance("R   Dz  unit");
  param = new G4UIparameter("R",'d',false);
  param->SetParameterRange("R>0");
  mCylinderSizeCmd->SetParameter(param);
  param = new G4UIparameter("Dz",'d',false);
  param->SetParameterRange("Dz>0");
  mCylinderSizeCmd->SetParameter(param);
  param = new G4UIparameter("unit",'s',true);
  param->SetDefaultValue("mm");
  mCylinderSizeCmd->SetParameter(param);
  //
  //   Division command
  mBinCmd = new G4UIcommand("/score/mesh/nBin",this);
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
  //mBinCmd->SetGuidance("  Axis:(int) Axis of division ");
//  mBinCmd->SetGuidance("  P1..Pn-1  :(double) \"paramter from P1 to Pn-1 for division.\"");
  param = new G4UIparameter("Ni",'i',false);
  param->SetDefaultValue("1");
  param->SetParameterRange("Ni>0");
  mBinCmd->SetParameter(param);
  param = new G4UIparameter("Nj",'i',false);
  param->SetDefaultValue("1");
  param->SetParameterRange("Nj>0");
  mBinCmd->SetParameter(param);
  param = new G4UIparameter("Nk",'i',false);
  param->SetDefaultValue("1");
  mBinCmd->SetParameter(param);
  param->SetParameterRange("Nk>0");
  //param = new G4UIparameter("Axis",'i',true);
  //param->SetDefaultValue("3");
  //mBinCmd->SetParameter(param);
  //
  //   Placement command
  mTransDir = new G4UIdirectory("/score/mesh/translate/");
  mTransDir->SetGuidance("Mesh translation commands.");
  //
  mTResetCmd = new G4UIcmdWithoutParameter("/score/mesh/translate/reset",this);
  mTResetCmd->SetGuidance("Reset translated position of the scoring mesh.");
  //
  mTXyzCmd = new G4UIcmdWith3VectorAndUnit("/score/mesh/translate/xyz",this);
  mTXyzCmd->SetGuidance("Translate the scoring mesh.");
  mTXyzCmd->SetParameterName("X","Y","Z",false,false);
  mTXyzCmd->SetDefaultUnit("mm");
  //
  mRotDir = new G4UIdirectory("/score/mesh/rotate/");
  mRotDir->SetGuidance("Mesh rotation commands.");
  //
  mRResetCmd = new G4UIcmdWithoutParameter("/score/mesh/rotate/reset",this);
  mRResetCmd->SetGuidance("Reset rotation angles of the scoring mesh.");
  //
  mRotXCmd = new G4UIcmdWithADoubleAndUnit("/score/mesh/rotate/rotateX",this);
  mRotXCmd->SetGuidance("Rotate the scoring mesh in X axis.");
  mRotXCmd->SetParameterName("Rx",false);
  mRotXCmd->SetDefaultUnit("deg");
  //
  mRotYCmd = new G4UIcmdWithADoubleAndUnit("/score/mesh/rotate/rotateY",this);
  mRotYCmd->SetGuidance("Rotate the scoring mesh in Y axis.");
  mRotYCmd->SetParameterName("Ry",false);
  mRotYCmd->SetDefaultUnit("deg");
  //
  mRotZCmd = new G4UIcmdWithADoubleAndUnit("/score/mesh/rotate/rotateZ",this);
  mRotZCmd->SetGuidance("Rotate the scoring mesh in Z axis.");
  mRotZCmd->SetParameterName("Rz",false);
  mRotZCmd->SetDefaultUnit("deg");
  //

  // Draw Scoring result
  drawCmd = new G4UIcommand("/score/drawProjection",this);
  drawCmd->SetGuidance("Draw projection(s) of scored quantities.");
  drawCmd->SetGuidance("Parameter <proj> specified which projection(s) to be drawn.");
  drawCmd->SetGuidance("  100 : xy-plane, 010 : yz-plane,    001 : zx-plane -- default 111");
  drawCmd->SetGuidance("  100 : N/A,      010 : z_phi-plane, 001 : r_phi-plane -- default 111");
  param = new G4UIparameter("meshName",'s',false);
  drawCmd->SetParameter(param);
  param = new G4UIparameter("psName",'s',false);
  drawCmd->SetParameter(param);
  param = new G4UIparameter("colorMapName",'s',true);
  param->SetDefaultValue("defaultLinearColorMap");
  drawCmd->SetParameter(param);
  param = new G4UIparameter("proj",'i',true);
  param->SetDefaultValue(111);
  drawCmd->SetParameter(param);
  drawCmd->SetToBeBroadcasted(false);

  // Draw column
  drawColumnCmd = new G4UIcommand("/score/drawColumn",this);
  drawColumnCmd->SetGuidance("Draw a cell column.");
  drawColumnCmd->SetGuidance(" plane = 0 : x-y, 1: y-z, 2: z-x  for box mesh");
  drawColumnCmd->SetGuidance("         0 : z-phi, 1: r-phi, 2: r-z  for cylinder mesh");
  param = new G4UIparameter("meshName",'s',false);
  drawColumnCmd->SetParameter(param);
  param = new G4UIparameter("psName",'s',false);
  drawColumnCmd->SetParameter(param);
  param = new G4UIparameter("plane",'i',false);
  param->SetParameterRange("plane>=0 && plane<=2");
  drawColumnCmd->SetParameter(param);
  param = new G4UIparameter("column",'i',false);
  drawColumnCmd->SetParameter(param);
  param = new G4UIparameter("colorMapName",'s',true);
  param->SetDefaultValue("defaultLinearColorMap");
  drawColumnCmd->SetParameter(param);
  drawColumnCmd->SetToBeBroadcasted(false);

  colorMapDir = new G4UIdirectory("/score/colorMap/");
  colorMapDir->SetGuidance("Color map commands.");

  listColorMapCmd = new G4UIcmdWithoutParameter("/score/colorMap/listScoreColorMaps",this);
  listColorMapCmd->SetGuidance("List registered score color maps.");
  listColorMapCmd->SetToBeBroadcasted(false);

  floatMinMaxCmd = new G4UIcmdWithAString("/score/colorMap/floatMinMax",this);
  floatMinMaxCmd->SetGuidance("Min/Max of the color map is calculated accorging to the actual scores.");
  floatMinMaxCmd->SetParameterName("colorMapName",true,false);
  floatMinMaxCmd->SetDefaultValue("defaultLinearColorMap");
  floatMinMaxCmd->SetToBeBroadcasted(false);

  colorMapMinMaxCmd = new G4UIcommand("/score/colorMap/setMinMax",this);
  colorMapMinMaxCmd->SetGuidance("Define min/max value of the color map.");
  param = new G4UIparameter("colorMapMame",'s',true);
  param->SetDefaultValue("defaultLinearColorMap");
  colorMapMinMaxCmd->SetParameter(param);
  param = new G4UIparameter("minValue",'d',false);
  colorMapMinMaxCmd->SetParameter(param);
  param = new G4UIparameter("maxValue",'d',false);
  colorMapMinMaxCmd->SetParameter(param);
  colorMapMinMaxCmd->SetToBeBroadcasted(false);

  /*
  chartCmd = new G4UIcommand("/score/drawChart",this);
  chartCmd->SetGuidance("Draw color chart on the screen.");
  chartCmd->SetGuidance("[usage] /score/drawChart");
  chartCmd->SetGuidance("  mesh    :(String) Mesh name.");
  chartCmd->SetGuidance("  psname  :(String) PS name.");
  chartCmd->SetGuidance("  On/Off  :(boolean) On or Off the color chart.");
  chartCmd->SetGuidance("  scale   :(String) default=linear, or log ");
  param = new G4UIparameter("meshName",'s',false);
  chartCmd->SetParameter(param);
  param = new G4UIparameter("psName",'s',false);
  chartCmd->SetParameter(param);
  param = new G4UIparameter("On",'s',true);
  param->SetDefaultValue("true");
  chartCmd->SetParameter(param);
  param = new G4UIparameter("scale",'s',true);
  param->SetDefaultValue("linear");
  chartCmd->SetParameter(param);
  */

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

}

G4ScoringMessenger::~G4ScoringMessenger()
{
    delete listCmd;
    delete verboseCmd;
    //
    delete           meshCreateDir;
    delete           meshBoxCreateCmd;
    delete           meshCylinderCreateCmd;
//    delete           meshSphereCreateCmd;
    //
    delete          meshOpnCmd;
    //
    delete    meshClsCmd;
//    delete    meshActCmd;
    delete          meshDir;
    //
    delete  mBoxSizeCmd;
    delete  mCylinderSizeCmd;
//    delete  mSphereSizeCmd;
    //
    delete      mBinCmd;
    //
    delete   mTResetCmd;
    delete   mTXyzCmd;
    delete   mTransDir;
    delete   mRResetCmd;
    delete   mRotXCmd;
    delete   mRotYCmd;
    delete   mRotZCmd;
    delete   mRotDir;
    //
    //delete     chartCmd;
    delete     dumpCmd;
    delete     drawCmd;
    delete     drawColumnCmd;
    delete     listColorMapCmd;
    delete     floatMinMaxCmd;
    delete     colorMapMinMaxCmd;
    delete     colorMapDir;
    delete     dumpQtyToFileCmd;
    delete     dumpAllQtsToFileCmd;
    //
    delete scoreDir;
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
      G4String colorMapName = next();
      G4int axflg = StoI(next());
      fSMan->DrawMesh(meshName,psName,colorMapName,axflg);
  } else if(command==drawColumnCmd) { 
      G4Tokenizer next(newVal);
      G4String meshName = next();
      G4String psName = next();
      G4int iPlane = StoI(next());
      G4int iColumn = StoI(next());
      G4String colorMapName = next();
      fSMan->DrawMesh(meshName,psName,iPlane,iColumn,colorMapName);
//  } else if(command==chartCmd ){
//      G4Tokenizer next(newVal);
//      G4String meshName = next();
//      G4String psName = next();
//      //G4bool   onFlag = StoB(next());
//      G4String scaleOption = next();
//      fSMan->DrawChart(meshName,psName,onFlag,scaleOption);

  } else if(command==dumpQtyToFileCmd) { 
      G4Tokenizer next(newVal);
      G4String meshName = next();
      G4String psName = next();
      G4String fileName = next();
      G4String option = next("\n");
      fSMan->DumpQuantityToFile(meshName, psName, fileName, option);
  } else if(command==dumpAllQtsToFileCmd) { 
      G4Tokenizer next(newVal);
      G4String meshName = next();
      G4String fileName = next();
      G4String option = next("\n");
      fSMan->DumpAllQuantitiesToFile(meshName, fileName, option);
  } else if(command==verboseCmd) { 
      fSMan->SetVerboseLevel(verboseCmd->GetNewIntValue(newVal)); 
  } else if(command==meshBoxCreateCmd) {
      G4VScoringMesh* currentmesh = fSMan->GetCurrentMesh(); 
      if ( currentmesh ){
        G4ExceptionDescription ed;
	ed << "ERROR[" << meshBoxCreateCmd->GetCommandPath()
	   << "] : Mesh <" << currentmesh->GetWorldName() 
	   << "> is still open. Close it first. Command ignored.";
        command->CommandFailed(ed);
      } else {

	G4VScoringMesh*  mesh = fSMan->FindMesh(newVal);
	if ( !mesh ){
	  mesh = new G4ScoringBox(newVal);
	  fSMan->RegisterScoringMesh(mesh);
	}else{
          G4ExceptionDescription ed;
	  ed << "ERROR[" << meshBoxCreateCmd->GetCommandPath()
	     << "] : Scoring mesh <" << newVal
	     << "> already exists. Command ignored.";
          command->CommandFailed(ed);
	}
      }
  } else if(command==meshCylinderCreateCmd) {
      G4VScoringMesh* currentmesh = fSMan->GetCurrentMesh(); 
      if ( currentmesh ){
        G4ExceptionDescription ed;
	ed << "ERROR[" << meshCylinderCreateCmd->GetCommandPath()
	   << "] : Mesh <" << currentmesh->GetWorldName() 
	   << "> is still open. Close it first. Command ignored.";
        command->CommandFailed(ed);
      } else {

	G4VScoringMesh*  mesh = fSMan->FindMesh(newVal);
	if ( !mesh ){
	  mesh = new G4ScoringCylinder(newVal);
	  fSMan->RegisterScoringMesh(mesh);
	}else{
          G4ExceptionDescription ed;
	  ed << "ERROR[" << meshCylinderCreateCmd->GetCommandPath()
	     << "] : Scoring mesh <" << newVal
	     << "> already exists. Command ignored.";
          command->CommandFailed(ed);
	}
      }
  } else if(command==listColorMapCmd) {
      fSMan->ListScoreColorMaps();
  } else if(command==floatMinMaxCmd) {
      G4VScoreColorMap* colorMap = fSMan->GetScoreColorMap(newVal);
      if(colorMap)
      { colorMap->SetFloatingMinMax(true); }
      else
      {
        G4ExceptionDescription ed;
        ed << "ERROR[" << floatMinMaxCmd->GetCommandPath()
	       << "] : color map <" << newVal << "> is not defined. Command ignored." ;
        command->CommandFailed(ed);
      }
  } else if(command==colorMapMinMaxCmd) {
      G4Tokenizer next(newVal);
      G4String mapName = next();
      G4double minVal = StoD(next());
      G4double maxVal = StoD(next());
      G4VScoreColorMap* colorMap = fSMan->GetScoreColorMap(mapName);
      if(colorMap)
      { colorMap->SetFloatingMinMax(false);
        colorMap->SetMinMax(minVal,maxVal); }
      else
      { 
        G4ExceptionDescription ed;
        ed << "ERROR[" << colorMapMinMaxCmd->GetCommandPath()
	   << "] : color map <" << newVal << "> is not defined. Command ignored." 
	   << G4endl; 
        command->CommandFailed(ed);
      }
  } else if(command==meshOpnCmd) {
      G4VScoringMesh* currentmesh = fSMan->GetCurrentMesh(); 
      if ( currentmesh ){
        G4ExceptionDescription ed;
	ed << "ERROR[" << meshOpnCmd->GetCommandPath()
	   << "] : Mesh <" << currentmesh->GetWorldName()
           << "> is still open. Close it first. Command ignored.";
        command->CommandFailed(ed);
      } else {
	G4VScoringMesh* mesh = fSMan->FindMesh(newVal); 
	if ( !mesh ){
          G4ExceptionDescription ed;
          ed << "ERROR[" << meshOpnCmd->GetCommandPath()
	     << "] : Scoring mesh <" << newVal << "> does not exist. Command ignored.";
          command->CommandFailed(ed);
	} else {
	  fSMan->SetCurrentMesh(mesh);
	}
      }
  } else if(command==meshClsCmd) {
      fSMan->CloseCurrentMesh();
  } else {
      //
      // Get Current Mesh
      //
      G4VScoringMesh* mesh = fSMan->GetCurrentMesh();
      //
      // Commands for Current Mesh
      if ( mesh ){
          // Tokens
          G4TokenVec token;
          FillTokenVec(newVal,token);
	  // 
	  // Mesh Geometry
	  //
//	  if(command==meshActCmd) {
//	      mesh->Activate(meshActCmd->GetNewBoolValue(newVal)); 
//	  } else
          if(command==mBoxSizeCmd) {
	      MeshShape shape = mesh->GetShape();
	      if ( shape == boxMesh ){
		  G4ThreeVector size = mBoxSizeCmd->GetNew3VectorValue(newVal);
		  G4double vsize[3];
		  vsize[0] = size.x();
		  vsize[1] = size.y();
		  vsize[2] = size.z();
		  mesh->SetSize(vsize);
	      } else {
                 G4ExceptionDescription ed;
                 ed << "ERROR[" << mBoxSizeCmd->GetCommandPath()
		    << "] : This mesh is not Box. Command ignored.";
                 command->CommandFailed(ed);
	      }
	  }else if(command==mCylinderSizeCmd) {
	      MeshShape shape = mesh->GetShape();
	      if ( shape == cylinderMesh ){
		  G4double vsize[3];
		  vsize[0]     = StoD(token[0]);
		  vsize[1]     = StoD(token[1]);
		  G4double unt = mCylinderSizeCmd->ValueOf(token[2]);
		  vsize[0] *= unt;
		  vsize[1] *= unt;
		  vsize[2]  = 0.0;
		  mesh->SetSize(vsize);
	      } else {
                 G4ExceptionDescription ed;
		 ed << "ERROR[" << mBoxSizeCmd->GetCommandPath()
		     << "] : This mesh is not Box. Command ignored.";
                 command->CommandFailed(ed);
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
	  }
      }else{
        G4ExceptionDescription ed;
        ed << "ERROR: No mesh is currently open. Open/create a mesh first. Command ignored.";
        command->CommandFailed(ed);
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
    while ( !(val = next()).isNull() ) { // Loop checking 12.18.2015 M.Asai
	token.push_back(val);
    }
}


void G4ScoringMessenger::MeshBinCommand(G4VScoringMesh* mesh,G4TokenVec& token){
    G4int Ni = StoI(token[0]);
    G4int Nj = StoI(token[1]);
    G4int Nk = StoI(token[2]);
    G4int nSegment[3];

    if(dynamic_cast<G4ScoringBox*>(mesh)) {
      G4cout << ".... G4ScoringMessenger::MeshBinCommand - G4ScoringBox" << G4endl;
      nSegment[0] = Ni;
      nSegment[1] = Nj;
      nSegment[2] = Nk;
    } else if(dynamic_cast<G4ScoringCylinder*>(mesh)) {
      G4cout << ".... G4ScoringMessenger::MeshBinCommand - G4ScoringCylinder" << G4endl;
      nSegment[0] = Nj;
      nSegment[1] = Nk;
      nSegment[2] = Ni;
    } else {
      G4Exception("G4ScoringMessenger::MeshBinCommand()", "001", FatalException, "invalid mesh type");
      return;
    }
    //
    mesh->SetNumberOfSegments(nSegment);
}

 
