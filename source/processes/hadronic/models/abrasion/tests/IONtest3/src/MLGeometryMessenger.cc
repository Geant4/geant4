////////////////////////////////////////////////////////////////////////////////
//
#include "MLGeometryMessenger.hh"

#include "MLGeometryConstruction.hh"
#include "MLMaterial.hh"
////////////////////////////////////////////////////////////////////////////////
//
MLGeometryMessenger::MLGeometryMessenger (MLGeometryConstruction * MLDet)
  :MLGeometry(MLDet)
{
  MLdetDir = new G4UIdirectory("/geometry/");
  MLdetDir->SetGuidance("MULASSIS geometry controls.");

  LayerDir = new G4UIdirectory("/geometry/layer/");
  LayerDir->SetGuidance("Controls for defining the layers.");

  ShapeCmd = new G4UIcmdWithAString("/geometry/layer/shape",this);
  ShapeCmd->SetGuidance("Select the geometry shape.");
  ShapeCmd->SetGuidance("slab or sphere");
  ShapeCmd->SetParameterName("Shape",false);
  ShapeCmd->SetCandidates("slab SLAB sphere SPHERE");
  ShapeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);


  DeleteCmd = new G4UIcmdWithAnInteger("/geometry/layer/delete",this);
  DeleteCmd->SetGuidance("Delete a layer or all layers.");
  DeleteCmd->SetGuidance("(If Nb = 0 : all layers are deleted)");
  DeleteCmd->SetParameterName("Delete",false);
  DeleteCmd->SetRange("Delete >= 0");
  DeleteCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ListCmd = new G4UIcmdWithAnInteger("/geometry/layer/list",this);
  ListCmd->SetGuidance("List a layer or all layers.");
  ListCmd->SetGuidance("(If Nb = 0: all layers are listed");
  ListCmd->SetParameterName("List",true);
  ListCmd->SetRange("List >= 0");
  ListCmd->SetDefaultValue(0);
  ListCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  AddCmd = new G4UIcommand("/geometry/layer/add",this);
  AddCmd->SetGuidance("Add a layer after layer NL, of material MA, colour index CI and thickness (unit).");
  AddCmd->SetGuidance(
    "  position is after layer NL: (NL>=0 && NL <= NbOfLayers) ");
  AddCmd->SetGuidance("  material: name");
  AddCmd->SetGuidance("  colour index: ci >0 && <= 12");
  AddCmd->SetGuidance("  thickness (with unit) : t>0.");
  //
  G4UIparameter* AbsNbPrm = new G4UIparameter("AddNb",'i',false);
  AbsNbPrm->SetGuidance("position number : from 0 to NbOfLayers");
  AbsNbPrm->SetParameterRange("AddNb >= 0");
  AddCmd->SetParameter(AbsNbPrm);
  //
  G4UIparameter* MatPrm = new G4UIparameter("material",'s',false);
  MatPrm->SetGuidance("material name : from the material list");
  AddCmd->SetParameter(MatPrm);
  //
  //
  G4UIparameter* ColPrm = new G4UIparameter("colour",'i',false);
  ColPrm->SetGuidance("colour index : from 1 to 12 (default max)");
  ColPrm->SetParameterRange("colour>0");
  AddCmd->SetParameter(ColPrm);
  //
  G4UIparameter* ThickPrm = new G4UIparameter("thickness",'d',false);
  ThickPrm->SetGuidance("thickness of absorber");
  ThickPrm->SetParameterRange("thickness>0.");
  AddCmd->SetParameter(ThickPrm);
  //
  G4UIparameter* unitPrm = new G4UIparameter("unit",'s',false);
  unitPrm->SetGuidance("unit of thickness");
  unitPrm->SetParameterCandidates("mum mm cm m km mg/cm2 g/cm2 kg/m2");
  AddCmd->SetParameter(unitPrm);
  //
  AddCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  UseDefaultCmd = new G4UIcmdWithoutParameter("/geometry/default",this);
  UseDefaultCmd->SetGuidance("use default geometry.");
  UseDefaultCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  UpdateCmd = new G4UIcmdWithoutParameter("/geometry/update",this);
  UpdateCmd->SetGuidance("Update geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you have changed the geometry.");
  UpdateCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  // Energy deposition layers
  //  PHSDir = new G4UIdirectory("/geometry/PHS/");
  // PHSDir->SetGuidance(" Controls for defining the PHS detectors ");

  AddELayerCmd = new G4UIcmdWithAnInteger("/analysis/phs/add",this);
  AddELayerCmd->SetGuidance("Specify a layer as energy deposition detector");
  AddELayerCmd->SetGuidance("0 <= Nb <= NbofLayers  ");
  AddELayerCmd->SetParameterName("AddELayer",false);
  AddELayerCmd->SetRange("AddELayer >= 0");
  AddELayerCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  DeleteELayerCmd = new G4UIcmdWithAnInteger("/analysis/phs/delete",this);
  DeleteELayerCmd->SetGuidance("Delete a PHS detector");
  DeleteELayerCmd->SetGuidance("Nb <=  Number of layers");
  DeleteELayerCmd->SetGuidance("If Nb=0 do not perform any PHS analysis");
  DeleteELayerCmd->SetParameterName("DeleteELayer",false);
  DeleteELayerCmd->SetRange("DeleteELayer >= 0");
  DeleteELayerCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ListELayerCmd = new G4UIcmdWithoutParameter("/analysis/phs/list",this);
  ListELayerCmd->SetGuidance("List the PHS detectors");
  ListELayerCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  //  Particle fluxex layers
  //  FluxDir = new G4UIdirectory("/geometry/fluence/");
  //FluxDir->SetGuidance(" Controls for defining the fluence detectors ");

  AddFLayerCmd = new G4UIcmdWithAnInteger("/analysis/fluence/add",this);
  AddFLayerCmd->SetGuidance("Select a layer as fluence detector");
  AddFLayerCmd->SetGuidance("0 =< Nb <= NbofLayers  ");
  AddFLayerCmd->SetParameterName("AddFLayer",false);
  AddFLayerCmd->SetRange("AddFLayer >= 0");
  AddFLayerCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  DeleteFLayerCmd = new G4UIcmdWithAnInteger("/analysis/fluence/delete",this);
  DeleteFLayerCmd->SetGuidance("Delete a fluence detector");
  DeleteFLayerCmd->SetGuidance("Nb <= Num of layers");
  DeleteFLayerCmd->SetGuidance("if Nb=0 do not perform any fluence analysis");
  DeleteFLayerCmd->SetParameterName("DeleteFLayer",false);
  DeleteFLayerCmd->SetRange("DeleteFLayer >= 0");
  DeleteFLayerCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ListFLayerCmd = new G4UIcmdWithoutParameter("/analysis/fluence/list",this);
  ListFLayerCmd->SetGuidance("List the layers selected as fluence detectors");
  ListFLayerCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  
}
////////////////////////////////////////////////////////////////////////////////
//
//
MLGeometryMessenger::~MLGeometryMessenger ()
{
  delete MLdetDir;
  
  delete LayerDir;
  delete ShapeCmd;
  delete DeleteCmd;
  delete ListCmd;
  delete AddCmd;
  delete UseDefaultCmd;
  delete UpdateCmd;

  delete AddELayerCmd;
  delete DeleteELayerCmd;
  delete ListELayerCmd;

  delete AddFLayerCmd;
  delete DeleteFLayerCmd;
  delete ListFLayerCmd;
}
////////////////////////////////////////////////////////////////////////////////
//
void MLGeometryMessenger::SetNewValue (G4UIcommand* command,G4String newValue)
{    
  if (command == ShapeCmd) {
    MLGeometry->SetShape(newValue);
    
  } else if (command == DeleteCmd) {
    MLGeometry->DeleteLayer(DeleteCmd->GetNewIntValue(newValue));

  } else if (command == ListCmd) {
    MLGeometry->ListLayer(ListCmd->GetNewIntValue(newValue));

  } else if (command == AddCmd) {
    G4int num      = 0;
    G4int col      = 0;
    G4double thick = 0.;
    char mate[40];
    char unts[40];
    const char* t = newValue;
    std::istrstream is((char*)t);
    is >>num >>mate >>col >>thick >>unts;
    G4String mat=mate;
    G4String unt=unts;
    if (unt == "mg/cm2") {
      thick /= (MLGeometry->GetMLMaterial()->GetMaterial(mat)->GetDensity()
        /(mg/cm3));
      unt = "cm";
    } else if (unt == "g/cm2") {
      thick /= (MLGeometry->GetMLMaterial()->GetMaterial(mat)->GetDensity()
        /(g/cm3));
      unt = "cm";
    } else if (unt == "kg/m2") {
      thick /= (MLGeometry->GetMLMaterial()->GetMaterial(mat)->GetDensity()
        /(kg/m3));
      unt = "m";
    }
    thick *= G4UIcommand::ValueOf(unt);
    MLGeometry->AddLayer(num,mat,col,thick);

  } else if (command == UseDefaultCmd) {
    MLGeometry->SetToDefault();

  } else if (command == UpdateCmd) {
    MLGeometry->UpdateGeometry();

  } else if (command == AddELayerCmd) {
    MLGeometry->AddEdepLayer(AddELayerCmd->GetNewIntValue(newValue));

  } else if (command == DeleteELayerCmd) {
    MLGeometry->DeleteEdepLayer(DeleteELayerCmd->GetNewIntValue(newValue));

  } else if (command == ListELayerCmd) {
    MLGeometry->ListEdepLayer();

  } else if (command == AddFLayerCmd) {
    MLGeometry->AddFluxLayer(AddFLayerCmd->GetNewIntValue(newValue));

  } else if (command == DeleteFLayerCmd) {
    MLGeometry->DeleteFluxLayer(DeleteFLayerCmd->GetNewIntValue(newValue));

  } else if (command == ListFLayerCmd) {
    MLGeometry->ListFluxLayer();
  }
}
////////////////////////////////////////////////////////////////////////////////
