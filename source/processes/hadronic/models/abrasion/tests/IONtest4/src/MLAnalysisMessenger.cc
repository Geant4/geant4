////////////////////////////////////////////////////////////////////////////////
//
#include "MLAnalysisMessenger.hh"

#include "MLAnalysisManager.hh"
#include "MLFluenceAnalyser.hh"
#include "MLDoseAnalyser.hh"
#include "MLNIELAnalyser.hh"
#include "MLPHSAnalyser.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithADouble.hh"

#include "G4UnitsTable.hh"
////////////////////////////////////////////////////////////////////////////////
//
MLAnalysisMessenger::MLAnalysisMessenger(MLAnalysisManager* analysis)
  :analysisManager(analysis)

{ 
  MLAnalysisDir = new G4UIdirectory("/analysis/");
  MLAnalysisDir->SetGuidance("MULASSIS analysis control.");
  
  FileNameCmd = new G4UIcmdWithAString("/analysis/file",this);
  FileNameCmd->SetGuidance(
    "Input the filename prefix for CVS and report file.");
  FileNameCmd->SetParameterName("filename",true,true);
  FileNameCmd->SetDefaultValue("mulassis");
  FileNameCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  NormalCmd = new G4UIcmdWithADoubleAndUnit("/analysis/normalise",this);
  NormalCmd->SetGuidance("Set the incident particle fluence.");
  NormalCmd->SetParameterName("normal",false);
  NormalCmd->SetRange("normal > 0.");
  NormalCmd->SetUnitCandidates("cm2 m2");
  NormalCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  MLAnaFluxDir = new G4UIdirectory ("/analysis/fluence/");
  MLAnaFluxDir->
    SetGuidance("Controls for fluence analysis.");

  FluenceTypeCmd = new G4UIcmdWithAString("/analysis/fluence/type",this);
  FluenceTypeCmd->SetGuidance("Selects whether the standard (omnidirectional) fluence or planar fluence is measured");
  FluenceTypeCmd->SetGuidance("Candidates are: OMNI or PLANAR");
  FluenceTypeCmd->SetParameterName("fluenceType",false);
  FluenceTypeCmd->SetCandidates("OMNI omni OMNIDIRECTIONAL omnidirectional PLANAR planar BOUNDARY boundary");
  FluenceTypeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  FluenceUnitCmd = new G4UIcmdWithAString("/analysis/fluence/unit",this);
  FluenceUnitCmd->SetGuidance("Set the area units for fluence calculation:");
  FluenceUnitCmd->SetGuidance("Candidates are: cm2 or m2");
  FluenceUnitCmd->SetParameterName("fluenceunit",false);
  FluenceUnitCmd->SetCandidates("cm2 m2 cm2s m2s");
  FluenceUnitCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  // control of the fluence particle energy spectra
  //
  MLAnaFluEngDir = new G4UIdirectory ("/analysis/fluence/energy/");
  MLAnaFluEngDir->
    SetGuidance("Control fluence energy spectra.");

  ModeEngCmd = new G4UIcmdWithAString("/analysis/fluence/energy/mode",this);
  ModeEngCmd->SetGuidance("Set the histogram energy binning scheme.");
  ModeEngCmd->SetGuidance("Candidates are: lin, log or arb");
  ModeEngCmd->SetParameterName("modeeng",false);
  ModeEngCmd->SetCandidates(
    "lin LIN linear LINEAR log LOG logarithmic LOGARITHMIC arb ARB arbitrary ARBITRARY");
  ModeEngCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  MaxEngCmd = new G4UIcmdWithADoubleAndUnit("/analysis/fluence/energy/max",this);
  MaxEngCmd->SetGuidance("Sets the maximum energy of fluence histogram.");
  MaxEngCmd->SetParameterName("maxeng",false);
  MaxEngCmd->SetRange("maxeng > 0.");
  MaxEngCmd->SetUnitCandidates("eV keV MeV GeV TeV PeV");
  MaxEngCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  MinEngCmd = new G4UIcmdWithADoubleAndUnit(
    "/analysis/fluence/energy/min",this);
  MinEngCmd->SetGuidance("Sets the minimum energy of fluence histogram.");
  MinEngCmd->SetParameterName("mineng",false);
  MinEngCmd->SetRange("mineng >= 0.");
  MinEngCmd->SetUnitCandidates("eV keV MeV GeV TeV PeV");
  MinEngCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  NbinEngCmd = new G4UIcmdWithAnInteger("/analysis/fluence/energy/nbin",this);
  NbinEngCmd->SetGuidance("Sets the number of bins in fluence histogram.");
  NbinEngCmd->SetParameterName("nbineng",true,true);
  NbinEngCmd->SetRange("nbineng > 0");
  NbinEngCmd->SetDefaultValue(1);
  NbinEngCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  AddEngCmd = new G4UIcmdWithADoubleAndUnit(
    "/analysis/fluence/energy/add",this);
  AddEngCmd->SetGuidance("Add an energy bin edge (arbitrary bin scheme).");
  AddEngCmd->SetParameterName("addeng",true,true);
  AddEngCmd->SetRange("addeng > 0.");
  AddEngCmd->SetUnitCandidates("eV keV MeV GeV TeV PeV");
  AddEngCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  DeleteEngCmd = new G4UIcmdWithADoubleAndUnit(
    "/analysis/fluence/energy/delete",this);
  DeleteEngCmd->SetGuidance(
    "Delete an energy bin edge (arbitrary bin scheme).");
  DeleteEngCmd->SetParameterName("deleteeng",false);
  DeleteEngCmd->SetRange("deleteeng > 0.");
  DeleteEngCmd->SetUnitCandidates("eV keV MeV GeV TeV PeV");
  DeleteEngCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ListEngCmd = new G4UIcmdWithoutParameter(
    "/analysis/fluence/energy/list",this);
  ListEngCmd->SetGuidance("List the energy binning scheme.");
  ListEngCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ClearEngCmd = new G4UIcmdWithoutParameter(
    "/analysis/fluence/energy/clear",this);
  ClearEngCmd->SetGuidance("Clear the arbitrary binning scheme.");
  ClearEngCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  DefaultEngCmd = new G4UIcmdWithoutParameter(
    "/analysis/fluence/energy/default",this);
  DefaultEngCmd->SetGuidance("Use the default energy binning scheme.");
  DefaultEngCmd->AvailableForStates(G4State_PreInit,G4State_Idle);


  // control of the fluence particle angular spectra
  //
  MLAnaFluAngDir = new G4UIdirectory("/analysis/fluence/angle/");
  MLAnaFluAngDir->
    SetGuidance("Control fluence angle spectra.");

  ModeAngCmd = new G4UIcmdWithAString("/analysis/fluence/angle/mode",this);
  ModeAngCmd->SetGuidance("Set the histogram angle binning scheme.");
  ModeAngCmd->SetGuidance("Candidates are: lin, arb");
  ModeAngCmd->SetParameterName("modeang",false);
  ModeAngCmd->SetCandidates(
    "lin LIN linear LINEAR arb ARB arbitrary ARBITRARY");
  ModeAngCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  MaxAngCmd = new G4UIcmdWithADoubleAndUnit("/analysis/fluence/angle/max",this);
  MaxAngCmd->SetGuidance("Sets the maximum angle of fluence histograms.");
  MaxAngCmd->SetParameterName("maxang",false);
  MaxAngCmd->SetUnitCandidates("deg rad");
  MaxAngCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  MinAngCmd = new G4UIcmdWithADoubleAndUnit("/analysis/fluence/angle/min",this);
  MinAngCmd->SetGuidance("Sets the minimum angle of fluence histograms.");
  MinAngCmd->SetParameterName("minang",false);
  MinAngCmd->SetUnitCandidates("deg rad");
  MinAngCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  NbinAngCmd = new G4UIcmdWithAnInteger("/analysis/fluence/angle/nbin",this);
  NbinAngCmd->SetGuidance("Sets the number of bins in fluence histograms.");
  NbinAngCmd->SetParameterName("nbinang",true,true);
  NbinAngCmd->SetRange("nbinang > 0");
  NbinAngCmd->SetDefaultValue(1);
  NbinAngCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  AddAngCmd = new G4UIcmdWithADoubleAndUnit("/analysis/fluence/angle/add",this);
  AddAngCmd->SetGuidance("Add an angle bin edge (arbitrary bin scheme).");
  AddAngCmd->SetParameterName("addang",false);
  AddAngCmd->SetUnitCandidates("deg rad");
  AddAngCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  DeleteAngCmd = new G4UIcmdWithADoubleAndUnit(
    "/analysis/fluence/angle/delete",this);
  DeleteAngCmd->SetGuidance("Delete an angle bin edge (arbitrary bin scheme).");
  DeleteAngCmd->SetParameterName("deleteang",false);
  DeleteAngCmd->SetUnitCandidates("deg rad");
  DeleteAngCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  ListAngCmd = new G4UIcmdWithoutParameter("/analysis/fluence/angle/list",this);
  ListAngCmd->SetGuidance("List the angle binning scheme.");
  ListAngCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  


  ClearAngCmd = new G4UIcmdWithoutParameter(
    "/analysis/fluence/angle/clear",this);
  ClearAngCmd->SetGuidance("Clear the arbitrary binning scheme.");
  ClearAngCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  


  DefaultAngCmd = new G4UIcmdWithoutParameter(
    "/analysis/fluence/angle/default",this);
  DefaultAngCmd->SetGuidance("Use the default angle binning scheme.");
  DefaultAngCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  // particle selection
  MLAnaFluPartDir = new G4UIdirectory ("/analysis/fluence/particle/");
  MLAnaFluPartDir->SetGuidance("Controls for fluence particle selection.");

  AddPartCmd = new G4UIcmdWithAString("/analysis/fluence/particle/add",this);
  AddPartCmd->SetGuidance("Add a particle to the list.");
  AddPartCmd->SetParameterName("addpart",false);
  AddPartCmd->SetCandidates("proton neutron gamma e- muon pion");
  AddPartCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  DeletePartCmd = new G4UIcmdWithAString(
    "/analysis/fluence/particle/delete",this);
  DeletePartCmd->SetGuidance("Remove a particle from fluence analysis list.");
  DeletePartCmd->SetParameterName("deletpart",false);
  DeletePartCmd->SetCandidates("proton neutron gamma e- muon pion");
  DeletePartCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  ListPartCmd = new G4UIcmdWithoutParameter(
    "/analysis/fluence/particle/list",this);
  ListPartCmd->SetGuidance("List the particles selected for fluence analysis.");
  ListPartCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  // control of the energy pulse height spectra
  //
  MLAnaPHSDir = new G4UIdirectory ("/analysis/phs/");
  MLAnaPHSDir->
    SetGuidance("Controls energy deposition PHS analysis.");

  MLAnaPHSEngDir = new G4UIdirectory ("/analysis/phs/energy/");
  MLAnaPHSEngDir->
    SetGuidance("Controls for pulse height spectra.");

  ModePHSCmd = new G4UIcmdWithAString("/analysis/phs/energy/mode",this);
  ModePHSCmd->SetGuidance("Set the histogram binning scheme.");
  ModePHSCmd->SetGuidance("Candidates are: lin, log, arb");
  ModePHSCmd->SetParameterName("modephs",false);
  ModePHSCmd->SetCandidates(
    "lin LIN linear LINEAR log LOG logarithmic LOGARITHMIC arb ARB arbitrary ARBITRARY");
  ModePHSCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  MaxPHSCmd = new G4UIcmdWithADoubleAndUnit("/analysis/phs/energy/max",this);
  MaxPHSCmd->SetGuidance("Sets the maximum energy of PHS histograms.");
  MaxPHSCmd->SetParameterName("maxphs",false);
  MaxPHSCmd->SetRange("maxphs > 0.");
  MaxPHSCmd->SetUnitCandidates("eV keV MeV GeV TeV PeV");
  MaxPHSCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  MinPHSCmd = new G4UIcmdWithADoubleAndUnit("/analysis/phs/energy/min",this);
  MinPHSCmd->SetGuidance("Sets the minimum energy of PHS histograms.");
  MinPHSCmd->SetParameterName("minphs",false);
  MinPHSCmd->SetRange("minphs >= 0.");
  MinPHSCmd->SetUnitCandidates("eV keV MeV GeV TeV PeV");
  MinPHSCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  NbinPHSCmd = new G4UIcmdWithAnInteger("/analysis/phs/energy/nbin",this);
  NbinPHSCmd->SetGuidance("Sets the number of bins in the PHS histograms.");
  NbinPHSCmd->SetParameterName("nbinphs",true,true);
  NbinPHSCmd->SetRange("nbinphs > 0");
  NbinPHSCmd->SetDefaultValue(1);
  NbinPHSCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  AddPHSCmd = new G4UIcmdWithADoubleAndUnit("/analysis/phs/energy/add",this);
  AddPHSCmd->SetGuidance("Add an energy bin edge (arbitrary bin scheme).");
  AddPHSCmd->SetParameterName("addphs",false);
  AddPHSCmd->SetRange("addphs > 0.");
  AddPHSCmd->SetUnitCandidates("eV keV MeV GeV TeV PeV");
  AddPHSCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  DeletePHSCmd = new G4UIcmdWithADoubleAndUnit(
    "/analysis/phs/energy/delete",this);
  DeletePHSCmd->SetGuidance(
    "Delete an energy bin edge (arbitrary bin scheme).");
  DeletePHSCmd->SetParameterName("deletephs",false);
  DeletePHSCmd->SetRange("deletephs > 0.");
  DeletePHSCmd->SetUnitCandidates("eV keV MeV GeV TeV PeV");
  DeletePHSCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  ListPHSCmd = new G4UIcmdWithoutParameter("/analysis/phs/energy/list",this);
  ListPHSCmd->SetGuidance("List the binning scheme.");
  ListPHSCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  ClearPHSCmd = new G4UIcmdWithoutParameter("/analysis/phs/energy/clear",this);
  ClearPHSCmd->SetGuidance("Clear the arbitrary binning scheme.");
  ClearPHSCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  DefaultPHSCmd = new G4UIcmdWithoutParameter(
    "/analysis/phs/energy/default",this);
  DefaultPHSCmd->SetGuidance("Use the default binning scheme.");
  DefaultPHSCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  // dose analysis
  MLAnaDoseDir = new G4UIdirectory ("/analysis/dose/");
  MLAnaDoseDir->
    SetGuidance("Controls for dose analysis.");

  AddDoseCmd = new G4UIcmdWithAnInteger("/analysis/dose/add",this);
  AddDoseCmd->SetGuidance("Add a layer to the dose analysis list.");
  AddDoseCmd->SetGuidance(
			  "  The layer has to be already specified as a PHS detector");
  AddDoseCmd->SetParameterName("adddose",false);
  AddDoseCmd->SetRange("adddose > 0");
  AddDoseCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  DeleteDoseCmd = new G4UIcmdWithAnInteger("/analysis/dose/delete",this);
  DeleteDoseCmd->SetGuidance("Remove a layer from dose analysis list.");
  DeleteDoseCmd->SetParameterName("deletedose",false);
  DeleteDoseCmd->SetRange("deletedose >= 0");
  DeleteDoseCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ListDoseCmd = new G4UIcmdWithoutParameter("/analysis/dose/list",this);
  ListDoseCmd->SetGuidance("List the layers selected for dose analysis.");
  ListDoseCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  DoseUnitCmd = new G4UIcmdWithAString("/analysis/dose/unit",this);
  DoseUnitCmd->SetGuidance("Set the dose unit. Candidates are:");
  DoseUnitCmd->SetGuidance(" eV keV MeV GeV TeV PeV rad Gy");
  DoseUnitCmd->SetParameterName("doseunit",false);
  DoseUnitCmd->SetCandidates("eV keV MeV GeV TeV PeV rad Gy");
  DoseUnitCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  // NIEL analysis
  MLAnaNIELDir = new G4UIdirectory ("/analysis/niel/");
  MLAnaNIELDir->SetGuidance("Controls for niel analysis.");

  AddNIELCmd = new G4UIcmdWithAnInteger("/analysis/niel/add",this);
  AddNIELCmd->SetGuidance("Add a layer to the niel analysis list.");
  AddNIELCmd->SetGuidance(
    "  The layer has to be already specified as a fluence detector.");
  AddNIELCmd->SetParameterName("addniel",false);
  AddNIELCmd->SetRange("addniel > 0");
  AddNIELCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  DeleteNIELCmd = new G4UIcmdWithAnInteger("/analysis/niel/delete",this);
  DeleteNIELCmd->SetGuidance("Remove a layer from niel analysis list.");
  DeleteNIELCmd->SetParameterName("deleteniel",false);
  DeleteNIELCmd->SetRange("deleteniel >= 0");
  DeleteNIELCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ListNIELCmd = new G4UIcmdWithoutParameter("/analysis/niel/list",this);
  ListNIELCmd->SetGuidance("List the layers selected for niel analysis.");
  ListNIELCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  NIELUnitCmd = new G4UIcmdWithAString("/analysis/niel/unit",this);
  NIELUnitCmd->SetGuidance("Set the dose unit. Candidates are:");
  NIELUnitCmd->SetGuidance("rad Gy");
  NIELUnitCmd->SetParameterName("doseunit",false);
  NIELUnitCmd->SetCandidates("rad Gy MeV/g");
  NIELUnitCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SelectNIELCmd = new G4UIcmdWithAString("/analysis/niel/function",this);
  SelectNIELCmd->SetGuidance("Select the niel function to be used.");
  SelectNIELCmd->SetGuidance("Candidates are: cern  or jpl ");
  SelectNIELCmd->SetParameterName("nielfunc",false);
  SelectNIELCmd->SetCandidates("cern CERN jpl JPL");
  SelectNIELCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}
////////////////////////////////////////////////////////////////////////////////
//
MLAnalysisMessenger::~MLAnalysisMessenger()
{
  delete MLAnalysisDir;
  delete FileNameCmd;
  delete NormalCmd;

  delete MLAnaFluxDir;

  delete MLAnaFluEngDir;
  delete ModeEngCmd;
  delete MaxEngCmd;
  delete MinEngCmd;
  delete NbinEngCmd;
  delete AddEngCmd;
  delete DeleteEngCmd;
  delete ListEngCmd;
  delete ClearEngCmd;
  delete DefaultEngCmd;

  delete MLAnaFluAngDir;
  delete ModeAngCmd;
  delete MaxAngCmd;
  delete MinAngCmd;
  delete NbinAngCmd;
  delete AddAngCmd;
  delete DeleteAngCmd;
  delete ListAngCmd;
  delete ClearAngCmd;
  delete DefaultAngCmd;
  delete FluenceTypeCmd;
  delete FluenceUnitCmd;

  delete MLAnaFluPartDir;
  delete AddPartCmd;
  delete DeletePartCmd;
  delete ListPartCmd;


  delete MLAnaPHSDir;
  delete MLAnaPHSEngDir;
  delete ModePHSCmd;
  delete MaxPHSCmd;
  delete MinPHSCmd;
  delete NbinPHSCmd;
  delete AddPHSCmd;
  delete DeletePHSCmd;
  delete ListPHSCmd;
  delete ClearPHSCmd;
  delete DefaultPHSCmd;

  delete MLAnaDoseDir;
  delete AddDoseCmd;
  delete DeleteDoseCmd;
  delete ListDoseCmd;
  delete DoseUnitCmd;

  delete MLAnaNIELDir;
  delete AddNIELCmd;
  delete DeleteNIELCmd;
  delete ListNIELCmd;
  delete NIELUnitCmd;
  delete SelectNIELCmd;

}
////////////////////////////////////////////////////////////////////////////////
//
//
void MLAnalysisMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{

  // 1D Histograms

  if (command == FileNameCmd) {
    analysisManager->SetFilename(newValue);

  } else if (command == NormalCmd) {
    G4double theValue = NormalCmd->GetNewDoubleRawValue(newValue);
    if (newValue.contains("cm2")) {
      theValue /= cm2;
    } else if (newValue.contains("m2")) {
      theValue /= m2;
    }
    analysisManager->SetNormalFactor(theValue);
  // energy
  } else if (command == ModeEngCmd) {
    analysisManager->GetFluenceAnalyser()->SetEType(newValue);
  
  } else if (command == MaxEngCmd) {
    analysisManager->GetFluenceAnalyser()
      ->SetEngMax(MaxEngCmd->GetNewDoubleValue(newValue));

  } else if (command == MinEngCmd) {
    analysisManager->GetFluenceAnalyser()
      ->SetEngMin(MinEngCmd->GetNewDoubleValue(newValue));

  } else if (command == NbinEngCmd) {
    analysisManager->GetFluenceAnalyser()
      ->SetENBin(NbinEngCmd->GetNewIntValue(newValue));

  } else if (command == AddEngCmd) {
    analysisManager->GetFluenceAnalyser()
      ->AddEEdg(AddEngCmd->GetNewDoubleValue(newValue));

  } else if (command == DeleteEngCmd) {
    analysisManager->GetFluenceAnalyser()
      ->DeleteEEdg(DeleteEngCmd->GetNewDoubleValue(newValue));

  } else if (command == ListEngCmd) {
    analysisManager->GetFluenceAnalyser()->ListEEdg();

  } else if (command == ClearEngCmd) {
    analysisManager->GetFluenceAnalyser()->ClearEEdg();

  } else if (command == DefaultEngCmd) {
    analysisManager->GetFluenceAnalyser()->DefaultEEdg();
  // angle
  } else if (command == ModeAngCmd) {
    analysisManager->GetFluenceAnalyser()->SetAType(newValue);

  } else if (command == MaxAngCmd) {
    analysisManager->GetFluenceAnalyser()
      ->SetAngMax(MaxAngCmd->GetNewDoubleValue(newValue));

  } else if (command == MinAngCmd) {
    analysisManager->GetFluenceAnalyser()
      ->SetAngMin(MinAngCmd->GetNewDoubleValue(newValue));

  } else if (command == NbinAngCmd) {
    analysisManager->GetFluenceAnalyser()
      ->SetANBin(NbinAngCmd->GetNewIntValue(newValue));

  } else if (command == AddAngCmd) {
    analysisManager->GetFluenceAnalyser()
      ->AddAEdg(AddAngCmd->GetNewDoubleValue(newValue));

  } else if (command == DeleteAngCmd) {
    analysisManager->GetFluenceAnalyser()
      ->DeleteAEdg(DeleteAngCmd->GetNewDoubleValue(newValue));

  } else if (command == ListAngCmd) {
    analysisManager->GetFluenceAnalyser()->ListAEdg();

  } else if (command == ClearAngCmd) {
    analysisManager->GetFluenceAnalyser()->ClearAEdg();

  } else if (command == DefaultAngCmd) {
    analysisManager->GetFluenceAnalyser()->DefaultAEdg();

  } else if (command == FluenceTypeCmd) {
    analysisManager->GetFluenceAnalyser()->SetFluenceType(newValue);

  } else if (command == FluenceUnitCmd) {
    analysisManager->GetFluenceAnalyser()->SetFluenceUnit(newValue);

  //particle type
  } else if (command == AddPartCmd) {
    analysisManager->GetFluenceAnalyser()->AddSPart(newValue);

  } else if (command == DeletePartCmd) {
    analysisManager->GetFluenceAnalyser()->DeleteSPart(newValue);

  } else if (command == ListPartCmd) {
    analysisManager->GetFluenceAnalyser()->ListSPart();

  // phs
  } else if (command == ModePHSCmd) {
    analysisManager->GetPHSAnalyser()->SetPType(newValue);

  } else if (command == MaxPHSCmd) {
    analysisManager->GetPHSAnalyser()
      ->SetPHSMax(MaxPHSCmd->GetNewDoubleValue(newValue));

  } else if (command == MinPHSCmd) {
    analysisManager->GetPHSAnalyser()
      ->SetPHSMin(MinPHSCmd->GetNewDoubleValue(newValue));

  } else if (command == NbinPHSCmd) {
    analysisManager->GetPHSAnalyser()
      ->SetPNBin(NbinPHSCmd->GetNewIntValue(newValue));

  } else if (command == AddPHSCmd) {
    analysisManager->GetPHSAnalyser()
      ->AddPEdg(AddPHSCmd->GetNewDoubleValue(newValue));

  } else if (command == DeletePHSCmd) {
    analysisManager->GetPHSAnalyser()
      ->DeletePEdg(DeletePHSCmd->GetNewDoubleValue(newValue));

  } else if (command == ListPHSCmd) {
    analysisManager->GetPHSAnalyser()->ListPEdg();

  } else if (command == ClearPHSCmd) {
    analysisManager->GetPHSAnalyser()->ClearPEdg();

  } else if (command == DefaultPHSCmd) {
    analysisManager->GetPHSAnalyser()->DefaultPEdg();

  // dose
  } else if (command == AddDoseCmd) {
    analysisManager->GetDoseAnalyser()
      ->AddDoseLayer(AddDoseCmd->GetNewIntValue(newValue));

  } else if (command == DeleteDoseCmd) {
    analysisManager->GetDoseAnalyser()
      ->DeleteDoseLayer(DeleteDoseCmd->GetNewIntValue(newValue));

  } else if (command == ListDoseCmd) {
    analysisManager->GetDoseAnalyser()->ListDoseLayer();

  } else if (command == DoseUnitCmd) {
    analysisManager->GetDoseAnalyser()->SetDoseUnit(newValue);

  // Niel
  } else if (command == AddNIELCmd) {
    analysisManager->GetNIELAnalyser()
      ->AddNielLayer(AddNIELCmd->GetNewIntValue(newValue));

  } else if (command == DeleteNIELCmd) {
    analysisManager->GetNIELAnalyser()
      ->DeleteNielLayer(DeleteNIELCmd->GetNewIntValue(newValue));

  } else if (command == ListNIELCmd) {
    analysisManager->GetNIELAnalyser()->ListNielLayer();

  } else if (command == NIELUnitCmd) {
    analysisManager->GetNIELAnalyser()->SetDoseUnit(newValue);

  } else if (command == SelectNIELCmd) {
    analysisManager->GetNIELAnalyser()->SelectNielFunction(newValue);

  }
}
////////////////////////////////////////////////////////////////////////////////
