////////////////////////////////////////////////////////////////////////////////
//
#include "MLPhysicsListMessenger.hh"

#include "G4ios.hh"
#include "G4Tokenizer.hh"
#include <iomanip>               
#include <strstream>
////////////////////////////////////////////////////////////////////////////////
//
MLPhysicsListMessenger::MLPhysicsListMessenger(MLPhysicsList* pPhys)
  :pPhysicsList(pPhys)
{  

  PhysDir = new G4UIdirectory("/phys/");
  PhysDir->SetGuidance(" Controls of the physics processes");
 
  scenarioCmd = new G4UIcmdWithAString("/phys/scenario",this);
  scenarioCmd->SetGuidance(" select the simulation scenario ");
  scenarioCmd->SetGuidance("            em: EM interactions using standard G4EM ");
  scenarioCmd->SetGuidance("          leem: EM interactions using le G4EM ");
  scenarioCmd->SetGuidance("  hadron-em-ln: hadronic + std-ionisation ");
  scenarioCmd->SetGuidance("hadron-leem-ln: hadronic + le-ionisation");
  scenarioCmd->SetGuidance("  hadron-em+ln: hadronic + std-ionisation + le neutron");
  scenarioCmd->SetGuidance("hadron-leem+ln: hadronic + le-ionisation + le neutron");
  scenarioCmd->SetGuidance("  hadron+em-ln: hadronic + std EM interactions");
  scenarioCmd->SetGuidance("  hadron+em+ln: hadronic + std EM + le neutron interactions");
  scenarioCmd->SetGuidance("hadron+leem-ln: hadronic + le EM interactions");
  scenarioCmd->SetGuidance("hadron+leem+ln: hadronic + le EM + le neutron interactions");

  scenarioCmd->SetCandidates(
    " tra em leem hadron-em-ln hadron-leem-ln hadron-em+ln hadron-leem+ln hadron+em-ln hadron+em+ln hadron+leem-ln hadron+leem+ln mars5gev mars5gev-ln mars5gev-pc binary classic");
  scenarioCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  listscenCmd = new G4UIcmdWithoutParameter("/phys/showscen",this);
  listscenCmd->SetGuidance("Show the selected simulation scenario.");
  listscenCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  verboseCmd = new G4UIcmdWithAnInteger("/phys/verboseLevel",this);
  verboseCmd->SetGuidance("Set the verbosity level for physics list.");
  verboseCmd->AvailableForStates(G4State_PreInit,G4State_Idle);


  RegionDir = new G4UIdirectory("/phys/region/");
  RegionDir->SetGuidance("Define the regions for cuts.");

  addRegCmd = new G4UIcmdWithAString("/phys/region/add",this);
  addRegCmd->SetGuidance(" Create a new region");
  addRegCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  delRegCmd = new G4UIcmdWithAString("/phys/region/delete",this);
  delRegCmd->SetGuidance(" Delete a region");
  delRegCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  listRegCmd = new G4UIcmdWithoutParameter("/phys/region/list",this);
  listRegCmd->SetGuidance("List the regions defined");
  listRegCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  addRegLayerCmd =  new G4UIcommand("/phys/region/addalayer",this);
  addRegLayerCmd->SetGuidance(" Add a shield layer to a region.");
  addRegLayerCmd->SetGuidance("region-name layer-number");
  G4UIparameter * matName1 = new G4UIparameter("regname",'s',false);
  addRegLayerCmd->SetParameter(matName1);
  G4UIparameter * cutVal1 = new G4UIparameter("layer",'i',false);
  cutVal1->SetParameterRange("layer > 0");
  addRegLayerCmd->SetParameter(cutVal1);
  addRegLayerCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  delRegLayerCmd =  new G4UIcommand("/phys/region/delalayer",this);
  delRegLayerCmd->SetGuidance(" Remove a shield layer from a region.");
  delRegLayerCmd->SetGuidance("region-name layer-number");
  G4UIparameter * matName2 = new G4UIparameter("regname",'s',false);
  delRegLayerCmd->SetParameter(matName2);
  G4UIparameter * cutVal2 = new G4UIparameter("layer",'i',false);
  cutVal2->SetParameterRange("layer > 0");
  delRegLayerCmd->SetParameter(cutVal2);
  delRegLayerCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  listRegLayerCmd = new G4UIcmdWithAString("/phys/region/listlayers",this);
  listRegLayerCmd->SetGuidance("list the shield layers in a region");
  listRegLayerCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  CutsDir = new G4UIdirectory("/phys/cuts/");
  CutsDir->SetGuidance("Controls to set the cuts.");

  CutsGlobalDir = new G4UIdirectory("/phys/cuts/global/");
  CutsGlobalDir->SetGuidance("Set the global cuts");

  gDefaultCmd = new G4UIcmdWithADoubleAndUnit("/phys/cuts/global/default",
    this);
  gDefaultCmd->SetGuidance("Set default cut in range.");
  gDefaultCmd->SetGuidance("defaultCut unit");
  gDefaultCmd->SetParameterName("defaultCut",false);
  gDefaultCmd->SetRange("defaultCut > 0.0");
  gDefaultCmd->SetUnitCandidates("mum mm cm m");
  gDefaultCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  gGammaCmd = new G4UIcmdWithADoubleAndUnit("/phys/cuts/global/gamma", this);
  gGammaCmd->SetGuidance("Set cut for gamma in range.");
  gGammaCmd->SetGuidance("gammaCut unit");
  gGammaCmd->SetParameterName("gammaCut",false);
  gGammaCmd->SetRange("gammaCut > 0.0");
  gGammaCmd->SetUnitCandidates("mum mm cm m");
  gGammaCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  gElectronCmd = new G4UIcmdWithADoubleAndUnit("/phys/cuts/global/electron",
    this);
  gElectronCmd->SetGuidance("Set cut for electron in range.");
  gElectronCmd->SetGuidance("electronCut unit");
  gElectronCmd->SetParameterName("electronCut",false);
  gElectronCmd->SetRange("electronCut > 0.0");
  gElectronCmd->SetUnitCandidates("mum mm cm m");
  gElectronCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  gPositronCmd = new G4UIcmdWithADoubleAndUnit("/phys/cuts/global/positron",
    this);
  gPositronCmd->SetGuidance("Set cut for positron in range.");
  gPositronCmd->SetGuidance("positronCut unit");
  gPositronCmd->SetParameterName("positronCut",false);
  gPositronCmd->SetRange("positronCut > 0.0");
  gPositronCmd->SetUnitCandidates("mum mm cm m");
  gPositronCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  CutsRegionDir = new G4UIdirectory("/phys/cuts/region/");
  CutsRegionDir->SetGuidance("Set the cuts in a defined region.");

  mSetCutCmd =  new G4UIcommand("/phys/cuts/region/setcut",this);
  mSetCutCmd->SetGuidance("Set cut-in-range in a region.");
  mSetCutCmd->SetGuidance("region-name value unit");
  G4UIparameter * matName3 = new G4UIparameter("matname",'s',false);
  mSetCutCmd->SetParameter(matName3);
  G4UIparameter * cutVal3 = new G4UIparameter("cut",'d',false);
  cutVal3->SetParameterRange("cut > 0.0");
  mSetCutCmd->SetParameter(cutVal3);
  G4UIparameter * unit3 = new G4UIparameter("Unit",'s',false);
  unit3->SetParameterCandidates("mum mm cm m");
  mSetCutCmd->SetParameter(unit3);
  mSetCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  mSetByPartiCmd =  new G4UIcommand("/phys/cuts/region/cutperparticle",this);
  mSetByPartiCmd->SetGuidance(
    "Set cut-in-range in a region for a particle type.");
  mSetByPartiCmd->SetGuidance("region-name particle-type value unit");
  G4UIparameter * matName4 = new G4UIparameter("matname",'s',false);
  mSetByPartiCmd->SetParameter(matName4);
  G4UIparameter * particle = new G4UIparameter("Particle",'s',false);
  particle->SetParameterCandidates("gamma e- e+");
  mSetByPartiCmd->SetParameter(particle);
  G4UIparameter * cutVal4 = new G4UIparameter("cut",'d',false);
  cutVal4->SetParameterRange("cut > 0.0");
  mSetByPartiCmd->SetParameter(cutVal4);
  G4UIparameter * unit4 = new G4UIparameter("Unit",'s',false);
  unit4->SetParameterCandidates("mum mm cm m");
  mSetByPartiCmd->SetParameter(unit4);
  mSetByPartiCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}
////////////////////////////////////////////////////////////////////////////////
//
MLPhysicsListMessenger::~MLPhysicsListMessenger()
{
  delete PhysDir;

  delete scenarioCmd;
  delete listscenCmd;
  delete verboseCmd;

  delete RegionDir;
  delete addRegCmd;
  delete delRegCmd;
  delete listRegCmd;
  delete addRegLayerCmd;  
  delete delRegLayerCmd;  
  delete listRegLayerCmd;

  delete CutsDir;
  delete CutsGlobalDir;

  delete gDefaultCmd;
  delete gGammaCmd;
  delete gElectronCmd;
  delete gPositronCmd;

  delete CutsRegionDir;
  
  delete mSetCutCmd;
  delete mSetByPartiCmd;

}
////////////////////////////////////////////////////////////////////////////////
//
void MLPhysicsListMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{       
  if (command == scenarioCmd) {
    pPhysicsList->SetScenario(newValue) ; 

  } else if (command == listscenCmd) {
    pPhysicsList->ShowScenario() ;
  
  } else if (command == verboseCmd) {
    pPhysicsList->SetVerboseLevel(verboseCmd->GetNewIntValue(newValue)) ;

  } else if (command == addRegCmd) {
    pPhysicsList->AddARegion(newValue) ;

  } else if (command == delRegCmd) {
    pPhysicsList->DeleteARegion(newValue) ;

  } else if (command == listRegCmd) {
    pPhysicsList->ListRegions() ;
    
  } else if (command == addRegLayerCmd) {
    G4Tokenizer next(newValue);
   
    // check 1st argument
    G4String regionName = G4String(next());
    
    // check 2nd argument
    G4int layer;
    const char* t2 = G4String(next());
    std::istrstream is2((char*)t2);
    is2 >> layer;

    pPhysicsList->AddALayerToARegion(regionName,layer) ;

  } else if (command == delRegLayerCmd) {
    G4Tokenizer next(newValue);
   
    // check 1st argument
    G4String regionName = G4String(next());
    
    // check 2nd argument
    G4int layer;
    const char* t2 = G4String(next());
    std::istrstream is2((char*)t2);
    is2 >> layer;

    pPhysicsList->DeleteALayerFromARegion(regionName,layer) ;

  } else if (command == listRegLayerCmd) {
    pPhysicsList->ListTheLayersInARegion(newValue) ;

  } else if (command == gDefaultCmd) {
    pPhysicsList->SetGlobalDefault(gDefaultCmd->GetNewDoubleValue(newValue)) ;
  
  } else if (command == gGammaCmd) {
    pPhysicsList->SetGammaDefault(gGammaCmd->GetNewDoubleValue(newValue)) ;
  
  } else if (command == gElectronCmd) {
    pPhysicsList->SetElectronDefault(gElectronCmd->GetNewDoubleValue(newValue)) ;
  
  } else if (command == gPositronCmd) {
    pPhysicsList->SetPositronDefault(gPositronCmd->GetNewDoubleValue(newValue)) ;

  } else if (command == mSetCutCmd) {
    G4Tokenizer next(newValue);
   
    // check 1st argument
    G4String regionName = G4String(next());
    
    // check 2nd argument
    G4double cut;
    const char* t2 = G4String(next());
    std::istrstream is2((char*)t2);
    is2 >> cut;

    // check 3rd argument
    G4String unit = G4String(next());
    cut          *= G4UIcommand::ValueOf(unit);
    
    pPhysicsList->SetCutInRangeForRegion(cut, regionName);
  
  } else if (command == mSetByPartiCmd) {
    G4Tokenizer next(newValue);
   
    // check 1st argument
    G4String regionName = G4String(next());

    // check 2nd argument
    G4String part = G4String(next());
    
    // check 3rd argument
    G4double cut;
    const char* t2 = G4String(next());
    std::istrstream is2((char*)t2);
    is2 >> cut;

    // check 4th argument
    G4String unit = G4String(next());
    cut          *= G4UIcommand::ValueOf(unit);
    
    pPhysicsList->SetCutInRangeForRegionForParticle(cut, regionName, part);
  } 
}
////////////////////////////////////////////////////////////////////////////////
