
//
#include "G4RadioactiveDecaymessenger.hh"

#include "g4std/iostream"
////////////////////////////////////////////////////////////////////////////////
//
G4RadioactiveDecaymessenger::G4RadioactiveDecaymessenger 
(G4RadioactiveDecay* theRadioactiveDecayContainer1)
:theRadioactiveDecayContainer(theRadioactiveDecayContainer1)
{
  //
  //
  // main directory for control of the RDM
  //
  //
  grdmDirectory = new G4UIdirectory("/grdm/");
  grdmDirectory->SetGuidance("Controls for the Radioactive Decay Module.");
  //
  //
  // Command to define the limits on nucleus the RDM will treat.
  //
  nucleuslimitsCmd = new
    G4UIcmdWithNucleusLimits("/grdm/nucleusLimits",this);
  nucleuslimitsCmd->SetGuidance 
    ("Set the amotic weight and number limits for the RDM.");
  nucleuslimitsCmd->SetParameterName("aMin","aMax","zMin","zMax",true);
  //

  //
  // The next command contols whether the decay will be treated analoguely or 
  // with variance reduction
  //
  analoguemcCmd = new G4UIcmdWithABool ("/grdm/analogueMC",this);
  analoguemcCmd->SetGuidance("false: variance reduction method; true: analogue method");
  analoguemcCmd->SetParameterName("AnalogueMC",true);
  analoguemcCmd->SetDefaultValue(true);

  //
  //
  // Command to selete a logical volume for RDM.
  //
  avolumeCmd = new
    G4UIcmdWithAString("/grdm/selectVolume",this);
  avolumeCmd->SetGuidance 
    ("Suppply a logical volumes name to add it to the RDM apply list");
  avolumeCmd->SetParameterName("aVolume",false);
  //
  //
  //
  // Command to de-selete a logical volume for RDM.
  //
  deavolumeCmd = new
    G4UIcmdWithAString("/grdm/deselectVolume",this);
  deavolumeCmd->SetGuidance 
    ("Suppply a logical volumes name to remove it from the RDM apply list");
  deavolumeCmd->SetParameterName("aVolume",false);
  //
  //
  // Command to selete all logical volumes for RDM.
  //
  allvolumesCmd = new
    G4UIcmdWithoutParameter("/grdm/allVolumes",this);
  allvolumesCmd->SetGuidance 
    (" apply RDM to all logical volumes. No parameter required.");
  //  allvolumeCmd->SetParameterName("AddAVolume",true);

  //
  // Command to de-selete a logical volume for RDM.
  //
  deallvolumesCmd = new
    G4UIcmdWithoutParameter("/grdm/noVolumes",this);
  deallvolumesCmd->SetGuidance 
    (" RDM is not applied to any logical volumes");

  //  deallvolumesCmd->SetParameterName("RemoveAVolume",true);
  //
  // The next command contols whether the branching ratio biasing will be applied or not
  //
  brbiasCmd = new G4UIcmdWithABool ("/grdm/BRbias",this);
  brbiasCmd->SetGuidance("false: no biasing; true: all branches are treated as equal");
  brbiasCmd->SetParameterName("BRBias",true);
  brbiasCmd->SetDefaultValue(true);

  //
  // Command to define the incident particle source time profile.
  //
  sourcetimeprofileCmd = new
    G4UIcmdWithAString("/grdm/sourceTimeProfile",this);
  sourcetimeprofileCmd->SetGuidance 
    ("Supply the name of the ascii file containing the source particle time profile");
  sourcetimeprofileCmd->SetParameterName("STimeProfile",true);
  sourcetimeprofileCmd->SetDefaultValue("source.data");
  //
  //
  // Command to define the incident particle source time profile.
  //
  decaybiasprofileCmd = new
    G4UIcmdWithAString("/grdm/decayBiasProfile",this);
  decaybiasprofileCmd->SetGuidance 
    ("Supply the name of the ascii file containing the decay bias time profile");
  decaybiasprofileCmd->SetParameterName("DBiasProfile",true);
  decaybiasprofileCmd->SetDefaultValue("bias.data");
  //

  //
  // This command setup the nuclei spliting parameter
  //
  splitnucleiCmd = new G4UIcmdWithAnInteger("/grdm/splitNuclei",this);
  splitnucleiCmd->SetGuidance("Set number of spliting for the isotopes.");
  splitnucleiCmd->SetParameterName("NSplit",true);
  splitnucleiCmd->SetDefaultValue(1);
  splitnucleiCmd->SetRange("NSplit>=1");

  //
  // This command setup the verbose level of radioactive decay
  //
  verboseCmd = new G4UIcmdWithAnInteger("/grdm/verbose",this);
  verboseCmd->SetGuidance("Set verbose level: 0, 1, 2 or 3");
  verboseCmd->SetParameterName("VerboseLevel",true);
  verboseCmd->SetDefaultValue(1);
  verboseCmd->SetRange("VerboseLevel>=0");

}
////////////////////////////////////////////////////////////////////////////////
//
G4RadioactiveDecaymessenger::~G4RadioactiveDecaymessenger ()
{
  delete grdmDirectory;
  delete nucleuslimitsCmd;
  delete sourcetimeprofileCmd;
  delete decaybiasprofileCmd;
  delete analoguemcCmd;
  delete brbiasCmd;
  delete splitnucleiCmd;
  delete verboseCmd;
  delete avolumeCmd;
  delete deavolumeCmd;
  delete allvolumesCmd;
  delete deallvolumesCmd;
}
////////////////////////////////////////////////////////////////////////////////
//
void G4RadioactiveDecaymessenger::SetNewValue (G4UIcommand *command, G4String newValues)
{
  if (command==nucleuslimitsCmd) {theRadioactiveDecayContainer->
				    SetNucleusLimits(nucleuslimitsCmd->GetNewNucleusLimitsValue(newValues));}
  else if  (command==analoguemcCmd) {
    G4int vl;
    const char* t = newValues;
    G4std::istrstream is((char*)t);
    is >> vl;
    theRadioactiveDecayContainer->SetAnalogueMonteCarlo(vl!=0);}
  else if  (command==avolumeCmd) {theRadioactiveDecayContainer->
				   SelectAVolume(newValues);}
  else if  (command==deavolumeCmd) {theRadioactiveDecayContainer->
				   DeselectAVolume(newValues);}
  else if  (command==allvolumesCmd) {theRadioactiveDecayContainer->
				   SelectAllVolumes();}
  else if  (command==deallvolumesCmd) {theRadioactiveDecayContainer->
				   DeselectAllVolumes();}
  else if  (command==brbiasCmd) {
    G4int vl;
    const char* t = newValues;
    G4std::istrstream is((char*)t);
    is >> vl;
    theRadioactiveDecayContainer->SetBRBias(vl!=0);}
  else if (command==sourcetimeprofileCmd) {theRadioactiveDecayContainer->
					     SetSourceTimeProfile(newValues);}
  else if (command==decaybiasprofileCmd) {theRadioactiveDecayContainer->
					    SetDecayBias(newValues);}
  else if (command==splitnucleiCmd) {theRadioactiveDecayContainer->
				       SetSplitNuclei(splitnucleiCmd->GetNewIntValue(newValues));}
  else if (command==verboseCmd) {theRadioactiveDecayContainer->
				       SetVerboseLevel(verboseCmd->GetNewIntValue(newValues));}
}






