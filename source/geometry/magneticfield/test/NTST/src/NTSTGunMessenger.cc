// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: NTSTGunMessenger.cc,v 1.2 2003-11-07 22:09:32 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "NTSTGunMessenger.hh"
#include "NTSTGunGenerator.hh"
#include "G4Geantino.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleTable.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4ios.hh"

NTSTGunMessenger::NTSTGunMessenger(NTSTGunGenerator * NTSTGun)
  : fNTSTGun(NTSTGun), gunDirectory(0), listCmd(0), particleCmd(0),
    plowCmd(0), phighCmd(0), t0Cmd(0), polCmd(0),  numberCmd(0),
    meanVertexCmd(0), rmsVertexCmd(0), 
    coslowCmd(0), coshighCmd(0), philowCmd(0), phihighCmd(0)
{
  particleTable = G4ParticleTable::GetParticleTable();

  // directory
  gunDirectory = new G4UIdirectory("/gen/gun/");
  gunDirectory->SetGuidance("NTST Particle Gun control commands.");

  // list command
  listCmd = new G4UIcmdWithoutParameter("/gen/gun/List",this);
  listCmd->SetGuidance("List available particles.");
  listCmd->SetGuidance(" Invokes G4ParticleTable.");

  // particle command
  particleCmd = new G4UIcmdWithAString("/gen/gun/particle",this);
  particleCmd->SetGuidance("Set particle to be generated.");
  particleCmd->SetGuidance(" (geantino is default)");
  particleCmd->SetParameterName("particleName",true);
  particleCmd->SetDefaultValue("geantino");
  G4String candidateList; 
  G4int nPtcl = particleTable->entries();
  for(G4int i=0;i<nPtcl;i++)
  {
    candidateList += particleTable->GetParticleName(i);
    candidateList += " ";
  }
  particleCmd->SetCandidates(candidateList);

  // meanVertex command

  meanVertexCmd = new G4UIcmdWith3VectorAndUnit("/gen/gun/vertex", this);
  meanVertexCmd->SetGuidance("Set vertex position.");
  meanVertexCmd->SetParameterName("x","y","z", true, true);
  meanVertexCmd->SetDefaultUnit("mm");
  meanVertexCmd->SetDefaultValue(G4ThreeVector(0,0,0));

  // rmsVertex command

  rmsVertexCmd = new G4UIcmdWith3VectorAndUnit("/gen/gun/spotsize", this);
  rmsVertexCmd->SetGuidance("Set vertex size.");
  rmsVertexCmd->SetParameterName("sx","sy","sz", true, true);
  rmsVertexCmd->SetDefaultUnit("mm");
  rmsVertexCmd->SetDefaultValue(G4ThreeVector(0,0,0));

  // plow command
  plowCmd = new G4UIcmdWithADoubleAndUnit("/gen/gun/plow",this);
  plowCmd->SetGuidance("Set low momentum limit.");
  plowCmd->SetParameterName("plow",true,true);
  plowCmd->SetDefaultValue(1.0);
  plowCmd->SetDefaultUnit("GeV");

  // phigh command
  phighCmd = new G4UIcmdWithADoubleAndUnit("/gen/gun/phigh",this);
  phighCmd->SetGuidance("Set high momentum limit.");
  phighCmd->SetParameterName("phigh",true,true);
  phighCmd->SetDefaultValue(1.0);
  phighCmd->SetDefaultUnit("GeV");

  // coslow command
  coslowCmd=new G4UIcmdWithADouble("/gen/gun/coslow",this);
  coslowCmd->SetGuidance("Set lower cosine limit.");
  coslowCmd->SetDefaultValue(-1.0);
  coslowCmd->SetParameterName("coslow",true,true);

  // coshigh command
  coshighCmd=new G4UIcmdWithADouble("/gen/gun/coshigh",this);
  coshighCmd->SetGuidance("Set upper cosine limit.");
  coshighCmd->SetDefaultValue(1.0);
  coshighCmd->SetParameterName("coshigh",true,true);

  // philow command
  philowCmd = new G4UIcmdWithADoubleAndUnit("/gen/gun/philow",this);
  philowCmd->SetGuidance("Set low phi limit.");
  philowCmd->SetParameterName("philow",true,true);
  philowCmd->SetDefaultValue(0.0);
  philowCmd->SetDefaultUnit("degree");

  // phihigh command
  phihighCmd = new G4UIcmdWithADoubleAndUnit("/gen/gun/phihigh",this);
  phihighCmd->SetGuidance("Set high phi limit.");
  phihighCmd->SetParameterName("phihigh",true,true);
  phihighCmd->SetDefaultValue(360.0);
  phihighCmd->SetDefaultUnit("degree");

  // t0 command
  t0Cmd = new G4UIcmdWithADoubleAndUnit("/gen/gun/t0",this);
  t0Cmd->SetGuidance("Set initial t0 of the particle.");
  t0Cmd->SetParameterName("t0",true,true);
  t0Cmd->SetDefaultValue(0.0);
  t0Cmd->SetDefaultUnit("ns");
  
  // polarization command
  polCmd = new G4UIcmdWith3Vector("/gen/gun/polarization",this);
  polCmd->SetGuidance("Set polarization.");
  polCmd->SetParameterName("Px","Py","Pz",true,true); 
  polCmd->SetRange("Px>=-1.&&Px<=1.&&Py>=-1.&&Py<=1.&&Pz>=-1.&&Pz<=1.");
  polCmd->SetDefaultValue(G4ThreeVector(0,0,0));

  // number command
  numberCmd = new G4UIcmdWithAnInteger("/gen/gun/number",this);
  numberCmd->SetGuidance("Set number of particles to be generated.");
  numberCmd->SetParameterName("N",true,true);
  numberCmd->SetRange("N>0");
  
  // set initial value to NTSTGunGenerator
  fNTSTGun->SetParticleDefinition( G4Geantino::Geantino() );
}

NTSTGunMessenger::~NTSTGunMessenger()
{
  delete listCmd;
  delete particleCmd;
  delete meanVertexCmd;
  delete rmsVertexCmd;
  delete plowCmd;
  delete phighCmd;
  delete coslowCmd;
  delete coshighCmd;
  delete philowCmd;
  delete phihighCmd;
  delete t0Cmd;
  delete polCmd;
  delete numberCmd;
  delete gunDirectory;
}

void NTSTGunMessenger::SetNewValue(G4UIcommand * command,G4String newValues)
{
  if( command==listCmd )
  { particleTable->DumpTable(); }
  else if( command==particleCmd )
  {
    G4ParticleDefinition* pd = particleTable->FindParticle(newValues);
    if(pd != NULL)
    { fNTSTGun->SetParticleDefinition( pd ); }
  }
  else if( command==meanVertexCmd )
  { fNTSTGun->SetMeanVertex(meanVertexCmd->GetNew3VectorValue(newValues)); }
  else if( command==rmsVertexCmd )
  { fNTSTGun->SetRmsVertex(rmsVertexCmd->GetNew3VectorValue(newValues)); }
  else if( command==plowCmd )
  { fNTSTGun->SetPlow(plowCmd->GetNewDoubleValue(newValues)); }
  else if( command==phighCmd )
  { fNTSTGun->SetPhigh(phighCmd->GetNewDoubleValue(newValues)); }
  else if( command==coslowCmd )
  { fNTSTGun->SetCoslow(coslowCmd->GetNewDoubleValue(newValues)); }
  else if( command==coshighCmd )
  { fNTSTGun->SetCoshigh(coshighCmd->GetNewDoubleValue(newValues)); }
  else if( command==philowCmd )
  { fNTSTGun->SetPhilow(philowCmd->GetNewDoubleValue(newValues)); }
  else if( command==phihighCmd )
  { fNTSTGun->SetPhihigh(phihighCmd->GetNewDoubleValue(newValues)); }
  else if( command==t0Cmd )
  { fNTSTGun->SetT0(t0Cmd->GetNewDoubleValue(newValues)); }
  else if( command==polCmd )
  { fNTSTGun->SetPolarization(polCmd->GetNew3VectorValue(newValues)); }
  else if( command==numberCmd )
  { fNTSTGun->SetNumberOfParticles(numberCmd->GetNewIntValue(newValues)); }
}

G4String NTSTGunMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String cv;
  
  if( command==particleCmd )
    { cv = fNTSTGun->GetParticleDefinition()->GetParticleName(); }
  else if( command==plowCmd )
    { cv = plowCmd->ConvertToString(fNTSTGun->GetPlow(),"GeV"); }
  else if( command==plowCmd )
    { cv = phighCmd->ConvertToString(fNTSTGun->GetPhigh(),"GeV"); }
  else if( command==coslowCmd )
    { cv = coslowCmd->ConvertToString(fNTSTGun->GetCoslow());}
  else if( command==coshighCmd )
    { cv = coshighCmd->ConvertToString(fNTSTGun->GetCoshigh()); }
  else if( command==philowCmd )
    { cv = philowCmd->ConvertToString(fNTSTGun->GetPhilow(),"degree"); }
  else if( command==plowCmd )
    { cv = phihighCmd->ConvertToString(fNTSTGun->GetPhihigh(),"degree"); }
  else if( command==t0Cmd )
    { cv = t0Cmd->ConvertToString(fNTSTGun->GetT0(),"ns"); }
  else if( command==polCmd )
    { cv = polCmd->ConvertToString(fNTSTGun->GetPolarization()); }
  else if( command==numberCmd )
    { cv = numberCmd->ConvertToString(fNTSTGun->GetNumberOfParticles()); }
  
  return cv;
}




