///////////////////////////////////////////////////////////////////////////////
// 
// MODULE:           G4GeneralParticleSourceMessenger.cc
//
// Version:          1.
// Date:             28/02/00
// Author:           C Ferguson, F Lei and P Truscott
// Organisation:     University of Southampton / DERA
// Customer:         ESA/ESTEC
//
///////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//
// Version 1.0, 28 February 2000, C Ferguson, Created.
//
// Version 1.1, 18 October 2000, Modified to inherit from G4VPrimaryGenerator.
// New name at the request of M. Asai.
//
///////////////////////////////////////////////////////////////////////////////
//
#include "G4Geantino.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleTable.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithABool.hh"
#include "G4ios.hh"

#include <iostream.h>

#include "g4std/iomanip"               
#include "g4rw/ctoken.h"               
#include "g4rw/rstream.h"               
#include "g4std/strstream"

#include "G4GeneralParticleSourceMessenger.hh"
#include "G4GeneralParticleSource.hh"
///////////////////////////////////////////////////////////////////////////////
//
G4GeneralParticleSourceMessenger::G4GeneralParticleSourceMessenger
  (G4GeneralParticleSource *fPtclGun) 
    : fParticleGun(fPtclGun),fShootIon(false)
{
  particleTable = G4ParticleTable::GetParticleTable();
  histtype = "biasx";

  gpsDirectory = new G4UIdirectory("/gps/");
  gpsDirectory->SetGuidance("General Paricle Source control commands.");

  typeCmd = new G4UIcmdWithAString("/gps/type",this);
  typeCmd->SetGuidance("Sets source distribution type.");
  typeCmd->SetGuidance("Either Point, Plane, Surface or Volume");
  typeCmd->SetParameterName("DisType",true,true);
  typeCmd->SetDefaultValue("Point");
  typeCmd->SetCandidates("Point Plane Surface Volume");

  shapeCmd = new G4UIcmdWithAString("/gps/shape",this);
  shapeCmd->SetGuidance("Sets source shape type.");
  shapeCmd->SetParameterName("Shape",true,true);
  shapeCmd->SetDefaultValue("NULL");
  shapeCmd->SetCandidates("Circle Annulus Ellipse Square Rectangle Sphere Ellipsoid Cylinder Para");

  centreCmd = new G4UIcmdWith3VectorAndUnit("/gps/centre",this);
  centreCmd->SetGuidance("Set centre coordinates of source.");
  centreCmd->SetParameterName("X","Y","Z",true,true);
  centreCmd->SetDefaultUnit("cm");
  centreCmd->SetUnitCandidates("micron mm cm m km");

  posrot1Cmd = new G4UIcmdWith3Vector("/gps/posrot1",this);
  posrot1Cmd->SetGuidance("Set rotation matrix of x'.");
  posrot1Cmd->SetGuidance("Posrot1 does not need to be a unit vector.");
  posrot1Cmd->SetParameterName("R1x","R1y","R1z",true,true); 
  posrot1Cmd->SetRange("R1x != 0 || R1y != 0 || R1z != 0");

  posrot2Cmd = new G4UIcmdWith3Vector("/gps/posrot2",this);
  posrot2Cmd->SetGuidance("Set rotation matrix of y'.");
  posrot2Cmd->SetGuidance("Posrot2 does not need to be a unit vector.");
  posrot2Cmd->SetParameterName("R2x","R2y","R2z",true,true); 
  posrot2Cmd->SetRange("R2x != 0 || R2y != 0 || R2z != 0");

  halfxCmd = new G4UIcmdWithADoubleAndUnit("/gps/halfx",this);
  halfxCmd->SetGuidance("Set x half length of source.");
  halfxCmd->SetParameterName("Halfx",true,true);
  halfxCmd->SetDefaultUnit("cm");
  halfxCmd->SetUnitCandidates("micron mm cm m km");

  halfyCmd = new G4UIcmdWithADoubleAndUnit("/gps/halfy",this);
  halfyCmd->SetGuidance("Set y half length of source.");
  halfyCmd->SetParameterName("Halfy",true,true);
  halfyCmd->SetDefaultUnit("cm");
  halfyCmd->SetUnitCandidates("micron mm cm m km");

  halfzCmd = new G4UIcmdWithADoubleAndUnit("/gps/halfz",this);
  halfzCmd->SetGuidance("Set z half length of source.");
  halfzCmd->SetParameterName("Halfz",true,true);
  halfzCmd->SetDefaultUnit("cm");
  halfzCmd->SetUnitCandidates("micron mm cm m km");

  radiusCmd = new G4UIcmdWithADoubleAndUnit("/gps/radius",this);
  radiusCmd->SetGuidance("Set radius of source.");
  radiusCmd->SetParameterName("Radius",true,true);
  radiusCmd->SetDefaultUnit("cm");
  radiusCmd->SetUnitCandidates("micron mm cm m km");

  radius0Cmd = new G4UIcmdWithADoubleAndUnit("/gps/radius0",this);
  radius0Cmd->SetGuidance("Set inner radius of source.");
  radius0Cmd->SetParameterName("Radius0",true,true);
  radius0Cmd->SetDefaultUnit("cm");
  radius0Cmd->SetUnitCandidates("micron mm cm m km");

  paralpCmd = new G4UIcmdWithADoubleAndUnit("/gps/paralp",this);
  paralpCmd->SetGuidance("Angle from y-axis of y' in Para");
  paralpCmd->SetParameterName("paralp",true,true);
  paralpCmd->SetDefaultUnit("rad");
  paralpCmd->SetUnitCandidates("rad deg");

  partheCmd = new G4UIcmdWithADoubleAndUnit("/gps/parthe",this);
  partheCmd->SetGuidance("Polar angle through centres of z faces");
  partheCmd->SetParameterName("parthe",true,true);
  partheCmd->SetDefaultUnit("rad");
  partheCmd->SetUnitCandidates("rad deg");

  parphiCmd = new G4UIcmdWithADoubleAndUnit("/gps/parphi",this);
  parphiCmd->SetGuidance("Azimuth angle through centres of z faces");
  parphiCmd->SetParameterName("parphi",true,true);
  parphiCmd->SetDefaultUnit("rad");
  parphiCmd->SetUnitCandidates("rad deg");

  confineCmd = new G4UIcmdWithAString("/gps/confine",this);
  confineCmd->SetGuidance("Confine source to volume (NULL to unset).");
  confineCmd->SetGuidance("usage: confine VolName");
  confineCmd->SetParameterName("VolName",true,true);
  confineCmd->SetDefaultValue("NULL");

  // Angular distribution commands
  angtypeCmd = new G4UIcmdWithAString("/gps/angtype",this);
  angtypeCmd->SetGuidance("Sets angular source distribution type");
  angtypeCmd->SetGuidance("Possible variables are: iso, cos or user");
  angtypeCmd->SetParameterName("AngDis",true,true);
  angtypeCmd->SetDefaultValue("iso");
  angtypeCmd->SetCandidates("iso cos user");

  angrot1Cmd = new G4UIcmdWith3Vector("/gps/angrot1",this);
  angrot1Cmd->SetGuidance("Sets the x' vector for angular distribution");
  angrot1Cmd->SetGuidance("Need not be a unit vector");
  angrot1Cmd->SetParameterName("AR1x","AR1y","AR1z",true,true);
  angrot1Cmd->SetRange("AR1x != 0 || AR1y != 0 || AR1z != 0");

  angrot2Cmd = new G4UIcmdWith3Vector("/gps/angrot2",this);
  angrot2Cmd->SetGuidance("Sets the y' vector for angular distribution");
  angrot2Cmd->SetGuidance("Need not be a unit vector");
  angrot2Cmd->SetParameterName("AR2x","AR2y","AR2z",true,true);
  angrot2Cmd->SetRange("AR2x != 0 || AR2y != 0 || AR2z != 0");

  minthetaCmd = new G4UIcmdWithADoubleAndUnit("/gps/mintheta",this);
  minthetaCmd->SetGuidance("Set minimum theta");
  minthetaCmd->SetParameterName("MinTheta",true,true);
  minthetaCmd->SetDefaultUnit("rad");
  minthetaCmd->SetUnitCandidates("rad deg");

  maxthetaCmd = new G4UIcmdWithADoubleAndUnit("/gps/maxtheta",this);
  maxthetaCmd->SetGuidance("Set maximum theta");
  maxthetaCmd->SetParameterName("MaxTheta",true,true);
  maxthetaCmd->SetDefaultValue(3.1416);
  maxthetaCmd->SetDefaultUnit("rad");
  maxthetaCmd->SetUnitCandidates("rad deg");

  minphiCmd = new G4UIcmdWithADoubleAndUnit("/gps/minphi",this);
  minphiCmd->SetGuidance("Set minimum phi");
  minphiCmd->SetParameterName("MinPhi",true,true);
  minphiCmd->SetDefaultUnit("rad");
  minphiCmd->SetUnitCandidates("rad deg");

  maxphiCmd = new G4UIcmdWithADoubleAndUnit("/gps/maxphi",this);
  maxphiCmd->SetGuidance("Set maximum phi");
  maxphiCmd->SetParameterName("MaxPhi",true,true);
  maxphiCmd->SetDefaultUnit("rad");
  maxphiCmd->SetUnitCandidates("rad deg");

  surfnormCmd = new G4UIcmdWithABool("/gps/surfnorm",this);
  surfnormCmd->SetGuidance("Makes a user-defined distribution with respect to surface normals rather than x,y,z axes.");
  surfnormCmd->SetGuidance("Default is false");
  surfnormCmd->SetParameterName("surfnorm",true);
  surfnormCmd->SetDefaultValue(false);

  // Energy
  energytypeCmd = new G4UIcmdWithAString("/gps/energytype",this);
  energytypeCmd->SetGuidance("Sets energy distribution type");
  energytypeCmd->SetParameterName("EnergyDis",true,true);
  energytypeCmd->SetDefaultValue("Mono");
  energytypeCmd->SetCandidates("Mono Lin Pow Exp Brem Bbody Cdg User Arb Epn");

  eminCmd = new G4UIcmdWithADoubleAndUnit("/gps/emin",this);
  eminCmd->SetGuidance("Sets Emin");
  eminCmd->SetParameterName("emin",true,true);
  eminCmd->SetDefaultUnit("keV");
  eminCmd->SetUnitCandidates("eV keV MeV GeV TeV PeV");

  emaxCmd = new G4UIcmdWithADoubleAndUnit("/gps/emax",this);
  emaxCmd->SetGuidance("Sets Emax");
  emaxCmd->SetParameterName("emax",true,true);
  emaxCmd->SetDefaultUnit("keV");
  emaxCmd->SetUnitCandidates("eV keV MeV GeV TeV PeV");

  monoenergyCmd = new G4UIcmdWithADoubleAndUnit("/gps/monoenergy",this);
  monoenergyCmd->SetGuidance("Sets Monoenergy");
  monoenergyCmd->SetParameterName("monoenergy",true,true);
  monoenergyCmd->SetDefaultUnit("keV");
  monoenergyCmd->SetUnitCandidates("eV keV MeV GeV TeV PeV");

  alphaCmd = new G4UIcmdWithADouble("/gps/alpha",this);
  alphaCmd->SetGuidance("Sets Alpha");
  alphaCmd->SetParameterName("alpha",true,true);
  
  tempCmd = new G4UIcmdWithADouble("/gps/temp",this);
  tempCmd->SetGuidance("Sets the temperature for Brem and BBody");
  tempCmd->SetParameterName("temp",true,true);

  ezeroCmd = new G4UIcmdWithADouble("/gps/ezero",this);
  ezeroCmd->SetGuidance("Sets ezero exponential distributions");
  ezeroCmd->SetParameterName("ezero",true,true);

  gradientCmd = new G4UIcmdWithADouble("/gps/gradient",this);
  gradientCmd->SetGuidance("Sets the gradient for Lin distributions");
  gradientCmd->SetParameterName("gradient",true,true);

  interceptCmd = new G4UIcmdWithADouble("/gps/intercept",this);
  interceptCmd->SetGuidance("Sets the intercept for Lin distributions");
  interceptCmd->SetParameterName("intercept",true,true);

  calculateCmd = new G4UIcmdWithoutParameter("/gps/calculate",this);
  calculateCmd->SetGuidance("Calculates distributions for Cdg and BBody");

  energyspecCmd = new G4UIcmdWithABool("/gps/energyspec",this);
  energyspecCmd->SetGuidance("Toggles between energy and momentum spectra");
  energyspecCmd->SetParameterName("energyspec",true);
  energyspecCmd->SetDefaultValue(true);

  diffspecCmd = new G4UIcmdWithABool("/gps/diffspec",this);
  diffspecCmd->SetGuidance("Toggles between differential and integral spectra");
  diffspecCmd->SetParameterName("diffspec",true,true);
  diffspecCmd->SetDefaultValue(true);

  // Biasing + histograms in general
  histnameCmd = new G4UIcmdWithAString("/gps/histname",this);
  histnameCmd->SetGuidance("Sets histogram type");
  histnameCmd->SetParameterName("HistType",true,true);
  histnameCmd->SetDefaultValue("biasx");
  histnameCmd->SetCandidates("biasx biasy biasz biast biasp biase theta phi energy arb epn");

  histpointCmd = new G4UIcmdWith3Vector("/gps/histpoint",this);
  histpointCmd->SetGuidance("Allows user to define a histogram");
  histpointCmd->SetGuidance("Enter: Ehi Weight");
  histpointCmd->SetParameterName("Ehi","Weight","Junk",true,true);
  histpointCmd->SetRange("Ehi >= 0. && Weight >= 0.");

  arbintCmd = new G4UIcmdWithAString("/gps/arbint",this);
  arbintCmd->SetGuidance("Sets Arbitrary Interpolation type.");
  arbintCmd->SetParameterName("int",true,true);
  arbintCmd->SetDefaultValue("NULL");
  arbintCmd->SetCandidates("Lin Log Exp Spline");

  // verbosity
  verbosityCmd = new G4UIcmdWithAnInteger("/gps/verbose",this);
  verbosityCmd->SetGuidance("Set Verbose level for GPS");
  verbosityCmd->SetGuidance(" 0 : Silent");
  verbosityCmd->SetGuidance(" 1 : Limited information");
  verbosityCmd->SetGuidance(" 2 : Detailed information");
  verbosityCmd->SetParameterName("level",false);
  verbosityCmd->SetRange("level>=0 && level <=2");

  particleCmd = new G4UIcmdWithAString("/gps/particle",this);
  particleCmd->SetGuidance("Set particle to be generated.");
  particleCmd->SetGuidance(" (geantino is default)");
  particleCmd->SetGuidance(" (ion can be specified for shooting ions)");
  particleCmd->SetParameterName("particleName",true);
  particleCmd->SetDefaultValue("geantino");
  G4String candidateList; 
  G4int nPtcl = particleTable->entries();
  for(G4int i=0;i<nPtcl;i++)
  {
    candidateList += particleTable->GetParticleName(i);
    candidateList += " ";
  }
  candidateList += "ion ";
  particleCmd->SetCandidates(candidateList);

  // SR1.3
  //  ionCmd = new UIcmdWithNucleusAndUnit("/gps/ion",this);
  //ionCmd->SetGuidance("define the primary ion (a,z,e)");
  //ionCmd->SetParameterName("A","Z","E",true);
  //ionCmd->SetDefaultUnit("keV");
  //ionCmd->SetUnitCandidates("keV MeV");

  ionCmd = new G4UIcommand("/gps/ion",this);
  ionCmd->SetGuidance("Set properties of ion to be generated.");
  ionCmd->SetGuidance("[usage] /gun/ion Z A Q E");
  ionCmd->SetGuidance("        Z:(int) AtomicNumber");
  ionCmd->SetGuidance("        A:(int) AtomicMass");
  ionCmd->SetGuidance("        Q:(int) Charge of Ion (in unit of e)");
  ionCmd->SetGuidance("        E:(double) Excitation energy (in keV)");
  
  G4UIparameter* param;
  param = new G4UIparameter("Z",'i',false);
  param->SetDefaultValue("1");
  ionCmd->SetParameter(param);
  param = new G4UIparameter("A",'i',false);
  param->SetDefaultValue("1");
  ionCmd->SetParameter(param);
  param = new G4UIparameter("Q",'i',true);
  param->SetDefaultValue("0");
  ionCmd->SetParameter(param);
  param = new G4UIparameter("E",'d',true);
  param->SetDefaultValue("0.0");
  ionCmd->SetParameter(param);


  timeCmd = new G4UIcmdWithADoubleAndUnit("/gps/time",this);
  timeCmd->SetGuidance("Set initial time of the particle.");
  timeCmd->SetParameterName("t0",true,true);
  timeCmd->SetDefaultUnit("ns");
  //timeCmd->SetUnitCategory("Time");
  //timeCmd->SetUnitCandidates("ns ms s");
  
  polCmd = new G4UIcmdWith3Vector("/gps/polarization",this);
  polCmd->SetGuidance("Set polarization.");
  polCmd->SetParameterName("Px","Py","Pz",true,true); 
  polCmd->SetRange("Px>=-1.&&Px<=1.&&Py>=-1.&&Py<=1.&&Pz>=-1.&&Pz<=1.");

  numberCmd = new G4UIcmdWithAnInteger("/gps/number",this);
  numberCmd->SetGuidance("Set number of particles to be generated.");
  numberCmd->SetParameterName("N",true,true);
  numberCmd->SetRange("N>0");
}

G4GeneralParticleSourceMessenger::~G4GeneralParticleSourceMessenger()
{
  delete typeCmd;
  delete shapeCmd;
  delete centreCmd;
  delete posrot1Cmd;
  delete posrot2Cmd;
  delete halfxCmd;
  delete halfyCmd;
  delete halfzCmd;
  delete radiusCmd;
  delete radius0Cmd;
  delete paralpCmd;
  delete partheCmd;
  delete parphiCmd;
  delete confineCmd;

  delete angtypeCmd;
  delete angrot1Cmd;
  delete angrot2Cmd;
  delete minthetaCmd;
  delete maxthetaCmd;
  delete minphiCmd;
  delete maxphiCmd;
  delete surfnormCmd;

  delete energytypeCmd;
  delete eminCmd;
  delete emaxCmd;
  delete monoenergyCmd;
  delete alphaCmd;
  delete tempCmd;
  delete ezeroCmd;
  delete gradientCmd;
  delete interceptCmd;
  delete calculateCmd;
  delete energyspecCmd;
  delete diffspecCmd;

  delete histnameCmd;
  delete histpointCmd;
  delete arbintCmd;

  delete verbosityCmd;
  delete ionCmd;
  delete particleCmd;
  delete timeCmd;
  delete polCmd;
  delete numberCmd;

  delete gpsDirectory;
}

void G4GeneralParticleSourceMessenger::SetNewValue(G4UIcommand *command, G4String newValues)
{
  if(command == typeCmd)
    {
      fParticleGun->SetPosDisType(newValues);
    }
  else if(command == shapeCmd)
    {
      fParticleGun->SetPosDisShape(newValues);
    }
  else if(command == centreCmd)
    {
      fParticleGun->SetCentreCoords(centreCmd->GetNew3VectorValue(newValues));
    }
  else if(command == posrot1Cmd)
    {
      fParticleGun->SetPosRot1(posrot1Cmd->GetNew3VectorValue(newValues));
    }
  else if(command == posrot2Cmd)
    {
      fParticleGun->SetPosRot2(posrot2Cmd->GetNew3VectorValue(newValues));
    }
  else if(command == halfxCmd)
    {
      fParticleGun->SetHalfX(halfxCmd->GetNewDoubleValue(newValues));
    }
  else if(command == halfyCmd)
    {
      fParticleGun->SetHalfY(halfyCmd->GetNewDoubleValue(newValues));
    }
  else if(command == halfzCmd)
    {
      fParticleGun->SetHalfZ(halfzCmd->GetNewDoubleValue(newValues));
    }
  else if(command == radiusCmd)
    {
      fParticleGun->SetRadius(radiusCmd->GetNewDoubleValue(newValues));
    }
  else if(command == radius0Cmd)
    {
      fParticleGun->SetRadius0(radius0Cmd->GetNewDoubleValue(newValues));
    }
  else if(command == paralpCmd)
    {
      fParticleGun->SetParAlpha(paralpCmd->GetNewDoubleValue(newValues));
    }
  else if(command == partheCmd)
    {
      fParticleGun->SetParTheta(partheCmd->GetNewDoubleValue(newValues));
    }
  else if(command == parphiCmd)
    {
      fParticleGun->SetParPhi(parphiCmd->GetNewDoubleValue(newValues));
    }
  else if(command == confineCmd)
    {
      fParticleGun->ConfineSourceToVolume(newValues);
    }
  else if(command == angtypeCmd)
    {
      fParticleGun->SetAngDistType(newValues);
    }
  else if(command == angrot1Cmd)
    {
      G4String a = "angref1";
      fParticleGun->DefineAngRefAxes(a,angrot1Cmd->GetNew3VectorValue(newValues));
    }
  else if(command == angrot2Cmd)
    {
      G4String a = "angref2";
      fParticleGun->DefineAngRefAxes(a,angrot2Cmd->GetNew3VectorValue(newValues));
    }
  else if(command == minthetaCmd)
    {
      fParticleGun->SetMinTheta(minthetaCmd->GetNewDoubleValue(newValues));
    }
  else if(command == minphiCmd)
    {
      fParticleGun->SetMinPhi(minphiCmd->GetNewDoubleValue(newValues));
    }
  else if(command == maxthetaCmd)
    {
      fParticleGun->SetMaxTheta(maxthetaCmd->GetNewDoubleValue(newValues));
    }
  else if(command == maxphiCmd)
    {
      fParticleGun->SetMaxPhi(maxphiCmd->GetNewDoubleValue(newValues));
    }
  else if(command == surfnormCmd)
    {
      fParticleGun->SetUserWRTSurface(surfnormCmd->GetNewBoolValue(newValues));
    }
  else if(command == energytypeCmd)
    {
      fParticleGun->SetEnergyDisType(newValues);
    }
  else if(command == eminCmd)
    {
      fParticleGun->SetEmin(eminCmd->GetNewDoubleValue(newValues));
    }
  else if(command == emaxCmd)
    {
      fParticleGun->SetEmax(emaxCmd->GetNewDoubleValue(newValues));
    }
  else if(command == monoenergyCmd)
    {
      fParticleGun->SetMonoEnergy(monoenergyCmd->GetNewDoubleValue(newValues));
    }
  else if(command == alphaCmd)
    {
      fParticleGun->SetAlpha(alphaCmd->GetNewDoubleValue(newValues));
    }
  else if(command == tempCmd)
    {
      fParticleGun->SetTemp(tempCmd->GetNewDoubleValue(newValues));
    }
  else if(command == ezeroCmd)
    {
      fParticleGun->SetEzero(ezeroCmd->GetNewDoubleValue(newValues));
    }
  else if(command == gradientCmd)
    {
      fParticleGun->SetGradient(gradientCmd->GetNewDoubleValue(newValues));
    }
  else if(command == interceptCmd)
    {
      fParticleGun->SetInterCept(interceptCmd->GetNewDoubleValue(newValues));
    }
  else if(command == calculateCmd)
    {
      fParticleGun->Calculate();
    }
  else if(command == energyspecCmd)
    {
      G4cout << "here in energyspecCmd" << G4endl;
      fParticleGun->InputEnergySpectra(energyspecCmd->GetNewBoolValue(newValues));
    }
  else if(command == diffspecCmd)
    {
      fParticleGun->InputDifferentialSpectra(diffspecCmd->GetNewBoolValue(newValues));
    }
  else if(command == histnameCmd)
    {
      histtype = newValues;
    }
  else if(command == histpointCmd)
    {
      if(histtype == "biasx")
	fParticleGun->SetXBias(histpointCmd->GetNew3VectorValue(newValues));
      if(histtype == "biasy")
	fParticleGun->SetYBias(histpointCmd->GetNew3VectorValue(newValues));
      if(histtype == "biasz")
	fParticleGun->SetZBias(histpointCmd->GetNew3VectorValue(newValues));
      if(histtype == "biast")
	fParticleGun->SetThetaBias(histpointCmd->GetNew3VectorValue(newValues));
      if(histtype == "biasp")
	fParticleGun->SetPhiBias(histpointCmd->GetNew3VectorValue(newValues));
      if(histtype == "biase")
	fParticleGun->SetEnergyBias(histpointCmd->GetNew3VectorValue(newValues));
      if(histtype == "theta")
	fParticleGun->UserDefAngTheta(histpointCmd->GetNew3VectorValue(newValues));
      if(histtype == "phi")
	fParticleGun->UserDefAngPhi(histpointCmd->GetNew3VectorValue(newValues));
      if(histtype == "energy")
	fParticleGun->UserEnergyHisto(histpointCmd->GetNew3VectorValue(newValues));
      if(histtype == "arb")
	fParticleGun->ArbEnergyHisto(histpointCmd->GetNew3VectorValue(newValues));
      if(histtype == "epn")
	fParticleGun->EpnEnergyHisto(histpointCmd->GetNew3VectorValue(newValues));
    }
  else if(command == arbintCmd)
    {
      fParticleGun->ArbInterpolate(newValues);
    }
  else if(command == verbosityCmd)
    {
      fParticleGun->SetVerbosity(verbosityCmd->GetNewIntValue(newValues));
    }
  else if( command==particleCmd )
    {
      if (newValues =="ion") {
	fShootIon = true;
      } else {
	fShootIon = false;
      G4ParticleDefinition* pd = particleTable->FindParticle(newValues);
      if(pd != NULL)
	{ fParticleGun->SetParticleDefinition( pd ); }
      }
    }
  //  else if(command == ionCmd)
  //    {
  //   fParticleGun->SetNucleus(ionCmd->GetNewNucleusValue(newValues));
  //  }
  else if( command==timeCmd )
    { fParticleGun->SetParticleTime(timeCmd->GetNewDoubleValue(newValues)); }
  else if( command==polCmd )
    { fParticleGun->SetParticlePolarization(polCmd->GetNew3VectorValue(newValues)); }
  else if( command==numberCmd )
    { fParticleGun->SetNumberOfParticles(numberCmd->GetNewIntValue(newValues)); }
  else if( command==ionCmd )
    { 
      if (fShootIon) {
	G4Tokenizer next( newValues );
	// check argument
	fAtomicNumber = StoI(next());
	fAtomicMass = StoI(next());
	G4String sQ = next();
	if (sQ.isNull()) {
	  fIonCharge = fAtomicNumber;
	} else {
	  fIonCharge = StoI(sQ);
	  sQ = next();
	  if (sQ.isNull()) {
	    fIonExciteEnergy = 0.0;
	  } else {
	    fIonExciteEnergy = StoD(sQ) * keV;
	  }
	}
	
	G4ParticleDefinition* ion;
	ion =  particleTable->GetIon( fAtomicNumber, fAtomicMass, fIonExciteEnergy);
	if (ion==0) {
	  G4cout << "Ion with Z=" << fAtomicNumber;
	  G4cout << " A=" << fAtomicMass << "is not be defined" << G4endl;    
	} else {
	  fParticleGun->SetParticleDefinition(ion);
	  fParticleGun->SetParticleCharge(fIonCharge*eplus);
	}
      } else {
	G4cout << "Set /gps/particle to ion before using /gps/ion command";
	G4cout << G4endl; 
      }
    }
  else
    {
      G4cout << "Error entering command" << G4endl;
    }
}

G4String G4GeneralParticleSourceMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String cv;
  
  //  if( command==directionCmd )
  //  { cv = directionCmd->ConvertToString(fParticleGun->GetParticleMomentumDirection()); }
  //  else if( command==energyCmd )
  //  { cv = energyCmd->ConvertToString(fParticleGun->GetParticleEnergy(),"GeV"); }
  //  else if( command==positionCmd )
  //  { cv = positionCmd->ConvertToString(fParticleGun->GetParticlePosition(),"cm"); }
  //  else if( command==timeCmd )
  //  { cv = timeCmd->ConvertToString(fParticleGun->GetParticleTime(),"ns"); }
  //  else if( command==polCmd )
  //  { cv = polCmd->ConvertToString(fParticleGun->GetParticlePolarization()); }
  //  else if( command==numberCmd )
  //  { cv = numberCmd->ConvertToString(fParticleGun->GetNumberOfParticles()); }
  
  cv = "Not implemented yet";

  return cv;
}
