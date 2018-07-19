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
///////////////////////////////////////////////////////////////////////////////
//
// MODULE:       G4GeneralParticleSourceMessenger.cc
//
// Version:      2.0
// Date:         5/02/04
// Author:       Fan Lei 
// Organisation: QinetiQ ltd.
// Customer:     ESA/ESTEC
//
///////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//
// Version 2.0, 05/02/2004, Fan Lei, Created.
//    After changes to version 1.1 as in Geant4 v6.0
//     - Mutilple particle source definition
//     - Re-structured commands
//     - old commands have been retained for backward compatibility, will be
//       removed in the future. 
//
//
// 20/03/2014, Andrew Green - Modifications for MT
//      - Added a check to force only one thread to parse the macro file
//          This information is fed into the GPS which now has a split mechanism for the large data (hence need to only read in 1 thread)
//      - Thread ID used is 0, so *should* work under interactive mode as well - may need checking, or the may be another way...
///////////////////////////////////////////////////////////////////////////////
//

#include "G4GeneralParticleSourceMessenger.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Geantino.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
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

#include "G4Tokenizer.hh"
#include "G4SingleParticleSource.hh"
#include "G4GeneralParticleSource.hh"

#include "G4AutoLock.hh"

namespace {
    G4Mutex creationM = G4MUTEX_INITIALIZER;
    G4GeneralParticleSourceMessenger* theInstance = 0;
}

G4GeneralParticleSourceMessenger* G4GeneralParticleSourceMessenger::GetInstance(G4GeneralParticleSource* psc)
{
    G4AutoLock l(&creationM);
    if ( theInstance == 0 ) theInstance = new G4GeneralParticleSourceMessenger(psc);
    return theInstance;
}

void G4GeneralParticleSourceMessenger::Destroy() {
    G4AutoLock l(&creationM);
    if ( theInstance != 0 ) {
        delete theInstance;
        theInstance = 0;
    }
}
///////////////////////////////////////////////////////////////////////////////
//
G4GeneralParticleSourceMessenger::G4GeneralParticleSourceMessenger
  (G4GeneralParticleSource *fPtclGun) 
    : fGPS(fPtclGun),fParticleGun(0),fShootIon(false),
      fAtomicNumber(0),fAtomicMass(0),fIonCharge(0),fIonExciteEnergy(0.),
      fAtomicNumberL(0),fAtomicMassL(0),fIonChargeL(0),fIonEnergyLevel(0)
      
{
  //A.Dotti - 10th October 2014
  //This messenger is special: it is instantiated in a user action but (e.g. in
  // a thread).
  //the UI commands it defines should be executed by the *master* thread because
  //they operate on shared resources and we want the UI commands to take
  //effect BEFORE the threads do some work (so all data are properly initialized)
  //To achieve this behavior we set to true a base calss protected
  //data member. Since it makes no sense to have more than one instance
  //of the messenger we check that we actually have only one.
  //Note that the logic of implementing, in a given worker thread only one
  //messenger is deletefated to the creator
  commandsShouldBeInMaster = true;

  
    
  particleTable = G4ParticleTable::GetParticleTable();
  histtype = "biasx";

  //UI Commands only for master
    G4bool broadcast = false;
  gpsDirectory = new G4UIdirectory("/gps/",broadcast);
    
  gpsDirectory->SetGuidance("General Paricle Source control commands.");
  //  gpsDirectory->SetGuidance(" The first 9 commands are the same as in G4ParticleGun ");

  // now the commands for mutiple sources
  sourceDirectory = new G4UIdirectory("/gps/source/");
  sourceDirectory->SetGuidance("Multiple source control sub-directory");

  addsourceCmd = new G4UIcmdWithADouble("/gps/source/add",this);
  addsourceCmd->SetGuidance("add a new source defintion to the particle gun with the specified intensity");
  addsourceCmd->SetParameterName("addsource",false,false);
  addsourceCmd->SetRange("addsource > 0.");

  listsourceCmd = new G4UIcmdWithoutParameter("/gps/source/list",this);
  listsourceCmd->SetGuidance("List the defined particle sources");

  clearsourceCmd = new G4UIcmdWithoutParameter("/gps/source/clear",this);
  clearsourceCmd->SetGuidance("Remove all the defined particle sources");

  getsourceCmd = new G4UIcmdWithoutParameter("/gps/source/show",this);
  getsourceCmd->SetGuidance("Show the current source index and intensity");

  setsourceCmd = new G4UIcmdWithAnInteger("/gps/source/set",this);
  setsourceCmd->SetGuidance("set the indexed source as the current one");
  setsourceCmd->SetGuidance("  so one can change its source definition");
  setsourceCmd->SetParameterName("setsource",false,false);
  setsourceCmd->SetRange("setsource >= 0");

  deletesourceCmd = new G4UIcmdWithAnInteger("/gps/source/delete",this);
  deletesourceCmd->SetGuidance("delete the indexed source from the list");
  deletesourceCmd->SetParameterName("deletesource",false,false);
  deletesourceCmd->SetRange("deletesource > 0");

  setintensityCmd = new G4UIcmdWithADouble("/gps/source/intensity",this);
  setintensityCmd->SetGuidance("reset the current source to the specified intensity");
  setintensityCmd->SetParameterName("setintensity",false,false);
  setintensityCmd->SetRange("setintensity > 0."); 

  multiplevertexCmd = new G4UIcmdWithABool("/gps/source/multiplevertex",this);
  multiplevertexCmd->SetGuidance("true for simulaneous generation mutiple vertex");
  multiplevertexCmd->SetGuidance("Default is false");
  multiplevertexCmd->SetParameterName("multiplevertex",true);
  multiplevertexCmd->SetDefaultValue(false);

  flatsamplingCmd = new G4UIcmdWithABool("/gps/source/flatsampling",this);
  flatsamplingCmd->SetGuidance("true for applying flat (biased) sampling among the sources");
  flatsamplingCmd->SetGuidance("Default is false");
  flatsamplingCmd->SetParameterName("flatsampling",true);
  flatsamplingCmd->SetDefaultValue(false);

  // below we reproduce commands awailable in G4Particle Gun

  listCmd = new G4UIcmdWithoutParameter("/gps/List",this);
  listCmd->SetGuidance("List available particles.");
  listCmd->SetGuidance(" Invoke G4ParticleTable.");

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

  directionCmd = new G4UIcmdWith3Vector("/gps/direction",this);
  directionCmd->SetGuidance("Set momentum direction.");
  directionCmd->SetGuidance("Direction needs not to be a unit vector.");
  directionCmd->SetGuidance("Angular distribution type is set to planar.");
  directionCmd->SetParameterName("Px","Py","Pz",false,false); 
  directionCmd->SetRange("Px != 0 || Py != 0 || Pz != 0");
  
  energyCmd = new G4UIcmdWithADoubleAndUnit("/gps/energy",this);
  energyCmd->SetGuidance("Set kinetic energy.");
  energyCmd->SetParameterName("Energy",false,false);
  energyCmd->SetDefaultUnit("GeV");
  //energyCmd->SetUnitCategory("Energy");
  //energyCmd->SetUnitCandidates("eV keV MeV GeV TeV");

  positionCmd = new G4UIcmdWith3VectorAndUnit("/gps/position",this);
  positionCmd->SetGuidance("Set starting position of the particle for a Point like source.");
  positionCmd->SetGuidance("Same effect as the two /gps/pos/type Point /gps/pos/centre commands.");
  positionCmd->SetParameterName("X","Y","Z",false,false);
  positionCmd->SetDefaultUnit("cm");
  //positionCmd->SetUnitCategory("Length");
  //positionCmd->SetUnitCandidates("microm mm cm m km");

  ionCmd = new G4UIcommand("/gps/ion",this);
  ionCmd->SetGuidance("Set properties of ion to be generated.");
  ionCmd->SetGuidance("[usage] /gps/ion Z A Q E");
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

  ionLvlCmd = new G4UIcommand("/gps/ionLvl",this);
  ionLvlCmd->SetGuidance("Set properties of ion to be generated.");
  ionLvlCmd->SetGuidance("[usage] /gps/ion Z A Q Lvl");
  ionLvlCmd->SetGuidance("        Z:(int) AtomicNumber");
  ionLvlCmd->SetGuidance("        A:(int) AtomicMass");
  ionLvlCmd->SetGuidance("        Q:(int) Charge of Ion (in unit of e)");
  ionLvlCmd->SetGuidance("        Lvl:(int) Number of metastable state excitation level (0-9)");

  G4UIparameter* paramL;
  paramL = new G4UIparameter("Z",'i',false);
  paramL->SetDefaultValue("1");
  ionLvlCmd->SetParameter(paramL);
  paramL = new G4UIparameter("A",'i',false);
  paramL->SetDefaultValue("1");
  ionLvlCmd->SetParameter(paramL);
  paramL = new G4UIparameter("Q",'i',true);
  paramL->SetDefaultValue("0");
  ionLvlCmd->SetParameter(paramL);
  paramL = new G4UIparameter("Lvl",'i',true);
  paramL->SetDefaultValue("0.0");
  ionLvlCmd->SetParameter(paramL);

  timeCmd = new G4UIcmdWithADoubleAndUnit("/gps/time",this);
  timeCmd->SetGuidance("Set initial time of the particle.");
  timeCmd->SetParameterName("t0",false,false);
  timeCmd->SetDefaultUnit("ns");
  //timeCmd->SetUnitCategory("Time");
  //timeCmd->SetUnitCandidates("ns ms s");
  
  polCmd = new G4UIcmdWith3Vector("/gps/polarization",this);
  polCmd->SetGuidance("Set polarization.");
  polCmd->SetParameterName("Px","Py","Pz",false,false); 
  polCmd->SetRange("Px>=-1.&&Px<=1.&&Py>=-1.&&Py<=1.&&Pz>=-1.&&Pz<=1.");

  numberCmd = new G4UIcmdWithAnInteger("/gps/number",this);
  numberCmd->SetGuidance("Set number of particles to be generated per vertex.");
  numberCmd->SetParameterName("N",false,false);
  numberCmd->SetRange("N>0");

  // verbosity
  verbosityCmd = new G4UIcmdWithAnInteger("/gps/verbose",this);
  verbosityCmd->SetGuidance("Set Verbose level for GPS");
  verbosityCmd->SetGuidance(" 0 : Silent");
  verbosityCmd->SetGuidance(" 1 : Limited information");
  verbosityCmd->SetGuidance(" 2 : Detailed information");
  verbosityCmd->SetParameterName("level",false);
  verbosityCmd->SetRange("level>=0 && level <=2");

  // now extended commands
  // Positional ones:
  positionDirectory = new G4UIdirectory("/gps/pos/");
  positionDirectory->SetGuidance("Positional commands sub-directory");

  typeCmd1 = new G4UIcmdWithAString("/gps/pos/type",this);
  typeCmd1->SetGuidance("Sets source distribution type.");
  typeCmd1->SetGuidance("Either Point, Beam, Plane, Surface or Volume");
  typeCmd1->SetParameterName("DisType",false,false);
  typeCmd1->SetDefaultValue("Point");
  typeCmd1->SetCandidates("Point Beam Plane Surface Volume");

  shapeCmd1 = new G4UIcmdWithAString("/gps/pos/shape",this);
  shapeCmd1->SetGuidance("Sets source shape for Plan, Surface or Volume type source.");
  shapeCmd1->SetParameterName("Shape",false,false);
  shapeCmd1->SetDefaultValue("NULL");
  shapeCmd1->SetCandidates("Circle Annulus Ellipse Square Rectangle Sphere Ellipsoid Cylinder EllipticCylinder Para");

  centreCmd1 = new G4UIcmdWith3VectorAndUnit("/gps/pos/centre",this);
  centreCmd1->SetGuidance("Set centre coordinates of source.");
  //centreCmd1->SetGuidance("   same effect as the /gps/position command");
  centreCmd1->SetParameterName("X","Y","Z",false,false);
  centreCmd1->SetDefaultUnit("cm");
  //  centreCmd1->SetUnitCandidates("micron mm cm m km");

  posrot1Cmd1 = new G4UIcmdWith3Vector("/gps/pos/rot1",this);
  posrot1Cmd1->SetGuidance("Set the 1st vector defining the rotation matrix'.");
  posrot1Cmd1->SetGuidance("It does not need to be a unit vector.");
  posrot1Cmd1->SetParameterName("R1x","R1y","R1z",false,false); 
  posrot1Cmd1->SetRange("R1x != 0 || R1y != 0 || R1z != 0");

  posrot2Cmd1 = new G4UIcmdWith3Vector("/gps/pos/rot2",this);
  posrot2Cmd1->SetGuidance("Set the 2nd vector defining the rotation matrix'.");
  posrot2Cmd1->SetGuidance("It does not need to be a unit vector.");
  posrot2Cmd1->SetParameterName("R2x","R2y","R2z",false,false); 
  posrot2Cmd1->SetRange("R2x != 0 || R2y != 0 || R2z != 0");

  halfxCmd1 = new G4UIcmdWithADoubleAndUnit("/gps/pos/halfx",this);
  halfxCmd1->SetGuidance("Set x half length of source.");
  halfxCmd1->SetParameterName("Halfx",false,false);
  halfxCmd1->SetDefaultUnit("cm");
  //  halfxCmd1->SetUnitCandidates("micron mm cm m km");

  halfyCmd1 = new G4UIcmdWithADoubleAndUnit("/gps/pos/halfy",this);
  halfyCmd1->SetGuidance("Set y half length of source.");
  halfyCmd1->SetParameterName("Halfy",false,false);
  halfyCmd1->SetDefaultUnit("cm");
  //  halfyCmd1->SetUnitCandidates("micron mm cm m km");

  halfzCmd1 = new G4UIcmdWithADoubleAndUnit("/gps/pos/halfz",this);
  halfzCmd1->SetGuidance("Set z half length of source.");
  halfzCmd1->SetParameterName("Halfz",false,false);
  halfzCmd1->SetDefaultUnit("cm");
  //  halfzCmd1->SetUnitCandidates("micron mm cm m km");

  radiusCmd1 = new G4UIcmdWithADoubleAndUnit("/gps/pos/radius",this);
  radiusCmd1->SetGuidance("Set radius of source.");
  radiusCmd1->SetParameterName("Radius",false,false);
  radiusCmd1->SetDefaultUnit("cm");
  //  radiusCmd1->SetUnitCandidates("micron mm cm m km");

  radius0Cmd1 = new G4UIcmdWithADoubleAndUnit("/gps/pos/inner_radius",this);
  radius0Cmd1->SetGuidance("Set inner radius of source when required.");
  radius0Cmd1->SetParameterName("Radius0",false,false);
  radius0Cmd1->SetDefaultUnit("cm");
  //  radius0Cmd1->SetUnitCandidates("micron mm cm m km");

  possigmarCmd1 = new G4UIcmdWithADoubleAndUnit("/gps/pos/sigma_r",this);
  possigmarCmd1->SetGuidance("Set standard deviation in radial of the beam positional profile");
  possigmarCmd1->SetGuidance(" applicable to Beam type source only");
  possigmarCmd1->SetParameterName("Sigmar",false,false);
  possigmarCmd1->SetDefaultUnit("cm");
  //  possigmarCmd1->SetUnitCandidates("micron mm cm m km");

  possigmaxCmd1 = new G4UIcmdWithADoubleAndUnit("/gps/pos/sigma_x",this);
  possigmaxCmd1->SetGuidance("Set standard deviation of beam positional profile in x-dir");
  possigmaxCmd1->SetGuidance(" applicable to Beam type source only");
  possigmaxCmd1->SetParameterName("Sigmax",false,false);
  possigmaxCmd1->SetDefaultUnit("cm");
  //  possigmaxCmd1->SetUnitCandidates("micron mm cm m km");

  possigmayCmd1 = new G4UIcmdWithADoubleAndUnit("/gps/pos/sigma_y",this);
  possigmayCmd1->SetGuidance("Set standard deviation of beam positional profile in y-dir");
  possigmayCmd1->SetGuidance(" applicable to Beam type source only");
  possigmayCmd1->SetParameterName("Sigmay",false,false);
  possigmayCmd1->SetDefaultUnit("cm");
  //  possigmayCmd1->SetUnitCandidates("micron mm cm m km");

  paralpCmd1 = new G4UIcmdWithADoubleAndUnit("/gps/pos/paralp",this);
  paralpCmd1->SetGuidance("Angle from y-axis of y' in Para");
  paralpCmd1->SetParameterName("paralp",false,false);
  paralpCmd1->SetDefaultUnit("rad");
  //  paralpCmd1->SetUnitCandidates("rad deg");

  partheCmd1 = new G4UIcmdWithADoubleAndUnit("/gps/pos/parthe",this);
  partheCmd1->SetGuidance("Polar angle through centres of z faces");
  partheCmd1->SetParameterName("parthe",false,false);
  partheCmd1->SetDefaultUnit("rad");
  //  partheCmd1->SetUnitCandidates("rad deg");

  parphiCmd1 = new G4UIcmdWithADoubleAndUnit("/gps/pos/parphi",this);
  parphiCmd1->SetGuidance("Azimuth angle through centres of z faces");
  parphiCmd1->SetParameterName("parphi",false,false);
  parphiCmd1->SetDefaultUnit("rad");
  //  parphiCmd1->SetUnitCandidates("rad deg");

  confineCmd1 = new G4UIcmdWithAString("/gps/pos/confine",this);
  confineCmd1->SetGuidance("Confine source to volume (NULL to unset).");
  confineCmd1->SetGuidance("usage: confine VolName");
  confineCmd1->SetParameterName("VolName",false,false);
  confineCmd1->SetDefaultValue("NULL");

  // old implementations
//  typeCmd = new G4UIcmdWithAString("/gps/type",this);
//  typeCmd->SetGuidance("Sets source distribution type. (obsolete!)");
//  typeCmd->SetGuidance("Either Point, Beam, Plane, Surface or Volume");
//  typeCmd->SetParameterName("DisType",false,false);
//  typeCmd->SetDefaultValue("Point");
//  typeCmd->SetCandidates("Point Beam Plane Surface Volume");

//  shapeCmd = new G4UIcmdWithAString("/gps/shape",this);
//  shapeCmd->SetGuidance("Sets source shape type.(obsolete!)");
//  shapeCmd->SetParameterName("Shape",false,false);
//  shapeCmd->SetDefaultValue("NULL");
//  shapeCmd->SetCandidates("Circle Annulus Ellipse Square Rectangle Sphere Ellipsoid Cylinder Para");

//  centreCmd = new G4UIcmdWith3VectorAndUnit("/gps/centre",this);
//  centreCmd->SetGuidance("Set centre coordinates of source.(obsolete!)");
//  centreCmd->SetParameterName("X","Y","Z",false,false);
//  centreCmd->SetDefaultUnit("cm");
  //  centreCmd->SetUnitCandidates("micron mm cm m km");

//  posrot1Cmd = new G4UIcmdWith3Vector("/gps/posrot1",this);
//  posrot1Cmd->SetGuidance("Set rotation matrix of x'.(obsolete!)");
//  posrot1Cmd->SetGuidance("Posrot1 does not need to be a unit vector.");
//  posrot1Cmd->SetParameterName("R1x","R1y","R1z",false,false); 
//  posrot1Cmd->SetRange("R1x != 0 || R1y != 0 || R1z != 0");
//
//  posrot2Cmd = new G4UIcmdWith3Vector("/gps/posrot2",this);
//  posrot2Cmd->SetGuidance("Set rotation matrix of y'.(obsolete!)");
//  posrot2Cmd->SetGuidance("Posrot2 does not need to be a unit vector.");
//  posrot2Cmd->SetParameterName("R2x","R2y","R2z",false,false); 
//  posrot2Cmd->SetRange("R2x != 0 || R2y != 0 || R2z != 0");

//  halfxCmd = new G4UIcmdWithADoubleAndUnit("/gps/halfx",this);
//  halfxCmd->SetGuidance("Set x half length of source.(obsolete!)");
//  halfxCmd->SetParameterName("Halfx",false,false);
//  halfxCmd->SetDefaultUnit("cm");
  //  halfxCmd->SetUnitCandidates("micron mm cm m km");

//  halfyCmd = new G4UIcmdWithADoubleAndUnit("/gps/halfy",this);
//  halfyCmd->SetGuidance("Set y half length of source.(obsolete!)");
//  halfyCmd->SetParameterName("Halfy",false,false);
//  halfyCmd->SetDefaultUnit("cm");
  //  halfyCmd->SetUnitCandidates("micron mm cm m km");

//  halfzCmd = new G4UIcmdWithADoubleAndUnit("/gps/halfz",this);
//  halfzCmd->SetGuidance("Set z half length of source.(obsolete!)");
//  halfzCmd->SetParameterName("Halfz",false,false);
//  halfzCmd->SetDefaultUnit("cm");
  //  halfzCmd->SetUnitCandidates("micron mm cm m km");

//  radiusCmd = new G4UIcmdWithADoubleAndUnit("/gps/radius",this);
//  radiusCmd->SetGuidance("Set radius of source.(obsolete!)");
//  radiusCmd->SetParameterName("Radius",false,false);
//  radiusCmd->SetDefaultUnit("cm");
  //  radiusCmd->SetUnitCandidates("micron mm cm m km");

//  radius0Cmd = new G4UIcmdWithADoubleAndUnit("/gps/radius0",this);
//  radius0Cmd->SetGuidance("Set inner radius of source.(obsolete!)");
//  radius0Cmd->SetParameterName("Radius0",false,false);
//  radius0Cmd->SetDefaultUnit("cm");
  //  radius0Cmd->SetUnitCandidates("micron mm cm m km");

//  possigmarCmd = new G4UIcmdWithADoubleAndUnit("/gps/sigmaposr",this);
//  possigmarCmd->SetGuidance("Set standard deviation of beam position in radial(obsolete!)");
//  possigmarCmd->SetParameterName("Sigmar",false,false);
//  possigmarCmd->SetDefaultUnit("cm");
  // possigmarCmd->SetUnitCandidates("micron mm cm m km");

//  possigmaxCmd = new G4UIcmdWithADoubleAndUnit("/gps/sigmaposx",this);
//  possigmaxCmd->SetGuidance("Set standard deviation of beam position in x-dir(obsolete!)");
//  possigmaxCmd->SetParameterName("Sigmax",false,false);
//  possigmaxCmd->SetDefaultUnit("cm");
  //  possigmaxCmd->SetUnitCandidates("micron mm cm m km");

//  possigmayCmd = new G4UIcmdWithADoubleAndUnit("/gps/sigmaposy",this);
//  possigmayCmd->SetGuidance("Set standard deviation of beam position in y-dir(obsolete!)");
//  possigmayCmd->SetParameterName("Sigmay",false,false);
//  possigmayCmd->SetDefaultUnit("cm");
  //  possigmayCmd->SetUnitCandidates("micron mm cm m km");

//  paralpCmd = new G4UIcmdWithADoubleAndUnit("/gps/paralp",this);
//  paralpCmd->SetGuidance("Angle from y-axis of y' in Para(obsolete!)");
//  paralpCmd->SetParameterName("paralp",false,false);
//  paralpCmd->SetDefaultUnit("rad");
  //  paralpCmd->SetUnitCandidates("rad deg");

//  partheCmd = new G4UIcmdWithADoubleAndUnit("/gps/parthe",this);
//  partheCmd->SetGuidance("Polar angle through centres of z faces(obsolete!)");
//  partheCmd->SetParameterName("parthe",false,false);
//  partheCmd->SetDefaultUnit("rad");
  //  partheCmd->SetUnitCandidates("rad deg");

//  parphiCmd = new G4UIcmdWithADoubleAndUnit("/gps/parphi",this);
//  parphiCmd->SetGuidance("Azimuth angle through centres of z faces(obsolete!)");
//  parphiCmd->SetParameterName("parphi",false,false);
//  parphiCmd->SetDefaultUnit("rad");
  //  parphiCmd->SetUnitCandidates("rad deg");

//  confineCmd = new G4UIcmdWithAString("/gps/confine",this);
//  confineCmd->SetGuidance("Confine source to volume (NULL to unset)(obsolete!) .");
//  confineCmd->SetGuidance("usage: confine VolName");
//  confineCmd->SetParameterName("VolName",false,false);
//  confineCmd->SetDefaultValue("NULL");

  // Angular distribution commands
  angularDirectory = new G4UIdirectory("/gps/ang/");
  angularDirectory->SetGuidance("Angular commands sub-directory");

  angtypeCmd1 = new G4UIcmdWithAString("/gps/ang/type",this);
  angtypeCmd1->SetGuidance("Sets angular source distribution type");
  angtypeCmd1->SetGuidance("Possible variables are: iso, cos, planar, beam1d, beam2d, focused or user");
  angtypeCmd1->SetParameterName("AngDis",false,false);
  angtypeCmd1->SetDefaultValue("iso");
  angtypeCmd1->SetCandidates("iso cos planar beam1d beam2d focused user");

  angrot1Cmd1 = new G4UIcmdWith3Vector("/gps/ang/rot1",this);
  angrot1Cmd1->SetGuidance("Sets the 1st vector for angular distribution rotation matrix");
  angrot1Cmd1->SetGuidance("Need not be a unit vector");
  angrot1Cmd1->SetParameterName("AR1x","AR1y","AR1z",false,false);
  angrot1Cmd1->SetRange("AR1x != 0 || AR1y != 0 || AR1z != 0");

  angrot2Cmd1 = new G4UIcmdWith3Vector("/gps/ang/rot2",this);
  angrot2Cmd1->SetGuidance("Sets the 2nd vector for angular distribution rotation matrix");
  angrot2Cmd1->SetGuidance("Need not be a unit vector");
  angrot2Cmd1->SetParameterName("AR2x","AR2y","AR2z",false,false);
  angrot2Cmd1->SetRange("AR2x != 0 || AR2y != 0 || AR2z != 0");

  minthetaCmd1 = new G4UIcmdWithADoubleAndUnit("/gps/ang/mintheta",this);
  minthetaCmd1->SetGuidance("Set minimum theta");
  minthetaCmd1->SetParameterName("MinTheta",false,false);
  minthetaCmd1->SetDefaultValue(0.);
  minthetaCmd1->SetDefaultUnit("rad");
  //  minthetaCmd1->SetUnitCandidates("rad deg");

  maxthetaCmd1 = new G4UIcmdWithADoubleAndUnit("/gps/ang/maxtheta",this);
  maxthetaCmd1->SetGuidance("Set maximum theta");
  maxthetaCmd1->SetParameterName("MaxTheta",false,false);
  maxthetaCmd1->SetDefaultValue(pi);
  maxthetaCmd1->SetDefaultUnit("rad");
  //  maxthetaCmd1->SetUnitCandidates("rad deg");

  minphiCmd1 = new G4UIcmdWithADoubleAndUnit("/gps/ang/minphi",this);
  minphiCmd1->SetGuidance("Set minimum phi");
  minphiCmd1->SetParameterName("MinPhi",false,false);
  minphiCmd1->SetDefaultUnit("rad");
  //  minphiCmd1->SetUnitCandidates("rad deg");

  maxphiCmd1 = new G4UIcmdWithADoubleAndUnit("/gps/ang/maxphi",this);
  maxphiCmd1->SetGuidance("Set maximum phi");
  maxphiCmd1->SetParameterName("MaxPhi",false,false);
  maxphiCmd1->SetDefaultValue(pi);
  maxphiCmd1->SetDefaultUnit("rad");
  //  maxphiCmd1->SetUnitCandidates("rad deg");

  angsigmarCmd1 = new G4UIcmdWithADoubleAndUnit("/gps/ang/sigma_r",this);
  angsigmarCmd1->SetGuidance("Set standard deviation in direction for 1D beam.");
  angsigmarCmd1->SetParameterName("Sigmara",false,false);
  angsigmarCmd1->SetDefaultUnit("rad");
  //  angsigmarCmd1->SetUnitCandidates("rad deg");
  
  angsigmaxCmd1 = new G4UIcmdWithADoubleAndUnit("/gps/ang/sigma_x",this);
  angsigmaxCmd1->SetGuidance("Set standard deviation in direction in x-direc. for 2D beam");
  angsigmaxCmd1->SetParameterName("Sigmaxa",false,false);
  angsigmaxCmd1->SetDefaultUnit("rad");
  //  angsigmaxCmd1->SetUnitCandidates("rad deg");

  angsigmayCmd1 = new G4UIcmdWithADoubleAndUnit("/gps/ang/sigma_y",this);
  angsigmayCmd1->SetGuidance("Set standard deviation in direction in y-direc. for 2D beam");
  angsigmayCmd1->SetParameterName("Sigmaya",false,false);
  angsigmayCmd1->SetDefaultUnit("rad");
  //  angsigmayCmd1->SetUnitCandidates("rad deg");

  angfocusCmd = new G4UIcmdWith3VectorAndUnit("/gps/ang/focuspoint",this);
  angfocusCmd->SetGuidance("Set the focusing point for the beam");
  angfocusCmd->SetParameterName("x","y","z",false,false);
  angfocusCmd->SetDefaultUnit("cm");
  //  angfocusCmd->SetUnitCandidates("micron mm cm m km");

  useuserangaxisCmd1 = new G4UIcmdWithABool("/gps/ang/user_coor",this);
  useuserangaxisCmd1->SetGuidance("true for using user defined angular co-ordinates");
  useuserangaxisCmd1->SetGuidance("Default is false");
  useuserangaxisCmd1->SetParameterName("useuserangaxis",true);
  useuserangaxisCmd1->SetDefaultValue(false);

  surfnormCmd1 = new G4UIcmdWithABool("/gps/ang/surfnorm",this);
  surfnormCmd1->SetGuidance("Makes a user-defined distribution with respect to surface normals rather than x,y,z axes.");
  surfnormCmd1->SetGuidance("Default is false");
  surfnormCmd1->SetParameterName("surfnorm",true);
  surfnormCmd1->SetDefaultValue(false);

  // old ones
//  angtypeCmd = new G4UIcmdWithAString("/gps/angtype",this);
//  angtypeCmd->SetGuidance("Sets angular source distribution type (obsolete!)");
//  angtypeCmd->SetGuidance("Possible variables are: iso, cos planar beam1d beam2d or user");
//  angtypeCmd->SetParameterName("AngDis",false,false);
//  angtypeCmd->SetDefaultValue("iso");
//  angtypeCmd->SetCandidates("iso cos planar beam1d beam2d user");

//  angrot1Cmd = new G4UIcmdWith3Vector("/gps/angrot1",this);
//  angrot1Cmd->SetGuidance("Sets the x' vector for angular distribution(obsolete!) ");
//  angrot1Cmd->SetGuidance("Need not be a unit vector");
//  angrot1Cmd->SetParameterName("AR1x","AR1y","AR1z",false,false);
//  angrot1Cmd->SetRange("AR1x != 0 || AR1y != 0 || AR1z != 0");

//  angrot2Cmd = new G4UIcmdWith3Vector("/gps/angrot2",this);
//  angrot2Cmd->SetGuidance("Sets the y' vector for angular distribution (obsolete!)");
//  angrot2Cmd->SetGuidance("Need not be a unit vector");
//  angrot2Cmd->SetParameterName("AR2x","AR2y","AR2z",false,false);
//  angrot2Cmd->SetRange("AR2x != 0 || AR2y != 0 || AR2z != 0");

//  minthetaCmd = new G4UIcmdWithADoubleAndUnit("/gps/mintheta",this);
//  minthetaCmd->SetGuidance("Set minimum theta (obsolete!)");
//  minthetaCmd->SetParameterName("MinTheta",false,false);
//  minthetaCmd->SetDefaultUnit("rad");
  //  minthetaCmd->SetUnitCandidates("rad deg");

//  maxthetaCmd = new G4UIcmdWithADoubleAndUnit("/gps/maxtheta",this);
//  maxthetaCmd->SetGuidance("Set maximum theta (obsolete!)");
//  maxthetaCmd->SetParameterName("MaxTheta",false,false);
//  maxthetaCmd->SetDefaultValue(3.1416);
//  maxthetaCmd->SetDefaultUnit("rad");
  //  maxthetaCmd->SetUnitCandidates("rad deg");

//  minphiCmd = new G4UIcmdWithADoubleAndUnit("/gps/minphi",this);
//  minphiCmd->SetGuidance("Set minimum phi (obsolete!)");
//  minphiCmd->SetParameterName("MinPhi",false,false);
//  minphiCmd->SetDefaultUnit("rad");
  //  minphiCmd->SetUnitCandidates("rad deg");

//  maxphiCmd = new G4UIcmdWithADoubleAndUnit("/gps/maxphi",this);
//  maxphiCmd->SetGuidance("Set maximum phi(obsolete!)");
//  maxphiCmd->SetParameterName("MaxPhi",false,false);
//  maxphiCmd->SetDefaultUnit("rad");
  //  maxphiCmd->SetUnitCandidates("rad deg");

//  angsigmarCmd = new G4UIcmdWithADoubleAndUnit("/gps/sigmaangr",this);
//  angsigmarCmd->SetGuidance("Set standard deviation of beam direction in radial(obsolete!).");
//  angsigmarCmd->SetParameterName("Sigmara",false,false);
//  angsigmarCmd->SetDefaultUnit("rad");
  //  angsigmarCmd->SetUnitCandidates("rad deg");

//  angsigmaxCmd = new G4UIcmdWithADoubleAndUnit("/gps/sigmaangx",this);
//  angsigmaxCmd->SetGuidance("Set standard deviation of beam direction in x-direc(obsolete!).");
//  angsigmaxCmd->SetParameterName("Sigmaxa",false,false);
//  angsigmaxCmd->SetDefaultUnit("rad");
  //  angsigmaxCmd->SetUnitCandidates("rad deg");

//  angsigmayCmd = new G4UIcmdWithADoubleAndUnit("/gps/sigmaangy",this);
//  angsigmayCmd->SetGuidance("Set standard deviation of beam direction in y-direc.(obsolete!)");
//  angsigmayCmd->SetParameterName("Sigmaya",false,false);
//  angsigmayCmd->SetDefaultUnit("rad");
  //  angsigmayCmd->SetUnitCandidates("rad deg");

//  useuserangaxisCmd = new G4UIcmdWithABool("/gps/useuserangaxis",this);
//  useuserangaxisCmd->SetGuidance("true for using user defined angular co-ordinates(obsolete!)");
//  useuserangaxisCmd->SetGuidance("Default is false");
//  useuserangaxisCmd->SetParameterName("useuserangaxis",true);
//  useuserangaxisCmd->SetDefaultValue(false);

//  surfnormCmd = new G4UIcmdWithABool("/gps/surfnorm",this);
//  surfnormCmd->SetGuidance("Makes a user-defined distribution with respect to surface normals rather than x,y,z axes (obsolete!).");
//  surfnormCmd->SetGuidance("Default is false");
//  surfnormCmd->SetParameterName("surfnorm",true);
//  surfnormCmd->SetDefaultValue(false);

  // Energy commands

  energyDirectory = new G4UIdirectory("/gps/ene/");
  energyDirectory->SetGuidance("Spectral commands sub-directory");

  energytypeCmd1 = new G4UIcmdWithAString("/gps/ene/type",this);
  energytypeCmd1->SetGuidance("Sets energy distribution type");
  energytypeCmd1->SetParameterName("EnergyDis",false,false);
  energytypeCmd1->SetDefaultValue("Mono");
  energytypeCmd1->SetCandidates("Mono Lin Pow Exp CPow Gauss Brem Bbody Cdg User Arb Epn");

  eminCmd1 = new G4UIcmdWithADoubleAndUnit("/gps/ene/min",this);
  eminCmd1->SetGuidance("Sets minimum energy");
  eminCmd1->SetParameterName("emin",false,false);
  eminCmd1->SetDefaultUnit("keV");
  //  eminCmd1->SetUnitCandidates("eV keV MeV GeV TeV PeV");

  emaxCmd1 = new G4UIcmdWithADoubleAndUnit("/gps/ene/max",this);
  emaxCmd1->SetGuidance("Sets maximum energy");
  emaxCmd1->SetParameterName("emax",false,false);
  emaxCmd1->SetDefaultUnit("keV");
  //  emaxCmd1->SetUnitCandidates("eV keV MeV GeV TeV PeV");

  monoenergyCmd1 = new G4UIcmdWithADoubleAndUnit("/gps/ene/mono",this);
  monoenergyCmd1->SetGuidance("Sets a monocromatic energy (same as  gps/energy)");
  monoenergyCmd1->SetParameterName("monoenergy",false,false);
  monoenergyCmd1->SetDefaultUnit("keV");
  //  monoenergyCmd1->SetUnitCandidates("eV keV MeV GeV TeV PeV");

  engsigmaCmd1 = new G4UIcmdWithADoubleAndUnit("/gps/ene/sigma",this);
  engsigmaCmd1->SetGuidance("Sets the standard deviation for Gaussian energy dist.");
  engsigmaCmd1->SetParameterName("Sigmae",false,false);
  engsigmaCmd1->SetDefaultUnit("keV");
  //  engsigmaCmd1->SetUnitCandidates("eV keV MeV GeV TeV PeV");

  alphaCmd1 = new G4UIcmdWithADouble("/gps/ene/alpha",this);
  alphaCmd1->SetGuidance("Sets Alpha (index) for power-law energy dist.");
  alphaCmd1->SetParameterName("alpha",false,false);
  
  tempCmd1 = new G4UIcmdWithADouble("/gps/ene/temp",this);
  tempCmd1->SetGuidance("Sets the temperature for Brem and BBody distributions (in Kelvin)");
  tempCmd1->SetParameterName("temp",false,false);

  ezeroCmd1 = new G4UIcmdWithADouble("/gps/ene/ezero",this);
  ezeroCmd1->SetGuidance("Sets E_0 for exponential distribution (in MeV)");
  ezeroCmd1->SetParameterName("ezero",false,false);

  gradientCmd1 = new G4UIcmdWithADouble("/gps/ene/gradient",this);
  gradientCmd1->SetGuidance("Sets the gradient for Lin distribution (in 1/MeV)");
  gradientCmd1->SetParameterName("gradient",false,false);

  interceptCmd1 = new G4UIcmdWithADouble("/gps/ene/intercept",this);
  interceptCmd1->SetGuidance("Sets the intercept for Lin distributions (in MeV)");
  interceptCmd1->SetParameterName("intercept",false,false);

  arbeintCmd1 = new G4UIcmdWithADouble("/gps/ene/biasAlpha",this);
  arbeintCmd1->SetGuidance("Set the power-law index for the energy sampling distri. )");
  arbeintCmd1->SetParameterName("arbeint",false,false);

  calculateCmd1 = new G4UIcmdWithoutParameter("/gps/ene/calculate",this);
  calculateCmd1->SetGuidance("Calculates the distributions for Cdg and BBody");

  energyspecCmd1 = new G4UIcmdWithABool("/gps/ene/emspec",this);
  energyspecCmd1->SetGuidance("True for energy and false for momentum spectra");
  energyspecCmd1->SetParameterName("energyspec",true);
  energyspecCmd1->SetDefaultValue(true);

  diffspecCmd1 = new G4UIcmdWithABool("/gps/ene/diffspec",this);
  diffspecCmd1->SetGuidance("True for differential and flase for integral spectra");
  diffspecCmd1->SetParameterName("diffspec",true);
  diffspecCmd1->SetDefaultValue(true);

  //old ones
//  energytypeCmd = new G4UIcmdWithAString("/gps/energytype",this);
//  energytypeCmd->SetGuidance("Sets energy distribution type (obsolete!)");
//  energytypeCmd->SetParameterName("EnergyDis",false,false);
//  energytypeCmd->SetDefaultValue("Mono");
//  energytypeCmd->SetCandidates("Mono Lin Pow Exp Gauss Brem Bbody Cdg User Arb Epn");

//  eminCmd = new G4UIcmdWithADoubleAndUnit("/gps/emin",this);
//  eminCmd->SetGuidance("Sets Emin (obsolete!)");
//  eminCmd->SetParameterName("emin",false,false);
//  eminCmd->SetDefaultUnit("keV");
  //  eminCmd->SetUnitCandidates("eV keV MeV GeV TeV PeV");

//  emaxCmd = new G4UIcmdWithADoubleAndUnit("/gps/emax",this);
//  emaxCmd->SetGuidance("Sets Emax (obsolete!)");
//  emaxCmd->SetParameterName("emax",false,false);
//  emaxCmd->SetDefaultUnit("keV");
  //  emaxCmd->SetUnitCandidates("eV keV MeV GeV TeV PeV");

//  monoenergyCmd = new G4UIcmdWithADoubleAndUnit("/gps/monoenergy",this);
//  monoenergyCmd->SetGuidance("Sets Monoenergy (obsolete, use gps/energy instead!)");
//  monoenergyCmd->SetParameterName("monoenergy",false,false);
//  monoenergyCmd->SetDefaultUnit("keV");
  //  monoenergyCmd->SetUnitCandidates("eV keV MeV GeV TeV PeV");

//  engsigmaCmd = new G4UIcmdWithADoubleAndUnit("/gps/sigmae",this);
//  engsigmaCmd->SetGuidance("Sets the standard deviation for Gaussian energy dist.(obsolete!)");
//  engsigmaCmd->SetParameterName("Sigmae",false,false);
//  engsigmaCmd->SetDefaultUnit("keV");
  //  engsigmaCmd->SetUnitCandidates("eV keV MeV GeV TeV PeV");

//  alphaCmd = new G4UIcmdWithADouble("/gps/alpha",this);
//  alphaCmd->SetGuidance("Sets Alpha (index) for power-law energy dist(obsolete!).");
//  alphaCmd->SetParameterName("alpha",false,false);
  
//  tempCmd = new G4UIcmdWithADouble("/gps/temp",this);
//  tempCmd->SetGuidance("Sets the temperature for Brem and BBody (in Kelvin)(obsolete!)");
//  tempCmd->SetParameterName("temp",false,false);

//  ezeroCmd = new G4UIcmdWithADouble("/gps/ezero",this);
//  ezeroCmd->SetGuidance("Sets ezero exponential distributions (in MeV)(obsolete!)");
//  ezeroCmd->SetParameterName("ezero",false,false);
//
//  gradientCmd = new G4UIcmdWithADouble("/gps/gradient",this);
//  gradientCmd->SetGuidance("Sets the gradient for Lin distributions (in 1/MeV)(obsolete!)");
//  gradientCmd->SetParameterName("gradient",false,false);
//
//  interceptCmd = new G4UIcmdWithADouble("/gps/intercept",this);
//  interceptCmd->SetGuidance("Sets the intercept for Lin distributions (in MeV)(obsolete!)");
//  interceptCmd->SetParameterName("intercept",false,false);
//
//  calculateCmd = new G4UIcmdWithoutParameter("/gps/calculate",this);
//  calculateCmd->SetGuidance("Calculates distributions for Cdg and BBody(obsolete!)");
//
//  energyspecCmd = new G4UIcmdWithABool("/gps/energyspec",this);
//  energyspecCmd->SetGuidance("True for energy and false for momentum spectra(obsolete!)");
//  energyspecCmd->SetParameterName("energyspec",true);
//  energyspecCmd->SetDefaultValue(true);
//
//  diffspecCmd = new G4UIcmdWithABool("/gps/diffspec",this);
//  diffspecCmd->SetGuidance("True for differential and flase for integral spectra(obsolete!)");
//  diffspecCmd->SetParameterName("diffspec",true);
//  diffspecCmd->SetDefaultValue(true);

  // Biasing + histograms in general
  histDirectory = new G4UIdirectory("/gps/hist/");
  histDirectory->SetGuidance("Histogram, biasing commands sub-directory");

  histnameCmd1 = new G4UIcmdWithAString("/gps/hist/type",this);
  histnameCmd1->SetGuidance("Sets histogram type");
  histnameCmd1->SetParameterName("HistType",false,false);
  histnameCmd1->SetDefaultValue("biasx");
  histnameCmd1->SetCandidates("biasx biasy biasz biast biasp biase biaspt biaspp theta phi energy arb epn");

  resethistCmd1 = new G4UIcmdWithAString("/gps/hist/reset",this);
  resethistCmd1->SetGuidance("Reset (clean) the histogram ");
  resethistCmd1->SetParameterName("HistType",false,false);
  resethistCmd1->SetDefaultValue("energy");
  resethistCmd1->SetCandidates("biasx biasy biasz biast biasp biase biaspt biaspp theta phi energy arb epn");

  histpointCmd1 = new G4UIcmdWith3Vector("/gps/hist/point",this);
  histpointCmd1->SetGuidance("Allows user to define a histogram");
  histpointCmd1->SetGuidance("Enter: Ehi Weight");
  histpointCmd1->SetParameterName("Ehi","Weight","Junk",true,true);
  histpointCmd1->SetRange("Ehi >= 0. && Weight >= 0.");

  histfileCmd1 = new G4UIcmdWithAString("/gps/hist/file",this);
  histfileCmd1->SetGuidance("import the arb energy hist in an ASCII file");
  histfileCmd1->SetParameterName("HistFile",false,false);

  arbintCmd1 = new G4UIcmdWithAString("/gps/hist/inter",this);
  arbintCmd1->SetGuidance("Sets the interpolation method for arbitrary distribution.");
  arbintCmd1->SetParameterName("int",false,false);
  arbintCmd1->SetDefaultValue("Lin");
  arbintCmd1->SetCandidates("Lin Log Exp Spline");

  // old ones
//  histnameCmd = new G4UIcmdWithAString("/gps/histname",this);
//  histnameCmd->SetGuidance("Sets histogram type (obsolete!)");
//  histnameCmd->SetParameterName("HistType",false,false);
//  histnameCmd->SetDefaultValue("biasx");
//  histnameCmd->SetCandidates("biasx biasy biasz biast biasp biase biaspt biaspp theta phi energy arb epn");
//
//  // re-set the histograms
//  resethistCmd = new G4UIcmdWithAString("/gps/resethist",this);
//  resethistCmd->SetGuidance("Re-Set the histogram (obsolete!)");
//  resethistCmd->SetParameterName("HistType",false,false);
//  resethistCmd->SetDefaultValue("energy");
//  resethistCmd->SetCandidates("biasx biasy biasz biast biasp biase biaspt biaspp theta phi energy arb epn");
//
//  histpointCmd = new G4UIcmdWith3Vector("/gps/histpoint",this);
//  histpointCmd->SetGuidance("Allows user to define a histogram (obsolete!)");
//  histpointCmd->SetGuidance("Enter: Ehi Weight");
//  histpointCmd->SetParameterName("Ehi","Weight","Junk",false,false);
//  histpointCmd->SetRange("Ehi >= 0. && Weight >= 0.");
//
//  arbintCmd = new G4UIcmdWithAString("/gps/arbint",this);
//  arbintCmd->SetGuidance("Sets Arbitrary Interpolation type.(obsolete!) ");
//  arbintCmd->SetParameterName("int",false,false);
//  arbintCmd->SetDefaultValue("NULL");
//  arbintCmd->SetCandidates("Lin Log Exp Spline");

}

G4GeneralParticleSourceMessenger::~G4GeneralParticleSourceMessenger()
{
  delete positionDirectory;
//  delete typeCmd;
//  delete shapeCmd;
//  delete centreCmd;
//  delete posrot1Cmd;
//  delete posrot2Cmd;
//  delete halfxCmd;
//  delete halfyCmd;
//  delete halfzCmd;
//  delete radiusCmd;
//  delete radius0Cmd;
//  delete possigmarCmd;
//  delete possigmaxCmd;
//  delete possigmayCmd;
//  delete paralpCmd;
//  delete partheCmd;
//  delete parphiCmd;
//  delete confineCmd;
  delete typeCmd1;
  delete shapeCmd1;
  delete centreCmd1;
  delete posrot1Cmd1;
  delete posrot2Cmd1;
  delete halfxCmd1;
  delete halfyCmd1;
  delete halfzCmd1;
  delete radiusCmd1;
  delete radius0Cmd1;
  delete possigmarCmd1;
  delete possigmaxCmd1;
  delete possigmayCmd1;
  delete paralpCmd1;
  delete partheCmd1;
  delete parphiCmd1;
  delete confineCmd1;

  delete angularDirectory;
//  delete angtypeCmd;
//  delete angrot1Cmd;
//  delete angrot2Cmd;
//  delete minthetaCmd;
//  delete maxthetaCmd;
//  delete minphiCmd;
//  delete maxphiCmd;
//  delete angsigmarCmd;
//  delete angsigmaxCmd;
//  delete angsigmayCmd;
//  delete useuserangaxisCmd;
//  delete surfnormCmd;
  delete angtypeCmd1;
  delete angrot1Cmd1;
  delete angrot2Cmd1;
  delete minthetaCmd1;
  delete maxthetaCmd1;
  delete minphiCmd1;
  delete maxphiCmd1;
  delete angsigmarCmd1;
  delete angsigmaxCmd1;
  delete angsigmayCmd1;
  delete angfocusCmd;
  delete useuserangaxisCmd1;
  delete surfnormCmd1;

  delete energyDirectory;
//  delete energytypeCmd;
//  delete eminCmd;
//  delete emaxCmd;
//  delete monoenergyCmd;
//  delete engsigmaCmd;
//  delete alphaCmd;
//  delete tempCmd;
//  delete ezeroCmd;
//  delete gradientCmd;
//  delete interceptCmd;
//  delete calculateCmd;
//  delete energyspecCmd;
//  delete diffspecCmd;
  delete energytypeCmd1;
  delete eminCmd1;
  delete emaxCmd1;
  delete monoenergyCmd1;
  delete engsigmaCmd1;
  delete alphaCmd1;
  delete tempCmd1;
  delete ezeroCmd1;
  delete gradientCmd1;
  delete interceptCmd1;
  delete arbeintCmd1;
  delete calculateCmd1;
  delete energyspecCmd1;
  delete diffspecCmd1;

  delete histDirectory;
  //delete histnameCmd;
  //delete resethistCmd;
  //delete histpointCmd;
  //delete arbintCmd;
  delete histnameCmd1;
  delete resethistCmd1;
  delete histpointCmd1;
  delete histfileCmd1;
  delete arbintCmd1;

  delete verbosityCmd;
  delete ionCmd;
  delete ionLvlCmd;
  delete particleCmd;
  delete timeCmd;
  delete polCmd;
  delete numberCmd;
  delete positionCmd;
  delete directionCmd;
  delete energyCmd;
  delete listCmd;

  delete sourceDirectory;
  delete addsourceCmd;
  delete listsourceCmd;
  delete clearsourceCmd;
  delete getsourceCmd;
  delete setsourceCmd;
  delete setintensityCmd;
  delete deletesourceCmd;
  delete multiplevertexCmd;
  delete flatsamplingCmd;

  delete gpsDirectory;
  theInstance = 0;
}

#define CHECKPG() { if (fParticleGun==NULL) { \
                 G4ExceptionDescription msg; \
                 msg << "Command "<< command->GetCommandPath()<<"/";\
                 msg<<command->GetCommandName() << " used but no particle sources are set.";\
                 msg <<" Add at least a source with: /gps/source/add.";\
		      G4Exception("G4GeneralParticleSourceMessenger::SetNewValue","G4GPS003",\
				  FatalException,msg); return;\
    } }

void G4GeneralParticleSourceMessenger::SetNewValue(G4UIcommand *command, G4String newValues)
{
//  if(command == typeCmd)
//    {
//      CHECKPG(); fParticleGun->GetPosDist()->SetPosDisType(newValues);
//      G4cout << " G4GeneralParticleSourceMessenger - Warning:" << G4endl
//             << " The command is obsolete and will be removed soon." << G4endl
//             << " Please try to use the new structured commands!" << G4endl;
//    }
//  else if(command == shapeCmd)
//    {
//      CHECKPG(); fParticleGun->GetPosDist()->SetPosDisShape(newValues);
//      G4cout << " G4GeneralParticleSourceMessenger - Warning:" << G4endl
//             << " The command is obsolete and will be removed soon." << G4endl
//             << " Please try to use the new structured commands!" << G4endl;
//    }
//  else if(command == centreCmd)
//    {
//      CHECKPG(); fParticleGun->GetPosDist()->SetCentreCoords(centreCmd->GetNew3VectorValue(newValues));
//      G4cout << " G4GeneralParticleSourceMessenger - Warning:" << G4endl
//             << " The command is obsolete and will be removed soon." << G4endl
//             << " Please try to use the new structured commands!" << G4endl;
//    }
//  else if(command == posrot1Cmd)
//    {
//      CHECKPG(); fParticleGun->GetPosDist()->SetPosRot1(posrot1Cmd->GetNew3VectorValue(newValues));
//      G4cout << " G4GeneralParticleSourceMessenger - Warning:" << G4endl
//             << " The command is obsolete and will be removed soon." << G4endl
//             << " Please try to use the new structured commands!" << G4endl;
//    }
//  else if(command == posrot2Cmd)
//    {
//      CHECKPG(); fParticleGun->GetPosDist()->SetPosRot2(posrot2Cmd->GetNew3VectorValue(newValues));
//      G4cout << " G4GeneralParticleSourceMessenger - Warning:" << G4endl
//             << " The command is obsolete and will be removed soon." << G4endl
//             << " Please try to use the new structured commands!" << G4endl;
//    }
//  else if(command == halfxCmd)
//    {
//      CHECKPG(); fParticleGun->GetPosDist()->SetHalfX(halfxCmd->GetNewDoubleValue(newValues));
//      G4cout << " G4GeneralParticleSourceMessenger - Warning:" << G4endl
//             << " The command is obsolete and will be removed soon." << G4endl
//             << " Please try to use the new structured commands!" << G4endl;
//    }
//  else if(command == halfyCmd)
//    {
//      CHECKPG(); fParticleGun->GetPosDist()->SetHalfY(halfyCmd->GetNewDoubleValue(newValues));
//      G4cout << " G4GeneralParticleSourceMessenger - Warning:" << G4endl
//             << " The command is obsolete and will be removed soon." << G4endl
//             << " Please try to use the new structured commands!" << G4endl;
//    }
//  else if(command == halfzCmd)
//    {
//      CHECKPG(); fParticleGun->GetPosDist()->SetHalfZ(halfzCmd->GetNewDoubleValue(newValues));
//      G4cout << " G4GeneralParticleSourceMessenger - Warning:" << G4endl
//             << " The command is obsolete and will be removed soon." << G4endl
//             << " Please try to use the new structured commands!" << G4endl;
//    }
//  else if(command == radiusCmd)
//    {
//      CHECKPG(); fParticleGun->GetPosDist()->SetRadius(radiusCmd->GetNewDoubleValue(newValues));
//      G4cout << " G4GeneralParticleSourceMessenger - Warning:" << G4endl
//             << " The command is obsolete and will be removed soon." << G4endl
//             << " Please try to use the new structured commands!" << G4endl;
//    }
//  else if(command == radius0Cmd)
//    {
//      CHECKPG(); fParticleGun->GetPosDist()->SetRadius0(radius0Cmd->GetNewDoubleValue(newValues));
//      G4cout << " G4GeneralParticleSourceMessenger - Warning:" << G4endl
//             << " The command is obsolete and will be removed soon." << G4endl
//             << " Please try to use the new structured commands!" << G4endl;
//    }
//  else if(command == possigmarCmd)
//    {
//      CHECKPG(); fParticleGun->GetPosDist()->SetBeamSigmaInR(possigmarCmd->GetNewDoubleValue(newValues));
//      G4cout << " G4GeneralParticleSourceMessenger - Warning:" << G4endl
//             << " The command is obsolete and will be removed soon." << G4endl
//             << " Please try to use the new structured commands!" << G4endl;
//    }
//  else if(command == possigmaxCmd)
//    {
//      CHECKPG(); fParticleGun->GetPosDist()->SetBeamSigmaInX(possigmaxCmd->GetNewDoubleValue(newValues));
//      G4cout << " G4GeneralParticleSourceMessenger - Warning:" << G4endl
//             << " The command is obsolete and will be removed soon." << G4endl
//             << " Please try to use the new structured commands!" << G4endl;
//    }
//  else if(command == possigmayCmd)
//    {
//      CHECKPG(); fParticleGun->GetPosDist()->SetBeamSigmaInY(possigmayCmd->GetNewDoubleValue(newValues));
//      G4cout << " G4GeneralParticleSourceMessenger - Warning:" << G4endl
//             << " The command is obsolete and will be removed soon." << G4endl
//             << " Please try to use the new structured commands!" << G4endl;
//    }
//  else if(command == paralpCmd)
//    {
//      CHECKPG(); fParticleGun->GetPosDist()->SetParAlpha(paralpCmd->GetNewDoubleValue(newValues));
//      G4cout << " G4GeneralParticleSourceMessenger - Warning:" << G4endl
//             << " The command is obsolete and will be removed soon." << G4endl
//             << " Please try to use the new structured commands!" << G4endl;
//    }
//  else if(command == partheCmd)
//    {
//      CHECKPG(); fParticleGun->GetPosDist()->SetParTheta(partheCmd->GetNewDoubleValue(newValues));
//      G4cout << " G4GeneralParticleSourceMessenger - Warning:" << G4endl
//             << " The command is obsolete and will be removed soon." << G4endl
//             << " Please try to use the new structured commands!" << G4endl;
//    }
//  else if(command == parphiCmd)
//    {
//      CHECKPG(); fParticleGun->GetPosDist()->SetParPhi(parphiCmd->GetNewDoubleValue(newValues));
//      G4cout << " G4GeneralParticleSourceMessenger - Warning:" << G4endl
//             << " The command is obsolete and will be removed soon." << G4endl
//             << " Please try to use the new structured commands!" << G4endl;
//    }
//  else if(command == confineCmd)
//    {
//      CHECKPG(); fParticleGun->GetPosDist()->ConfineSourceToVolume(newValues);
//      G4cout << " G4GeneralParticleSourceMessenger - Warning:" << G4endl
//             << " The command is obsolete and will be removed soon." << G4endl
//             << " Please try to use the new structured commands!" << G4endl;
//    }
//  else if(command == angtypeCmd)
//    {
//      CHECKPG(); fParticleGun->GetAngDist()->SetAngDistType(newValues);
//      G4cout << " G4GeneralParticleSourceMessenger - Warning:" << G4endl
//             << " The command is obsolete and will be removed soon." << G4endl
//             << " Please try to use the new structured commands!" << G4endl;
//    }
//  else if(command == angrot1Cmd)
//    {
//      CHECKPG();
//      G4String a = "angref1";
//      fParticleGun->GetAngDist()->DefineAngRefAxes(a,angrot1Cmd->GetNew3VectorValue(newValues));
//      G4cout << " G4GeneralParticleSourceMessenger - Warning:" << G4endl
//             << " The command is obsolete and will be removed soon." << G4endl
//             << " Please try to use the new structured commands!" << G4endl;
//    }
//  else if(command == angrot2Cmd)
//    {
//      CHECKPG();
//      G4String a = "angref2";
//      fParticleGun->GetAngDist()->DefineAngRefAxes(a,angrot2Cmd->GetNew3VectorValue(newValues));
//      G4cout << " G4GeneralParticleSourceMessenger - Warning:" << G4endl
//             << " The command is obsolete and will be removed soon." << G4endl
//             << " Please try to use the new structured commands!" << G4endl;
//    }
//  else if(command == minthetaCmd)
//    {
//      CHECKPG(); fParticleGun->GetAngDist()->SetMinTheta(minthetaCmd->GetNewDoubleValue(newValues));
//      G4cout << " G4GeneralParticleSourceMessenger - Warning:" << G4endl
//             << " The command is obsolete and will be removed soon." << G4endl
//             << " Please try to use the new structured commands!" << G4endl;
//    }
//  else if(command == minphiCmd)
//    {
//      CHECKPG(); fParticleGun->GetAngDist()->SetMinPhi(minphiCmd->GetNewDoubleValue(newValues));
//      G4cout << " G4GeneralParticleSourceMessenger - Warning:" << G4endl
//             << " The command is obsolete and will be removed soon." << G4endl
//             << " Please try to use the new structured commands!" << G4endl;
//    }
//  else if(command == maxthetaCmd)
//    {
//      CHECKPG(); fParticleGun->GetAngDist()->SetMaxTheta(maxthetaCmd->GetNewDoubleValue(newValues));
//      G4cout << " G4GeneralParticleSourceMessenger - Warning:" << G4endl
//             << " The command is obsolete and will be removed soon." << G4endl
//             << " Please try to use the new structured commands!" << G4endl;
//    }
//  else if(command == maxphiCmd)
//    {
//      CHECKPG(); fParticleGun->GetAngDist()->SetMaxPhi(maxphiCmd->GetNewDoubleValue(newValues));
//      G4cout << " G4GeneralParticleSourceMessenger - Warning:" << G4endl
//             << " The command is obsolete and will be removed soon." << G4endl
//             << " Please try to use the new structured commands!" << G4endl;
//    }
//  else if(command == angsigmarCmd)
//    {
//      CHECKPG(); fParticleGun->GetAngDist()->SetBeamSigmaInAngR(angsigmarCmd->GetNewDoubleValue(newValues));
//      G4cout << " G4GeneralParticleSourceMessenger - Warning:" << G4endl
//             << " The command is obsolete and will be removed soon." << G4endl
//             << " Please try to use the new structured commands!" << G4endl;
//    }
//  else if(command == angsigmaxCmd)
//    {
//      CHECKPG(); fParticleGun->GetAngDist()->SetBeamSigmaInAngX(angsigmaxCmd->GetNewDoubleValue(newValues));
//      G4cout << " G4GeneralParticleSourceMessenger - Warning:" << G4endl
//             << " The command is obsolete and will be removed soon." << G4endl
//             << " Please try to use the new structured commands!" << G4endl;
//    }
//  else if(command == angsigmayCmd)
//    {
//      CHECKPG(); fParticleGun->GetAngDist()->SetBeamSigmaInAngY(angsigmayCmd->GetNewDoubleValue(newValues));
//      G4cout << " G4GeneralParticleSourceMessenger - Warning:" << G4endl
//             << " The command is obsolete and will be removed soon." << G4endl
//             << " Please try to use the new structured commands!" << G4endl;
//    }
//  else if(command == useuserangaxisCmd)
//    {
//      CHECKPG(); fParticleGun->GetAngDist()->SetUseUserAngAxis(useuserangaxisCmd->GetNewBoolValue(newValues));
//      G4cout << " G4GeneralParticleSourceMessenger - Warning:" << G4endl
//             << " The command is obsolete and will be removed soon." << G4endl
//             << " Please try to use the new structured commands!" << G4endl;
//    }
//  else if(command == surfnormCmd)
//    {
//      CHECKPG(); fParticleGun->GetAngDist()->SetUserWRTSurface(surfnormCmd->GetNewBoolValue(newValues));
//      G4cout << " G4GeneralParticleSourceMessenger - Warning:" << G4endl
//             << " The command is obsolete and will be removed soon." << G4endl
//             << " Please try to use the new structured commands!" << G4endl;
//    }
//  else if(command == energytypeCmd)
//    {
//      CHECKPG(); fParticleGun->GetEneDist()->SetEnergyDisType(newValues);
//      G4cout << " G4GeneralParticleSourceMessenger - Warning:" << G4endl
//             << " The command is obsolete and will be removed soon." << G4endl
//             << " Please try to use the new structured commands!" << G4endl;
//    }
//  else if(command == eminCmd)
//    {
//      CHECKPG(); fParticleGun->GetEneDist()->SetEmin(eminCmd->GetNewDoubleValue(newValues));
//      G4cout << " G4GeneralParticleSourceMessenger - Warning:" << G4endl
//             << " The command is obsolete and will be removed soon." << G4endl
//             << " Please try to use the new structured commands!" << G4endl;
//    }
//  else if(command == emaxCmd)
//    {
//      CHECKPG(); fParticleGun->GetEneDist()->SetEmax(emaxCmd->GetNewDoubleValue(newValues));
//      G4cout << " G4GeneralParticleSourceMessenger - Warning:" << G4endl
//             << " The command is obsolete and will be removed soon." << G4endl
//             << " Please try to use the new structured commands!" << G4endl;
//    }
//  else if(command == monoenergyCmd)
//    {
//      CHECKPG(); fParticleGun->GetEneDist()->SetMonoEnergy(monoenergyCmd->GetNewDoubleValue(newValues));
//      G4cout << " G4GeneralParticleSourceMessenger - Warning:" << G4endl
//             << " The command is obsolete and will be removed soon." << G4endl
//             << " Please try to use the new structured commands!" << G4endl;
//    }
//  else if(command == engsigmaCmd)
//    {
//      CHECKPG(); fParticleGun->GetEneDist()->SetBeamSigmaInE(engsigmaCmd->GetNewDoubleValue(newValues));
//      G4cout << " G4GeneralParticleSourceMessenger - Warning:" << G4endl
//             << " The command is obsolete and will be removed soon." << G4endl
//             << " Please try to use the new structured commands!" << G4endl;
//    }
//  else if(command == alphaCmd)
//    {
//      CHECKPG(); fParticleGun->GetEneDist()->SetAlpha(alphaCmd->GetNewDoubleValue(newValues));
//      G4cout << " G4GeneralParticleSourceMessenger - Warning:" << G4endl
//             << " The command is obsolete and will be removed soon." << G4endl
//             << " Please try to use the new structured commands!" << G4endl;
//    }
//  else if(command == tempCmd)
//    {
//      CHECKPG(); fParticleGun->GetEneDist()->SetTemp(tempCmd->GetNewDoubleValue(newValues));
//      G4cout << " G4GeneralParticleSourceMessenger - Warning:" << G4endl
//             << " The command is obsolete and will be removed soon." << G4endl
//             << " Please try to use the new structured commands!" << G4endl;
//    }
//  else if(command == ezeroCmd)
//    {
//      CHECKPG(); fParticleGun->GetEneDist()->SetEzero(ezeroCmd->GetNewDoubleValue(newValues));
//      G4cout << " G4GeneralParticleSourceMessenger - Warning:" << G4endl
//             << " The command is obsolete and will be removed soon." << G4endl
//             << " Please try to use the new structured commands!" << G4endl;
//    }
//  else if(command == gradientCmd)
//    {
//      CHECKPG(); fParticleGun->GetEneDist()->SetGradient(gradientCmd->GetNewDoubleValue(newValues));
//      G4cout << " G4GeneralParticleSourceMessenger - Warning:" << G4endl
//             << " The command is obsolete and will be removed soon." << G4endl
//             << " Please try to use the new structured commands!" << G4endl;
//    }
//  else if(command == interceptCmd)
//    {
//      CHECKPG(); fParticleGun->GetEneDist()->SetInterCept(interceptCmd->GetNewDoubleValue(newValues));
//      G4cout << " G4GeneralParticleSourceMessenger - Warning:" << G4endl
//             << " The command is obsolete and will be removed soon." << G4endl
//             << " Please try to use the new structured commands!" << G4endl;
//    }
//  else if(command == calculateCmd)
//    {
//      CHECKPG(); fParticleGun->GetEneDist()->Calculate();
//      G4cout << " G4GeneralParticleSourceMessenger - Warning:" << G4endl
//             << " The command is obsolete and will be removed soon." << G4endl
//             << " Please try to use the new structured commands!" << G4endl;
//    }
//  else if(command == energyspecCmd)
//    {
//      CHECKPG(); fParticleGun->GetEneDist()->InputEnergySpectra(energyspecCmd->GetNewBoolValue(newValues));
//      G4cout << " G4GeneralParticleSourceMessenger - Warning:" << G4endl
//             << " The command is obsolete and will be removed soon." << G4endl
//             << " Please try to use the new structured commands!" << G4endl;
//    }
//  else if(command == diffspecCmd)
//    {
//      CHECKPG(); fParticleGun->GetEneDist()->InputDifferentialSpectra(diffspecCmd->GetNewBoolValue(newValues));
//      G4cout << " G4GeneralParticleSourceMessenger - Warning:" << G4endl
//             << " The command is obsolete and will be removed soon." << G4endl
//             << " Please try to use the new structured commands!" << G4endl;
//    }
//  else if(command == histnameCmd)
//    {
//      histtype = newValues;
//      G4cout << " G4GeneralParticleSourceMessenger - Warning:" << G4endl
//             << " The command is obsolete and will be removed soon." << G4endl
//             << " Please try to use the new structured commands!" << G4endl;
//    }
//  else
//      if(command == histpointCmd)
//    {
//      CHECKPG();
//      if(histtype == "biasx")
//	fParticleGun->GetBiasRndm()->SetXBias(histpointCmd->GetNew3VectorValue(newValues));
//      if(histtype == "biasy")
//	fParticleGun->GetBiasRndm()->SetYBias(histpointCmd->GetNew3VectorValue(newValues));
//      if(histtype == "biasz")
//	fParticleGun->GetBiasRndm()->SetZBias(histpointCmd->GetNew3VectorValue(newValues));
//      if(histtype == "biast")
//	fParticleGun->GetBiasRndm()->SetThetaBias(histpointCmd->GetNew3VectorValue(newValues));
//      if(histtype == "biasp")
//	fParticleGun->GetBiasRndm()->SetPhiBias(histpointCmd->GetNew3VectorValue(newValues));
//      if(histtype == "biase")
//	fParticleGun->GetBiasRndm()->SetEnergyBias(histpointCmd->GetNew3VectorValue(newValues));
//      if(histtype == "theta")
//	fParticleGun->GetAngDist()->UserDefAngTheta(histpointCmd->GetNew3VectorValue(newValues));
//      if(histtype == "phi")
//	fParticleGun->GetAngDist()->UserDefAngPhi(histpointCmd->GetNew3VectorValue(newValues));
//      if(histtype == "energy")
//	fParticleGun->GetEneDist()->UserEnergyHisto(histpointCmd->GetNew3VectorValue(newValues));
//      if(histtype == "arb")
//	fParticleGun->GetEneDist()->ArbEnergyHisto(histpointCmd->GetNew3VectorValue(newValues));
//      if(histtype == "epn")
//	fParticleGun->GetEneDist()->EpnEnergyHisto(histpointCmd->GetNew3VectorValue(newValues));
//      G4cout << " G4GeneralParticleSourceMessenger - Warning: The command is obsolete and will be removed soon. Please try to use the new structured commands!" << G4endl;
//    }
//  else if(command == resethistCmd)
//    {
//      CHECKPG();
//      if(newValues == "theta" || newValues == "phi") {
//	fParticleGun->GetAngDist()->ReSetHist(newValues);
//      } else if (newValues == "energy" || newValues == "arb" || newValues == "epn") {
//	fParticleGun->GetEneDist()->ReSetHist(newValues);
//      } else {
//	fParticleGun->GetBiasRndm()->ReSetHist(newValues);
//      }
//      G4cout << " G4GeneralParticleSourceMessenger - Warning:" << G4endl
//             << " The command is obsolete and will be removed soon." << G4endl
//             << " Please try to use the new structured commands!" << G4endl;
//    }
//  else if(command == arbintCmd)
//    {
//      CHECKPG();
//      fParticleGun->GetEneDist()->ArbInterpolate(newValues);
//      G4cout << " G4GeneralParticleSourceMessenger - Warning:" << G4endl
//             << " The command is obsolete and will be removed soon." << G4endl
//             << " Please try to use the new structured commands!" << G4endl;
//    }
// else
	if( command==directionCmd )
    {
      CHECKPG();
      fParticleGun->GetAngDist()->SetAngDistType("planar");
      fParticleGun->GetAngDist()->SetParticleMomentumDirection(directionCmd->GetNew3VectorValue(newValues));
    }
  else if( command==energyCmd )
    {
      CHECKPG();
      fParticleGun->GetEneDist()->SetEnergyDisType("Mono");
      fParticleGun->GetEneDist()->SetMonoEnergy(energyCmd->GetNewDoubleValue(newValues));
    }
  else if( command==positionCmd )
    {
      CHECKPG();
      fParticleGun->GetPosDist()->SetPosDisType("Point");    
      fParticleGun->GetPosDist()->SetCentreCoords(positionCmd->GetNew3VectorValue(newValues));
    }
  else if(command == verbosityCmd)
    {
	  fGPS->SetVerbosity(verbosityCmd->GetNewIntValue(newValues));
      //CHECKPG();
      //fParticleGun->SetVerbosity(verbosityCmd->GetNewIntValue(newValues));
    }
  else if( command==particleCmd )
    {
      if (newValues =="ion") {
	fShootIon = true;
      } else {
	fShootIon = false;
	G4ParticleDefinition* pd = particleTable->FindParticle(newValues);
	if(pd != NULL)
	  { CHECKPG(); fParticleGun->SetParticleDefinition( pd ); }
      }
    }
  else if( command==timeCmd )
    { CHECKPG(); fParticleGun->SetParticleTime(timeCmd->GetNewDoubleValue(newValues)); }
  else if( command==polCmd )
    { CHECKPG(); fParticleGun->SetParticlePolarization(polCmd->GetNew3VectorValue(newValues)); }
  else if( command==numberCmd )
    { CHECKPG(); fParticleGun->SetNumberOfParticles(numberCmd->GetNewIntValue(newValues)); }
  else if( command==ionCmd )
    { IonCommand(newValues); }
  else if( command==ionLvlCmd )
    { IonLvlCommand(newValues); }
  else if( command==listCmd ){ 
    particleTable->DumpTable(); 
  }
  else if( command==addsourceCmd )
    {
      fGPS->AddaSource(addsourceCmd->GetNewDoubleValue(newValues));
    }
  else if( command==listsourceCmd )
    { 
      fGPS->ListSource();
    }
  else if( command==clearsourceCmd )
    { 
      fGPS->ClearAll();
      fParticleGun = 0;      
    }
  else if( command==getsourceCmd )
    { 
      G4cout << " Current source index:" << fGPS->GetCurrentSourceIndex() 
	     << " ; Intensity:" << fGPS->GetCurrentSourceIntensity() << G4endl;
    }
  else if( command==setsourceCmd )
    { 
	  //NOTE: This will also sets fParticleGun to the courrent source
	  //      Not very clean, the GPS::SetCurrentSourceto will call:
	  //      this::SetParticleSource( G4ParticleSource* )
	  //      The point is that GPS has no public API to get a source by
	  //      index
	  //TODO: Can we add this API?
	  const G4int sn = setsourceCmd->GetNewIntValue(newValues);
	  if ( sn >= fGPS->GetNumberofSource() )
	  {
		G4ExceptionDescription msg;
		msg << "Using command "<<setsourceCmd->GetCommandPath()<<"/"<<setsourceCmd->GetCommandName()<<" "<<sn;
		msg << " is invalid "<<fGPS->GetNumberofSource()<<" source(s) are defined.";
		G4Exception("G4GeneralParticleSourceMessenger::SetNewValue","G4GPS005",FatalException,msg);
	  }
	  fGPS->SetCurrentSourceto(setsourceCmd->GetNewIntValue(newValues));

    }
  else if( command==setintensityCmd )
    { 
      fGPS->SetCurrentSourceIntensity(setintensityCmd->GetNewDoubleValue(newValues));
    }
  else if( command==deletesourceCmd )
    { 
      fGPS->DeleteaSource(deletesourceCmd->GetNewIntValue(newValues));
    }
  else if(command == multiplevertexCmd)
    {
      fGPS->SetMultipleVertex(multiplevertexCmd->GetNewBoolValue(newValues));
    }
  else if(command == flatsamplingCmd)
    {
      fGPS->SetFlatSampling(flatsamplingCmd->GetNewBoolValue(newValues));
    }
  //
  // new implementations
  //
  //
  else if(command == typeCmd1)
    {
      CHECKPG();
      fParticleGun->GetPosDist()->SetPosDisType(newValues);
    }
  else if(command == shapeCmd1)
    {
      CHECKPG();
      fParticleGun->GetPosDist()->SetPosDisShape(newValues);
    }
  else if(command == centreCmd1)
    {
      CHECKPG();
      fParticleGun->GetPosDist()->SetCentreCoords(centreCmd1->GetNew3VectorValue(newValues));
    }
  else if(command == posrot1Cmd1)
    {
      CHECKPG();
      fParticleGun->GetPosDist()->SetPosRot1(posrot1Cmd1->GetNew3VectorValue(newValues));
    }
  else if(command == posrot2Cmd1)
    {
      CHECKPG();
      fParticleGun->GetPosDist()->SetPosRot2(posrot2Cmd1->GetNew3VectorValue(newValues));
    }
  else if(command == halfxCmd1)
    {
      CHECKPG();
      fParticleGun->GetPosDist()->SetHalfX(halfxCmd1->GetNewDoubleValue(newValues));
    }
  else if(command == halfyCmd1)
    {
      CHECKPG();
      fParticleGun->GetPosDist()->SetHalfY(halfyCmd1->GetNewDoubleValue(newValues));
    }
  else if(command == halfzCmd1)
    {
      CHECKPG();
      fParticleGun->GetPosDist()->SetHalfZ(halfzCmd1->GetNewDoubleValue(newValues));
    }
  else if(command == radiusCmd1)
    {
      CHECKPG();
      fParticleGun->GetPosDist()->SetRadius(radiusCmd1->GetNewDoubleValue(newValues));
    }
  else if(command == radius0Cmd1)
    {
      CHECKPG();
      fParticleGun->GetPosDist()->SetRadius0(radius0Cmd1->GetNewDoubleValue(newValues));
    }
  else if(command == possigmarCmd1)
    {
      CHECKPG();
      fParticleGun->GetPosDist()->SetBeamSigmaInR(possigmarCmd1->GetNewDoubleValue(newValues));
    }
  else if(command == possigmaxCmd1)
    {
      CHECKPG(); fParticleGun->GetPosDist()->SetBeamSigmaInX(possigmaxCmd1->GetNewDoubleValue(newValues));
    }
  else if(command == possigmayCmd1)
    {
      CHECKPG(); fParticleGun->GetPosDist()->SetBeamSigmaInY(possigmayCmd1->GetNewDoubleValue(newValues));
    }
  else if(command == paralpCmd1)
    {
      CHECKPG(); fParticleGun->GetPosDist()->SetParAlpha(paralpCmd1->GetNewDoubleValue(newValues));
    }
  else if(command == partheCmd1)
    {
      CHECKPG(); fParticleGun->GetPosDist()->SetParTheta(partheCmd1->GetNewDoubleValue(newValues));
    }
  else if(command == parphiCmd1)
    {
      CHECKPG(); fParticleGun->GetPosDist()->SetParPhi(parphiCmd1->GetNewDoubleValue(newValues));
    }
  else if(command == confineCmd1)
    {
      CHECKPG(); fParticleGun->GetPosDist()->ConfineSourceToVolume(newValues);
    }
  else if(command == angtypeCmd1)
    {
      CHECKPG(); fParticleGun->GetAngDist()->SetAngDistType(newValues);
    }
  else if(command == angrot1Cmd1)
    {
      CHECKPG(); 
      G4String a = "angref1";
      fParticleGun->GetAngDist()->DefineAngRefAxes(a,angrot1Cmd1->GetNew3VectorValue(newValues));
    }
  else if(command == angrot2Cmd1)
    {
      CHECKPG(); 
      G4String a = "angref2";
      fParticleGun->GetAngDist()->DefineAngRefAxes(a,angrot2Cmd1->GetNew3VectorValue(newValues));
    }
  else if(command == minthetaCmd1)
    {
      CHECKPG(); 
      fParticleGun->GetAngDist()->SetMinTheta(minthetaCmd1->GetNewDoubleValue(newValues));
    }
  else if(command == minphiCmd1)
    {
      CHECKPG(); fParticleGun->GetAngDist()->SetMinPhi(minphiCmd1->GetNewDoubleValue(newValues));
    }
  else if(command == maxthetaCmd1)
    {
      CHECKPG(); fParticleGun->GetAngDist()->SetMaxTheta(maxthetaCmd1->GetNewDoubleValue(newValues));
    }
  else if(command == maxphiCmd1)
    {
      CHECKPG(); fParticleGun->GetAngDist()->SetMaxPhi(maxphiCmd1->GetNewDoubleValue(newValues));
    }
  else if(command == angsigmarCmd1)
    {
      CHECKPG(); fParticleGun->GetAngDist()->SetBeamSigmaInAngR(angsigmarCmd1->GetNewDoubleValue(newValues));
    }
  else if(command == angsigmaxCmd1)
    {
      CHECKPG(); fParticleGun->GetAngDist()->SetBeamSigmaInAngX(angsigmaxCmd1->GetNewDoubleValue(newValues));
    }
  else if(command == angsigmayCmd1)
    {
      CHECKPG(); fParticleGun->GetAngDist()->SetBeamSigmaInAngY(angsigmayCmd1->GetNewDoubleValue(newValues));
    }
  else if(command == angfocusCmd)
    {
      CHECKPG(); fParticleGun->GetAngDist()->SetFocusPoint(angfocusCmd->GetNew3VectorValue(newValues));
    }
  else if(command == useuserangaxisCmd1)
    {
      CHECKPG(); fParticleGun->GetAngDist()->SetUseUserAngAxis(useuserangaxisCmd1->GetNewBoolValue(newValues));
    }
  else if(command == surfnormCmd1)
    {
      CHECKPG(); fParticleGun->GetAngDist()->SetUserWRTSurface(surfnormCmd1->GetNewBoolValue(newValues));
    }
  else if(command == energytypeCmd1)
    {
      CHECKPG(); fParticleGun->GetEneDist()->SetEnergyDisType(newValues);
    }
  else if(command == eminCmd1)
    {
      CHECKPG(); fParticleGun->GetEneDist()->SetEmin(eminCmd1->GetNewDoubleValue(newValues));
    }
  else if(command == emaxCmd1)
    {
      CHECKPG(); fParticleGun->GetEneDist()->SetEmax(emaxCmd1->GetNewDoubleValue(newValues));
    }
  else if(command == monoenergyCmd1)
    {
      CHECKPG(); fParticleGun->GetEneDist()->SetMonoEnergy(monoenergyCmd1->GetNewDoubleValue(newValues));
    }
  else if(command == engsigmaCmd1)
    {
      CHECKPG(); fParticleGun->GetEneDist()->SetBeamSigmaInE(engsigmaCmd1->GetNewDoubleValue(newValues));
    }
  else if(command == alphaCmd1)
    {
      CHECKPG(); fParticleGun->GetEneDist()->SetAlpha(alphaCmd1->GetNewDoubleValue(newValues));
    }
  else if(command == tempCmd1)
    {
      CHECKPG(); fParticleGun->GetEneDist()->SetTemp(tempCmd1->GetNewDoubleValue(newValues));
    }
  else if(command == ezeroCmd1)
    {
      CHECKPG(); fParticleGun->GetEneDist()->SetEzero(ezeroCmd1->GetNewDoubleValue(newValues));
    }
  else if(command == gradientCmd1)
    {
      CHECKPG(); fParticleGun->GetEneDist()->SetGradient(gradientCmd1->GetNewDoubleValue(newValues));
    }
  else if(command == interceptCmd1)
    {
      CHECKPG(); fParticleGun->GetEneDist()->SetInterCept(interceptCmd1->GetNewDoubleValue(newValues));
    }
  else if(command == arbeintCmd1)
    {
      CHECKPG(); fParticleGun->GetEneDist()->SetBiasAlpha(arbeintCmd1->GetNewDoubleValue(newValues));
    }
  else if(command == calculateCmd1)
    {
      CHECKPG(); fParticleGun->GetEneDist()->Calculate();
    }
  else if(command == energyspecCmd1)
    {
      CHECKPG(); fParticleGun->GetEneDist()->InputEnergySpectra(energyspecCmd1->GetNewBoolValue(newValues));
    }
  else if(command == diffspecCmd1)
    {
      CHECKPG(); fParticleGun->GetEneDist()->InputDifferentialSpectra(diffspecCmd1->GetNewBoolValue(newValues));
    }
  else if(command == histnameCmd1)
    {
      histtype = newValues;
    }
  else if(command == histfileCmd1)
    {
      histtype = "arb";
      CHECKPG(); fParticleGun->GetEneDist()->ArbEnergyHistoFile(newValues);
    }
  else if(command == histpointCmd1)
    {
      CHECKPG(); 
      if(histtype == "biasx")
	fParticleGun->GetBiasRndm()->SetXBias(histpointCmd1->GetNew3VectorValue(newValues));
      if(histtype == "biasy")
	fParticleGun->GetBiasRndm()->SetYBias(histpointCmd1->GetNew3VectorValue(newValues));
      if(histtype == "biasz")
	fParticleGun->GetBiasRndm()->SetZBias(histpointCmd1->GetNew3VectorValue(newValues));
      if(histtype == "biast")
	fParticleGun->GetBiasRndm()->SetThetaBias(histpointCmd1->GetNew3VectorValue(newValues));
      if(histtype == "biasp")
	fParticleGun->GetBiasRndm()->SetPhiBias(histpointCmd1->GetNew3VectorValue(newValues));
      if(histtype == "biaspt")
	fParticleGun->GetBiasRndm()->SetPosThetaBias(histpointCmd1->GetNew3VectorValue(newValues));
      if(histtype == "biaspp")
	fParticleGun->GetBiasRndm()->SetPosPhiBias(histpointCmd1->GetNew3VectorValue(newValues));
      if(histtype == "biase")
	fParticleGun->GetBiasRndm()->SetEnergyBias(histpointCmd1->GetNew3VectorValue(newValues));
      if(histtype == "theta")
	fParticleGun->GetAngDist()->UserDefAngTheta(histpointCmd1->GetNew3VectorValue(newValues));
      if(histtype == "phi")
	fParticleGun->GetAngDist()->UserDefAngPhi(histpointCmd1->GetNew3VectorValue(newValues));
      if(histtype == "energy")
	fParticleGun->GetEneDist()->UserEnergyHisto(histpointCmd1->GetNew3VectorValue(newValues));
      if(histtype == "arb")
	fParticleGun->GetEneDist()->ArbEnergyHisto(histpointCmd1->GetNew3VectorValue(newValues));
      if(histtype == "epn")
	fParticleGun->GetEneDist()->EpnEnergyHisto(histpointCmd1->GetNew3VectorValue(newValues));
    }
  else if(command == resethistCmd1)
    {
      CHECKPG(); 
      if(newValues == "theta" || newValues == "phi") {
	fParticleGun->GetAngDist()->ReSetHist(newValues);
      } else if (newValues == "energy" || newValues == "arb" || newValues == "epn") {
	fParticleGun->GetEneDist()->ReSetHist(newValues);
      } else {
	fParticleGun->GetBiasRndm()->ReSetHist(newValues);
      }
    }
  else if(command == arbintCmd1)
    {
      CHECKPG(); fParticleGun->GetEneDist()->ArbInterpolate(newValues);
    }
  else 
    {
      G4cout << "Error entering command" << G4endl;
    }
}

G4String G4GeneralParticleSourceMessenger::GetCurrentValue(G4UIcommand *)
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

void G4GeneralParticleSourceMessenger::IonCommand(G4String newValues)
{
  if (fShootIon)
  {
    G4Tokenizer next( newValues );
    // check argument
    fAtomicNumber = StoI(next());
    fAtomicMass = StoI(next());
    G4String sQ = next();
    if (sQ.isNull())
    {
	fIonCharge = fAtomicNumber;
    }
    else
    {
	fIonCharge = StoI(sQ);
	sQ = next();
	if (sQ.isNull())
        {
	  fIonExciteEnergy = 0.0;
        }
        else
        {
	  fIonExciteEnergy = StoD(sQ) * keV;
        }
    }
    G4ParticleDefinition* ion;
    ion =  G4IonTable::GetIonTable()->GetIon( fAtomicNumber, fAtomicMass, fIonExciteEnergy);
    if (ion==0)
    {
      G4ExceptionDescription ed;
      ed << "Ion with Z=" << fAtomicNumber;
      ed << " A=" << fAtomicMass << " is not defined";    
      ionCmd->CommandFailed(ed);
    }
    else
    {
      fParticleGun->SetParticleDefinition(ion);
      fParticleGun->SetParticleCharge(fIonCharge*eplus);
    }
  }
  else
  {
    G4ExceptionDescription ed;
    ed << "Set /gps/particle to ion before using /gps/ion command";
    ionCmd->CommandFailed(ed);
  }
}

void G4GeneralParticleSourceMessenger::IonLvlCommand(G4String newValues)
{
  if (fShootIon)
  {
    G4Tokenizer next(newValues);
    // check argument
    fAtomicNumberL = StoI(next());
    fAtomicMassL = StoI(next());
    G4String sQ = next();
    if (sQ.isNull()) {
      fIonChargeL = fAtomicNumberL;
    } else {
      fIonChargeL = StoI(sQ);
      sQ = next();
      if (sQ.isNull()) {
        fIonEnergyLevel = 0;
      } else {
        fIonEnergyLevel = StoI(sQ);
      }
    }

    G4ParticleDefinition* ion;
    ion =  G4IonTable::GetIonTable()->GetIon(fAtomicNumberL, fAtomicMassL, fIonEnergyLevel);
    if (ion == 0) {
      G4ExceptionDescription ed;
      ed << "Ion with Z=" << fAtomicNumberL;
      ed << " A=" << fAtomicMassL << " is not defined";
      ionLvlCmd->CommandFailed(ed);
    } else {
      fParticleGun->SetParticleDefinition(ion);
      fParticleGun->SetParticleCharge(fIonChargeL*eplus);
    }

  } else {
    G4ExceptionDescription ed;
    ed << "Set /gps/particle to ion before using /gps/ionLvl command";
    ionLvlCmd->CommandFailed(ed);
  }
}

