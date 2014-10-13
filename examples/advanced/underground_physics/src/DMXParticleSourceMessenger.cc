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
// --------------------------------------------------------------
//   GEANT 4 - Underground Dark Matter Detector Advanced Example
//
//      For information related to this code contact: Alex Howard
//      e-mail: alexander.howard@cern.ch
// --------------------------------------------------------------
// Comments
//
//                  Underground Advanced
//               by A. Howard and H. Araujo 
//                    (27th November 2001)
//
// ParticleSourceMessenger program
// --------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////////
// This particle source is a shortened version of G4GeneralParticleSource by
// C Ferguson, F Lei & P Truscott (University of Southampton / DERA), with
// some minor modifications.
//////////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <iomanip>               

#include "DMXParticleSourceMessenger.hh"
#include "DMXParticleSource.hh"

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

///////////////////////////////////////////////////////////////////////////////
DMXParticleSourceMessenger::DMXParticleSourceMessenger
(DMXParticleSource *fPtclGun) : fParticleGun(fPtclGun),fShootIon(false) {

  particleTable = G4ParticleTable::GetParticleTable();

  // create directory
  gunDirectory = new G4UIdirectory("/dmx/gun/");
  gunDirectory->SetGuidance("Particle Source control commands.");

  // list available particles
  listCmd = new G4UIcmdWithoutParameter("/dmx/gun/List",this);
  listCmd->SetGuidance("List available particles.");
  listCmd->SetGuidance(" Invoke G4ParticleTable.");

  // set particle  
  particleCmd = new G4UIcmdWithAString("/dmx/gun/particle",this);
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
  

  // particle direction
  directionCmd = new G4UIcmdWith3Vector("/dmx/gun/direction",this);
  directionCmd->SetGuidance("Set momentum direction.");
  directionCmd->SetGuidance("Direction needs not to be a unit vector.");
  directionCmd->SetParameterName("Px","Py","Pz",true,true); 
  directionCmd->SetRange("Px != 0 || Py != 0 || Pz != 0");
  
  // particle energy
  energyCmd = new G4UIcmdWithADoubleAndUnit("/dmx/gun/energy",this);
  energyCmd->SetGuidance("Set kinetic energy.");
  energyCmd->SetParameterName("Energy",true,true);
  energyCmd->SetDefaultUnit("GeV");
  //energyCmd->SetUnitCategory("Energy");
  //energyCmd->SetUnitCandidates("eV keV MeV GeV TeV");

  positionCmd = new G4UIcmdWith3VectorAndUnit("/dmx/gun/position",this);
  positionCmd->SetGuidance("Set starting position of the particle.");
  positionCmd->SetParameterName("X","Y","Z",true,true);
  positionCmd->SetDefaultUnit("cm");
  //positionCmd->SetUnitCategory("Length");
  //positionCmd->SetUnitCandidates("microm mm cm m km");

 
  // ion 
  ionCmd = new G4UIcommand("/dmx/gun/ion",this);
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
  

  // source distribution type
  typeCmd = new G4UIcmdWithAString("/dmx/gun/type",this);
  typeCmd->SetGuidance("Sets source distribution type.");
  typeCmd->SetGuidance("Either Point or Volume");
  typeCmd->SetParameterName("DisType",true,true);
  typeCmd->SetDefaultValue("Point");
  typeCmd->SetCandidates("Point Volume");
  
  // source shape
  shapeCmd = new G4UIcmdWithAString("/dmx/gun/shape",this);
  shapeCmd->SetGuidance("Sets source shape type.");
  shapeCmd->SetParameterName("Shape",true,true);
  shapeCmd->SetDefaultValue("NULL");
  shapeCmd->SetCandidates("Sphere Cylinder");
  
  // centre coordinates
  centreCmd = new G4UIcmdWith3VectorAndUnit("/dmx/gun/centre",this);
  centreCmd->SetGuidance("Set centre coordinates of source.");
  centreCmd->SetParameterName("X","Y","Z",true,true);
  centreCmd->SetDefaultUnit("cm");
  centreCmd->SetUnitCandidates("nm um mm cm m km");

  // half height of source
  halfzCmd = new G4UIcmdWithADoubleAndUnit("/dmx/gun/halfz",this);
  halfzCmd->SetGuidance("Set z half length of source.");
  halfzCmd->SetParameterName("Halfz",true,true);
  halfzCmd->SetDefaultUnit("cm");
  halfzCmd->SetUnitCandidates("nm um mm cm m km");

  // radius of source  
  radiusCmd = new G4UIcmdWithADoubleAndUnit("/dmx/gun/radius",this);
  radiusCmd->SetGuidance("Set radius of source.");
  radiusCmd->SetParameterName("Radius",true,true);
  radiusCmd->SetDefaultUnit("cm");
  radiusCmd->SetUnitCandidates("nm um mm cm m km");
  
  // confine to volume
  confineCmd = new G4UIcmdWithAString("/dmx/gun/confine",this);
  confineCmd->SetGuidance("Confine source to volume (NULL to unset).");
  confineCmd->SetGuidance("usage: confine VolName");
  confineCmd->SetParameterName("VolName",true,true);
  confineCmd->SetDefaultValue("NULL");
  
  // angular distribution
  angtypeCmd = new G4UIcmdWithAString("/dmx/gun/angtype",this);
  angtypeCmd->SetGuidance("Sets angular source distribution type");
  angtypeCmd->SetGuidance("Possible variables are: iso direction");
  angtypeCmd->SetParameterName("AngDis",true,true);
  angtypeCmd->SetDefaultValue("iso");
  angtypeCmd->SetCandidates("iso direction");
  
  // energy distribution
  energytypeCmd = new G4UIcmdWithAString("/dmx/gun/energytype",this);
  energytypeCmd->SetGuidance("Sets energy distribution type");
  energytypeCmd->SetGuidance("Possible variables are: Mono");
  energytypeCmd->SetParameterName("EnergyDis",true,true);
  energytypeCmd->SetDefaultValue("Mono");
  energytypeCmd->SetCandidates("Mono");

  // verbosity
  verbosityCmd = new G4UIcmdWithAnInteger("/dmx/gun/verbose",this);
  verbosityCmd->SetGuidance("Set Verbose level for gun");
  verbosityCmd->SetGuidance(" 0 : Silent");
  verbosityCmd->SetGuidance(" 1 : Limited information");
  verbosityCmd->SetGuidance(" 2 : Detailed information");
  verbosityCmd->SetParameterName("level",false);
  verbosityCmd->SetRange("level>=0 && level <=2");
  
}


DMXParticleSourceMessenger::~DMXParticleSourceMessenger() {

  delete typeCmd;
  delete shapeCmd;
  delete centreCmd;
  delete halfzCmd;
  delete radiusCmd;
  delete confineCmd;
  delete angtypeCmd;
  delete energytypeCmd;
  delete verbosityCmd;
  delete ionCmd;
  delete particleCmd;
  delete positionCmd;
  delete directionCmd;
  delete energyCmd;
  delete listCmd;

  delete gunDirectory;
}

void DMXParticleSourceMessenger::SetNewValue
   (G4UIcommand *command, G4String newValues) {

  if(command == typeCmd)
    fParticleGun->SetPosDisType(newValues);

  else if(command == shapeCmd)
    fParticleGun->SetPosDisShape(newValues);

  else if(command == centreCmd)
    fParticleGun->SetCentreCoords(centreCmd->GetNew3VectorValue(newValues));

  else if(command == halfzCmd)
    fParticleGun->SetHalfZ(halfzCmd->GetNewDoubleValue(newValues));

  else if(command == radiusCmd)
    fParticleGun->SetRadius(radiusCmd->GetNewDoubleValue(newValues));

  else if(command == angtypeCmd)
      fParticleGun->SetAngDistType(newValues);

  else if(command == confineCmd)
    fParticleGun->ConfineSourceToVolume(newValues);

  else if(command == energytypeCmd)
    fParticleGun->SetEnergyDisType(newValues);
  
  else if(command == verbosityCmd)
    fParticleGun->SetVerbosity(verbosityCmd->GetNewIntValue(newValues));

  else if( command==particleCmd ) {
    if (newValues =="ion") {
      fShootIon = true;
    } else {
      fShootIon = false;
      G4ParticleDefinition* pd = particleTable->FindParticle(newValues);
      if(pd != NULL)
  	{ fParticleGun->SetParticleDefinition( pd ); }
    }
  }

  else if( command==ionCmd ) {
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
      ion =  G4IonTable::GetIonTable()->GetIon(fAtomicNumber,fAtomicMass,fIonExciteEnergy);
      if (ion==0) {
	G4cout << "Ion with Z=" << fAtomicNumber;
	G4cout << " A=" << fAtomicMass << "is not be defined" << G4endl;    
      } else {
	fParticleGun->SetParticleDefinition(ion);
	fParticleGun->SetParticleCharge(fIonCharge*eplus);
      }
    } else {
      G4cout<<"Set /dmx/gun/particle to ion before using /dmx/gun/ion command";
      G4cout<<G4endl; 
    }
  }

  else if( command==listCmd )
    particleTable->DumpTable();

  else if( command==directionCmd ) { 
    fParticleGun->SetAngDistType("direction");
    fParticleGun->SetParticleMomentumDirection
      (directionCmd->GetNew3VectorValue(newValues));
  }

  else if( command==energyCmd ) {
    fParticleGun->SetEnergyDisType("Mono");
    fParticleGun->SetMonoEnergy(energyCmd->GetNewDoubleValue(newValues));
  }  

  else if( command==positionCmd ) { 
    fParticleGun->SetPosDisType("Point");    
    fParticleGun->SetCentreCoords(positionCmd->GetNew3VectorValue(newValues));
  }
  else
    G4cout << "Error entering command" << G4endl;
}




