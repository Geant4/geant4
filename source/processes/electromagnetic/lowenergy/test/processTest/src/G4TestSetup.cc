//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4TestSetup.cc,v 1.1 2001-10-15 13:03:12 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 07 Oct 2001   MGP        Created
//
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4TestSetup.hh"
#include "G4VProcess.hh"
#include "G4ios.hh"
#include "g4std/fstream"
#include "g4std/iomanip"

#include "G4VContinuousDiscreteProcess.hh"
#include "G4ProcessManager.hh"
#include "G4VProcess.hh"

#include "G4LowEnergyPhotoElectric.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4LowEnergyGammaConversion.hh"
#include "G4GammaConversion.hh"

#include "G4LowEnergyRayleigh.hh"

#include "G4LowEnergyCompton.hh"
#include "G4LowEnergyPolarizedCompton.hh"
#include "G4ComptonScattering.hh"

#include "G4LowEnergyBremsstrahlung.hh"
#include "G4eBremsstrahlung.hh"

#include "G4LowEnergyIonisation.hh"
#include "G4eIonisation.hh"

#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4GRSVolume.hh"
#include "G4UnitsTable.hh"

#include "G4EnergyLossTables.hh"

#include "G4VParticleChange.hh"
#include "G4ParticleChange.hh"
#include "G4DynamicParticle.hh"

#include "Randomize.hh"


G4TestSetup::G4TestSetup()
{
  ioniProcess = 0;
  bremProcess = 0;
  ioniSel = false;
  bremSel = false;
  eProcessManager = 0;
  positronProcessManager = 0;
  processManager = 0;
  track = 0;
  step = new G4Step;
}

G4TestSetup:: ~G4TestSetup()
{
  if (! ioniSel) delete ioniProcess;
  if (! bremSel) delete bremProcess;
  delete eProcessManager;
  delete positronProcessManager;
  delete processManager;
  delete physicalFrame;
  physicalFrame = 0;
  //  delete step;
  step = 0;
  // delete track;
  track = 0;
}
 
G4VProcess* G4TestSetup::createTestProcess() 
{
  G4VProcess* process = 0;

  if (processType == 1)
    {
      if (selection == 1) process = new G4LowEnergyCompton;
      if (selection == 2) process = new G4ComptonScattering;
      if (selection == 3) process = new G4LowEnergyPolarizedCompton;
    }
     if (processType == 2)
    {
      if (selection == 1) process = new G4LowEnergyGammaConversion;
      if (selection == 2) process = new G4GammaConversion;
    }
   if (processType == 3)
    {
      if (selection == 1) process = new G4LowEnergyPhotoElectric;
      if (selection == 2) process = new G4PhotoElectricEffect;
    }
   if (processType == 4)
    {
      if (selection == 1) process = new G4LowEnergyRayleigh;
    }

   if (processType == 5)
    {
      bremSel = true;
      process = bremProcess;
    }

  if (processType == 6)
    {
      ioniSel = true;
      process = ioniProcess;
    }
    
  if (process == 0) G4Exception("The selected process is not available");

  if (processType < 5)
    {
      G4ParticleDefinition* gamma = G4Gamma::GammaDefinition();
      G4ProcessManager* processManager = new G4ProcessManager(gamma);
      gamma->SetProcessManager(processManager);
      processManager->AddProcess(process);
      process->BuildPhysicsTable(*gamma);
    }

  G4cout << "The selected process is " << process->GetProcessName() << G4endl;
  
  return process;
}

void G4TestSetup::createElectronProcesses()
{
  G4ParticleDefinition* gamma = G4Gamma::GammaDefinition();
  G4ParticleDefinition* electron = G4Electron::ElectronDefinition();
  G4ParticleDefinition* positron = G4Positron::PositronDefinition();

  gamma->SetCuts(1e-3*mm);
  electron->SetCuts(1e-3*mm);
  positron->SetCuts(1e-3*mm);
  
  if (selection == 1 || selection == 3)
    {
      bremProcess = new G4LowEnergyBremsstrahlung;
      ioniProcess = new G4LowEnergyIonisation;
    }
  else
    {
      bremProcess = new G4eBremsstrahlung;
      ioniProcess = new G4eIonisation;
    }
  // Initialize the physics tables 


  G4ProcessManager* eProcessManager = new G4ProcessManager(electron);
  electron->SetProcessManager(eProcessManager);
  eProcessManager->AddProcess(bremProcess);  
  eProcessManager->AddProcess(ioniProcess);  

  G4ProcessManager* positronProcessManager = new G4ProcessManager(positron);
  positron->SetProcessManager(positronProcessManager);
  positronProcessManager->AddProcess(bremProcess);
  
  bremProcess->BuildPhysicsTable(*electron);
  ioniProcess->BuildPhysicsTable(*electron);
}

  void G4TestSetup::makeGeometry()
{
  G4double dimX = 10.*cm;
  G4double dimY = 10.*cm;
  G4double dimZ = 10.*cm;
  
  G4Box box("Frame",dimX, dimY, dimZ);
  
  G4LogicalVolume logicalVol(&box,material,"LFrame", 0, 0, 0);
  logicalVol.SetMaterial(material); 
  
  physicalFrame = new G4PVPlacement(0,G4ThreeVector(),
				    "PFrame",&logicalVol,0,false,0);

}

  void G4TestSetup::makeMaterials()
{
  G4Material* Be = new G4Material("Beryllium",    4.,  9.01*g/mole, 1.848*g/cm3);
  G4Material* Graphite = new G4Material("Graphite",6., 12.00*g/mole, 2.265*g/cm3 );
  G4Material* Al  = new G4Material("Aluminium", 13., 26.98*g/mole, 2.7 *g/cm3);
  G4Material* Si  = new G4Material("Silicon",   14., 28.055*g/mole, 2.33*g/cm3);
  G4Material* LAr = new G4Material("LArgon",   18., 39.95*g/mole, 1.393*g/cm3);
  G4Material* Fe  = new G4Material("Iron",      26., 55.85*g/mole, 7.87*g/cm3);
  G4Material* Cu  = new G4Material("Copper",    29., 63.55*g/mole, 8.96*g/cm3);
  G4Material*  W  = new G4Material("Tungsten", 74., 183.85*g/mole, 19.30*g/cm3);
  G4Material* Pb  = new G4Material("Lead",      82., 207.19*g/mole, 11.35*g/cm3);
  G4Material*  U  = new G4Material("Uranium", 92., 238.03*g/mole, 18.95*g/cm3);

  G4Element*   H  = new G4Element ("Hydrogen", "H", 1. ,  1.01*g/mole);
  G4Element*   O  = new G4Element ("Oxygen"  , "O", 8. , 16.00*g/mole);
  G4Element*   C  = new G4Element ("Carbon"  , "C", 6. , 12.00*g/mole);
  G4Element*  Cs  = new G4Element ("Cesium"  , "Cs", 55. , 132.905*g/mole);
  G4Element*   N  = new G4Element("Nitrogen",   "N" , 7., 14.01*g/mole);
  G4Element*   I  = new G4Element ("Iodide"  , "I", 53. , 126.9044*g/mole);

  G4Material*  maO = new G4Material("Oxygen", 8., 16.00*g/mole, 1.1*g/cm3);

  G4Material* water = new G4Material ("Water" , 1.*g/cm3, 2);
  water->AddElement(H,2);
  water->AddElement(O,1);

  G4Material* ethane = new G4Material ("Ethane" , 0.4241*g/cm3, 2);
  ethane->AddElement(H,6);
  ethane->AddElement(C,2);
  
  G4Material* csi = new G4Material ("CsI" , 4.53*g/cm3, 2);
  csi->AddElement(Cs,1);
  csi->AddElement(I,1);

  G4Material* Air = new G4Material("Air"  ,  1.290*mg/cm3, 2);
  Air->AddElement(N,0.7);
  Air->AddElement(O,0.3);

  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

 G4int nMaterials = G4Material::GetNumberOfMaterials();

  G4cout << "Select the material among the available ones: " << G4endl;
  for (G4int mat = 0; mat < nMaterials; mat++)
    {
      G4cout << mat << ") "
	     << (*theMaterialTable)[mat]->GetName()
	     << G4endl;
    }
  G4int materialId;
  G4cin >> materialId;

  material = (*theMaterialTable)[materialId] ;

  G4cout << "The selected material is: " << material->GetName() << G4endl;
}

const G4Track* G4TestSetup::makeTrack()
{
  G4double energy = eMin + (eMax - eMin) * G4UniformRand();

  if (track == 0)
    {   
      G4double initX = 0.; 
      G4double initY = 0.; 
      G4double initZ = 1.;
      G4ParticleMomentum direction(initX,initY,initZ);
      G4DynamicParticle dynamicPart(part,direction,energy);      
      G4ThreeVector position(0.,0.,0.);
      G4double time = 0. ;
      
      track = new G4Track(&dynamicPart,time,position);
      
      // do I really need this?     
      G4GRSVolume* touche = new G4GRSVolume(physicalFrame,0,position);   
      track->SetTouchable(touche);
    }
  else
    {
      track->SetKineticEnergy(energy);
    }

  if (track == 0) G4cout << "Track = 0" << G4endl;
  G4double e = track->GetKineticEnergy();
  G4cout << "Track energy = " << e << G4endl;

  return track;
}

const G4Step* G4TestSetup::makeStep()
{
  step->SetTrack(track);

  G4ThreeVector aPosition(0.,0.,0.);
  G4ThreeVector newPosition(0.,0.,1.*mm);
  
  G4StepPoint* aPoint = new G4StepPoint();
  aPoint->SetPosition(aPosition);
  aPoint->SetMaterial(material);
  G4double safety = 10000.*cm;
  aPoint->SetSafety(safety);
  step->SetPreStepPoint(aPoint);
  G4StepPoint* newPoint = new G4StepPoint();
  newPoint->SetPosition(newPosition);
  newPoint->SetMaterial(material);
  step->SetPostStepPoint(newPoint);
  step->SetStepLength(1*mm);

  track->SetStep(step); 
  
return step;
}

void G4TestSetup::init()
{
  G4int processSelection;
  G4cout << "LowEnergy [1] or Standard [2] or Polarised? [3]" << G4endl;
  G4cin >> processSelection;
  if (processSelection < 1 || processSelection > 3) G4Exception("Wrong input");

  selection = processSelection;
  if (selection == 1) selName = "lowe";
  if (selection == 2) selName = "std";
  if (selection == 3) selName = "pol";

  G4cout << "Process to be tested: " << G4endl
	 << "Compton [1], GammaConversion [2], Photoelectric [3], Rayleigh [4]" << G4endl
	 << "Bremsstrahlung [5], eIonisation [6]" << G4endl;
  G4cin >> processType;
  if (processType < 1 || processType > 6) G4Exception("Wrong input");

  if (processType == 1) pName = "compton";
  if (processType == 2) pName = "conversion";
  if (processType == 3) pName = "photoel";
  if (processType == 4) pName = "rayleigh";
  if (processType == 5) pName = "brem";
  if (processType == 6) pName = "ionisation";

  if (processType < 5) part = G4Gamma::GammaDefinition();
  else part = G4Electron::ElectronDefinition();

  G4cout << "Min and max energy (MeV) of the incident particle" << G4endl;
  G4cin >> eMin >> eMax;
  eMin = eMin * MeV;
  eMax = eMax * MeV;

  makeMaterials();
  makeGeometry();
  createElectronProcesses();
}
