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
// $Id: G4mfpTest.cc,v 1.2 2001-10-12 13:10:56 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// File name:     G4mfpTest
//
// Author:        Maria Grazia Pia
//
// Description:   Tests MeanFreePath in electromagnetic processes
//                Output: ntuple with MeanFreePath and cross sections for 
//                lowenergy and standard processes
//
// Modifications: 
// 2 August 2001  MGP            Created
//
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4ios.hh"
#include "g4std/fstream"
#include "g4std/iomanip"

#include "G4Material.hh"
#include "G4VContinuousDiscreteProcess.hh"
#include "G4ProcessManager.hh"
#include "G4VeLowEnergyLoss.hh"
#include "G4EnergyLossTables.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleChange.hh"
#include "G4DynamicParticle.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4Step.hh"
#include "G4GRSVolume.hh"
#include "Randomize.hh"
#include "G4UnitsTable.hh"

#include "G4LowEnergyCompton.hh"
#include "G4LowEnergyPolarizedCompton.hh"
#include "G4ComptonScattering.hh"

#include "G4LowEnergyPhotoElectric.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4LowEnergyGammaConversion.hh"
#include "G4GammaConversion.hh"

#include "G4LowEnergyRayleigh.hh"
#include "G4LowEnergyBremsstrahlung.hh"
#include "G4eBremsstrahlung.hh"

#include "G4LowEnergyIonisation.hh"
#include "G4eIonisation.hh"

#include "CLHEP/Hist/TupleManager.h"
#include "CLHEP/Hist/HBookFile.h"
#include "CLHEP/Hist/Histogram.h"
#include "CLHEP/Hist/Tuple.h"

#include "G4LowEnergyGammaConversionMG.hh"
#include "G4LowEnergyRayleighMG.hh"
#include "G4LowEnergyComptonMG.hh"
#include "G4LowEnergyPhotoElectricMG.hh"


HepTupleManager* hbookManager;

int main()
{
  // Setup

  G4int nIterations = 100000;
  G4int materialId = 3;
  G4int test = 0;

  G4cout.setf( ios::scientific, ios::floatfield );

  // -------------------------------------------------------------------

  // ---- HBOOK initialization


  hbookManager = new HBookFile("mfp.hbook", 58);
  assert (hbookManager != 0);
  
  // ---- Book a histogram and ntuples
  G4cout<<"Hbook file name: "<<((HBookFile*) hbookManager)->filename()<<endl;
  
  // ---- primary ntuple ------
  HepTuple* ntuple1 = hbookManager->ntuple("MeanFreePath Ntuple");
  assert (ntuple1 != 0);
  
  //--------- Material definitions ---------

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
  G4Element*   I  = new G4Element ("Iodide"  , "I", 53. , 126.9044*g/mole);
  G4Element*   N  = new G4Element("Nitrogen" , "N" , 7., 14.01*g/mole);

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

  G4double fractionMass;
  G4Material* air = new G4Material("Air", 1.290*mg/cm3, 2);
  air->AddElement(N, fractionMass=0.7);
  air->AddElement(O, fractionMass=0.3);

  // Interactive set-up

  G4cout << "How many interactions? " << G4endl;
  G4cin >> nIterations;

  if (nIterations <= 0) G4Exception("Wrong input");

  G4cout << "1) Compton 2) Photoelectric 3) GammaConversion 4) Rayleigh" << G4endl;
  G4cout << "5) Bremsstrahlung 6) Ionisation" << G4endl;
  G4int type;
  G4cin >> type;

  if (nIterations <= 0) G4Exception("Wrong input");
  if (type < 1 || type > 6) G4Exception("Wrong input");

  G4double initEnergy = 1*MeV; 
  G4double initX = 0.; 
  G4double initY = 0.; 
  G4double initZ = 1.;

  static const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

 G4int nMaterials = G4Material::GetNumberOfMaterials();

  G4cout << "Available materials are: " << G4endl;
  for (G4int mat = 0; mat < nMaterials; mat++)
    {
      G4cout << mat << ") "
	     << (*theMaterialTable)[mat]->GetName()
	     << G4endl;
    }
  G4cout << "Which material? " << G4endl;
  G4cin >> materialId;

  G4Material* material = (*theMaterialTable)[materialId] ;

  G4cout << "The selected material is: "
	 << material->GetName()
	 << G4endl;

  G4double dimX = 1*mm;
  G4double dimY = 1*mm;
  G4double dimZ = 1*mm;
  
  // Geometry 

  G4Box* theFrame = new G4Box ("Frame",dimX, dimY, dimZ);
  
  G4LogicalVolume* logicalFrame = new G4LogicalVolume(theFrame,
						      (*theMaterialTable)[materialId],
						      "LFrame", 0, 0, 0);
  logicalFrame->SetMaterial(material); 
  
  G4PVPlacement* physicalFrame = new G4PVPlacement(0,G4ThreeVector(),
						   "PFrame",logicalFrame,0,false,0);
  
  // Particle definitions

  G4ParticleDefinition* gamma = G4Gamma::GammaDefinition();
  G4ParticleDefinition* electron = G4Electron::ElectronDefinition();
  G4ParticleDefinition* positron = G4Positron::PositronDefinition();
  
  gamma->SetCuts(1e-3*mm);
  electron->SetCuts(1e-3*mm);
  positron->SetCuts(1e-3*mm);
  
  // Electrons
  G4ProcessManager* eProcessManager = new G4ProcessManager(electron);
  electron->SetProcessManager(eProcessManager);

 // Bremsstrahlung
  G4LowEnergyBremsstrahlung* bremsstrahlung = new G4LowEnergyBremsstrahlung;
  eProcessManager->AddProcess(bremsstrahlung);
  bremsstrahlung->BuildPhysicsTable(*electron);
  G4eBremsstrahlung* bremsstrahlungStd = new G4eBremsstrahlung;
  eProcessManager->AddProcess(bremsstrahlungStd);
  bremsstrahlungStd->BuildPhysicsTable(*electron);

  // Ionisation
  G4LowEnergyIonisation* ionisation = new G4LowEnergyIonisation;
  eProcessManager->AddProcess(ionisation);
  ionisation->BuildPhysicsTable(*electron);
  G4eIonisation* ionisationStd = new G4eIonisation;
  eProcessManager->AddProcess(ionisationStd);
  ionisationStd->BuildPhysicsTable(*electron);

  // Positrons    
  G4ProcessManager* positronProcessManager = new G4ProcessManager(positron);
  positron->SetProcessManager(positronProcessManager);
  positronProcessManager->AddProcess(bremsstrahlung);
  
  // Photons
  G4ProcessManager* gProcessManager = new G4ProcessManager(gamma);
  gamma->SetProcessManager(gProcessManager);

  // Compton
  G4LowEnergyCompton* compton = new G4LowEnergyCompton;
  G4LowEnergyComptonMG* comptonMG = new G4LowEnergyComptonMG;
  G4ComptonScattering* comptonStd = new G4ComptonScattering;
  gProcessManager->AddProcess(compton);
  compton->BuildPhysicsTable(*gamma);
  gProcessManager->AddProcess(comptonMG);
  comptonMG->BuildPhysicsTable(*gamma);
  gProcessManager->AddProcess(comptonStd);
  comptonStd->BuildPhysicsTable(*gamma);

  // Photoelectric
  G4LowEnergyPhotoElectric* photoelectric = new G4LowEnergyPhotoElectric;
  gProcessManager->AddProcess(photoelectric);
  photoelectric->BuildPhysicsTable(*gamma);
  G4PhotoElectricEffect* photoelectricStd = new G4PhotoElectricEffect;
  gProcessManager->AddProcess(photoelectricStd);
  photoelectricStd->BuildPhysicsTable(*gamma);
  G4LowEnergyPhotoElectricMG* photoelectricMG = new G4LowEnergyPhotoElectricMG;
  gProcessManager->AddProcess(photoelectricMG);
  photoelectricMG->BuildPhysicsTable(*gamma);

  // GammaConversion
  G4LowEnergyGammaConversion* gammaConversion = new G4LowEnergyGammaConversion;
  gProcessManager->AddProcess(gammaConversion);
  gammaConversion->BuildPhysicsTable(*gamma);
  G4GammaConversion* gammaConversionStd = new G4GammaConversion;
  gProcessManager->AddProcess(gammaConversionStd);
  gammaConversionStd->BuildPhysicsTable(*gamma);
  G4LowEnergyGammaConversionMG* gammaConversionMG = new G4LowEnergyGammaConversionMG;
  gProcessManager->AddProcess(gammaConversionMG);
  gammaConversionMG->BuildPhysicsTable(*gamma);

  // Rayleigh
  G4LowEnergyRayleigh* rayleigh = new G4LowEnergyRayleigh;
  gProcessManager->AddProcess(rayleigh);
  rayleigh->BuildPhysicsTable(*gamma);
  G4LowEnergyRayleighMG* rayleighMG = new G4LowEnergyRayleighMG;
  gProcessManager->AddProcess(rayleighMG);
  rayleighMG->BuildPhysicsTable(*gamma);

  // --------- Test GetMeanFreePath

  G4cout << "Mean Free Path in material " << material->GetName() << G4endl;

  G4double eMin = 250. * eV;
  G4double eMax = 200. * GeV;
  G4double diff = eMax - eMin;

  for (G4int iter=0; iter<nIterations; iter++)
    {
      G4double range = G4UniformRand();
      G4double gEnergy;
      if (range < 0.1)
	{
	  eMin = 250. * eV;
          eMax = 10. * keV;
	  diff = eMax - eMin;
	  G4double random =  G4UniformRand();
	  gEnergy = eMin + diff * random;
	}
      else if (range < 0.2)
	{
	  eMin = 10. * keV;
          eMax = 100. * keV;
	  diff = eMax - eMin;
	  G4double random =  G4UniformRand();
	  gEnergy = eMin + diff * random;
	}
      else if (range < 0.3)
	{
	  eMin = 100. * keV;
          eMax = 1. * MeV;
	  diff = eMax - eMin;
	  G4double random =  G4UniformRand();
	  gEnergy = eMin + diff * random;
	}
      else if (range < 0.4)
	{
	  eMin = 1. * MeV;
          eMax = 10. * MeV;
	  diff = eMax - eMin;
	  G4double random =  G4UniformRand();
	  gEnergy = eMin + diff * random;
	}
      else if (range < 0.5)
	{
	  eMin = 10. * MeV;
          eMax = 100. * MeV;
	  diff = eMax - eMin;
	  G4double random =  G4UniformRand();
	  gEnergy = eMin + diff * random;
	}
      else if (range < 0.6)
	{
	  eMin = 100. * MeV;
          eMax = 1. * GeV;
	  diff = eMax - eMin;
	  G4double random =  G4UniformRand();
	  gEnergy = eMin + diff * random;
	}
      else if (range < 0.8)
	{
	  eMin = 1. * GeV;
          eMax = 100. * GeV;
	  diff = eMax - eMin;
	  G4double random =  G4UniformRand();
	  gEnergy = eMin + diff * random;
	}
         else
	{
	  gEnergy = eMin + G4UniformRand();
	}

      if (type == 3 && gEnergy < 1.0220*MeV) gEnergy = 1.02200 * MeV;

      // Create a DynamicParticle  
      G4ParticleMomentum gDirection(initX,initY,initZ);
      G4DynamicParticle dynamicPhoton(G4Gamma::Gamma(),gDirection,gEnergy);
      
      // Track 
      
      G4ThreeVector aPosition(0.,0.,0.);
      G4ThreeVector newPosition(0.,0.,1.*mm);
      G4double aTime = 0. ;
      
      G4Track* gTrack = new G4Track(&dynamicPhoton,aTime,aPosition);
      
      // do I really need this?
      
      G4GRSVolume* touche = new G4GRSVolume(physicalFrame, 0, aPosition);   
      gTrack->SetTouchable(touche);
      
      // Step 
      
      G4Step* step = new G4Step();  
      step->SetTrack(gTrack);      
      G4StepPoint* aPoint = new G4StepPoint();
      aPoint->SetPosition(aPosition);
      aPoint->SetMaterial(material);
      G4double safety = 100.*cm;
      aPoint->SetSafety(safety);
      step->SetPreStepPoint(aPoint);
      G4StepPoint* newPoint = new G4StepPoint();
      newPoint->SetPosition(newPosition);
      newPoint->SetMaterial(material);
      newPoint->SetSafety(safety);
      step->SetPostStepPoint(newPoint);
      step->SetStepLength(1*micrometer);
      gTrack->SetStep(step); 

      G4ForceCondition* force = new G4ForceCondition;

      G4double mfp;
      G4double mfpStd;
      G4double mfpMG = 0.;

      if (type == 1)
	{
	  mfp = compton->GetMeanFreePath(*gTrack,0.,force);
	  mfpStd = comptonStd->GetMeanFreePath(*gTrack,0.,force);
	  mfpMG = comptonMG->DumpMeanFreePath(*gTrack,0.,force);
	}
      else if (type == 2)
	{
	  mfp = photoelectric->GetMeanFreePath(*gTrack,0.,force);
	  mfpStd = photoelectricStd->GetMeanFreePath(*gTrack,0.,force);
	  mfpMG = photoelectricMG->DumpMeanFreePath(*gTrack,0.,force);
	}
	else if (type == 3)
	{
	  mfp = gammaConversion->GetMeanFreePath(*gTrack,0.,force);
	  mfpStd = gammaConversionStd->GetMeanFreePath(*gTrack,0.,force);
	  mfpMG = gammaConversionMG->DumpMeanFreePath(*gTrack,0.,force);
	}
      else if (type == 4)
	{
	  mfp = rayleigh->GetMeanFreePath(*gTrack,0.,force);
	  //  mfpStd = rayleighStd->GetMeanFreePath(*gTrack,0.,force);
	  mfpMG = rayleighMG->DumpMeanFreePath(*gTrack,0.,force);
	}
      else if (type == 5)
	{
	  mfp = bremsstrahlung->GetMeanFreePath(*gTrack,0.,force);
	  mfpStd = bremsstrahlungStd->GetMeanFreePath(*gTrack,0.,force);
	  //  mfpMG = bremsstrahlungMG->DumpMeanFreePath(*gTrack,0.,force);
	}
      else 
	{
	  mfp = ionisation->GetMeanFreePath(*gTrack,0.,force);
	  mfpStd = ionisationStd->GetMeanFreePath(*gTrack,0.,force);
	  //  mfpMG = ionisationMG->DumpMeanFreePath(*gTrack,0.,force);
	}

      G4double sigma = 0.;
      if (mfp > 0.) sigma = 1./mfp;
      G4double sigmaMG = 0.;
      if (mfpMG > 0.) sigmaMG = 1./mfpMG;
      G4double sigmaStd = 0.;
      if (mfpStd > 0.) sigmaStd = 1./mfpStd;

      G4double kineticEnergy = gTrack->GetKineticEnergy();

      G4cout << iter << ") e = " << gEnergy 
	     << " - mfp = " << mfp 
	     << ";  mfpMG = " << mfpMG << G4endl;

      ntuple1->column("e", gEnergy);
      ntuple1->column("ekin", kineticEnergy);
      ntuple1->column("mfp", mfp);
      ntuple1->column("mfpmg", mfpMG); 
      ntuple1->column("mfpstd", mfpStd); 
      ntuple1->column("sigma", sigma);
      ntuple1->column("sigmamg", sigmaMG); 
      ntuple1->column("sigmastd", sigmaStd);
      ntuple1->dumpData();      
    } 
  

  cout  << "End of iteration "  <<  G4endl;
  hbookManager->write();
  delete hbookManager;
  
  // delete materials and elements
  //  delete Be;
  //  delete Graphite;
  //  delete Al;
  //  delete Si;
  //  delete LAr;
  //  delete Fe;
  //  delete Cu;
  //  delete W;
  //  delete Pb;
  //  delete U;
  //  delete H;
  //  delete maO;
  //  delete C;
  //  delete Cs;
  //  delete I;
  //  delete O;
  //  delete water;
  //  delete ethane;
  //  delete csi;
  //  delete step;
  //  delete touche;
  //  delete Be;
  //  delete Graphite;
  //  delete Al;
  //  delete Si;
  //  delete LAr;
  //  delete Fe;
  //  delete Cu;
  //  delete W;
  //  delete Pb;
  //  delete U;
  //  delete H;
  //  delete maO;
  //  delete C;
  //  delete Cs;
  //  delete I;
  //  delete O;
  //  delete water;
  //  delete ethane;
  //  delete csi;
  // delete step;
  //  delete touche;

  cout << "END OF THE MAIN PROGRAM" << G4endl;
}












