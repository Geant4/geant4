// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4hLETestLossTable.cc,v 1.1 2000-05-03 16:32:20 chauvie Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      File name:     CremonesiUT.cc
//
//      Author:        Stephane Chauvie
// 
//      Creation date: 5 January 2000
//
//      Modifications: 
//
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4ios.hh"
#include <fstream>
#include <iomanip>

#include "G4Material.hh"
#include "G4ProcessManager.hh"
#include "G4hLowEnergyIonisation.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleChange.hh"
#include "G4DynamicParticle.hh"

#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4AntiProton.hh"

#include "CLHEP/Hist/TupleManager.h"
#include "CLHEP/Hist/HBookFile.h"
#include "CLHEP/Hist/Histogram.h"
#include "CLHEP/Hist/Tuple.h"

HepTupleManager* hbookManager;

main()
{

  // Setup

  G4int particleID = 1;
  G4cout <<"Which particle?"<<G4endl<<setw(30)<<"[1] Proton, [2] AntiProton,"<<G4endl;
  cin>>particleID;
  
  G4cout.setf( ios::scientific, ios::floatfield );
  // -------------------------------------------------------------------
  ofstream out("g4hletestlosstable.dat");
  // ---- HBOOK initialization


  hbookManager = new HBookFile("g4hletestlosstable.hbook", 38);
  assert (hbookManager != 0);
  
  // ---- Book a histogram and ntuples
  G4cout<<"Hbook file name: "<<((HBookFile*) hbookManager)->filename()<<G4endl;
  
  // ---- primary ntuple ------
  HepTuple* ntuple1 = hbookManager->ntuple("Primary Ntuple");
  assert (ntuple1 != 0);
  
  //--------- Materials definition ---------

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

  static const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

  //--------- Particle definition ---------
  G4Electron* theElectron = G4Electron::Electron();
  
  G4ParticleDefinition* electron = G4Electron::ElectronDefinition();
  G4ParticleDefinition* proton = G4Proton::ProtonDefinition();
  G4ParticleDefinition* antiproton = G4AntiProton::AntiProtonDefinition();
  
  electron->SetCuts(1e-6*mm);
  proton->SetCuts(1e-6*mm);
  
  //--------- Processes definition ---------
  
  G4hLowEnergyIonisation* hIonisationProcess = new G4hLowEnergyIonisation("ionLowEIoni");
  hIonisationProcess->SetNuclearStoppingOn();
  hIonisationProcess->SetStoppingPowerTableName("ICRU_R49p"); 
  
  // Ionisation loss with/without Barkas effect
  hIonisationProcess->SetAntiProtonStoppingOn();
  //hIonisationProcess->SetAntiProtonStoppingOff();
      
    G4ProcessManager* theProtonProcessManager = new G4ProcessManager(proton);
      proton->SetProcessManager(theProtonProcessManager);
      theProtonProcessManager->AddProcess(hIonisationProcess);
      
    G4ProcessManager* theAntiProtonProcessManager = new G4ProcessManager(antiproton);
      antiproton->SetProcessManager(theAntiProtonProcessManager);
      theAntiProtonProcessManager->AddProcess(hIonisationProcess);
      
    // -------- create 1 Dynamic Particle  ----

  G4double partEnergy = 200*MeV;

  G4ParticleMomentum partDirection(1,0,0);
  
  G4DynamicParticle p(proton,partDirection,partEnergy);
  G4DynamicParticle pbar(antiproton,partDirection,partEnergy);
  
  // --------- check applicability
  
  G4ParticleDefinition* ProtonDefinition = p.GetDefinition();
  G4ParticleDefinition* AntiProtonDefinition = pbar.GetDefinition();
   if(! (hIonisationProcess->IsApplicable(*ProtonDefinition)
	|| hIonisationProcess->IsApplicable(*AntiProtonDefinition) )  )
    {
      G4Exception("FAIL: *** Not Applicable ***\n");
    }

  // Initialize the physics tables for ALL processes

  	hIonisationProcess->BuildPhysicsTable(*ProtonDefinition);
	
  	hIonisationProcess->BuildPhysicsTable(*AntiProtonDefinition);

  //------------------------------- Loss Table Test------------------------
  
  G4Material* apttoMaterial ;
  G4String MaterialName ;
 
  G4double minArg = 1*eV, maxArg = 200*MeV, argStp;
  const G4int pntNum = 1000;
  G4double Tkin[pntNum+1];
  argStp = (log10(maxArg)-log10(minArg))/pntNum;
  for(G4int d = 0; d < pntNum+1; d++){ 
    Tkin[d] = pow(10,(log10(minArg) + d*argStp));
  }
  
  //____________________LOSS TABLE TEST________________________________________________
    
  for ( G4int J = 0 ; J < theMaterialTable->length() ; J++ ){

    apttoMaterial = (*theMaterialTable)[ J ] ;
    MaterialName  = apttoMaterial->GetName() ;

     //G4cout<<"Material: "<<MaterialName<<G4endl;
       
    for (G4int ipnt=0 ; ipnt<pntNum; ipnt++){

	p.SetKineticEnergy(Tkin[ipnt]);
        pbar.SetKineticEnergy(Tkin[ipnt]);
         
    G4double dedxnow = 0 , cf = 0 , dedx =0 ;
    G4double* deltaCut = theElectron->GetCutsInEnergy();
    
    if( particleID == 1){
      dedxnow = hIonisationProcess->
        		GetPreciseDEDX(	apttoMaterial,
    		        		p.GetKineticEnergy(),
		   			p.GetDefinition()) ;
      if(Tkin[ipnt]<=2*MeV) dedxnow+= hIonisationProcess->
                            GetDeltaRaysEnergy( apttoMaterial,
                       	      			p.GetKineticEnergy(),
		       	      			deltaCut[J]);
    }
    if( particleID == 2){
      dedxnow = hIonisationProcess->
        		GetPreciseDEDX(	apttoMaterial,
    		        		pbar.GetKineticEnergy(),
		   			pbar.GetDefinition()) ;
      if(Tkin[ipnt]<=2*MeV) dedxnow+= hIonisationProcess->
        		    GetDeltaRaysEnergy( apttoMaterial,
                       	      			pbar.GetKineticEnergy(),
		       	      			deltaCut[J]);
    }
    
       
  
  //-----------ntuple-------------------------------------      
       
    if(MaterialName=="Silicon") out<<Tkin[ipnt]/MeV<<"   "<<dedxnow/MeV/mm<<G4endl;
    ntuple1->column("matind",J);
    ntuple1->column("kinen",Tkin[ipnt]/MeV);
    ntuple1->column("dedx",dedxnow/MeV/mm);
    ntuple1->dumpData();

    }
  }// for loop on materials
   
  
  hbookManager->write();
  delete hbookManager;

  // delete materials and elements
  delete Be;
  delete Graphite;
  delete Al;
  delete Si;
  delete LAr;
  delete Fe;
  delete Cu;
  delete W;
  delete Pb;
  delete U;
  delete H;
  delete maO;
  delete C;
  delete Cs;
  delete I;
  delete O;
  delete water;
  delete ethane;
  delete csi;

  cout<<"END OF THE MAIN PROGRAM"<<G4endl;
}  
