// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PiMinusAbsorptionAtRestTest.cc,v 1.2 1999-12-15 14:53:39 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      File name:     G4PiMinusAbsorptionAtRestTest.cc 
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it), 
//                     (from Ch. Voelcker, from M. Maire)
// 
//      Creation date: 8 May 1998
//
//      Modifications: 
//      MG Pia         6 Jul 1998 Modified handling of ProcessManager
//                                to be consistent with changes in 
//                                ParticleDefinition
//
// -------------------------------------------------------------------

#include "G4ios.hh"
#include "g4std/fstream"
#include "g4std/iomanip"

#include "G4Material.hh"
#include "G4Element.hh"

#include "G4ProcessManager.hh"
#include "G4PiMinusAbsorptionAtRest.hh"
#include "G4PionMinusAbsorptionAtRest.hh"

#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include "G4PionMinus.hh"
#include "G4ParticleMomentum.hh"

#include "G4Box.hh"
#include "G4PVPlacement.hh"

#include "G4GRSVolume.hh"
#include "G4LogicalVolume.hh"

#include "G4ProcessManager.hh"
#include "G4ForceCondition.hh"

#include "G4Step.hh"
#include "G4Track.hh"

#include "CLHEP/Hist/TupleManager.h"
#include "CLHEP/Hist/HBookFile.h"
#include "CLHEP/Hist/Histogram.h"
#include "CLHEP/Hist/Tuple.h"

#include "G4IonTable.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4ParticleTypes.hh"
#include "G4IonTable.hh"
#include "G4ParticleTable.hh"
#include "G4GenericIon.hh"

int main()
{
  G4ParticleDefinition* electron = G4Electron::ElectronDefinition();
  G4ParticleDefinition* ion  = G4GenericIon::GenericIonDefinition();

  G4ParticleDefinition* proton = G4Proton::ProtonDefinition();
  proton->SetCuts(0.01);
  G4ParticleDefinition* neutron  = G4Neutron::NeutronDefinition();
  neutron->SetCuts(0.01);


  G4ParticleTable* theTableOfParticles;
  theTableOfParticles = G4ParticleTable::GetParticleTable();
  G4IonTable* theTable = new G4IonTable();

  // -------------------------------------------------------------------
  // MGP ---- HBOOK initialization
  HepTupleManager* hbookManager;
  hbookManager = new HBookFile("mg.hbook", 58);
  assert (hbookManager != 0);

  // MGP ---- Book a histogram
  HepHistogram* hEKin;
  hEKin = hbookManager->histogram("Kinetic Energy", 100,0.,150.);
  assert (hEKin != 0);  

  HepHistogram* hP;
  hP = hbookManager->histogram("Momentum", 100,0.,1000.);
  assert (hP != 0);  

  HepHistogram* hNSec;
  hNSec = hbookManager->histogram("Number of secondaries", 25,0.,25.);
  assert (hNSec != 0);  

  HepHistogram* hDebug;
  hDebug = hbookManager->histogram("Debug", 100,0.,140.);
  assert (hDebug != 0);  

  // MGP ---- Book a ntuple
  HepTuple* ntuple;
  ntuple = hbookManager->ntuple("Pion absorption at rest Ntuple");
  assert (ntuple != 0);

  G4String name, symbol;
  G4double a, iz, z, density;
  G4int nEl;

  G4int nIter = 10;
  G4int imat = 3;
  G4int verboseLevel = 1;
  G4int processId = 1;
  G4int deexcitationIndex =0;

  cout << " 0) Copper      1) Lead         2) Iron         3) Carbon" << G4endl;
  cout << " 4) LArgon      5) Polystyrene  6) Tungsten     7) Oxygen" << G4endl;
  cout << " 8) Beryllium   9) Aluminium   10) Uranium     11) BGO" << G4endl;
  cout << "12) NaI        13) CsI         14) Kapton" << G4endl;

  G4cout 
    << "Enter number of absorptions [10], material [3], Verbose level [1]"
    << G4endl;
  G4cin >> nIter >> imat >> verboseLevel;

  G4cout 
    << "Enter process: 1 G4PiMinusAbsorptionAtRest, 2 G4GheishaAbsAtRest" 
    << G4endl;
  G4cin >> processId;

  if (processId < 1 || processId >2)
    {
      G4cout << "Wrong process ID, set to default" << G4endl;
      processId = 1;
    }

  G4cout << "Enter deexcitation algorithm: 0 Theo, 1 Dummy" << G4endl;
  G4cin >> deexcitationIndex;
  
  G4int hnt;
  G4cout << "Enter histo (0) or ntuple (1) " << G4endl;
  G4cin >> hnt;

  // Materials definition 

  // Materials definition 

    G4Material* Cu = new G4Material(name="Copper", density=8.96*g/cm3, nEl=1);
    G4Element* elCu = new G4Element(name="Copper", symbol="Cu", iz=29., a=63.55*g/mole);
    Cu->AddElement( elCu, 1 );
    
    G4Material* Pb = new G4Material(name="Lead", density=11.35*g/cm3, nEl=1);
    G4Element* elPb = new G4Element(name="Lead", symbol="Pb", iz=82., a=207.19*g/mole);
    Pb->AddElement( elPb, 1 );
    
    G4Material*  Fe = new G4Material(name="Iron", density=7.87*g/cm3, nEl=1);
    G4Element* elFe = new G4Element(name="Iron", symbol="Fe", iz=26., a=55.85*g/mole);
    Fe->AddElement( elFe, 1 );
    
    G4Material*  Graphite = new G4Material(name="Graphite", density=2.265*g/cm3, nEl=1);
    G4Element* elC  = new G4Element(name="Carbon", symbol="C", iz=6., a=12.0107*g/mole);
    Graphite->AddElement( elC, 1 );

    G4Material*  LAr= new G4Material(name="LArgon", density=1.393*g/cm3, nEl=1);
    G4Element* elAr  = new G4Element(name="Argon", symbol="Ar", iz=18., a=39.95*g/mole);
    LAr->AddElement( elAr, 1 );
    
    G4Material*  PS = new G4Material(name="PolyStyrene", density=1.032*g/cm3, nEl=2);
    //    G4Element* elC = new G4Element(name="Carbon", symbol="C", iz=6., a=12.01*g/mole);
    G4Element* elH = new G4Element(name="Hydrogen", symbol="H", iz=1., a=1.01*g/mole);
    PS->AddElement( elC, 8 );
    PS->AddElement( elH, 8 );
    
    G4Material*  W  = new G4Material(name="Tungsten", density=19.30*g/cm3, nEl=1);
    G4Element* elW  = new G4Element(name="Tungsten", symbol="W", iz=74., a=183.85*g/mole);
    W->AddElement( elW, 1 );
    
    // approximate numbers for O
    G4Material*  O = new G4Material(name="Oxygen", density=1.1*g/cm3, nEl=1);
    G4Element* elO  = new G4Element(name="Oxygen", symbol="O", iz=8., a=15.9994*g/mole);
    O->AddElement( elO,  1 );
    
    G4Material*  Be = new G4Material(name="Beryllium", density=1.848*g/cm3, nEl=1);
    G4Element* elBe  = new G4Element(name="Beryllium", symbol="Be", iz=4., a=9.01*g/mole);
    Be->AddElement( elBe, 1 );
    
    G4Material*  Al = new G4Material(name="Aluminium", density=2.70*g/cm3, nEl=1);
    G4Element* elAl  = new G4Element(name="Aluminium", symbol="Al", iz=13., a=26.98*g/mole);
    Al->AddElement( elAl, 1 );
    
    
    G4Material*  U = new G4Material(name="Uranium", density=18.95*g/cm3, nEl=1);
    G4Element* elU  = new G4Element(name="Uranium", symbol="U", iz=92., a=238.03*g/mole);
    U->AddElement( elU, 1 );
    
    G4Material*  BGO = new G4Material(name="BGO", density=2.15*g/cm3, nEl=3);
    G4Element* elBi = new G4Element(name="Bismuth", symbol="Bi", iz=83., a=208.98*g/mole);
    G4Element* elGe = new G4Element(name="Germanium", symbol="Ge", iz=32., a=72.59*g/mole);
    BGO->AddElement( elBi, 4 );
    BGO->AddElement( elGe, 3 );
    BGO->AddElement( elO, 12 );
    
    G4Material*  NaI = new G4Material(name="NaI", density=3.67*g/cm3, nEl=2);
    G4Element* elNa = new G4Element(name="Sodium", symbol="Na", iz=11., a=22.990*g/mole);
    G4Element* elI = new G4Element(name="Iodine", symbol="I", iz=53., a=126.904*g/mole);
    NaI->AddElement( elNa, 1 );
    NaI->AddElement( elI, 1 );
    
    G4Material*  CsI = new G4Material(name="CsI", density=4.53*g/cm3, nEl=2);
    G4Element* elCs = new G4Element(name="Cesium", symbol="Cs", iz=55., a=132.905*g/mole);
    CsI->AddElement( elCs, 1 );
    CsI->AddElement( elI, 1 );

    G4Material*  Kapton = new G4Material(name="Kapton", density=1.53*g/cm3, nEl=4); 
    // formula: private communications, see mail.
    Kapton->AddElement( elC, 22 );
    Kapton->AddElement( elH, 10 );
    Kapton->AddElement( elO,  5 );
    G4Element* elN = new G4Element(name="Nitrogen", symbol="N", iz=7., a=14.007*g/mole);
    Kapton->AddElement( elN, 2 );


    //  G4Material* Al = new G4Material("Aluminium", 13., 26.98*g/mole, 2.7 *g/cm3 );
    //  G4Material* Fe = new G4Material("Iron",      26., 55.85*g/mole, 7.87*g/cm3 );
    //  G4Material* Pb = new G4Material("Lead",      82., 207.19*g/mole, 11.35*g/cm3 );
    //  G4Material* Graphite = new G4Material("Graphite", 6., 12.00*g/mole, 2.265*g/cm3 );
  
    //  G4Element*   H = new G4Element ("Hydrogen", "H", 1. ,  1.01*g/mole);
    //  G4Element*   O = new G4Element ("Oxygen"  , "O", 8. , 16.00*g/mole);
    //  G4Element*   C = new G4Element ("Carbon"  , "C", 6. , 12.00*g/mole);
    //  G4Element*  Cs = new G4Element ("Cesium"  , "Cs", 55. , 132.905*g/mole);
    //  G4Element*   I = new G4Element ("Iodide"  , "I", 53. , 126.9044*g/mole);
  
    //  G4Material* water = new G4Material ("Water" , 1.*g/cm3, 2);
    //  water->AddElement(H,2);
    //  water->AddElement(O,1);
    //  
    //  G4Material* ethane = new G4Material ("Ethane" , 0.4241*g/cm3, 2);
    //  ethane->AddElement(H,6);
    //  ethane->AddElement(C,2);
  
    //  G4Material* csi = new G4Material ("CsI" , 4.53*g/cm3, 2);
    //  csi->AddElement(Cs,1);
    //  csi->AddElement(I,1);
  
  //  G4Element::DumpInfo(); 
  //  G4Material::DumpInfo();
  
  static const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  
  // Fix a material
  // G4int imat = 0;   // graphite
  
  // Geometry definitions

  G4Box* theFrame = new G4Box ("Frame",1*m, 1*m, 1*m);
  G4LogicalVolume* LogicalFrame = 
    new G4LogicalVolume (theFrame,(*theMaterialTable)(imat),"LFrame", 0,0,0);
  
  G4PVPlacement* PhysicalFrame = 
    new G4PVPlacement(0,G4ThreeVector(),
		      "PFrame",LogicalFrame,0,false,0);
  
  // The center-of-mass of the cube should be located at the origin!

  // Particle 
  G4ParticleDefinition* pionMinus = G4PionMinus::PionMinusDefinition();

  // Processes
  G4ProcessManager* thePionProcessManager = new G4ProcessManager(pionMinus);
  pionMinus->SetProcessManager(thePionProcessManager);
  
  //  thePionProcessManager->SetParticleType(pionMinus);

  G4PiMinusAbsorptionAtRest* thePiAbsorptionProcess = 
    new G4PiMinusAbsorptionAtRest("G4PiMinusAbsorptionAtRest");
  G4PionMinusAbsorptionAtRest* gheishaAbsorptionProcess =
    new G4PionMinusAbsorptionAtRest();

  switch (processId)
    {
    case 1:
      {
	thePionProcessManager->AddProcess(thePiAbsorptionProcess,0,-1,0);
	thePiAbsorptionProcess->SetDeexcitationAlgorithm(deexcitationIndex);  
	break;
      }
    case 2:
      {
	thePionProcessManager->AddProcess(gheishaAbsorptionProcess,0,-1,0);
        break;
      }
    default:
      {
	thePionProcessManager->AddProcess(thePiAbsorptionProcess,0,-1,0);
	thePiAbsorptionProcess->SetDeexcitationAlgorithm(deexcitationIndex);  
	break;
      }
    }

  G4ForceCondition* condition;
  
  // Set cut and Build CrossSection Tables 
  pionMinus->SetCuts(1.*mm);
  
  // Create one Dynamic Particle  
  G4double pionEnergy = 0.*MeV;
  G4ParticleMomentum pionDirection(0.,0.,1.);
  G4DynamicParticle aPionMinus(G4PionMinus::PionMinus(),pionDirection,pionEnergy);
  
  // Track definition (for this test ONLY!)
  G4ThreeVector aPosition(0.,0.,0.);
  G4double aTime = 0. ;
  G4Track* ptrack = new G4Track(&aPionMinus,aTime,aPosition) ;
  G4Track& aTrack = (*ptrack) ;
  
  //ptrack->SetVolume(PhysicalFrame);
  
  // Do I really need this?
  G4GRSVolume* touche = new G4GRSVolume(PhysicalFrame, NULL, aPosition);   
  ptrack->SetTouchable(touche);
  
  // Create 1 Step (for this test only)
  G4Step* Step = new G4Step();  
  G4Step& aStep = (*Step);
  Step->SetTrack(ptrack);
  
  // Check applicability
  G4ParticleDefinition* PionMinusDefinition = aPionMinus.GetDefinition();

  if (! thePiAbsorptionProcess->IsApplicable(*PionMinusDefinition)) 
    {
      G4cout 
	<< PionMinusDefinition->GetParticleName() << " is not a PionMinus!" << G4endl;
      G4Exception("FAIL: *** exit program ***\n");
      //     return ;
    }

  // Test the DoIt for the Pion Absorption

  G4Material* apttoMaterial ;
  apttoMaterial = (*theMaterialTable)(imat) ;
  LogicalFrame->SetMaterial(apttoMaterial); 
  aPionMinus.SetKineticEnergy(0.*MeV);
  aPionMinus.SetMomentumDirection(0., 0., 1.);
  G4VParticleChange* aParticleChange;
  G4Track* aFinalParticle;
  G4String aParticleName;  
  
  thePiAbsorptionProcess->SetVerboseLevel(verboseLevel);
  
  G4int iteration = 0;   

  do {
    
    if (processId == 1)
      { aParticleChange = thePiAbsorptionProcess->AtRestDoIt(aTrack, aStep); }
    else 
      {
	if (processId == 2) 
	  { aParticleChange = gheishaAbsorptionProcess->AtRestDoIt(aTrack, aStep); }
      }
    
    // Loop over final particle List

    G4double e = 0;
    G4double eKin = 0;
    G4double Px = 0;
    G4double Py = 0;
    G4double Pz = 0;


    hNSec->accumulate(aParticleChange->GetNumberOfSecondaries());
    hDebug->accumulate(aParticleChange->GetLocalEnergyDeposit());

    G4int i;
    for (i=0; i<(aParticleChange->GetNumberOfSecondaries()); i++) 
      {
	// The following two items should be filled per event, not
	// per secondary; filled here just for convenience, to avoid
	// complicated logic 
	ntuple->column("nsec",aParticleChange->GetNumberOfSecondaries());
	ntuple->column("excitation",aParticleChange->GetLocalEnergyDeposit());

	aFinalParticle = aParticleChange->GetSecondary(i) ;
      
	e    =  aFinalParticle->GetTotalEnergy();
	eKin = aFinalParticle->GetKineticEnergy();
	Px   = (aFinalParticle->GetMomentum()).x() ;
	Py   = (aFinalParticle->GetMomentum()).y() ;
	Pz   = (aFinalParticle->GetMomentum()).z() ;
      
	aParticleName = aFinalParticle->GetDefinition()->GetParticleName();

	if (aFinalParticle->GetDefinition() == G4Proton::ProtonDefinition() ||
	    aFinalParticle->GetDefinition() == G4Neutron::NeutronDefinition())
	  {
	    hEKin->accumulate(eKin);
	    hP->accumulate(sqrt(Px*Px+Py*Py+Pz*Pz));

	    //	ntuple->column("px", Px);
	    //	ntuple->column("py", Py);
	    //	ntuple->column("pz", Pz);
	    //        ntuple->column("p", sqrt(Px*Px+Py*Py+Pz*Pz));
	    //        ntuple->column("e", e);
	    if (hnt == 1)
	      {
		ntuple->column("ekin", eKin);
		ntuple->dumpData(); 
	      }

	    delete aParticleChange->GetSecondary(i);
	    
	  }
      }

    G4cout << "******* Iteration = " << iteration << G4endl;
    iteration++;
    aParticleChange->Clear();
    
  }  while (iteration < nIter) ; // end of do-while
  
  hbookManager->write();

  // Clean up
  //  delete aFinalParticle;
  delete Al;
  delete Fe;
  delete Pb;
  delete Graphite;
  //  delete H;
  delete O;
  //  delete C;
  //  delete Cs;
  //  delete I;
  //  delete water;
  //  delete ethane;
  //  delete csi;

  delete touche;
  delete Step;

  
  return EXIT_SUCCESS;
}
