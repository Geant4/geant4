// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4KaonMinusAbsorptionAtRestTest.cc,v 1.2 1999-12-15 14:53:39 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// KaonMinusAtRestTest.cc 
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      File name:     G4KaonMinusAbsorptionAtRestTest.cc 
//
//      Author:        Christian Voelcker (from M. Maire)
// 
//      Creation date: ?
//
//      Modifications: 
//
//      Maria Grazia Pia 30 Jun 1998 Added histograms (with CLHEP)
//      Maria Grazia Pia  6 Jul 1998 Modified handling of ProcessManager
//                                   to be consistent with changes in 
//                                   ParticleDefinition
//
// -------------------------------------------------------------------

#include "G4ios.hh"
#include "g4std/fstream"
#include "g4std/iomanip"

#include "G4Material.hh"

#include "G4ProcessManager.hh"
#include "G4KaonMinusAbsorptionAtRest.hh"

#include "G4DynamicParticle.hh"
#include "G4KaonMinus.hh"

#include "G4Box.hh"
#include "G4PVPlacement.hh"

#include "G4Step.hh"
#include "G4GRSVolume.hh"

#include "CLHEP/Hist/TupleManager.h"
#include "CLHEP/Hist/HBookFile.h"
#include "CLHEP/Hist/Histogram.h"
#include "CLHEP/Hist/Tuple.h"


//    it tests the G4ComptonScattering process ----------------
//    ---------- M.Maire on 19/04/96 --------------------------
//
//  Modifs:
//  19-02-97, Modifs in DoIt for new physic scheme
//  25-02-97, 'simplified' version for the process manager
//  12-03-97, GetMeanFreePath and PostStepDoIt
//  21-03-97, calling sequence of AddProcess modified

int main()
{
  //-------- set output format-------
   G4cout.setf( G4std::ios::scientific, G4std::ios::floatfield );
  //-------- write results onto a file --------
   G4std::ofstream outFile( "KaonMinusAbsorptionAtRest.out", G4std::ios::out);
   outFile.setf( G4std::ios::scientific, G4std::ios::floatfield );

  // -------------------------------------------------------------------
  // MGP ---- HBOOK initialization
  HepTupleManager* hbookManager;
  hbookManager = new HBookFile("mgkaon.hbook", 58);
  assert (hbookManager != 0);

  // MGP ---- Book a histogram
  HepHistogram* hEKin;
  hEKin = hbookManager->histogram("Kinetic Energy", 100,0.,200.);
  assert (hEKin != 0);  

  //  HepHistogram* hE;
  //  hE = hbookManager->histogram("Energy", 100,0.,2000.);
  //  assert (hE != 0);  

  HepHistogram* hP;
  hP = hbookManager->histogram("Momentum", 100,0.,1000.);
  assert (hP != 0);  

  HepHistogram* hNSec;
  hNSec = hbookManager->histogram("Number of secondaries", 40,0.,40.);
  assert (hNSec != 0);  

  HepHistogram* hDebug;
  hDebug = hbookManager->histogram("Debug", 100,0.,200.);
  assert (hDebug != 0);  

  // MGP ---- Book a ntuple
  HepTuple* ntuple;
  ntuple = hbookManager->ntuple("Pion absorption at rest Ntuple");
  assert (ntuple != 0);

   G4int niter=10;
   G4int imat=3;
   G4int verboseLevel=1;
   G4int processID=1;
   G4cout << " How many absorptions? [10], Which material? [3], which Verbose Level? [1]" << G4endl;
   G4cin >> niter >> imat >> verboseLevel;
   G4cout << " which process? (1)=G4KaonMinusAbsorptionAtRest, (2)=G4HadronAtRest  [1]" << G4endl;
   G4cin >> processID;
   
  G4int hnt;
  G4cout << "Enter histo (0) or ntuple (1) " << G4endl;
  G4cin >> hnt;
  

  //
  //--------- Materials definition ---------
  //
  G4Material* Al = new G4Material("Aluminium", 13., 26.98*g/mole, 2.7 *g/cm3 );
  G4Material* Fe = new G4Material("Iron",      26., 55.85*g/mole, 7.87*g/cm3 );
  G4Material* Pb = new G4Material("Lead",      82., 207.19*g/mole, 11.35*g/cm3 );
  G4Material* Graphite = new G4Material("Graphite",6., 12.00*g/mole, 2.265*g/cm3 );

  G4Element*   H = new G4Element ("Hydrogen", "H", 1. ,  1.01*g/mole);
  G4Element*   O = new G4Element ("Oxygen"  , "O", 8. , 16.00*g/mole);
  G4Element*   C = new G4Element ("Carbon"  , "C", 6. , 12.00*g/mole);
  G4Element*  Cs = new G4Element ("Cesium"  , "Cs", 55. , 132.905*g/mole);
  G4Element*   I = new G4Element ("Iodide"  , "I", 53. , 126.9044*g/mole);


  G4Material* water = new G4Material ("Water" , 1.*g/cm3, 2);
  water->AddElement(H,2);
  water->AddElement(O,1);

  G4Material* ethane = new G4Material ("Ethane" , 0.4241*g/cm3, 2);
  ethane->AddElement(H,6);
  ethane->AddElement(C,2);
  
  G4Material* csi = new G4Material ("CsI" , 4.53*g/cm3, 2);
  csi->AddElement(Cs,1);
  csi->AddElement(I,1);
  

//  G4Element::DumpInfo(); 
//  G4Material::DumpInfo();

  static const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

  // fix a material
  // G4int imat = 0;   // graphite

  // Geometry definitions
  //
    G4Box* theFrame = new G4Box ("Frame",1*m, 1*m, 1*m);

    G4LogicalVolume* LogicalFrame = new G4LogicalVolume(theFrame,
                                                  (*theMaterialTable)(imat),
						  "LFrame", 0, 0, 0);

    G4PVPlacement* PhysicalFrame = new G4PVPlacement(0,G4ThreeVector(),
                                          "PFrame",LogicalFrame,0,false,0);

// the center-of-mass of the cube should be located at the origin!

  //--------- Particle definition ---------
  //
  G4ParticleDefinition* kaonMinus = G4KaonMinus::KaonMinusDefinition();

  //--------- Processes definition ---------
  //

  G4ProcessManager* theKaonProcessManager = new G4ProcessManager(kaonMinus);
  kaonMinus->SetProcessManager(theKaonProcessManager);
  
  //  G4ProcessManager* theKaonProcessManager = kaonMinus->GetProcessManager();

  G4KaonMinusAbsorptionAtRest  kaonMinusAbsorptionAtRest;
  //  G4HadronAtRest hadronAtRest;
  
  G4KaonMinusAbsorptionAtRest theKaonAbsorptionProcess;
  theKaonProcessManager->AddProcess(&theKaonAbsorptionProcess,-1,-1,0);
  //  G4HadronAtRest theHadronAbsorptionProcess;
  //  theKaonProcessManager->AddProcess(&theHadronAbsorptionProcess,-1,-1,0);

  G4ForceCondition* condition;

  // ------- set cut and Build CrossSection Tables -------
  //
  kaonMinus->SetCuts(1.*mm);

  // -------- create 1 Dynamic Particle  ----

  G4double kaonEnergy = 0.*MeV;
  G4ParticleMomentum kaonDirection(0.,0.,1.);
  G4DynamicParticle aKaonMinus(G4KaonMinus::KaonMinus(),kaonDirection,kaonEnergy);

  //--------- track definition (for this test ONLY!)------------
  //
    G4ThreeVector aPosition(0.,0.,0.);
    G4double aTime = 0. ;

    G4Track* ptrack = new G4Track(&aKaonMinus,aTime,aPosition) ;
    G4Track& aTrack = (*ptrack) ;

  //ptrack->SetVolume(PhysicalFrame);

  // do I really need this?
    G4GRSVolume* touche = new G4GRSVolume(PhysicalFrame, NULL, aPosition);   
    ptrack->SetTouchable(touche);

  // -------- create 1 Step (for this test only)----  

    G4Step* Step = new G4Step();  G4Step& aStep = (*Step);
    Step->SetTrack(ptrack);
  
 //
  // --------- check applicability
  //
    G4ParticleDefinition* KaonMinusDefinition=aKaonMinus.GetDefinition();
    if(!theKaonAbsorptionProcess.IsApplicable(*KaonMinusDefinition)) {
       G4cout 
            << KaonMinusDefinition->GetParticleName()
            << " is not a KaonMinus!" << G4endl;
       G4Exception("FAIL: *** exit program ***\n");
  //     return ;
     }
  //
  // --------- Test the DoIt for the Kaon Absorption
  //
    G4Material* apttoMaterial ;
    apttoMaterial = (*theMaterialTable)(imat) ;

    LogicalFrame->SetMaterial(apttoMaterial); 
    aKaonMinus.SetKineticEnergy(0.*MeV);
    aKaonMinus.SetMomentumDirection(0., 0., 1.);
    G4VParticleChange* aParticleChange;
    G4Track* aFinalParticle;
    G4String aParticleName;  

    theKaonAbsorptionProcess.SetVerboseLevel(verboseLevel);

    G4int iteration = 0;   
    do {

        if(processID==1){
           aParticleChange = theKaonAbsorptionProcess.AtRestDoIt(aTrack, aStep);
        } else if(processID==2) {
	  //           aParticleChange = theHadronAbsorptionProcess.AtRestDoIt(aTrack, aStep);
        } else {
           G4Exception("No such process!\n");
        }


	//        outFile << "! ---------------------------------------------------------------" << G4endl;   
	//        outFile << "! Energy deposit = " << aParticleChange->GetLocalEnergyDeposit() << G4endl;
	//
	//        G4double E = 0., Px = 0., Py = 0., Pz = 0., Edep = 0. ;
	//        
	//        outFile << "! Name         Px             Py            Pz          Etot" << G4endl;
	//        outFile << "! ---------------------------------------------------------------" << G4endl
        // loop over final particle List.

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
	// complicated logic to dump ntuple when there are no secondaries
	ntuple->column("nsec",aParticleChange->GetNumberOfSecondaries());
	ntuple->column("deposit",aParticleChange->GetLocalEnergyDeposit());

	aFinalParticle = aParticleChange->GetSecondary(i) ;
      
	e    =  aFinalParticle->GetTotalEnergy();
	eKin = aFinalParticle->GetKineticEnergy();
	Px   = (aFinalParticle->GetMomentum()).x() ;
	Py   = (aFinalParticle->GetMomentum()).y() ;
	Pz   = (aFinalParticle->GetMomentum()).z() ;
      
	aParticleName = aFinalParticle->GetDefinition()->GetParticleName();

        hEKin->accumulate(eKin);
        hP->accumulate(sqrt(Px*Px+Py*Py+Pz*Pz));

	ntuple->column("px", Px);
	ntuple->column("py", Py);
	ntuple->column("pz", Pz);
        ntuple->column("p", sqrt(Px*Px+Py*Py+Pz*Pz));
        ntuple->column("e", e);
        ntuple->column("ekin", eKin);
        
	// MGP ---- Dump the ntuple content at the end of current iteration
	if (hnt == 1) { ntuple->dumpData(); }
        delete aParticleChange->GetSecondary(i);

	//	outFile << "  " << aParticleName 
	//		<< "  " << Px 
	//		<< "  " << Py 
	//		<< "  " << Pz
	//		<< "  " << E 
	//		<< G4endl ;
      }

    aParticleChange->Clear();
    iteration++;
    G4cout << "******* Iteration = " << iteration << G4endl;
        
    }  while (iteration < niter) ;

  hbookManager->write();

  delete aFinalParticle;

 // delete materials and elements
  delete Al;
  delete Fe;
  delete Pb;
  delete Graphite;
  delete H;
  delete O;
  delete C;
  delete Cs;
  delete I;
  delete water;
  delete ethane;
  delete csi;
  
  return EXIT_SUCCESS;
}
