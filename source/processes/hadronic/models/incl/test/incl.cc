#include <stdio.h>

#include <iostream>
#include <fstream>

#include <string>

#include <stdlib.h>
#include <math.h>
#include <time.h>

// For Geant4:
#include "globals.hh"
#include "Randomize.hh"

#include "G4StateManager.hh"
#include "G4ParticleTable.hh"

#include "G4ios.hh"
//#include "G4BertiniData.hh"
#include "G4IonTable.hh"
#include "G4Nucleus.hh"
#include "G4NucleiModel.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4Proton.hh"
#include "G4KineticTrack.hh"
#include "G4KineticTrackVector.hh"
#include "G4Fancy3DNucleus.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleChange.hh"
#include "G4Track.hh"
#include "G4HadFinalState.hh"
#include "G4DynamicParticle.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"

#include "G4InclCascadeInterface.hh"
#include "G4InclAblaCascadeInterface.hh"

#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPhysicsList.hh"
#include "G4RunManager.hh"
#include "G4DecayPhysics.hh"
// End of G4 simulation declarations

#include "G4Incl.hh"

//#include "G4InclAblaLightIonInterface.hh"
#include "G4KineticTrackVector.hh"
#include "G4FragmentVector.hh"
#include "G4ParticleChange.hh"
#include "G4ReactionProductVector.hh"
#include "G4ReactionProduct.hh"
#include "G4FermiBreakUp.hh"
#include "math.h"
#include "G4GenericIon.hh"
#include "CLHEP/Random/Random.h"

// Needed for Dresner
const int volant1Size = 200;
struct volant1_ {
  int nombre1;
  float acv1[volant1Size], zcv1[volant1Size],ecv1[volant1Size],tcv1[volant1Size],pcv1[volant1Size];

  double exc;
};

struct cv_ {
  float tabener[2001];
};

//#define DEBUG 1
#ifdef USEROOT // Use ROOT for data analysis

#include "TROOT.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"

#include "G4VInclLogger.hh"
#include "G4InclRootLogger.hh"

#endif // USEROOT

using namespace std;

enum type {fullrun = 0, inclrun = 1};

class pList: public G4VUserPhysicsList
{
public:
  pList(){;};
  ~pList(){;};

protected:
  void ConstructParticle(){  G4Geantino::GeantinoDefinition();};
  void ConstructProcess(){  AddTransportation();};
  void SetCuts(){  SetCutsWithDefault();   }
};


class G4LogicalVolume;
class G4VPhysicalVolume;

#include "G4VUserDetectorConstruction.hh"

class det : public G4VUserDetectorConstruction
{
public:
  det():  experimentalHall_log(0), tracker_log(0),
	  calorimeterBlock_log(0), calorimeterLayer_log(0),
	  experimentalHall_phys(0), calorimeterLayer_phys(0),
	  calorimeterBlock_phys(0), tracker_phys(0)
  {;};
  ~det(){;};

  G4VPhysicalVolume* Construct(){

    G4double a;  // atomic mass
    G4double z;  // atomic number
    G4double density;

    G4Material* Ar = 
      new G4Material("ArgonGas", z= 18., a= 39.95*g/mole, density= 1.782*mg/cm3);

    G4Material* Al = 
      new G4Material("Aluminum", z= 13., a= 26.98*g/mole, density= 2.7*g/cm3);

    G4Material* Pb = 
      new G4Material("Lead", z= 82., a= 207.19*g/mole, density= 11.35*g/cm3);

    // volumes: experimental hall (world volume). Beam line along x axis

    G4double expHall_x = 3.0*m;
    G4double expHall_y = 1.0*m;
    G4double expHall_z = 1.0*m;
    G4Box* experimentalHall_box
      = new G4Box("expHall_box",expHall_x,expHall_y,expHall_z);
    experimentalHall_log = new G4LogicalVolume(experimentalHall_box,
					       Ar,"expHall_log",0,0,0);
    experimentalHall_phys = new G4PVPlacement(0,G4ThreeVector(),
					      experimentalHall_log,"expHall",0,false,0);

    G4double innerRadiusOfTheTube = 0.*cm;   // a tracker tube
    G4double outerRadiusOfTheTube = 60.*cm;
    G4double hightOfTheTube = 50.*cm;
    G4double startAngleOfTheTube = 0.*deg;
    G4double spanningAngleOfTheTube = 360.*deg;
    G4Tubs* tracker_tube = new G4Tubs("tracker_tube",innerRadiusOfTheTube,
				      outerRadiusOfTheTube,hightOfTheTube,
				      startAngleOfTheTube,spanningAngleOfTheTube);
    tracker_log = new G4LogicalVolume(tracker_tube,Al,"tracker_log",0,0,0);
    G4double trackerPos_x = -1.0*m;
    G4double trackerPos_y = 0.*m;
    G4double trackerPos_z = 0.*m;
    tracker_phys = new G4PVPlacement(0,
				     G4ThreeVector(trackerPos_x,trackerPos_y,trackerPos_z),
				     tracker_log,"tracker",experimentalHall_log,false,0);


    G4double block_x = 1.0*m;   // a calorimeter block
    G4double block_y = 50.0*cm;
    G4double block_z = 50.0*cm;
    G4Box* calorimeterBlock_box = new G4Box("calBlock_box",block_x,
					    block_y,block_z);
    calorimeterBlock_log = new G4LogicalVolume(calorimeterBlock_box,
					       Pb,"caloBlock_log",0,0,0);
    G4double blockPos_x = 1.0*m;
    G4double blockPos_y = 0.0*m;
    G4double blockPos_z = 0.0*m;
    calorimeterBlock_phys = new G4PVPlacement(0,
					      G4ThreeVector(blockPos_x,blockPos_y,blockPos_z),
					      calorimeterBlock_log,"caloBlock",experimentalHall_log,false,0);

    G4double calo_x = 1.*cm;   // calorimeter layers
    G4double calo_y = 40.*cm;
    G4double calo_z = 40.*cm;
    G4Box* calorimeterLayer_box = new G4Box("caloLayer_box",
					    calo_x,calo_y,calo_z);
    calorimeterLayer_log = new G4LogicalVolume(calorimeterLayer_box,
					       Al,"caloLayer_log",0,0,0);
    for(G4int i=0;i<19;i++) // loop for 19 layers
      {
	G4double caloPos_x = (i-9)*10.*cm;
	G4double caloPos_y = 0.0*m;
	G4double caloPos_z = 0.0*m;
	calorimeterLayer_phys = new G4PVPlacement(0,
						  G4ThreeVector(caloPos_x,caloPos_y,caloPos_z),
						  calorimeterLayer_log,"caloLayer",calorimeterBlock_log,false,i);
      }

    return experimentalHall_phys;
  };

private:
    
  // Logical volumes
  //
  G4LogicalVolume* experimentalHall_log;
  G4LogicalVolume* tracker_log;
  G4LogicalVolume* calorimeterBlock_log;
  G4LogicalVolume* calorimeterLayer_log;

  // Physical volumes
  //
  G4VPhysicalVolume* experimentalHall_phys;
  G4VPhysicalVolume* calorimeterLayer_phys;
  G4VPhysicalVolume* calorimeterBlock_phys;
  G4VPhysicalVolume* tracker_phys;
};

int main(int argc, char *argv[])
{
  G4RunManager* runManager = new G4RunManager;

  G4VUserDetectorConstruction* d = new det;
  runManager->SetUserInitialization(d);

  G4VUserPhysicsList* physics = new pList();

  runManager->SetUserInitialization(physics);

  G4Gamma::Gamma();
  const G4ParticleDefinition* proton = G4Proton::Proton();
  const G4ParticleDefinition* neutron = G4Neutron::Neutron();
  const G4ParticleDefinition* pin = G4PionMinus::PionMinus();
  const G4ParticleDefinition* pip = G4PionPlus::PionPlus();
  const G4ParticleDefinition* pi0 = G4PionZero::PionZero();

  const G4ParticleDefinition* deu = G4Deuteron::DeuteronDefinition();
  const G4ParticleDefinition* tri = G4Triton::TritonDefinition();
  const G4ParticleDefinition* he3 = G4He3::He3Definition();
  const G4ParticleDefinition* alp = G4Alpha::AlphaDefinition();
  const G4ParticleDefinition* ion = G4GenericIon::GenericIon();

  G4DecayPhysics decays;
  decays.ConstructParticle();

  runManager->Initialize();
  G4ParticleTable::GetParticleTable()->SetReadiness();

  int verboseLevel = 0;
  G4ParticleTable *theTableOfParticles = G4ParticleTable::GetParticleTable();
  G4FragmentVector *theSpectatorFermiBreakupResult = 0;
  G4FragmentVector *theFermiBreakupResult = 0;

  bool useProjSpect = false;
  bool useFermiBreakup = false;

  char *filename = new char [2000];
  char *summaryFilename = new char [2000];

  // For CPU time measurement:
  long starttime = 0;
  long stoptime = 0;

  // Summary and diagnostics data:
  double fCrossSection;
  double geomCrossSection;

  int transparent = 0;

  int doinit = 1;
  int type = fullrun;
  bool usingRoot = false;

  // Run type: full or cascade only
  if(strcmp("cascade", argv[1]) == 0) {
    type = inclrun;
  }

  // For ASCII output file
  ofstream out;

  G4Hazard *hazard = new G4Hazard();
  G4VarNtp *varntp = new G4VarNtp();
  G4Calincl *calincl = new G4Calincl();
  G4Ws *ws = new G4Ws();
  G4Mat *mat = new G4Mat();

  // set initial values:
  // First random seed:
  //  hazard->ial = 38035;
  hazard->ial = 979678188;

  // other seeds:
  hazard->igraine[0] = 3997;
  hazard->igraine[1] = 15573;
  hazard->igraine[2] = 9971;
  hazard->igraine[3] = 9821; 
  hazard->igraine[4] = 99233; 
  hazard->igraine[5] = 11167; 
  hazard->igraine[6] = 12399;
  hazard->igraine[7] = 11321; 
  hazard->igraine[8] = 9825;
  hazard->igraine[9] = 2587; 
  hazard->igraine[10] = 1775;
  hazard->igraine[11] = 56799; 
  hazard->igraine[12] = 1156;
  //  hazard->igraine[13] = 11207;
  hazard->igraine[13] = 38957; 
  hazard->igraine[14] = 35779; 
  hazard->igraine[15] = 10055; 
  hazard->igraine[16] = 76533; 
  hazard->igraine[17] = 33759;
  hazard->igraine[18] = 13227;

  G4FermiBreakUp *fermiBreakUp = new G4FermiBreakUp();

  G4Incl *incl = new G4Incl(hazard, calincl, ws, mat, varntp);
  incl->setVerboseLevel(1);

#ifdef USEROOT 
  TFile *dataFile = NULL;

  TString rootfilename(argv[7]);
  if(rootfilename.Contains(".root") == 1) {
    dataFile = new TFile(rootfilename, "RECREATE");
    usingRoot = true;
  }

#ifdef G4INCLDEBUG
  // Create a logger
  G4VInclLogger *theLogger = new G4InclRootLogger();

  // INCL variables
  theLogger->bookHistogram1D("bimpact", 100, -10.0, 10.0);
  theLogger->bookHistogram1D("mzini", 210, 0.0, 210.0);
  theLogger->bookHistogram1D("interpolationResult", 600, -1.0, 30.0);
  theLogger->bookHistogram2D("interpolationPoints", 600, 0.0, 20.0, 600, 0.0, 1.0);
  theLogger->bookHistogram1D("r_distrib", 900, 0.0, 15.0);   // CALL HBOOK1(9,'r_space',200,0.,20.,0.)
  theLogger->bookHistogram1D("p_distrib", 900, 0.0, 300.0);  // CALL HBOOK1(10,'p_space',200,0.,400.,0.)
  theLogger->bookHistogram2D("r-p_correl", 50, 0.0, 12.0, 50, 0.0,500.0);
  theLogger->bookHistogram1D("refInterpolXV", 800, -2.0, 6.0);
  theLogger->bookHistogram1D("refInterpolRes", 600, -1.0, 30.0);
  theLogger->bookHistogram1D("refInterpolResWhenXVgt1", 600, -1.0, 30.0);
  theLogger->bookHistogram1D("totalCXResult", 1000, -10.0, 200.0);
  theLogger->bookHistogram1D("totalLowECX", 1000, -10.0, 200.0);
  theLogger->bookHistogram2D("ppCX", 1000, -10.0, 10000.0, 1000, -10.0, 200.0);
  theLogger->bookHistogram2D("pnCX", 1000, -10.0, 10000.0, 1000, -10.0, 200.0);
  theLogger->bookHistogram2D("nnCX", 1000, -10.0, 10000.0, 1000, -10.0, 200.0);
  theLogger->bookHistogram2D("NDCX", 1000, -10.0, 10000.0, 1000, -10.0, 200.0);
  theLogger->bookHistogram2D("DDCX", 1000, -10.0, 10000.0, 1000, -10.0, 200.0);

  // ABLA variables
  theLogger->bookHistogram1D("pcorem", 1000, 0.0, 2000.0);

  incl->registerLogger(theLogger);
#endif

  // A tree into which we will put the output data
  TTree *h101 = new TTree("h101", "Data from INCL4+ABLA");

  const Int_t maxpart = 255;
  // Variables for Ntuple (TTree)
  Int_t Event;
  Int_t BulletType;
  Double_t BulletE;

  Int_t Massini, Mzini; 
  Double_t Exini;
  Double_t Pcorem, Mcorem, Pxrem, Pyrem, Pzrem;
  Int_t Mulncasc, Mulnevap,Mulntot;

  Int_t Masp, Mzsp;
  Double_t Exsp;

  Double_t Bimpact;
  Int_t Jremn, Kfis;
  Double_t Estfis;
  Int_t Izfis, Iafis;
  Int_t Ntrack;
  Int_t baryonNumber;
  Int_t Itypcasc[maxpart], Avv[maxpart], Zvv[maxpart];
  Double_t Enerj[maxpart], Plab[maxpart], Tetlab[maxpart], Philab[maxpart];
  Double_t momX[maxpart], momY[maxpart], momZ[maxpart];
  Double_t Tsp, Pxsp, Pysp, Pzsp;

  // With these branches (variables)
  h101->Branch("Event", &Event, "Event/I");
  h101->Branch("BulletType", &BulletType, "BulletType/I");
  h101->Branch("BulletE", &BulletE, "BulletE/D");

  h101->Branch("Masp", &Masp, "Masp/I");
  h101->Branch("Mzsp", &Mzsp, "Mzsp/I"); 
  h101->Branch("Exsp", &Exsp, "Exsp/D");
  h101->Branch("Tsp", &Tsp, "Tsp/D");
  h101->Branch("Pxsp", &Pxsp, "Pxsp/D");
  h101->Branch("Pysp", &Pysp, "Pysp/D");
  h101->Branch("Pzsp", &Pzsp, "Pzsp/D");

  h101->Branch("Massini", &Massini, "Massini/I");
  h101->Branch("Mzini", &Mzini, "Mzini/I"); 
  h101->Branch("Exini", &Exini, "Exini/D");
  h101->Branch("Pcorem", &Pcorem, "Pcorem/D");
  h101->Branch("Mcorem", &Mcorem, "Mcorem/D");
  h101->Branch("Pxrem", &Pxrem, "Pxrem/D");
  h101->Branch("Pyrem", &Pyrem, "Pyrem/D");
  h101->Branch("Pzrem", &Pzrem, "Pzrem/D");
  h101->Branch("Mulncasc", &Mulncasc, "Mulncasc/I");
  h101->Branch("Mulnevap", &Mulnevap, "Mulnevap/I");
  h101->Branch("Mulntot", &Mulntot, "Mulntot/I");

  h101->Branch("Bimpact", &Bimpact, "Bimpact/D");
  h101->Branch("Jremn", &Jremn, "Jremn/I");
  h101->Branch("Kfis", &Kfis, "Kfis/I");
  h101->Branch("Estfis", &Estfis, "Estfis/D");
  h101->Branch("Izfis", &Izfis, "Izfis/I");
  h101->Branch("Iafis", &Iafis, "Iafis/I");

  h101->Branch("baryonNumber", &baryonNumber, "baryonNumber/I");
  h101->Branch("Ntrack", &Ntrack, "Ntrack/I");
  h101->Branch("Ityp", Itypcasc, "Ityp[Ntrack]/I");
  h101->Branch("Avv", Avv, "Avv[Ntrack]/I"); 
  h101->Branch("Zvv", Zvv, "Zvv[Ntrack]/I");
  h101->Branch("Enerj", Enerj, "Enerj[Ntrack]/D");
  h101->Branch("Plab", Plab, "Plab[Ntrack]/D");
  h101->Branch("Tetlab", Tetlab, "Tetlab[Ntrack]/D");
  h101->Branch("Philab", Philab, "Philab[Ntrack]/D");
  h101->Branch("momX", momX, "momX[Ntrack]/D");
  h101->Branch("momY", momY, "momY[Ntrack]/D");
  h101->Branch("momZ", momZ, "momZ[Ntrack]/D");
  
#endif // USEROOT

  if (filename == NULL)
  {
    cout << "Memory allocation error." << endl;
    return -1;
  }

  if(argc < 8) {
    cout <<"Usage: incl <cascade|full> targetA targetZ bulletType bulletEnergy events outputfile" << endl;
    return -1;
  }

  if(!usingRoot) {
    strcpy (filename, argv[7]);
    // Open the output file:
    out.open(filename);
  }

  // Initialize FINPUT:
  for(int i = 0; i < 15; i++) {
    calincl->f[i] = 0.0;
  }
  // End of Initialization

  // Input parameters:
  // targetA, targetZ, bulletType, bulletE, filename
  // Usage: cppinterfacetest targetA targetZ bulletType bulletE filename

  // Set input parameters:

  // Target nucleus:
  // Mass number:
  //FINPUT(1)
  calincl->f[0] = atof(argv[2]);
  // Charge number:
  // FINPUT(2)
  calincl->f[1] = atof(argv[3]);

  int baryonNumberBalanceInINCL = calincl->f[0];
  int chargeNumberBalanceInINCL = calincl->f[1];

  if(calincl->f[0] < 17) useFermiBreakup = true;

  // Bullet:
  // Bullet type (1.00 = proton):
  // FINPUT(7)
  calincl->f[6] = atof(argv[4]);
//       IF(IA_BE.EQ.1.AND.IZ_BE.EQ.1) finput(7)=1 
//       IF(IA_BE.EQ.1.AND.IZ_BE.EQ.0) finput(7)=2 
//       IF(IA_BE.EQ.-1.AND.IZ_BE.EQ.1) finput(7)=3 
//       IF(IA_BE.EQ.-1.AND.IZ_BE.EQ.0) finput(7)=4 
//       IF(IA_BE.EQ.-1.AND.IZ_BE.EQ.-1) finput(7)=5 
//       IF(IA_BE.EQ.2.AND.IZ_BE.EQ.1) finput(7)=6 
//       IF(IA_BE.EQ.3.AND.IZ_BE.EQ.1) finput(7)=7 
//       IF(IA_BE.EQ.3.AND.IZ_BE.EQ.2) finput(7)=8 
//       IF(IA_BE.EQ.4.AND.IZ_BE.EQ.2) finput(7)=9 
//       IF(IA_BE.GT.4) finput(7)= -IA_BE
  if(calincl->f[6] == 1) { // Proton
    baryonNumberBalanceInINCL++;
    chargeNumberBalanceInINCL++;
  } else if (calincl->f[6] == 2) { // Neutron
    baryonNumberBalanceInINCL++;
  } else if (calincl->f[6] == 3) { // Pi+
    chargeNumberBalanceInINCL++;
  } else if (calincl->f[6] == 5) { // Pi-
    chargeNumberBalanceInINCL--;
  } else if (calincl->f[6] == 6) { // Deuteron
    baryonNumberBalanceInINCL = baryonNumberBalanceInINCL + 2;
    chargeNumberBalanceInINCL++;
  } else if (calincl->f[6] == 7) { // Triton
    baryonNumberBalanceInINCL = baryonNumberBalanceInINCL + 3;
    chargeNumberBalanceInINCL++;
  } else if (calincl->f[6] == 8) { // He3
    baryonNumberBalanceInINCL = baryonNumberBalanceInINCL + 3;
    chargeNumberBalanceInINCL = chargeNumberBalanceInINCL + 2;
  } else if (calincl->f[6] == 9) { // Alpha
    baryonNumberBalanceInINCL = baryonNumberBalanceInINCL + 4;
    chargeNumberBalanceInINCL = chargeNumberBalanceInINCL + 2;
  } else if(calincl->f[6] < 0) { // Carbon-12
    baryonNumberBalanceInINCL = baryonNumberBalanceInINCL + 12;
    chargeNumberBalanceInINCL = chargeNumberBalanceInINCL + 6;
  }    

  if(calincl->f[6] < 0 || calincl->f[6] > 5) useProjSpect = true;

  // Bullet energy:
  // FINPUT(3)
  calincl->f[2] = atof(argv[5]);

  // Run settings:
  // Time scaling:
  // FINPUT(6)
  calincl->f[5] = 1.0;

  // Nuclear potential:
  // FINPUT(5)
  calincl->f[4] = 45.0;

  // NOSURF:
  ws->nosurf = -2;

  // XFOISA
  ws->xfoisa = 10;

  // NPAULSTR
  ws->npaulstr = 0;

  // Events (only for compatibility)
  // Deprecated
  calincl->icoup = 1;

  // Events:
  int totalevents = atoi(argv[6]);

  // Dresner initialization
  int iflag_initdres = 0;
  // CALL init_dresner(RACINE)
  // iflag_initdres = 1;
  int particleI = 0;
  
  // End of input parameters

  cout << ";; Outputfile: " << argv[7] << endl;
  cout << ";; Projectile: " << endl;
  cout << ";; Type: " << calincl->f[6] << endl;
  cout << ";; energy: " << calincl->f[2] << " mev " << endl;
  cout << ";; target: " << endl;
  cout << ";; a: " << calincl->f[0] << endl;
  cout << ";; z: " << calincl->f[1] << endl;
  cout << ";; events: " << argv[6] << endl;
  cout << endl;

  cout << ";; Running..." << endl;

  if(useProjSpect || useFermiBreakup) {
    incl->setUseProjectileSpectators(true);
    incl->setUseFermiBreakUp(true);
  }

  mat->nbmat = 1;
  mat->amat[0] = int(calincl->f[0]);
  mat->zmat[0] = int(calincl->f[1]);

  cout <<";; Initializing INCL" << endl;
  incl->initIncl(false);
  cout <<";; Initialization complete." << endl;

  starttime = clock();
  for(int n = 1; n <= totalevents; n++) {
   // std::cout <<"Starting event " << n << std::endl;
    if(n == 1) {
      doinit = 1;
    } 
    else {
      doinit = 0;
    }

    if(type == fullrun) {
      incl->processEventInclAbla(n);
    } else if(type == inclrun) {
      incl->processEventIncl();
    } else {
      cout <<";; Unknown run type." << endl;
      exit(1);
    }
    // The coulomb transparency:
    if(n == 1) {
      //      totalevents = totalevents - debugval_.ntranscoul;
      totalevents = totalevents - 1;
    }

    //   *      1   * I*4  *             * VARNTP  * massini	     A of the remnant
    // *      2   * I*4  *             * VARNTP  * MZINI	     Z    "        "
    // *      3   * R*4  *             * VARNTP  * EXINI	     Excit energy " "
    // *      4   * I*4  *             * VARNTP  * MULNCASC	     Cascade n multip.
    // *      5   * I*4  *             * VARNTP  * MULNEVAP	     Evapo   "      "
    // *      6   * I*4  *             * VARNTP  * MULNTOT	     Total   "      "
    // *      7   * R*4  *             * VARNTP  * BIMPACT	     Impact parameter
    // *      8   * I*4  *             * VARNTP  * JREMN	     Remnant Intrinsic Spin
    // *      9   * I*4  *             * VARNTP  * KFIS	     Fission 1/0=Y/N
    // *     10   * R*4  *             * VARNTP  * ESTFIS		Excit energy at fis
    // *     11   * I*4  *             * VARNTP  * IZFIS		Z of fiss nucleus
    // *     12   * I*4  *             * VARNTP  * IAFIS		A of "          "
    // *     13   * I*4  *[0,250]      * VARNTP  * NTRACK		Number of particles
    // *     14   * I*4  *             * VARNTP  * ITYP(NTRACK)	emitted in cascade (0)
    // *                                                                       or evapo   (1)
    // *     15   * I*4  *             * VARNTP  * AVV(NTRACK)		A (-1 for pions)
    // *     16   * I*4  *             * VARNTP  * ZVV(NTRACK)		Z
    // *     17   * R*4  *             * VARNTP  * ENERJ(NTRACK)		kinetic energy
    // *     18   * R*4  *             * VARNTP  * PLAB(NTRACK)		momentum
    // *     19   * R*4  *             * VARNTP  * TETLAB(NTRACK)		Theta (deg)
    // *     20   * R*4  *             * VARNTP  * PHILAB(NTRACK)		Phi   (deg)

// C *************************  DEBUT DRESNER ************************ 
// 	iia=iarem
// 	iiz=izrem
// 	exc=esrem
// 	IF(i.EQ.ievtest)  write(6,*)'dresner remnant',iia,iiz,exc,i
//         nb_dresner1=0
// 	call dresner(iia,iiz,exc,multen,multep,i)
// 	IF(i.EQ.ievtest)  write(6,*)'sortie dresner remn',nombre1
//         nb_dresner1=nombre1
//         NTRACKdeb1_sav=NTRACK
	
// 	IF(nb_dresner1.EQ.1) THEN
//       	NTRACK=NTRACK+1		! on recopie le remnant dans le ntuple
//         ITYPCASC(NTRACK)=1
//       	AVV(NTRACK)=IAREM
//       	ZVV(NTRACK)=IZREM
//       	PLAB(NTRACK)=PCOREM
//       	ENERJ(NTRACK)=SQRT(PCOREM**2+MCOREM**2)-MCOREM
// 	TETLAB(NTRACK)=180.*ACOS(GAREM)/3.141592654
// 	PHILAB(NTRACK)=180.*ATAN2(BEREM,ALREM)/3.141592654
//        ELSE
// C 
// C Dresner donne nombre1 particules (<200) dans les tableaux:
// C    acv1    A
// C    zcv1    Z
// C    ecv1    Kinetic energy
// C    tcv1    polar angle theta
// C    pcv1    polar angle phi (added 3/2008 AB)
// C  (No phi value... no chance to check momentum conservation in Dresner!) 
// C ATTENTION, il y a des particules avec que des ZERO (separateurs de DRESNER)

// cDRES  npart0(j,1) is the particle counter for the j-th type, for any            
// cDRES  evaporation of particles from the original nucleus before fissioni        
// cDRES  a zero value is inserted in energy tables after evap. products            
// cDRES  from original nucleus                                                     
// cDRES  npart0(j,2) is the particle counter for the j-th type for                 
// cDRES  evaporation of the first fragment                                         
// cDRES  a zero value will be inserted between evaporation product                 
// cDRES  energies of first and second fragments                                    
// cDRES  consequently npart(j) will all be a value of 2 greater than true          
// cDRES  value since zero also inserted after original nucleus evaporation  
// C
//       if(nombre1.gt.0) then 	
// c      WRITE(6,*) 'run: ',i

//       PXD=0.
//       PYD=0.
//       PZD=0.
//       EDRES=0.
//       ADRES=0.
//       ZDRES=0.
//       IFLAGRAP=0
            
//       nbpart_ss=0
//       NTRACKdeb=NTRACK
//  	do kcv=1,nombre1
// 	IF(i.EQ.ievtest) WRITE(6,*) 'k,A,Z:',kcv,acv1(kcv),zcv1(kcv)
// C Here we don't copy separators of Dresner in the ntuple:
// 	IF(acv1(kcv).NE.0.OR.zcv1(kcv).NE.0) THEN
// 	nbpart_ss = nbpart_ss+1		!(nombre1 includes separator=null particles)	
//  	NTRACK=NTRACK+1
	
// c      WRITE(6,*) 'Primary values from DRESNER:'
// c      WRITE(6,*) kcv,acv1(kcv),zcv1(kcv),ecv1(kcv),tcv1(kcv),pcv1(kcv)
//         ITYPCASC(NTRACK)=0
//         AVV(NTRACK)=acv1(kcv)
//         ZVV(NTRACK)=zcv1(kcv)
// 	ENERJ(NTRACK)=ecv1(kcv)
// 	TETLAB(NTRACK)=tcv1(kcv)
// 	PHILAB(NTRACK)=pcv1(kcv)
// C Phi not given by dresner!

// C On ajoute au moins l'impulsion! (AB 3/2008):
//       	IF (AVV(NTRACK).EQ.1.AND.ZVV(NTRACK).EQ.1) THEN
// 	RMASS=938.2723
// 	ELSE IF(AVV(NTRACK).EQ.1.AND.ZVV(NTRACK).EQ.0) THEN
// 	RMASS=939.5657
// 	ELSE IF(AVV(NTRACK).EQ.2.AND.ZVV(NTRACK).EQ.1) THEN
// 	RMASS=1874.34
// 	ELSE IF(AVV(NTRACK).EQ.3.AND.ZVV(NTRACK).EQ.1) THEN
// 	RMASS=2806.359
// 	ELSE IF(AVV(NTRACK).EQ.3.AND.ZVV(NTRACK).EQ.2) THEN
// 	RMASS=2807.119
// 	ELSE IF(AVV(NTRACK).EQ.4.AND.ZVV(NTRACK).EQ.2) THEN
// 	RMASS=3724.818
//       	ELSE
// 	APRF=AVV(NTRACK)
// 	ZPRF=ZVV(NTRACK)
// C Essais avec la masse de KHS (9/2002):
// 	CALL MGLMS(APRF,ZPRF,0,EL)
//         RMASS = ZPRF*FMP + (APRF-ZPRF)*FMN + EL
//       	ENDIF
// 	PLAB(NTRACK)=SQRT(ENERJ(NTRACK)*(ENERJ(NTRACK)+2.*RMASS))   
//         UCT=PLAB(NTRACK)*SIN(TETLAB(NTRACK)*3.14159/180.)
	
// C For pure DRESNER momentum conservation (not satisfied!):
// c	PX=UCT
// c     s     *COS(PHILAB(NTRACK)*3.14159/180.)	
// c	PY=UCT
// c     s     *SIN(PHILAB(NTRACK)*3.14159/180.)	
// c	PZ=PLAB(NTRACK)*COS(TETLAB(NTRACK)*3.14159/180.)
// c	PXD=PXD+PX
// c	PYD=PYD+PY
// c	PZD=PZD+PZ
// c	WRITE(6,*) 'Pure Dresner:',PX,PY,PZ
	
// C Phi random... (not given by dresner) or exact for 2 and 3 particles!
// 	IF(nbpart_ss.EQ.1) THEN
// 	U1=UCT
// 	ikeep1=NTRACK
// 	CALL RIBM(rndm,IY(18))
// 	phikeep=rndm*360.
// 	PHILAB(NTRACK)=phikeep - 180.
// c	WRITE(6,*) 'First p',ikeep1,phikeep,PHILAB(ikeep1)
// 	ENDIF
	
// 	IF(nbpart_ss.EQ.2) THEN
// 	U2=UCT
// 	ikeep2=NTRACK
// 	PHILAB(NTRACK)= AMOD(phikeep + 180.,360.) -180.
// c	WRITE(6,*) phikeep-180.,PHILAB(NTRACK)
// 	ENDIF 
	
// 	IF(nbpart_ss.EQ.3) THEN
// 	   U3=UCT
// 	   val_loc3=(-u1**2+u2**2-u3**2)/(2.*u1*u3)
// 	   val_loc2=(-u1**2-u2**2+u3**2)/(2.*u1*u2)
// 	   IF(ABS(val_loc2).LE.1.AND.ABS(val_loc2).LE.1) THEN
// 	   phi_2=ACOS(val_loc2)
// 	   phi_3=ACOS(val_loc3)
// c	   WRITE(6,*) 'OK sol 3'
// c	   WRITE(6,*) 'u1,u2,u3,phi2,phi3:',U1,U2,U3,phi_2,phi_3
// c	   WRITE(6,*) U1,U2*COS(phi_2),U3*COS(phi_3),
// c     s                U2*SIN(phi_2),U3*SIN(phi_3)
// 	   phi_2=phi_2*180./3.14159+phikeep
// 	   phi_3=360.-phi_3*180./3.14159+phikeep
// 	   PHILAB(ikeep2)=AMOD(phi_2,360.)-180.
// 	   PHILAB(NTRACK)=AMOD(phi_3,360.)-180.
	   
// C Verifs:
// c	   WRITE(6,*)phikeep,phi_2,phi_3
// c	   WRITE(6,*) 'ikeep1,ikeep2,NTRACK:',ikeep1,ikeep2,NTRACK
// c	   WRITE(6,*) 'philab',PHILAB(ikeep1),PHILAB(ikeep2),PHILAB(NTRACK)
// c	   WRITE(6,*) 'BILAN:'
// c	   SV1=U1*COS(PHILAB(ikeep1)*3.14159/180.)
// c	   SV2=U2*COS(PHILAB(ikeep2)*3.14159/180.)
// c	   SV3=U3*COS(PHILAB(NTRACK)*3.14159/180.)
// c	   ST=SV1+SV2+SV3
// c	   WRITE(6,*) 'PX:',SV1,SV2,SV3,ST
// c	   SV1=U1*SIN(PHILAB(ikeep1)*3.14159/180.)
// c	   SV2=U2*SIN(PHILAB(ikeep2)*3.14159/180.)
// c	   SV3=U3*SIN(PHILAB(NTRACK)*3.14159/180.)
// c	   ST=SV1+SV2+SV3	    
// c	   WRITE(6,*) 'PY:',SV1,SV2,SV3,ST
// c	   SV1=PLAB(ikeep1)*COS(TETLAB(ikeep1)*3.14159/180.)
// c	   SV2=PLAB(ikeep2)*COS(TETLAB(ikeep2)*3.14159/180.)
// c	   SV3=PLAB(NTRACK)*COS(TETLAB(NTRACK)*3.14159/180.)
// c	   ST=SV1+SV2+SV3	    	   
// c	   WRITE(6,*) 'PZ:',SV1,SV2,SV3,ST 

// 	   ELSE
// 	   CALL RIBM(rndm,IY(17))
// 	   phikeep=rndm*360.
// 	   PHILAB(ikeep2)=phikeep - 180.
// 	   CALL RIBM(rndm,IY(19))
// 	   phikeep=rndm*360.
// 	   PHILAB(NTRACK)=phikeep - 180.	   
// 	   ENDIF	   
// 	ENDIF
	
// 	IF(nbpart_ss.GT.3) THEN
// 	CALL RIBM(rndm,IY(18))
// 	phikeep=rndm*360.
// 	PHILAB(NTRACK)=phikeep - 180.
// 	ENDIF
	
// C Check conservations:
//         ADRES=ADRES+AVV(NTRACK)
//         ZDRES=ZDRES+ZVV(NTRACK)
// 	EDRES=EDRES+RMASS+ENERJ(NTRACK)
// 	TETDRES=TETLAB(NTRACK)*3.14159/180.
// 	PHIDRES=PHILAB(NTRACK)*3.14159/180.
//         UCT=PLAB(NTRACK)*SIN(TETDRES)
// 	PX=UCT*COS(PHIDRES)	
// 	PY=UCT*SIN(PHIDRES)	
// 	PZ=PLAB(NTRACK)*COS(TETDRES)
// 	PXD=PXD+PX
// 	PYD=PYD+PY
// 	PZD=PZD+PZ
	
	
// 	ENDIF             
//         enddo

// C Print for conservation checks:
//       IF(IFLAGRAP.EQ.1) THEN
//       WRITE(6,*) 'AREM,SA,ZREM,SZ:',IAREM,ADRES,IZREM,ZDRES
//       Diff=ESREM+MCOREM-EDRES
//       WRITE(6,*) 'EREM,SE,Diff:',ESREM+MCOREM,EDRES,Diff
//       WRITE(6,*) 'nb,SPX,SPY,SPZ:',nbpart_ss,PXD,PYD,PZD
//       ENDIF
      
// 	nb_dresnerv1=nbpart_ss
//       endif

// C Lorentz boost:
//       PXS=0.
//       PYS=0.
//       PZS=0.
//       PXSB=0.
//       PYSB=0.
//       PZSB=0.
      
// C Velocities of transfo
//       E=ERECREM+MCOREM
//       B1=PXREM/E
//       B2=PYREM/E
//       B3=PZREM/E
      
// 	DO iloc=NTRACKdeb+1,NTRACK
// 	RMASS=(PLAB(iloc)**2-ENERJ(iloc)**2)/(2.*ENERJ(iloc))
// 	E=ENERJ(iloc)+RMASS
//         UCT=PLAB(iloc)*SIN(TETLAB(iloc)*3.14159/180.)
// 	PX=UCT
//      s     *COS(PHILAB(iloc)*3.14159/180.)	
// 	PY=UCT
//      s     *SIN(PHILAB(iloc)*3.14159/180.)	
// 	PZ=PLAB(iloc)*COS(TETLAB(iloc)*3.14159/180.)
// 	PXS=PXS+PX
// 	PYS=PYS+PY
// 	PZS=PZS+PZ
// 	CALL LOREN(PX,PY,PZ,B1,B2,B3,E)
// 	PLAB(iloc)=SQRT(PX**2+PY**2+PZ**2)
// 	IF(PLAB(iloc).LE.0.0001) THEN
// 	TETLAB(iloc)=0.
// 	PHILAB(iloc)=0.
// 	ELSE
// 	TETLAB(iloc)=ACOS(PZ/PLAB(iloc))*180./3.14159
// 	PHILAB(iloc)=ATAN2(PY,PX)*180./3.14159
// 	ENDIF
// 	ENERJ(iloc)=E-RMASS
// 	PXSB=PXSB+PX
// 	PYSB=PYSB+PY
// 	PZSB=PZSB+PZ
	
// 	ENDDO
	
// c	WRITE(6,*) 'Bilan avant boost:',PXS,PYS,PZS
// c	WRITE(6,*) 'Bilan apres boost:',PXSB,PYSB,PZSB
// c	WRITE(6,*) 'Impulsion spectat:',
// c     s                P1_PROJSPEC,P2_PROJSPEC,P3_PROJSPEC
// c        WRITE(6,*) 'Bilan Dresner (sans Boost):',PXD,PYD,PZD

// C endif no fractioning of the remnant     
//       ENDIF

      
//       multeen=multn+multen
//       MULNCASC=multn
//       MULNEVAP=multen
//       MULNTOT=multeen
// c      write(1,*)MULNCASC,MULNEVAP,MULNTOT
//       multeep=multp+multep
// call hfill(4051,multeen,0.,1.)
// call hfill(4052,multeep,0.,1.)
// call hfill(308,multeen,exxc,1.)
// C 21/07/08 AB: Ne semble pas utilise ni ici ni dans dresner.f
// c      do ijk=0,2000
// c      if(tabener(ijk).ne.0.)then
// c      write(6,*)ijk,multeen,tabener(ijk)
// c      aijk=ijk
// c      write(6,*)aijk,multeen,tabener(ijk)
// ccall hfill(310,aijk,multeen,1.)
// c      endif
// c      enddo
	
// C *************************  FIN DRESNER ************************ 
		

    if(varntp->ntrack > 0) {
      //      std::cout <<"Filling event number: " << n << std::endl;
      Masp = varntp->masp;
      double momXsum = 0.0, momYsum = 0.0, momZsum = 0.0;
      Ntrack = 0;
      for(particleI = 1; particleI <= varntp->ntrack; particleI++) {
#ifdef DEBUG
	  momX[particleI] = varntp->plab[particleI]*TMath::Sin(varntp->tetlab[particleI]*TMath::Pi()/180.0)*TMath::Cos(varntp->philab[particleI]*TMath::Pi()/180.0);
	  momY[particleI] = varntp->plab[particleI]*TMath::Sin(varntp->tetlab[particleI]*TMath::Pi()/180.0)*TMath::Sin(varntp->philab[particleI]*TMath::Pi()/180.0);
	  momZ[particleI] = varntp->plab[particleI]*TMath::Cos(varntp->tetlab[particleI]*TMath::Pi()/180.0);
	  momXsum += momX[particleI];
	  momYsum += momY[particleI];
	  momZsum += momZ[particleI];
	  G4cout <<"A = " << varntp->avv[particleI] << " Z = " << varntp->zvv[particleI] << " mom = (" << momX[particleI] << ", " << momY[particleI] << ", " << momZ[particleI] <<")" << G4endl;
#endif
	if(!usingRoot) {
	  out << n << " " << calincl->f[6]  << " " << calincl->f[2] << " ";
	  out << varntp->massini << " " << varntp->mzini << " ";
	  out << varntp->exini << " " << varntp->masp << " " << varntp->mzsp << " " << varntp->mulncasc << " " << varntp->mulnevap << " " << varntp->mulntot << " ";
	  out << varntp->bimpact << " " << varntp->jremn << " " << varntp->kfis << " " << varntp->estfis << " ";
	  out << varntp->izfis << " " << varntp->iafis << " " << varntp->ntrack << " ";
	  out << varntp->itypcasc[particleI] << " ";
	  out << varntp->avv[particleI] << " " << varntp->zvv[particleI] << " " << varntp->enerj[particleI] << " ";
	  out << varntp->plab[particleI] << " " << varntp->tetlab[particleI] << " " << varntp->philab[particleI] << endl;
	}
#ifdef USEROOT
	// Do not store non-physical particles with A = 0 and Z = 0 to the ROOT tree:
	if(varntp->avv[particleI] == 0 && varntp->zvv[particleI] == 0) continue;
	if(usingRoot) {
	  //	  std::cout <<"Filling particle number: " << particleI << std::endl;
	  Event = n;
	  BulletType = (Int_t)calincl->f[6]; 
	  BulletE = calincl->f[2];

	  Masp = int(varntp->masp);
	  Mzsp = int(varntp->mzsp);
	  Exsp = varntp->exsp;

	  Tsp = varntp->spectatorT;
	  Pxsp = varntp->spectatorP1;
	  Pysp = varntp->spectatorP2;
	  Pzsp = varntp->spectatorP3;

	  Massini = int(varntp->massini);
	  Mzini = int(varntp->mzini);
	  Exini = varntp->exini;
	  Pcorem = varntp->pcorem;
	  Mcorem = varntp->mcorem;
	  Pxrem = varntp->pxrem;
	  Pyrem = varntp->pyrem;
	  Pzrem = varntp->pzrem;
	  Mulncasc = varntp->mulncasc;
	  Mulnevap = varntp->mulnevap;
	  Mulntot = varntp->mulntot;
	  Bimpact = varntp->bimpact; 
	  Jremn = varntp->jremn; 
	  Kfis = varntp->kfis;
	  Estfis = varntp->estfis;
	  Izfis = varntp->izfis; 
	  Iafis = varntp->iafis; 
	  baryonNumber = varntp->getTotalBaryonNumber();
	  //	  Ntrack = varntp->ntrack;
	  Itypcasc[Ntrack] = varntp->itypcasc[particleI];
	  Avv[Ntrack] = varntp->avv[particleI];
	  Zvv[Ntrack] = varntp->zvv[particleI];
	  Enerj[Ntrack] = varntp->enerj[particleI];
	  Plab[Ntrack] = varntp->plab[particleI];
	  Tetlab[Ntrack] = varntp->tetlab[particleI];
	  Philab[Ntrack] = varntp->philab[particleI];
	  momX[Ntrack] = varntp->plab[particleI]*TMath::Sin(varntp->tetlab[particleI]*TMath::Pi()/180.0)*TMath::Cos(varntp->philab[particleI]*TMath::Pi()/180.0);
	  momY[Ntrack] = varntp->plab[particleI]*TMath::Sin(varntp->tetlab[particleI]*TMath::Pi()/180.0)*TMath::Sin(varntp->philab[particleI]*TMath::Pi()/180.0);
	  momZ[Ntrack] = varntp->plab[particleI]*TMath::Cos(varntp->tetlab[particleI]*TMath::Pi()/180.0);
	  Ntrack++;
	  //	  std::cout <<"Particle " << particleI << " of " << Ntrack - 1<<" filled." << std::endl;
	}
#endif // USEROOT
      }

      // Use G4 Fermi break-up here:

      if(useProjSpect && varntp->masp > 1) { // Fermi break-up for the spectator nucleus
	baryonNumberBalanceInINCL -= varntp->masp;
	G4double nuclearMass = G4NucleiProperties::GetNuclearMass(G4int(varntp->masp), G4int(varntp->mzsp)) + varntp->exsp * MeV;
	// Use momentum scaling to compensate for different masses in G4 and INCL:
	G4double momentumScaling = (std::sqrt(std::pow(varntp->spectatorT, 2) 
					      + 2.0 * varntp->spectatorT * nuclearMass)) /
	  (std::sqrt(std::pow(varntp->spectatorP1, 2) + std::pow(varntp->spectatorP2, 2)
		     + std::pow(varntp->spectatorP3, 2)));
	G4LorentzVector p4(momentumScaling * varntp->spectatorP1 * MeV, momentumScaling * varntp->spectatorP2 * MeV,
			   momentumScaling * varntp->spectatorP3 * MeV,
			   varntp->spectatorT * MeV + nuclearMass);
	// Four-momentum, baryon number and charge balance:
	G4LorentzVector fourMomentumBalance = p4;
	G4int baryonNumberBalance = varntp->masp;
	chargeNumberBalanceInINCL -= varntp->mzsp;
	G4int chargeBalance = varntp->mzsp;

	G4LorentzRotation toFragmentZ;
	// Assume that Fermi breakup uses Z as the direction of the projectile
	toFragmentZ.rotateZ(-p4.theta());
	toFragmentZ.rotateY(-p4.phi());
	G4LorentzRotation toFragmentLab = toFragmentZ.inverse();
	//      p4 *= toFragmentZ;
      
	G4LorentzVector p4rest = p4;
	p4rest.boost(-p4.boostVector());
	if(verboseLevel > 0) {
	  G4cout <<"Spectator nucleus:" << G4endl;
	  G4cout <<"p4: " << G4endl;
	  G4cout <<" px: " << p4.px() <<" py: " << p4.py() <<" pz: " << p4.pz() << G4endl;
	  G4cout <<" E = " << p4.e() << G4endl;
	  G4cout <<"p4rest: " << G4endl;
	  G4cout <<" px: " << p4rest.px() <<" py: " << p4rest.py() <<" pz: " << p4rest.pz() << G4endl;
	  G4cout <<" E = " << p4rest.e() << G4endl;
	}
	G4Fragment theSpectatorNucleus(G4int(varntp->masp), G4int(varntp->mzsp), p4rest);
	theSpectatorFermiBreakupResult = fermiBreakUp->BreakItUp(theSpectatorNucleus);
	if(theSpectatorFermiBreakupResult != 0) {
	  G4FragmentVector::iterator fragment;
	  for(fragment = theSpectatorFermiBreakupResult->begin(); fragment != theSpectatorFermiBreakupResult->end(); fragment++) {
	    Itypcasc[Ntrack] = -2; // Fragment of the spectator nucleus
	    G4ParticleDefinition *theFragmentDefinition = 0;
	    if((*fragment)->GetA_asInt() == 1.0 && (*fragment)->GetZ_asInt() == 0.0) { // Neutron
	      theFragmentDefinition = G4Neutron::NeutronDefinition();
	      Avv[Ntrack] = 1;
	      Zvv[Ntrack] = 0;
	    } else if ((*fragment)->GetA_asInt() == 1.0 && (*fragment)->GetZ_asInt() == 1.0) {
	      theFragmentDefinition = G4Proton::ProtonDefinition();
	      Avv[Ntrack] = 1;
	      Zvv[Ntrack] = 1;
	    } else {
	      theFragmentDefinition = theTableOfParticles->GetIon(G4int((*fragment)->GetZ_asInt()), G4int((*fragment)->GetA_asInt()), (*fragment)->GetExcitationEnergy());
	      Avv[Ntrack] = (*fragment)->GetA_asInt();
	      Zvv[Ntrack] = (*fragment)->GetZ_asInt();
	    }
	    if(theFragmentDefinition != 0) {
	      G4DynamicParticle *theFragment = new G4DynamicParticle(theFragmentDefinition, (*fragment)->GetMomentum());
	      G4LorentzVector labMomentum = theFragment->Get4Momentum();
	      labMomentum.boost(p4.boostVector());
	      //	  labMomentum *= toFragmentLab;
	      //	  labMomentum *= toLabFrame;
	      Enerj[Ntrack] = (labMomentum.e() - theFragmentDefinition->GetPDGMass()) / MeV;
	      Plab[Ntrack] = std::sqrt(std::pow(labMomentum.px() / MeV, 2)
				       + std::pow(labMomentum.py() / MeV, 2)
				       + std::pow(labMomentum.pz() / MeV, 2));
	      Tetlab[Ntrack] = labMomentum.theta() / degree;
	      Philab[Ntrack] = labMomentum.phi() / degree;
	      theFragment->Set4Momentum(labMomentum);
	      fourMomentumBalance -= theFragment->Get4Momentum();
	      baryonNumberBalance -= theFragmentDefinition->GetAtomicMass();
	      chargeBalance -= theFragmentDefinition->GetAtomicNumber();
	      if(verboseLevel > 0) {
		G4cout <<"Resulting fragment: " << G4endl;
		G4cout <<" kinetic energy = " << theFragment->GetKineticEnergy() / MeV << " MeV" << G4endl;
		G4cout <<" momentum = " << theFragment->GetMomentum().mag() / MeV << " MeV" << G4endl;
	      }
	      Ntrack++; // Increment the particle counter
	    } else {
	      G4cout <<"G4InclAblaCascadeInterface: Error. Fragment produced by Fermi break-up does not exist." 
		     << G4endl;
	      G4cout <<"Resulting fragment: " << G4endl;
	      G4cout <<" Z = " << (*fragment)->GetZ_asInt() << G4endl;
	      G4cout <<" A = " << (*fragment)->GetA_asInt() << G4endl;
	      G4cout <<" Excitation : " << (*fragment)->GetExcitationEnergy() / MeV << " MeV" << G4endl;
	      G4cout <<" momentum = " << (*fragment)->GetMomentum().mag() / MeV << " MeV" << G4endl;
	    }
	  }
	  if(std::abs(fourMomentumBalance.mag() / MeV) > 0.1 * MeV) { 
	    G4cout <<"Four-momentum balance after spectator nucleus Fermi break-up:" << G4endl;
	    G4cout <<"Magnitude: " << fourMomentumBalance.mag() / MeV << " MeV" << G4endl;
	    G4cout <<"Vector components (px, py, pz, E) = ("
		   << fourMomentumBalance.px() << ", "
		   << fourMomentumBalance.py() << ", "
		   << fourMomentumBalance.pz() << ", "
		   << fourMomentumBalance.e() << ")" << G4endl;
	  }
	  if(baryonNumberBalance != 0) {
	    G4cout <<"Baryon number balance after spectator nucleus Fermi break-up: " << baryonNumberBalance << G4endl;
	  }
	  if(chargeBalance != 0) {
	    G4cout <<"Charge balance after spectator nucleus Fermi break-up: " << chargeBalance << G4endl;
	  }
	}
      }

      if(varntp->needsFermiBreakup) { // Fermi break-up for the remnant
	baryonNumberBalanceInINCL -= varntp->massini;
	chargeNumberBalanceInINCL -= varntp->mzini;
	// Call Fermi Break-up
	G4double nuclearMass = G4NucleiProperties::GetNuclearMass(G4int(varntp->massini), G4int(varntp->mzini)) + varntp->exini * MeV;
	G4LorentzVector fragmentMomentum(varntp->pxrem * MeV, varntp->pyrem * MeV, varntp->pzrem * MeV,
					 varntp->erecrem * MeV + nuclearMass);
	G4double momentumScaling = (std::sqrt(std::pow(varntp->erecrem, 2) 
					      + 2.0 * varntp->erecrem * nuclearMass)) /
	  (std::sqrt(std::pow(varntp->pxrem, 2) + std::pow(varntp->pyrem, 2)
		     + std::pow(varntp->pzrem, 2)));

	G4LorentzVector p4(momentumScaling * varntp->pxrem * MeV, momentumScaling * varntp->pyrem * MeV,
			   momentumScaling * varntp->pzrem * MeV,
			   varntp->erecrem + nuclearMass);

	// For four-momentum, baryon number and charge conservation check:
	G4LorentzVector fourMomentumBalance = p4;
	G4int baryonNumberBalance = varntp->massini;
	G4int chargeBalance = varntp->mzini;

	G4LorentzRotation toFragmentZ;
	toFragmentZ.rotateZ(-p4.theta());
	toFragmentZ.rotateY(-p4.phi());
	G4LorentzRotation toFragmentLab = toFragmentZ.inverse();
	//      p4 *= toFragmentZ;

	G4LorentzVector p4rest = p4;
	p4rest.boost(-p4.boostVector());
	if(verboseLevel > 0) {
	  G4cout <<"Cascade remnant nucleus:" << G4endl;
	  G4cout <<"p4: " << G4endl;
	  G4cout <<" px: " << p4.px() <<" py: " << p4.py() <<" pz: " << p4.pz() << G4endl;
	  G4cout <<" E = " << p4.e() << G4endl;

	  G4cout <<"p4rest: " << G4endl;
	  G4cout <<" px: " << p4rest.px() <<" py: " << p4rest.py() <<" pz: " << p4rest.pz() << G4endl;
	  G4cout <<" E = " << p4rest.e() << G4endl;
	}

	G4Fragment theCascadeRemnant(G4int(varntp->massini), G4int(varntp->mzini), p4rest);
	theFermiBreakupResult = fermiBreakUp->BreakItUp(theCascadeRemnant);
	if(theFermiBreakupResult != 0) {
	  G4FragmentVector::iterator fragment;
	  for(fragment = theFermiBreakupResult->begin(); fragment != theFermiBreakupResult->end(); fragment++) {
	    G4ParticleDefinition *theFragmentDefinition = 0;
	    Itypcasc[Ntrack] = 0; // Fragment of the remnant
	    if((*fragment)->GetA_asInt() == 1.0 && (*fragment)->GetZ_asInt() == 0.0) { // Neutron
	      theFragmentDefinition = G4Neutron::NeutronDefinition();
	      Avv[Ntrack] = 1;
	      Zvv[Ntrack] = 0;
	    } else if ((*fragment)->GetA_asInt() == 1.0 && (*fragment)->GetZ_asInt() == 1.0) {
	      theFragmentDefinition = G4Proton::ProtonDefinition();
	      Avv[Ntrack] = 1;
	      Zvv[Ntrack] = 1;
	    } else {
	      theFragmentDefinition = theTableOfParticles->GetIon(G4int((*fragment)->GetZ_asInt()), G4int((*fragment)->GetA_asInt()), (*fragment)->GetExcitationEnergy());
	      Avv[Ntrack] = (*fragment)->GetA_asInt();
	      Zvv[Ntrack] = (*fragment)->GetZ_asInt();
	    }

	    if(theFragmentDefinition != 0) {
	      G4DynamicParticle *theFragment = new G4DynamicParticle(theFragmentDefinition, (*fragment)->GetMomentum());
	      G4LorentzVector labMomentum = theFragment->Get4Momentum();
	      labMomentum.boost(p4.boostVector());
	      //	  labMomentum *= toFragmentLab;
	      //	  labMomentum *= toLabFrame;
	      Enerj[Ntrack] = (labMomentum.e() - theFragmentDefinition->GetPDGMass()) / MeV;
	      Plab[Ntrack] = std::sqrt(std::pow(labMomentum.px() / MeV, 2)
				       + std::pow(labMomentum.py() / MeV, 2)
				       + std::pow(labMomentum.pz() / MeV, 2));
	      Tetlab[Ntrack] = labMomentum.theta() / degree;
	      Philab[Ntrack] = labMomentum.phi() / degree;
	      theFragment->Set4Momentum(labMomentum);
	      fourMomentumBalance -= theFragment->Get4Momentum();
	      baryonNumberBalance -= theFragmentDefinition->GetAtomicMass();
	      chargeBalance -= theFragmentDefinition->GetAtomicNumber();
	      if(verboseLevel > 0) {
		G4cout <<"Resulting fragment: " << G4endl;
		G4cout <<" kinetic energy = " << theFragment->GetKineticEnergy() / MeV << " MeV" << G4endl;
		G4cout <<" momentum = " << theFragment->GetMomentum().mag() / MeV << " MeV" << G4endl;
	      }
	      Ntrack++; // Increment the particle counter
	    } else {
	      G4cout <<"G4InclAblaCascadeInterface: Error. Fragment produced by Fermi break-up does not exist." << G4endl;
	      G4cout <<"Resulting fragment: " << G4endl;
	      G4cout <<" Z = " << (*fragment)->GetZ_asInt() << G4endl;
	      G4cout <<" A = " << (*fragment)->GetA_asInt() << G4endl;
	      G4cout <<" Excitation : " << (*fragment)->GetExcitationEnergy() / MeV << " MeV" << G4endl;
	      G4cout <<" momentum = " << (*fragment)->GetMomentum().mag() / MeV << " MeV" << G4endl;
	    }
	  }
	  if(std::abs(fourMomentumBalance.mag() / MeV) > 0.1 * MeV) { 
	    G4cout <<"Four-momentum balance after remnant nucleus Fermi break-up:" << G4endl;
	    G4cout <<"Magnitude: " << fourMomentumBalance.mag() / MeV << " MeV" << G4endl;
	    G4cout <<"Vector components (px, py, pz, E) = ("
		   << fourMomentumBalance.px() << ", "
		   << fourMomentumBalance.py() << ", "
		   << fourMomentumBalance.pz() << ", "
		   << fourMomentumBalance.e() << ")" << G4endl;
	  }
	  if(baryonNumberBalance != 0) {
	    G4cout <<"Baryon number balance after remnant nucleus Fermi break-up: " << baryonNumberBalance << G4endl;
	  }
	  if(chargeBalance != 0) {
	    G4cout <<"Charge balance after remnant nucleus Fermi break-up: " << chargeBalance << G4endl;
	  }
	}
      }

// C Deexcitation of the projectile spectators if any:
//       IF(finput(7).LT.0.AND.A_PROJSPEC.NE.0) THEN                                    
// C Evapo or Breakup of the projectile spectators through Dresner deexcitation:
// C *************************  DEBUT DRESNER ************************ 
// C In the NTUPLE:
//         MASP=A_PROJSPEC
// 	MZSP=Z_PROJSPEC
// 	EXSP=EX_PROJSPEC
// 	iia=A_PROJSPEC
// 	iiz=Z_PROJSPEC
// 	exc=EX_PROJSPEC
//         nb_dresner2=0
// c	WRITE(6,*) 'Fermi-breakup,b,a,z,Ex/A:',
// c     s                  INUM,bimpact,iia,iiz,exc/iia
// 	call dresner(iia,iiz,exc,multen,multep,i)
// 	nb_dresner2=nombre1
//         NTRACKdeb2_sav=NTRACK
	
// 	IF(i.EQ.ievtest) THEN
// 	WRITE(6,*) 'DRESNER2 output nb_part:',nombre1
// 	DO iloc=1,nombre1
// 	WRITE(6,*) 'A,Z,T,Tet,Phi:',acv1(iloc),zcv1(iloc),ecv1(iloc),
//      s   tcv1(iloc),pcv1(iloc)  
// 	ENDDO
// 	ENDIF

//         IFLAGRAP=0
		
// 	IF(nb_dresner2.EQ.1.AND.ecv1(1).EQ.0.) THEN
//       	NTRACK=NTRACK+1		! on recopie l'unique spect-projo dans le ntuple
//         ITYPCASC(NTRACK)=-2
//         AVV(NTRACK)=A_PROJSPEC
//         ZVV(NTRACK)=Z_PROJSPEC
// 	ENERJ(NTRACK)=T_PROJSPEC
// 	P_PROJSPEC=SQRT(P1_PROJSPEC**2+P2_PROJSPEC**2+P3_PROJSPEC**2)
// 	PLAB(NTRACK)=P_PROJSPEC
// 	TETLAB(NTRACK)=180.*ACOS(P3_PROJSPEC/P_PROJSPEC)/3.141592654
// 	PHILAB(NTRACK)=180.*ATAN2(P2_PROJSPEC,P1_PROJSPEC)/3.141592654
		
// C (endif after the lorentz boost not needed for only one particle=spect-projo)	
//         ELSE
	  
// C 
// C Dresner donne nombre1 particules (<200) dans les tableaux:
// C    acv1    A
// C    zcv1    Z
// C    ecv1    Kinetic energy
// C    tcv1    polar angle theta
// C    pcv1    polar angle phi (added 3/2008 AB)
// C  (No phi value... no chance to check momentum conservation in Dresner!) 
// C ATTENTION, il y a des particules avec que des ZERO (separateurs de DRESNER)

// cDRES  npart0(j,1) is the particle counter for the j-th type, for any            
// cDRES  evaporation of particles from the original nucleus before fissioni        
// cDRES  a zero value is inserted in energy tables after evap. products            
// cDRES  from original nucleus                                                     
// cDRES  npart0(j,2) is the particle counter for the j-th type for                 
// cDRES  evaporation of the first fragment                                         
// cDRES  a zero value will be inserted between evaporation product                 
// cDRES  energies of first and second fragments                                    
// cDRES  consequently npart(j) will all be a value of 2 greater than true          
// cDRES  value since zero also inserted after original nucleus evaporation         

// C
//       IF(nombre1.gt.0) THEN
// c      WRITE(6,*) 'run: ',i

//       PXD=0.
//       PYD=0.
//       PZD=0.
            
//       nbpart_ss=0

//       NTRACKdeb=NTRACK
//  	DO kcv=1,nombre1
// C Here we don't copy separators of Dresner in the ntuple:
// 	IF(acv1(kcv).NE.0.OR.zcv1(kcv).NE.0) THEN
// 	nbpart_ss = nbpart_ss+1		!(nombre1 includes separator=null particles)	
//  	NTRACK=NTRACK+1
	
// c      WRITE(6,*) 'Primary values from DRESNER:'
// c      WRITE(6,*) kcv,acv1(kcv),zcv1(kcv),ecv1(kcv),tcv1(kcv),pcv1(kcv)
//         ITYPCASC(NTRACK)=-2
//         AVV(NTRACK)=acv1(kcv)
//         ZVV(NTRACK)=zcv1(kcv)
// 	ENERJ(NTRACK)=ecv1(kcv)
// 	TETLAB(NTRACK)=tcv1(kcv)
// 	PHILAB(NTRACK)=pcv1(kcv)
// C Phi not given by dresner!
	
// C On ajoute au moins l'impulsion! (AB 3/2008):
//       	IF (AVV(NTRACK).EQ.1.AND.ZVV(NTRACK).EQ.1) THEN
// 	RMASS=938.2723
// 	ELSE IF(AVV(NTRACK).EQ.1.AND.ZVV(NTRACK).EQ.0) THEN
// 	RMASS=939.5657
// 	ELSE IF(AVV(NTRACK).EQ.2.AND.ZVV(NTRACK).EQ.1) THEN
// 	RMASS=1874.34
// 	ELSE IF(AVV(NTRACK).EQ.3.AND.ZVV(NTRACK).EQ.1) THEN
// 	RMASS=2806.359
// 	ELSE IF(AVV(NTRACK).EQ.3.AND.ZVV(NTRACK).EQ.2) THEN
// 	RMASS=2807.119
// 	ELSE IF(AVV(NTRACK).EQ.4.AND.ZVV(NTRACK).EQ.2) THEN
// 	RMASS=3724.818
//       	ELSE
// 	APRF=AVV(NTRACK)
// 	ZPRF=ZVV(NTRACK)
// C Essais avec la masse de KHS (9/2002):
// 	CALL MGLMS(APRF,ZPRF,0,EL)
//         RMASS = ZPRF*FMP + (APRF-ZPRF)*FMN + EL
//       	ENDIF
// 	if(i.EQ.ievtest) WRITE(6,*) 'Proj Spect M,A,Z:',
//      s                     RMASS,AVV(NTRACK),ZVV(NTRACK)
// 	PLAB(NTRACK)=SQRT(ENERJ(NTRACK)*(ENERJ(NTRACK)+2.*RMASS))   
	
//         UCT=PLAB(NTRACK)*SIN(TETLAB(NTRACK)*3.14159/180.)
	
// C For pure DRESNER momentum conservation (not satisfied!):
// 	PX=UCT
//      s     *COS(PHILAB(NTRACK)*3.14159/180.)	
// 	PY=UCT
//      s     *SIN(PHILAB(NTRACK)*3.14159/180.)	
// 	PZ=PLAB(NTRACK)*COS(TETLAB(NTRACK)*3.14159/180.)
// 	PXD=PXD+PX
// 	PYD=PYD+PY
// 	PZD=PZD+PZ
// c	WRITE(6,*) 'Pure Dresner:',PX,PY,PZ
	
// C Phi random... (not given by dresner) or exact for 2 and 3 particles!
// 	IF(nbpart_ss.EQ.1) THEN
// 	U1=UCT
// 	ikeep1=NTRACK
// 	CALL RIBM(rndm,IY(18))
// 	phikeep=rndm*360.
// 	PHILAB(NTRACK)=phikeep - 180.
// c	WRITE(6,*) 'First p',ikeep1,phikeep,PHILAB(ikeep1)
// 	ENDIF
	
// 	IF(nbpart_ss.EQ.2) THEN
// 	U2=UCT
// 	ikeep2=NTRACK
// 	PHILAB(NTRACK)= AMOD(phikeep + 180.,360.) -180.
// c	WRITE(6,*) phikeep-180.,PHILAB(NTRACK)
// 	ENDIF 
	
// 	IF(nbpart_ss.EQ.3) THEN
// 	   U3=UCT
// 	   val_loc3=(-u1**2+u2**2-u3**2)/(2.*u1*u3)
// 	   val_loc2=(-u1**2-u2**2+u3**2)/(2.*u1*u2)
// 	   IF(ABS(val_loc2).LE.1.AND.ABS(val_loc2).LE.1) THEN
// 	   phi_2=ACOS(val_loc2)
// 	   phi_3=ACOS(val_loc3)
// c	   WRITE(6,*) 'OK sol 3'
// c	   WRITE(6,*) 'u1,u2,u3,phi2,phi3:',U1,U2,U3,phi_2,phi_3
// c	   WRITE(6,*) U1,U2*COS(phi_2),U3*COS(phi_3),
// c     s                U2*SIN(phi_2),U3*SIN(phi_3)
// 	   phi_2=phi_2*180./3.14159+phikeep
// 	   phi_3=360.-phi_3*180./3.14159+phikeep
// 	   PHILAB(ikeep2)=AMOD(phi_2,360.)-180.
// 	   PHILAB(NTRACK)=AMOD(phi_3,360.)-180.
// c	   WRITE(6,*)phikeep,phi_2,phi_3
// c	   WRITE(6,*) 'ikeep1,ikeep2,NTRACK:',ikeep1,ikeep2,NTRACK
// c	   WRITE(6,*) 'philab',PHILAB(ikeep1),PHILAB(ikeep2),PHILAB(NTRACK)
// c	   WRITE(6,*) 'BILAN:'
// 	   SV1=U1*COS(PHILAB(ikeep1)*3.14159/180.)
// 	   SV2=U2*COS(PHILAB(ikeep2)*3.14159/180.)
// 	   SV3=U3*COS(PHILAB(NTRACK)*3.14159/180.)
// 	   ST=SV1+SV2+SV3
// c	   WRITE(6,*) 'PX:',SV1,SV2,SV3,ST
// 	   SV1=U1*SIN(PHILAB(ikeep1)*3.14159/180.)
// 	   SV2=U2*SIN(PHILAB(ikeep2)*3.14159/180.)
// 	   SV3=U3*SIN(PHILAB(NTRACK)*3.14159/180.)
// 	   ST=SV1+SV2+SV3	    
// c	   WRITE(6,*) 'PY:',SV1,SV2,SV3,ST
// 	   SV1=PLAB(ikeep1)*COS(TETLAB(ikeep1)*3.14159/180.)
// 	   SV2=PLAB(ikeep2)*COS(TETLAB(ikeep2)*3.14159/180.)
// 	   SV3=PLAB(NTRACK)*COS(TETLAB(NTRACK)*3.14159/180.)
// 	   ST=SV1+SV2+SV3	    	   
// c	   WRITE(6,*) 'PZ:',SV1,SV2,SV3,ST 
// 	   ELSE
// 	   CALL RIBM(rndm,IY(17))
// 	   phikeep=rndm*360.
// 	   PHILAB(ikeep2)=phikeep - 180.
// 	   CALL RIBM(rndm,IY(19))
// 	   phikeep=rndm*360.
// 	   PHILAB(NTRACK)=phikeep - 180.	   
// 	   ENDIF	   
// 	ENDIF
	
// 	IF(nbpart_ss.GT.3) THEN
// 	CALL RIBM(rndm,IY(18))
// 	phikeep=rndm*360.
// 	PHILAB(NTRACK)=phikeep - 180.
// 	ENDIF
			
// 	ENDIF     
//       	ENDDO
// 	nb_dresnerv2=nbpart_ss
// C Lorentz boost:
//       PXS=0.
//       PYS=0.
//       PZS=0.
//       PXSB=0.
//       PYSB=0.
//       PZSB=0.
      
// C Velocities of transfo
//       E=T_PROJSPEC+M_PROJSPEC+EX_PROJSPEC
//       B1=P1_PROJSPEC/E
//       B2=P2_PROJSPEC/E
//       B3=P3_PROJSPEC/E
      
// 	DO iloc=NTRACKdeb+1,NTRACK
// 	RMASS=(PLAB(iloc)**2-ENERJ(iloc)**2)/(2.*ENERJ(iloc))
// 	E=ENERJ(iloc)+RMASS
//         UCT=PLAB(iloc)*SIN(TETLAB(iloc)*3.14159/180.)
// 	PX=UCT
//      s     *COS(PHILAB(iloc)*3.14159/180.)	
// 	PY=UCT
//      s     *SIN(PHILAB(iloc)*3.14159/180.)	
// 	PZ=PLAB(iloc)*COS(TETLAB(iloc)*3.14159/180.)
// 	PXS=PXS+PX
// 	PYS=PYS+PY
// 	PZS=PZS+PZ
// 	CALL LOREN(PX,PY,PZ,B1,B2,B3,E)
// 	PLAB(iloc)=SQRT(PX**2+PY**2+PZ**2)
// 	IF(PLAB(iloc).LE.0.0001) THEN
// 	TETLAB(iloc)=0.
// 	PHILAB(iloc)=0.
// 	ELSE
// 	TETLAB(iloc)=ACOS(PZ/PLAB(iloc))*180./3.14159
// 	PHILAB(iloc)=ATAN2(PY,PX)*180./3.14159
// 	ENDIF
// 	ENERJ(iloc)=E-RMASS
// 	PXSB=PXSB+PX
// 	PYSB=PYSB+PY
// 	PZSB=PZSB+PZ
	
// 	ENDDO

//         IF(IFLAGRAP.EQ.1) THEN	
// 	WRITE(6,*) 'Bilan avant boost:',PXS,PYS,PZS
// 	WRITE(6,*) 'Bilan apres boost:',PXSB,PYSB,PZSB
// 	WRITE(6,*) 'Impulsion spectat:',
//      s                P1_PROJSPEC,P2_PROJSPEC,P3_PROJSPEC
//         WRITE(6,*) 'Bilan Dresner (sans Boost):',PXD,PYD,PZD
//         ENDIF
// C endif no fractioning of the spect-projo     
//       ENDIF
// C endif break-up of the spect-projo      
//       ENDIF
	
// C *************************  FIN DRESNER ************************ 
//       ENDIF

#ifdef DEBUG
      G4cout <<"-------------------------------------------------" << G4endl;
      G4cout <<" mom = (" << momXsum << ", " << momYsum << ", " << momZsum << ")" << G4endl;
      G4cout <<"Total momentum: " << varntp->getMomentumSum() << G4endl;
#endif
#ifdef USEROOT 
      if(usingRoot) {
	//	std::cout <<"Filling..." << std::endl;
	h101->Fill();
	//	std::cout <<"Filling complete." << std::endl;
      }
#endif // USEROOT
      varntp->ntrack = 0;
    }
    else {
      if(varntp->ntrack == -2) {
	n = n - 1;
      }
      else {
	transparent++;
      }
    }
  }
  stoptime = clock();
  
  if(!usingRoot) {
    out.close();
  }

  sprintf(summaryFilename, "%s.runSummary", (char*)argv[7]);

  ofstream summaryFile;
  summaryFile.open(summaryFilename);
  
  summaryFile << "INCL4/ABLA C++ thin-target calculation:" << endl;
  summaryFile << endl;
  summaryFile << "Run setup:" << endl;
  summaryFile << "Bullet: " << endl;
  summaryFile << "\t Type: " << calincl->f[6] << endl;
  summaryFile << "\t Energy: " << calincl->f[2] << " MeV " << endl;
  summaryFile << "Target: " << endl;
  summaryFile << "\t A: " << calincl->f[0] << endl;
  summaryFile << "\t Z: " << calincl->f[1] << endl;
  summaryFile << "Events: " << argv[6] << endl;
  summaryFile << "CPU time: " << (stoptime - starttime)/1000 << " milliseconds";
  summaryFile << endl;
  if(!usingRoot) {
    summaryFile << "Calculation output in ASCII file: " << filename << endl;
  }
#ifdef USEROOT
  else {
    summaryFile << "Calculation output in ROOT file: " << rootfilename << endl;
  }
#endif // USEROOT
  summaryFile << endl;

//         f_cross_sect=(icoup-ntrans)
// 	f_cross_sect=f_cross_sect/(icoup+ntrans_coul)
  fCrossSection = (((double)totalevents) - ((double)transparent))/((double) totalevents + 0); //(double) debugval_.ntranscoul);


//  geomCrossSection = fCrossSection * 31.4159 * ws_.BMAX * ws_.BMAX;
  geomCrossSection = 31.4159*(ws->bmax)*(ws->bmax);

  summaryFile << "Output information: " << endl;
  summaryFile << "Total number of events: " << endl;
  summaryFile << "\t Asked by the user: " << (totalevents + 0 ) << endl; //debugval_.ntranscoul) << endl;
  summaryFile << "\t Transparent: " << transparent << endl;
  summaryFile << "\t Proper cascades: " << (totalevents - transparent) << endl;
  summaryFile << "Maximum impact parameter: " << ws->bmax << endl;
  summaryFile << "Geometrical cross-section: " << geomCrossSection << " mb" << endl;
  summaryFile << "Total reaction cross section: " << fCrossSection*31.4159*pow(ws->bmax,2) << " mb" << endl;
  summaryFile << "---------------------" << G4endl;
  summaryFile << "Normalization factor = " << geomCrossSection/(double(totalevents) - double(transparent)) << endl;
  summaryFile.close();

#ifdef USEROOT
  if(dataFile != NULL) {
    std::cout <<";; Writing data..." << std::endl;
    h101->Write();
#ifdef G4INCLDEBUG
    theLogger->saveHistograms();
#endif
    std::cout <<";; Closing ROOT file..." << std::endl;
    dataFile->Close();
  }
#endif // USEROOT
  
  cout << ";; Simulation run done." << endl;
  return 0;
}

