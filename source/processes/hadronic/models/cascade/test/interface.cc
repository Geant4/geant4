#include <vector>

#include <iomanip>
#include <iostream>
#include <string>
#include <stdio.h>
#include <time.h>

#include "globals.hh"
#include "Randomize.hh"

#include "G4StateManager.hh"
#include "G4ParticleTable.hh"

#include "G4Collider.hh"
#include "G4InuclCollider.hh"
#include "G4IntraNucleiCascader.hh"
#include "G4NonEquilibriumEvaporator.hh"
#include "G4EquilibriumEvaporator.hh"
#include "G4Fissioner.hh"
#include "G4BigBanger.hh"
#include "G4ElementaryParticleCollider.hh"
#include "G4InuclParticle.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclNuclei.hh"
#include "G4CollisionOutput.hh"
#include "G4Analyser.hh"
#include "G4WatcherGun.hh"
#include "G4ios.hh"
#include "G4BertiniData.hh"
#include "G4CascadSpecialFunctions.hh"
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

#include "G4CascadeInterface.hh"
#include "G4IBertini.hh"
#include "G4ElasticCascadeInterface.hh"
//#include "G4InuclEvaporation.hh"

#include "G4InclCascadeInterface.hh"
#include "G4InclAblaCascadeInterface.hh"

#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPhysicsList.hh"
#include "G4RunManager.hh"

enum particleType { nuclei = 0, proton = 1, neutron = 2, pionPlus = 3, pionMinus = 5, pionZero = 7, photon = 10 };
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

void test(std::string, int);

G4int tCascadeInterface(int runId, int nCascades, double eKin, int A, int Z);

int main(int argc, char **argv ) {
  G4int verboseLevel = 1;
  if (verboseLevel > 3) {
    G4cout << "Geant4 cascade region benchmarks" << G4endl;
  }
  if (argc < 2) {
      printf("usage: benchmarks <test ID> <parameters>\n");
      return(1);
    }

  enum particleType { nuclei = 0, proton = 1, neutron = 2, pionPlus = 3, pionMinus = 5, pionZero = 7, foton = 10 };

  // defaults
  G4int runId      = 1;
  G4int nCollisions = 10;      // collisions to be generated
  G4int  bulletType = proton;    // bullet particle
  G4double     eKin = 160;      // movement on in z-direction
  G4int        A = 27;      // target atomic weight Al
  G4int        Z = 13;      // target atomic number
  runId       =           (argc > 1) ? atoi(argv[1]) : runId;
  nCollisions =           (argc > 2) ? atoi(argv[2]) : nCollisions;
  bulletType  =           (argc > 3) ? atoi(argv[3]) : proton;
  eKin        = G4double(((argc > 4) ? atoi(argv[4]) : eKin));
  A           =           (argc > 5) ? atoi(argv[5]) : A;
  Z           =           (argc > 6) ? atoi(argv[6]) : Z;

  if (verboseLevel > 3) {
    G4cout << " runId        " << runId << G4endl;
    G4cout << " nCollisions " << nCollisions << G4endl;
    G4cout << "  bulletType " << bulletType  << G4endl;
    G4cout << "        eKin " << eKin        << G4endl;
    G4cout << "           A " << A           << G4endl;
    G4cout << "           Z " << Z           << G4endl;
  }

  G4RunManager* runManager = new G4RunManager;

  G4VUserDetectorConstruction* d = new det;
  runManager->SetUserInitialization(d);

  G4VUserPhysicsList* physics = new pList();

  runManager->SetUserInitialization(physics);
  runManager->Initialize();

  test("Cascade interface", tCascadeInterface(runId, nCollisions, eKin, A, Z));

  return 0;       
}

void test(std::string txt, int testStatus) {
  G4int verboseLevel=1;

  if (verboseLevel > 3) {
    G4cout << txt << ": ";
    if (testStatus){ 
      G4cout << "OK";
    } else {
      G4cout << "Fail" << G4endl;
    }; 

    G4cout << G4endl;  // test timing 
  }
}

int tCascadeInterface(int runId, int nCollisions, double eKin, int A, int Z) {
  G4int verboseLevel = 1;                          
  if (verboseLevel > 3) {
    G4cout << ">>> tCascadeInterface start" << G4endl;
  }

  //G4ParticleDefinition *particle = G4PionMinus::PionMinus();  
  //G4ParticleDefinition *particle = G4Neutron::Neutron();  
  //G4ParticleDefinition *particle = G4Proton::Proton();  

  G4Nucleus targetNucleus;                                        

  targetNucleus.SetParameters(A, Z);

  if (verboseLevel > 3) {
    G4cout << "target" << G4endl;
    G4cout << " a              : " << A                              << G4endl;
    G4cout << " z              : " << Z                              << G4endl;
    G4cout << " atomic mass    : " << targetNucleus.AtomicMass(A, Z) << G4endl;
  }

  G4ThreeVector       outVector, aPosition(0., 0., 0.);

  G4Proton          * aProton   = G4Proton::Proton();
  G4DynamicParticle aParticle;
  aParticle.SetDefinition(aProton);
  aParticle.SetKineticEnergy(eKin);
  aParticle.SetMomentumDirection(0.0, 0.0, 1.0);   
  aParticle.SetMomentum(sqrt(2*eKin*938.27));
 
  const  G4HadProjectile hadProj = G4HadProjectile(aParticle);

  G4HadFinalState   * o    = new G4HadFinalState();   

  // Set the particle table singleton ready:
  //  G4ParticleTable::GetParticleTable()->SetReadiness();

  G4CascadeInterface *theCascade        = new G4CascadeInterface();
  G4IBertini *theCascadeDev = new G4IBertini();
 
  for (G4int cascadeID =1 ; cascadeID <= nCollisions; cascadeID++) { 

    if (runId==1) {
      o = theCascade->ApplyYourself(hadProj, targetNucleus);  // released standard interface
    } else {
      o = theCascadeDev->ApplyYourself(hadProj, targetNucleus); // Interface to Development version  
    }

    G4int nPart = o->GetNumberOfSecondaries();

    G4double kE=1.0;
    G4int type;    
    for (G4int i =1 ; i < nPart; i++) { 
      G4HadSecondary    * NuclSecond = o->GetSecondary(i);
      G4DynamicParticle * op = NuclSecond->GetParticle();
      type=0; // default

      if (op->GetDefinition() ==  G4Proton::Proton()  ) type = proton;
      if (op->GetDefinition() ==  G4Neutron::Neutron()  ) type = neutron;
      if (op->GetDefinition() ==  G4PionPlus::PionPlus()  ) type = pionPlus;
      if (op->GetDefinition() ==  G4PionMinus::PionMinus()  ) type = pionMinus;
      if (op->GetDefinition() ==  G4PionZero::PionZero()  ) type = pionZero;
      if (op->GetDefinition() ==  G4Gamma::Gamma()  ) type = photon;

      G4ThreeVector mom =  op->GetMomentum();
      //      G4double outEtot =  op->GetTotalEnergy();
      G4double eKin =  op->GetKineticEnergy();

      G4cout <<
	std::setw(8)  << runId        << 
	std::setw(8)  << cascadeID    << 
	std::setw(8)  << type         << 
	std::setw(8)  << 0            << // dummy model
	std::setw(13) << eKin *kE     << 
	std::setw(13) << mom[0] *kE   << 
	std::setw(13) << mom[1] *kE   << 
	std::setw(13) << mom[2] *kE   << 
	std::setw(13) << 0            << // dummy A 
	std::setw(13) << 0            << //  dummy Z
	std::setw(13) << 0            <<// dummy E exitation
	std::setw(13) << 0            << // dummy coulumbOK
	// 	  " "  << cok   
        G4endl; 
    }

    if (verboseLevel > 3) {
      G4cout << "inc " << cascadeID << G4endl;
      G4cout << "  # secondaries " << nPart << G4endl;
      outVector  =  o->GetMomentumChange();
      G4cout << "  momentum change " << outVector << G4endl;
      G4double outE = o->GetEnergyChange();
      G4cout << "  energy change " << outE << G4endl;
    
      for (G4int iSecondary =1 ; iSecondary < nPart; iSecondary++) { 
	G4HadSecondary    * NuclSecond = o->GetSecondary(iSecondary);
	G4cout << "    secondary         " << iSecondary << G4endl;  
	G4DynamicParticle * secPart = NuclSecond->GetParticle();

	G4cout<<"      nucleus name      " << secPart->GetDefinition()->GetParticleName() << G4endl;

	G4ThreeVector outVectorN =  secPart->GetMomentum();
	G4cout << "      out vector      " << outVectorN << G4endl;
	G4double outEtot =  secPart->GetTotalEnergy();
	G4cout << "      particle  tot E " << outEtot << G4endl;  
	G4double outEkin =  secPart->GetKineticEnergy();
	G4cout << "      particle  kin E " << outEkin << G4endl;  
      }
    }
  }
  delete theCascade;
  return 1;   
}
