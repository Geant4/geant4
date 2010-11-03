
#include "globals.hh"
#include "G4ios.hh"

#include "G4Material.hh"
#include "G4ElementVector.hh"

#include "MaterialsList.hh"
#include "TestStoppingPhysics.hh"
#include "TestStoppingHisto.hh"

#include "G4VContinuousDiscreteProcess.hh"

#include "G4VRestProcess.hh"

#include "G4ProcessManager.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleChange.hh"
#include "G4HadronCrossSections.hh"
#include "G4VCrossSectionDataSet.hh"

#include "G4ParticleTable.hh"
#include "G4ParticleChange.hh"
#include "G4DynamicParticle.hh"

#include "G4AntiProton.hh"
#include "G4AntiNeutron.hh"
#include "G4PionMinus.hh"
#include "G4KaonMinus.hh"
#include "G4ForceCondition.hh"

#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4Step.hh"
#include "G4GRSVolume.hh"
#include "G4GRSSolid.hh"

#include "G4UnitsTable.hh"

#include "G4StateManager.hh"
#include "G4Navigator.hh"
#include "G4Timer.hh"

#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>

using namespace std;

int main(int argc, char** argv) {

  G4cout <<"========================================================" <<G4endl;
  G4cout <<"======      Stopping Particles Test Start       ========" <<G4endl;
  G4cout <<"========================================================" <<G4endl;

  // G4String  namePart = "anti_proton";
  G4String  namePart = "pi-";
  G4String  nameMat  = "Cu";
  G4String  nameGen  = "stopping";
  G4double  energy   = 0.;
  G4double  m_p      = 125.*MeV;
  G4int     nevt     = 1000;
  G4double theStep   = 0.01*micrometer;
  G4Material* material = 0;
  G4int  verbose  = 0; 

  // Track
  G4ThreeVector aPosition  = G4ThreeVector(0.,0.,0.);
  G4double      aTime      = 0. ;
  G4ThreeVector aDirection = G4ThreeVector(0.0,0.0,1.0);
  G4double nx = 0.0, ny = 0.0, nz = 0.0;

  // Control on input

  if (argc < 2) {
    G4cout << "Input file is not specified! Exit" << G4endl;
    exit(1);
  }

  std::ifstream* fin = new std::ifstream();
  G4String fname = argv[1];
  fin->open(fname.c_str());
  if( !fin->is_open()) {
    G4cout << "Input file <" << fname << "> does not exist! Exit" << G4endl;
    exit(1);
  }

  //--------- Materials definition ---------
  //
  MaterialsList*  mate = new MaterialsList();
  
  // physics needs to be initialized before the 1st use of particle table,
  // because it constructs particles - otherwise the table is just empty
  //
  TestStoppingPhysics*   phys = new TestStoppingPhysics();

  G4ParticleTable* partTable = G4ParticleTable::GetParticleTable();
  partTable->SetReadiness();

  // Geometry
  //
  G4Box* sFrame = new G4Box ("Box", 100.0*cm, 100.0*cm, 100.0*cm);
  G4LogicalVolume* lFrame = new G4LogicalVolume(sFrame,material,"Box",0,0,0);
  G4PVPlacement* pFrame = new G4PVPlacement(0,G4ThreeVector(),"Box",
                                            lFrame,0,false,0);

  assert(pFrame);

  // ---- Read input file
  G4cout << "Available commands are: " << G4endl;
  G4cout << "#events" << G4endl;
  G4cout << "#particle" << G4endl;
  G4cout << "#energy(MeV)" << G4endl;
  G4cout << "#momentum(MeV/c)" << G4endl;
  G4cout << "#step(mm)" << G4endl;
  G4cout << "#material" << G4endl;
  G4cout << "#generator" << G4endl;
  G4cout << "#verbose" << G4endl;
  G4cout << "#position(mm)" << G4endl;
  G4cout << "#direction" << G4endl;
  G4cout << "#time(ns)" << G4endl; // why would I need this ?
  G4cout << "#run" << G4endl;
  G4cout << "#exit" << G4endl;

  G4String line, line1;
  G4bool end = true;
  for(G4int run=0; run<100; run++) {
    do {
      (*fin) >> line;
      G4cout << "Next line " << line << G4endl;
      if(line == "#particle") {
        (*fin) >> namePart;
      } else if(line == "#momentum(MeV/c)") {
        (*fin) >> m_p;
        m_p *= MeV;
      } else if(line == "#events") {
        (*fin) >> nevt;
      } else if(line == "#step(mm)") {
        (*fin) >> theStep;
        theStep *= mm;
      } else if(line == "#material") {
        (*fin) >> nameMat;
      } else if(line == "#generator") {
        (*fin) >> nameGen;
      } else if(line == "#run") {
        break;
      } else if(line == "#verbose") {
        (*fin) >> verbose;
      } else if(line == "#position(mm)") {
        (*fin) >> nx >> ny >> nz;
        aPosition = G4ThreeVector(nx*mm, ny*mm, nz*mm);
      } else if(line == "#direction") {
        (*fin) >> nx >> ny >> nz;
        if(nx*nx + ny*ny + nz*nz > 0.0) {
          aDirection = G4ThreeVector(nx, ny, nz);
          aDirection = aDirection.unit();
	}
      } else if(line == "#time(ns)") {
        (*fin) >> aTime;
        aTime *= ns;
      } else if(line == "#exit") {
        end = false;
        break;
      }
    } while(end);

    if(!end) break;

    G4cout << "###### Start new run # " << run << "     #####" << G4endl;

    material = mate->GetMaterial(nameMat);
    if (!material) {
      G4cout << "Material <" << nameMat << "> is not found out" << G4endl;
      exit(1);
    }

    G4ParticleDefinition* part = 
      (G4ParticleTable::GetParticleTable())->FindParticle(namePart);

    G4VProcess* proc = phys->GetProcess(nameGen, namePart, material);
    
    
    G4double amass = phys->GetNucleusMass();

    if (!proc) { 
      G4cout << "For particle: " << part->GetParticleName()
	     << " generator " << nameGen << " is unavailable"
	     << G4endl;
      exit(1);
    } else {
      G4cout << "Nucleus mass " << amass << G4endl;
    }

    const G4Element* elm = material->GetElement(0); 


    G4int A = (G4int)(elm->GetN()+0.5);
    G4int Z = (G4int)(elm->GetZ()+0.5);

    G4cout << "The particle:  " << part->GetParticleName() << G4endl;
    G4cout << "The material:  " << material->GetName() 
	   << "  Z= " << Z << "  A= " << A << G4endl;
    G4cout << "The step:      " << theStep/mm << " mm" << G4endl;
    G4cout << "The position:  " << aPosition/mm << " mm" << G4endl;
    G4cout << "The direction: " << aDirection << G4endl;
    G4cout << "The time:      " << aTime/ns << " ns" << G4endl;

    G4double mass = part->GetPDGMass();
    energy = sqrt(m_p*m_p + mass*mass);

    G4cout << "energy = " << energy/GeV << " GeV" << G4endl;

    // Create a DynamicParticle
    G4DynamicParticle dParticle(part,aDirection,energy);

    // G4VCrossSectionDataSet* cs = 0;
    G4double cross_sec = 0.0;
    // for stopping code, it's a nonsense because it's only there for part.code=16
    cross_sec = (G4HadronCrossSections::Instance())->
                GetCaptureCrossSection(&dParticle, elm);
    cross_sec /= millibarn;

    G4Track* gTrack;
    gTrack = new G4Track(&dParticle,aTime,aPosition);

    // Step

    G4Step* step;
    step = new G4Step();
    step->SetTrack(gTrack);

    G4StepPoint *aPoint, *bPoint;
    aPoint = new G4StepPoint();
    aPoint->SetPosition(aPosition);
    aPoint->SetMaterial(material);
    G4double safety = 10000.*cm;
    aPoint->SetSafety(safety);
    step->SetPreStepPoint(aPoint);

    bPoint = aPoint;
    G4ThreeVector bPosition = aDirection*theStep;
    bPosition += aPosition;
    bPoint->SetPosition(bPosition);
    step->SetPostStepPoint(bPoint);
    step->SetStepLength(theStep);

    G4Navigator* nav = new G4Navigator();
    nav->SetWorldVolume(pFrame);
    G4TouchableHandle touch(nav->CreateTouchableHistory());
    
    if(!G4StateManager::GetStateManager()->SetNewState(G4State_Idle))
      G4cout << "G4StateManager PROBLEM! " << G4endl;

    
    // OK, we're almost done with prep work
    // final step: book histo's
    //

    TestStoppingHisto histo(namePart, nameMat, nameGen); 
                       // ?? better object than a pointer because it's done per run;
                       // so make the life cycle be run only, rather than delete ??

    
    G4Timer* timer = new G4Timer();
    timer->Start();
    G4ThreeVector   mom;
    G4LorentzVector labv;
    G4VParticleChange* aChange = 0;
    
    for (G4int iter=0; iter<nevt; iter++) {

      if(verbose>=1) 
        G4cout << "### " << iter << "-th event start " << G4endl;

      G4double e0 = energy-mass;

      dParticle.SetKineticEnergy(e0);

      gTrack->SetStep(step);
      gTrack->SetKineticEnergy(e0);
      gTrack->SetTouchableHandle(touch);// Set Box touchable history 
                                        // - this darn stuff is needed by CHIPS
      
      labv = G4LorentzVector(0., 0., sqrt(e0*(e0 + 2.0*mass))/GeV,
                             (e0 + mass + amass)/GeV);
      
      aChange = proc->AtRestDoIt(*gTrack,*step);
      
      if(verbose>=1 and iter == 1000*(iter/1000)) 
        G4cout << "##### " << iter << "-th event  #####" << G4endl;
      
      // everything will be picked up there
      //
      histo.FillEvt(aChange);

      // loop over secondaries and cleanup
      //
      for (G4int i=0; i<aChange->GetNumberOfSecondaries(); i++) 
      {
	delete aChange->GetSecondary(i);
      }
      aChange->Clear();
    }
 
    timer->Stop();
    G4cout << "  "  << *timer << G4endl;
    delete timer;
    
    // write out histo's (per run/test)
    //    
    histo.Write(nevt); 

  }

  delete mate;
  delete fin;
  delete phys;

  G4cout << "###### End of test #####" << G4endl;
}

