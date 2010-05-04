
#include "globals.hh"
#include "G4ios.hh"

#include "G4Material.hh"
#include "G4ElementVector.hh"
#include "MaterialTest47.hh"
#include "Test30Physics.hh"
#include "HistoBNLTest47.hh"
#include "HistoITEPTest47.hh"
#include "G4VContinuousDiscreteProcess.hh"
#include "G4ProcessManager.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleChange.hh"
#include "G4HadronCrossSections.hh"
#include "G4VCrossSectionDataSet.hh"
#include "G4ProtonInelasticCrossSection.hh"
#include "G4NeutronInelasticCrossSection.hh"
#include "G4HadronInelasticDataSet.hh"
#include "G4PiNuclearCrossSection.hh"
#include "G4HadronElasticDataSet.hh"

#include "G4ParticleTable.hh"
#include "G4ParticleChange.hh"
#include "G4DynamicParticle.hh"
#include "G4AntiProton.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "G4PionZero.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4Alpha.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4GenericIon.hh"
#include "G4ForceCondition.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4Step.hh"
#include "G4GRSVolume.hh"

#include "G4BGGNucleonInelasticXS.hh"
#include "G4BGGPionInelasticXS.hh"

#include "G4UnitsTable.hh"
#include "G4ExcitationHandler.hh"
#include "G4PreCompoundModel.hh"
#include "G4Evaporation.hh"

#include "G4StateManager.hh"
#include "G4DecayPhysics.hh"
#include "G4Timer.hh"
#include "Randomize.hh"

#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>

using namespace std;

int main(int argc, char** argv) {

  G4cout <<"========================================================" <<G4endl;
  G4cout <<"======      Intermediate Energy Test Start      ========" <<G4endl;
  G4cout <<"========================================================" <<G4endl;

  // Set default parameter values
  G4String  namePart = "proton";
  G4String  nameMat  = "Be";
  G4String  nameGen  = "Binary";
  G4double  energy   = 0;
  G4double  m_p      = 1.40*GeV;
  G4int     nevt     = 1000;
  G4double theStep   = 0.01*micrometer;
  G4Material* material = 0;
  G4bool isBNL = false;
  G4bool gtran = false;
  G4bool gemis = false;
  G4int  verbose  = 0; 
  G4bool xsbgg    = true;
  G4bool nevap    = true;
  G4long myseed   = 135799753;

  // Choose the Random engine
  //
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  CLHEP::HepRandom::setTheSeed(myseed);

  G4DecayPhysics decays;
  decays.ConstructParticle();  

  const G4ParticleDefinition* proton = G4Proton::Proton();
  const G4ParticleDefinition* neutron = G4Neutron::Neutron();
  const G4ParticleDefinition* pin = G4PionMinus::PionMinus();
  const G4ParticleDefinition* pip = G4PionPlus::PionPlus();

  G4ParticleTable* partTable = G4ParticleTable::GetParticleTable();
  partTable->SetReadiness();

  // Track
  G4ThreeVector aPosition  = G4ThreeVector(0.,0.,0.);
  G4double      aTime      = 0. ;
  G4ThreeVector aDirection = G4ThreeVector(0.0,0.0,1.0);
  G4double nx = 0.0, ny = 0.0, nz = 0.0;
  G4int jobid = -1;
  G4int clusterid = -1;

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

  MaterialTest47*  mate = new MaterialTest47();
  Test30Physics*   phys = new Test30Physics();
  HistoITEPTest47  histoITEP(namePart, nameMat, m_p/GeV, nameGen);
  HistoBNLTest47   histoBNL(namePart, nameMat, m_p/GeV, nameGen);

  // Geometry

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
  G4cout << "#time(ns)" << G4endl;
  G4cout << "#run" << G4endl;
  G4cout << "#exit" << G4endl;
  G4cout << "#HETCEmission" << G4endl;
  G4cout << "#GNASHTransition" << G4endl;
  G4cout << "#GEMEvaporation" << G4endl;
  G4cout << "#randomSeed" << G4endl;
  G4cout << "#jobID" << G4endl;
  G4cout << "#clusterID" << G4endl;

  G4String line, line1;
  G4bool end = true;
  for(G4int run=0; run<100; run++) {
    do {
      (*fin) >> line;
      G4cout << "Next line " << line << G4endl;
      if(line == "#particle") {
        (*fin) >> namePart;
	histoITEP.setParticle(namePart);
	histoBNL.setParticle(namePart);
      } else if(line == "#momentum(MeV/c)") {
        (*fin) >> m_p;
        m_p *= MeV;
	histoITEP.setMomentum(m_p/GeV);
	histoBNL.setMomentum(m_p/GeV);
      } else if(line == "#events") {
        (*fin) >> nevt;
      } else if(line == "#isBNL") {
        isBNL = true;
      } else if(line == "#isITEP") {
        isBNL = false;
      } else if(line == "#step(mm)") {
        (*fin) >> theStep;
        theStep *= mm;
      } else if(line == "#material") {
        (*fin) >> nameMat;
	histoITEP.setTarget(nameMat);
	histoBNL.setTarget(nameMat);
      } else if(line == "#generator") {
        (*fin) >> nameGen;
	histoITEP.setGenerator(nameGen);
	histoBNL.setGenerator(nameGen);
      } else if(line == "#xs_ghad") {
	xsbgg = false;
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
      } else if(line == "#HETCEmission") {
        gemis = true;
      } else if(line == "#GNASHTransition") {
        gtran = true;
      } else if(line == "#GEMEvaporation") {
        nevap = true;
      } else if(line == "#randomSeed") {
        (*fin) >> myseed;
	CLHEP::HepRandom::setTheSeed(myseed);
	G4cout << "###### Set Random Seed to " << myseed << "     #####" << G4endl;
      } else if ( line == "#jobID") {
	(*fin) >> jobid ;
	histoITEP.setJobID(jobid);
	histoBNL.setJobID(jobid);
      } else if ( line == "#clusterID") {
	(*fin) >> clusterid ;
	histoITEP.setClusterID(clusterid);
	histoBNL.setClusterID(clusterid);
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
    G4double eTot  = energy+amass;

    G4cout << "energy = " << energy/GeV << " GeV" << G4endl;
    G4cout << "Target Mass = " << amass/GeV << " GeV and Initial total energy = " << eTot/GeV << " GeV" << G4endl;

    // Create a DynamicParticle
    G4DynamicParticle dParticle(part,aDirection,energy);
    G4VCrossSectionDataSet* cs = 0;
    G4double cross_sec = 0.0;

    if(nameGen == "LElastic" || nameGen == "elastic") {
      cs = new G4HadronElasticDataSet();
    } else if(part == proton && Z > 1 && nameGen != "lepar") {
      if(xsbgg) cs = new G4BGGNucleonInelasticXS(part);
      else      cs = new G4ProtonInelasticCrossSection();
    } else if(part == neutron && Z > 1 && nameGen != "lepar") {
      if(xsbgg) cs = new G4BGGNucleonInelasticXS(part);
      else      cs = new G4NeutronInelasticCrossSection();
    } else if((part == pin || part == pip) && Z > 1 && nameGen != "lepar") {
      if(xsbgg) cs = new G4BGGPionInelasticXS(part);
      else cs = new G4PiNuclearCrossSection();
    } else { 
      cs = new G4HadronInelasticDataSet();
    }

    if(cs) {
      cs->BuildPhysicsTable(*part);
      cross_sec = cs->GetCrossSection(&dParticle, elm);
    } else {
      cross_sec = (G4HadronCrossSections::Instance())->
        GetInelasticCrossSection(&dParticle, elm);
    }
    cross_sec /= millibarn;
    G4cout << "    cross(b)= " << cross_sec << " mb" << G4endl;

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

    if(!G4StateManager::GetStateManager()->SetNewState(G4State_Idle))
      G4cout << "G4StateManager PROBLEM! " << G4endl;

    G4Timer* timer = new G4Timer();
    timer->Start();
    G4ThreeVector   mom;
    G4LorentzVector labv;
    G4VParticleChange* aChange = 0;
    G4int nsec2 = 0;

    for (G4int iter=0; iter<nevt; iter++) {

      if(verbose>=1) 
        G4cout << "### " << iter << "-th event start " << G4endl;

      G4double e0 = energy-mass;

      dParticle.SetKineticEnergy(e0);

      gTrack->SetStep(step);
      gTrack->SetKineticEnergy(e0);
      labv = G4LorentzVector(0., 0., sqrt(e0*(e0 + 2.0*mass))/GeV,
                             (e0 + mass + amass)/GeV);
      aChange = proc->PostStepDoIt(*gTrack,*step);

      G4int n = aChange->GetNumberOfSecondaries();
      if (verbose>1) {
	if (iter < 10 || (n <=2 && nsec2 < 10)) {
	  G4cout << "Event " << iter << " Process " << proc->GetProcessName() <<  " Type " << proc->GetProcessType() << "/" << proc->GetProcessSubType() << " Secondaries " << n << G4endl;
	}
      }
      if (n <= 2) nsec2++;

      if(verbose>=1 and iter == 1000*(iter/1000)) 
        G4cout << "##### " << iter << "-th event  #####" << G4endl;

      if (isBNL) histoBNL.fill(aChange, labv);
      else       histoITEP.fill(aChange, labv);

      for (G4int i=0; i<n; i++) {
        delete aChange->GetSecondary(i);
      }
      aChange->Clear();
    }
 
    timer->Stop();
    G4cout << "  "  << *timer << G4endl;
    delete timer;

    // Committing the transaction with the tree
    if (isBNL) histoBNL.write(cross_sec, nevt);
    else       histoITEP.write(cross_sec, nevt);
  }

  delete mate;
  delete fin;
  delete phys;

  G4cout << "###### End of test #####" << G4endl;
}

