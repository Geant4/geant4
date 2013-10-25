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

#include "globals.hh"
#include "G4ios.hh"

#include "G4Material.hh"
#include "G4ElementVector.hh"
#include "G4NistManager.hh"

#include "Test30Physics.hh"
#include "HistoBNLTest47.hh"
#include "HistoITEPTest47.hh"
#include "HistoMIPSTest47.hh"
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
#include "G4ForceCondition.hh"
#include "G4TouchableHistory.hh"

#include "G4UnitsTable.hh"
#include "G4ExcitationHandler.hh"
#include "G4PreCompoundModel.hh"
#include "G4Evaporation.hh"

#include "G4StateManager.hh"
#include "G4DecayPhysics.hh"
#include "G4Timer.hh"
#include "Randomize.hh"
#include "CLHEP/Units/PhysicalConstants.h"

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
  ostringstream osMat(ios_base::out|ios_base::app);// to enable appending in output operations
  G4String  nameGen  = "Binary";
  G4double  energy   = 0;
  G4double  m_p      = 1.40*CLHEP::GeV;
  G4int     nevt     = 1000;
  G4double theStep   = 0.01*CLHEP::micrometer;
  G4Material* material = 0;
  G4bool isBNL       = false;
  G4bool isMIPS      = false;
  //G4bool gtran     = false;
  //G4bool gemis     = false;
  G4int  verbose     = 0; 
  G4bool xsbgg       = true;
  //G4bool nevap     = true;
  G4long myseed      = 135799753;
  G4double dpByp     = 0.;
  G4double anglX     = 0.;
  G4double danglX    = 0.;
  G4double anglY     = 0.;
  G4double danglY    = 0.;
  G4double targetL   = 0.;
  G4double beamX     = 0.;
  G4double beamDX    = 0.;
  G4double beamY     = 0.;
  G4double beamDY    = 0.;
  G4double targetR   = 0.;

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

  G4GenericIon* gion = G4GenericIon::GenericIon();
  gion->SetProcessManager(new G4ProcessManager(gion));  

  G4ParticleTable* partTable = G4ParticleTable::GetParticleTable();
  partTable->SetReadiness();

  G4IonTable* ionTable = partTable->GetIonTable();
  ionTable->CreateAllIon();
  ionTable->CreateAllIsomer();

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
  // NOTE (JVY): local material definition replaced by the standard G4NistManager - see later in the code
  
  Test30Physics*   phys = new Test30Physics();
  HistoITEPTest47  histoITEP(namePart, nameMat, m_p/CLHEP::GeV, nameGen);
  HistoBNLTest47   histoBNL(namePart, nameMat, m_p/CLHEP::GeV, nameGen);
  HistoMIPSTest47  histoMIPS(namePart, nameMat, m_p/CLHEP::GeV, nameGen);

  // Geometry

  G4Box* sFrame = new G4Box ("Box", 100.0*CLHEP::cm, 100.0*CLHEP::cm, 100.0*CLHEP::cm);
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
  G4cout << "#isITEP" << G4endl;
  G4cout << "#isBNL" << G4endl;
  G4cout << "#isMIPS" << G4endl;
  G4cout << "#dpByp" << G4endl;
  G4cout << "#angleX(mrad)" << G4endl;
  G4cout << "#angleY(mrad)" << G4endl;
  G4cout << "#targetSize(mm)" << G4endl;

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
	histoMIPS.setParticle(namePart);
      } else if(line == "#momentum(MeV/c)") {
        (*fin) >> m_p;
        m_p *= CLHEP::MeV;
	histoITEP.setMomentum(m_p/CLHEP::GeV);
	histoBNL.setMomentum(m_p/CLHEP::GeV);
	histoMIPS.setMomentum(m_p/CLHEP::GeV);
      } else if(line == "#events") {
        (*fin) >> nevt;
      } else if(line == "#isBNL") {
        isBNL  = true;
	isMIPS = false;
      } else if(line == "#isITEP") {
        isBNL  = false;
	isMIPS = false;
      } else if(line == "#isMIPS") {
        isBNL  = false;
	isMIPS = true;
      } else if(line == "#step(mm)") {
        (*fin) >> theStep;
        theStep *= CLHEP::mm;
      } else if(line == "#material") {
        (*fin) >> nameMat;
	histoITEP.setTarget(nameMat);
	histoBNL.setTarget(nameMat);
	histoMIPS.setTarget(nameMat);
      } else if(line == "#generator") {
        (*fin) >> nameGen;
	histoITEP.setGenerator(nameGen);
	histoBNL.setGenerator(nameGen);
	histoMIPS.setGenerator(nameGen);
      } else if(line == "#xs_ghad") {
	xsbgg = false;
      } else if(line == "#run") {
        break;
      } else if(line == "#verbose") {
        (*fin) >> verbose;
	histoITEP.setDebug(verbose>2);
	histoBNL.setDebug(verbose>2);
	histoMIPS.setDebug(verbose>2);
      } else if(line == "#position(mm)") {
        (*fin) >> nx >> ny >> nz;
        aPosition = G4ThreeVector(nx*CLHEP::mm, ny*CLHEP::mm, nz*CLHEP::mm);
      } else if(line == "#direction") {
        (*fin) >> nx >> ny >> nz;
        if(nx*nx + ny*ny + nz*nz > 0.0) {
          aDirection = G4ThreeVector(nx, ny, nz);
          aDirection = aDirection.unit();
	}
      } else if(line == "#time(ns)") {
        (*fin) >> aTime;
        aTime *= CLHEP::ns;
      } else if(line == "#exit") {
        end = false;
        break;
      } else if(line == "#HETCEmission") {
        //gemis = true;
      } else if(line == "#GNASHTransition") {
        //gtran = true;
      } else if(line == "#GEMEvaporation") {
        //nevap = true;
      } else if(line == "#randomSeed") {
        (*fin) >> myseed;
	CLHEP::HepRandom::setTheSeed(myseed);
	G4cout << "###### Set Random Seed to " << myseed << "     #####" << G4endl;
      } else if ( line == "#jobID") {
        (*fin) >> jobid ;
        histoITEP.setJobID(jobid);
        histoBNL.setJobID(jobid);
        histoMIPS.setJobID(jobid);
      } else if ( line == "#clusterID") {
        (*fin) >> clusterid ;
        histoITEP.setClusterID(clusterid);
        histoBNL.setClusterID(clusterid);
        histoMIPS.setClusterID(clusterid);
      } else if ( line == "#dpByp") {
        (*fin) >> dpByp;
      } else if ( line == "#angleX") {
        (*fin) >> anglX >> danglX;
      } else if ( line == "#angleY") {
        (*fin) >> anglY >> danglY;
      } else if ( line == "#targetSize") {
        (*fin) >> targetL >> targetR >> beamX >> beamDX >> beamY >> beamDY;
      }
    } while(end);

    if(!end) break;

    G4cout << "###### Start new run # " << run << "     #####" << G4endl;

    osMat.clear();
    osMat.str("G4_");
    osMat << nameMat;
    G4String nameMatG4 = osMat.str();   
     
    G4cout << "###### Material: " << nameMatG4 << " derived from " << nameMat << G4endl;
    
    material = G4NistManager::Instance()->FindOrBuildMaterial(nameMatG4);

    if (!material) {
      G4cout << "Material <" << nameMat << "> is not found out" << G4endl;
      exit(1);
    }

    G4ParticleDefinition* part = 
      (G4ParticleTable::GetParticleTable())->FindParticle(namePart);

    G4VProcess* proc = phys->GetProcess(nameGen, namePart, material);

    const G4Element* elm = material->GetElement(0); 
    G4int A = (G4int)(elm->GetN()+0.5);
    G4int Z = (G4int)(elm->GetZ()+0.5);

    // G4double amass = phys->GetNucleusMass();
    G4double amass = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(Z, A);
    G4double pmass = (G4Proton::Proton())->GetPDGMass();

    if (!proc) { 
      G4cout << "For particle: " << part->GetParticleName()
	     << " generator " << nameGen << " is unavailable"
	     << G4endl;
      exit(1);
    } else {
      G4cout << "Nucleus mass " << amass << " proton mass " << pmass << G4endl;
    }

    G4cout << "The particle:  " << part->GetParticleName() << G4endl;
    G4cout << "The material:  " << material->GetName() 
	   << "  Z= " << Z << "  A= " << A << G4endl;
    G4cout << "The step:      " << theStep/CLHEP::mm << " mm" << G4endl;
    G4cout << "The position:  " << aPosition/CLHEP::mm << " mm" << G4endl;
    G4cout << "The direction: " << aDirection << G4endl;
    G4cout << "The time:      " << aTime/CLHEP::ns << " ns" << G4endl;
    if (isMIPS) {
      G4cout << "Delta(p)/p  " << dpByp << " Angle (X) " << anglX << " +- "
	     << danglX << " mrad; (Y) " << anglY << " +- " << danglY 
	     << " mrad; Beam size (X) " << beamX << " +- " << beamDX << " (Y) "
	     << beamY << " +- " << beamDY << " Target (Radius) "  << targetR
	     << "  (Length) " << targetL << " mm " << G4endl;
    }

    G4double mass = part->GetPDGMass();
    energy = sqrt(m_p*m_p + mass*mass);
    // amass defined earlier
    G4double eTot  = energy+amass;

    G4cout << "energy = " << energy/CLHEP::GeV << " GeV" << G4endl;
    G4cout << "Target Mass = " << amass/CLHEP::GeV << " GeV and Initial total energy = " << eTot/CLHEP::GeV << " GeV" << G4endl;

    // Create a DynamicParticle
    G4DynamicParticle dParticle(part,aDirection,energy);
    G4VCrossSectionDataSet* cs = 0;
    G4double cross_sec = 0.0;
    
    if( nameGen == "BertiniElastic" || nameGen == "elastic" ) {
      cs = new G4HadronElasticDataSet();
    } else if(part == proton && Z > 1) {
      if(xsbgg) cs = new G4BGGNucleonInelasticXS(part);
      else      cs = new G4ProtonInelasticCrossSection();
    } else if(part == neutron && Z > 1 && nameGen ) {
      if(xsbgg) cs = new G4BGGNucleonInelasticXS(part);
      else      cs = new G4NeutronInelasticCrossSection();
    } else if((part == pin || part == pip) && Z > 1 && nameGen ) {
      if(xsbgg) cs = new G4BGGPionInelasticXS(part);
      else cs = new G4PiNuclearCrossSection();
    } else { 
      cs = new G4HadronInelasticDataSet();
    }

    G4Track* gTrack;
    gTrack = new G4Track(&dParticle,aTime,aPosition);
    G4TouchableHandle fpTouchable(new G4TouchableHistory());
    gTrack->SetTouchableHandle(fpTouchable);

    // -------- Step

    G4Step* step;
    step = new G4Step();
    step->SetTrack(gTrack);
    gTrack->SetStep(step);
    
    G4StepPoint *aPoint, *bPoint;
    aPoint = new G4StepPoint();
    aPoint->SetPosition(aPosition);
    aPoint->SetMaterial(material);
    G4double safety = 10000.*CLHEP::cm;
    aPoint->SetSafety(safety);
    step->SetPreStepPoint(aPoint);

    bPoint = aPoint;
    G4ThreeVector bPosition = aDirection*theStep;
    bPosition += aPosition;
    bPoint->SetPosition(bPosition);
    step->SetPostStepPoint(bPoint);
    step->SetStepLength(theStep);

    if(cs) {
      cs->BuildPhysicsTable(*part);
      cross_sec = cs->GetCrossSection(&dParticle, elm);
    } else {
      cross_sec = (G4HadronCrossSections::Instance())->
        GetInelasticCrossSection(&dParticle, Z, A );
    }
    cross_sec /= CLHEP::millibarn;
    G4cout << "    cross(b)= " << cross_sec << " mb" << G4endl;

    if(!G4StateManager::GetStateManager()->SetNewState(G4State_Idle))
      G4cout << "G4StateManager PROBLEM! " << G4endl;

    G4Timer* timer = new G4Timer();
    timer->Start();
    G4ThreeVector   mom;
    G4LorentzVector labv, labp;
    G4VParticleChange* aChange = 0;
    G4int nsec2 = 0;
    G4ThreeVector newPosition(aPosition);

    for (G4int iter=0; iter<nevt; iter++) {

      if(verbose>=1) 
        G4cout << "### " << iter << "-th event start " << G4endl;

      if (isMIPS) {
	G4double angleX = G4RandGauss::shoot(anglX,danglX)*CLHEP::mrad;
	G4double angleY = G4RandGauss::shoot(anglY,danglY)*CLHEP::mrad;
	G4double zv     = targetL*(G4UniformRand() - 0.5)*CLHEP::mm;
	G4double xv     = targetR*CLHEP::mm;
	G4double yv     = targetR*CLHEP::mm;
	if (targetR > 0) {
	  while (std::sqrt(xv*xv+yv*yv) > targetR) {
	    xv          = G4RandGauss::shoot(beamX,beamDX)*CLHEP::mm;
	    yv          = G4RandGauss::shoot(beamY,beamDY)*CLHEP::mm;
	  }
	}
	G4double pmom   = m_p*G4RandGauss::shoot(1.0,dpByp);
	newPosition = aPosition + G4ThreeVector(xv,yv,zv);
	aDirection  = G4ThreeVector(-std::cos(angleX)*std::sin(angleY),
				    std::sin(angleX),
				    std::cos(angleX)*std::cos(angleY));
	bPosition   = newPosition+aDirection*theStep;
	energy = sqrt(pmom*pmom + mass*mass);
	mom    = pmom*aDirection;
	aPoint->SetPosition(newPosition);
	bPoint->SetPosition(bPosition);
	step->SetPreStepPoint(aPoint);
	step->SetPostStepPoint(bPoint);
	dParticle.SetMomentumDirection(aDirection.x(),aDirection.y(),aDirection.z());
	gTrack->SetPosition(newPosition);
	if (verbose>1) G4cout << "New Position " << newPosition 
			      << ", direction " << aDirection
			      << " and p " << mom << " " << pmom << G4endl;
      }
      G4double e0 = energy-mass;

      dParticle.SetKineticEnergy(e0);

      gTrack->SetStep(step);
      gTrack->SetKineticEnergy(e0);
      labv = G4LorentzVector(mom.x()/CLHEP::GeV, mom.y()/CLHEP::GeV, 
			     mom.z()/CLHEP::GeV, (e0+mass+amass)/CLHEP::GeV);
      labp = G4LorentzVector(mom.x()/CLHEP::GeV, mom.y()/CLHEP::GeV,
			     mom.z()/CLHEP::GeV, (e0+mass+pmass)/CLHEP::GeV);
      aChange = proc->PostStepDoIt(*gTrack,*step); 

      G4int n = aChange->GetNumberOfSecondaries();
      if (verbose>1) 
      {
	if (iter < 10 || (n <=2 && nsec2 < 10)) 
	{
	  G4cout << "Event " << iter << " Process " << proc->GetProcessName() <<  " Type " << proc->GetProcessType() << "/" << proc->GetProcessSubType() << " Secondaries " << n << G4endl;
	}
      }
      if (n <= 2) nsec2++;

      if(verbose>=1 and iter == 1000*(iter/1000)) 
        G4cout << "##### " << iter << "-th event  #####" << G4endl;

      if     (isMIPS) histoMIPS.fill(aChange, labv, newPosition, labp);
      else if (isBNL) histoBNL.fill (aChange, labv, newPosition, labp);
      else            histoITEP.fill(aChange, labv, newPosition, labp);

      for (G4int i=0; i<n; i++) {
        delete aChange->GetSecondary(i);
      }
      aChange->Clear();
    }
 
    timer->Stop();
    G4cout << "  "  << *timer << G4endl;
    delete timer;

    // Committing the transaction with the tree
    if      (isMIPS) histoMIPS.write(cross_sec, nevt);
    else if (isBNL)  histoBNL.write(cross_sec, nevt);
    else             histoITEP.write(cross_sec, nevt);
  }

  delete fin;
  delete phys;

  G4cout << "###### End of test #####" << G4endl;
}

