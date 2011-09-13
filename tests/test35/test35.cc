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
// $Id: test35.cc,v 1.30 2010-12-23 15:16:12 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//
//      File name:     Harp
//
//      Author:        V.Ivanchenko
//
//      Creation date: 12 March 2002
//
//      Modifications:
//      14.11.03 Renamed to cascade
//      24.11.05 Use binning corresponding to HARP 
//
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4ios.hh"

#include "G4Material.hh"
#include "G4ElementVector.hh"
#include "Test30Material.hh"
#include "Test30Physics.hh"
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

#include "Histo.hh"
#include "G4Timer.hh"
#include "G4NistManager.hh"

#include "G4QInelastic.hh"
#include "G4ForceCondition.hh"
#include "G4TouchableHistory.hh"

#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>

using namespace std;

int main(int argc, char** argv)
{
  cout << "========================================================" << endl;
  cout << "======             HARP Test Start           ========" << endl;
  cout << "========================================================" << endl;

  Histo  histo;
 
  G4String  namePart = "proton";
  G4String  nameMat  = "G4_Al";
  G4String  nameGen  = "Binary";
  G4bool    logx     = false;
  G4bool    usepaw   = false;
  G4bool    isInitH  = false;
  G4bool    inclusive= true;
  G4int     verbose  = 0;
  G4double  energy   = 100.*MeV;
  G4double  sigmae   = 0.0;
  G4double  m_pmax   = 0.0;
  G4double  m_p      = 0.0;
  G4double  elim     = 30.*MeV;
  G4double  dangl    = 5.0;
  G4int     nevt     = 1000;
  G4int     modu     = 10000;
  G4int     nbins    = 100;
  G4int     nbinsa   = 40;
  G4int     nbinse   = 80;
  G4int     nbinsd   = 20;
  G4int     nbinspi  = 20;
  G4int     nangl    = 0;
  G4int     nanglpi  = 0;
  G4int     nmompi   = 0;
  G4String hFile     = "test35";
  G4double theStep   = 0.01*micrometer;
  G4double range     = 1.0*micrometer;
  G4double  emax     = 160.*MeV;
  G4double  emaxpi   = 200.*MeV;
  G4double ebinlog   = 2.0*MeV;
  G4double eBound    = 70.*MeV;
  G4double kBound    = 0.2;
  G4Material* material = 0;
  G4bool nevap       = false;
  G4bool gtran       = false;
  G4bool gemis       = false;
  G4bool xssolang    = true;
  G4bool xsbgg       = true;
  G4bool saverand    = false;

  G4DecayPhysics decays;
  decays.ConstructParticle();  

  const G4ParticleDefinition* proton = G4Proton::Proton();
  const G4ParticleDefinition* neutron = G4Neutron::Neutron();
  const G4ParticleDefinition* pin = G4PionMinus::PionMinus();
  const G4ParticleDefinition* pip = G4PionPlus::PionPlus();
  const G4ParticleDefinition* pi0 = G4PionZero::PionZero();

  G4ParticleTable* partTable = G4ParticleTable::GetParticleTable();
  partTable->SetReadiness();

  G4double ang[20] = {0.0};
  G4double angpi[20] = {0.0};
  G4double mompi[20] = {0.0};

  // Track
  G4ThreeVector aPosition = G4ThreeVector(0.,0.,0.);
  G4double      aTime     = 0. ;
  G4ThreeVector aDirection      = G4ThreeVector(0.0,0.0,1.0);
  G4double nx = 0.0, ny = 0.0, nz = 0.0;

  cout.setf( std::ios::scientific, std::ios::floatfield );

  // Control on input

  if(argc < 2) {
    cout << "Input file is not specified! Exit" << endl;
    exit(1);
  }

  std::ifstream* fin = new std::ifstream();
  G4String fname = argv[1];
  fin->open(fname.c_str());
  if( !fin->is_open()) {
    cout << "Input file <" << fname << "> does not exist! Exit" << endl;
    exit(1);
  }

  //--------- Materials definition ---------

  Test30Material*  mate = new Test30Material();
  Test30Physics*   phys = new Test30Physics();
  G4NistManager::Instance()->SetVerbose(verbose);

  // Geometry

  G4double dimX = 100.0*cm;
  G4double dimY = 100.0*cm;
  G4double dimZ = 100.0*cm;

  G4Box* sFrame = new G4Box ("Box",dimX, dimY, dimZ);
  G4LogicalVolume* lFrame = new G4LogicalVolume(sFrame,material,"Box",0,0,0);
  G4PVPlacement* pFrame = new G4PVPlacement(0,G4ThreeVector(),"Box",
                                            lFrame,0,false,0);

  G4String line, line1;
  G4bool end = true;
  for(G4int run=0; run<20; ++run) {
    nameGen = "";
    do {
      line = "";
      if( fin->eof() ) {
        end = false;
        break;
      }
      (*fin) >> line;
      if(verbose > 0) { cout << "Next line " << line << endl; }
      if(line == "#particle") {
        (*fin) >> namePart;
      } else if(line == "#sigmae(MeV)") {
        (*fin) >> sigmae;
        sigmae *= MeV;
      } else if(line == "#emax(MeV)") {
        (*fin) >> emax;
        emax *= MeV;
      } else if(line == "#emaxpi(MeV)") {
        (*fin) >> emaxpi;
        emaxpi *= MeV;
      } else if(line == "#elim(MeV)") {
        (*fin) >> elim;
        elim *= MeV;
      } else if(line == "#ebinlog(MeV)") {
        (*fin) >> ebinlog;
	if (ebinlog < 1.1) { ebinlog = 1.1; }
        ebinlog *= MeV;
      } else if(line == "#pmax(MeV/c)") {
        (*fin) >> m_pmax;
        m_pmax *= MeV;
      } else if(line == "#momentum(MeV/c)") {
        (*fin) >> m_p;
        m_p *= MeV;
      } else if(line == "#events") {
	nevt = 1;
        (*fin) >> nevt;
	G4cout<<"Number of events according to macro "<< nevt << G4endl;
	char* st = getenv("statistic");
        if(st) {
	  G4int st0 = atol(st);
	  if (st0 > 0) {  
	    nevt = st0;
	    G4cout<<"Number of events is forced to "<< nevt << G4endl;
	  }
	}
      } else if(line == "#exclusive") {
        inclusive = false;
      } else if(line == "#inclusive") {
        inclusive = true;
      } else if(line == "#nbins") {
        (*fin) >> nbins;
      } else if(line == "#nbinse") {
        (*fin) >> nbinse;
      } else if(line == "#nbinsa") {
        (*fin) >> nbinsa;
      } else if(line == "#nbinsd") {
        (*fin) >> nbinsd;
      } else if(line == "#nbinspi") {
        (*fin) >> nbinspi;
      } else if(line == "#paw") {
        G4String sss;
        (*fin) >> sss;
      } else if(line == "#nangle") {
        (*fin) >> nangl;
      } else if(line == "#nanglepi") {
        (*fin) >> nanglpi;
      } else if(line == "#anglespi") {
        for(int k=0; k<nanglpi; k++) {(*fin) >> angpi[k];}
      } else if(line == "#nmompi") {
        (*fin) >> nmompi;
      } else if(line == "#mompi(MeV/c)") {
        for(int k=0; k<nmompi; k++) {(*fin) >> mompi[k];}
      } else if(line == "#dangle") {
        (*fin) >> dangl;
      } else if(line == "#angles") {
        for(int k=0; k<nangl; k++) {(*fin) >> ang[k];}
      } else if(line == "#range(mm)") {
        (*fin) >> range;
        range *= mm;
      } else if(line == "#step(mm)") {
        (*fin) >> theStep;
        theStep *= mm;
      } else if(line == "#eBound(MeV)") {
        (*fin) >> eBound;
        eBound *= MeV;
      } else if(line == "#kBound") {
        (*fin) >> kBound;
      } else if(line == "#material") {
        (*fin) >> nameMat;
      } else if(line == "#particle") {
        (*fin) >> namePart;
      } else if(line == "#generator") {
	nameGen = "";
        (*fin) >> nameGen;
	if (nameGen == "") {
	  G4cout << "Generator name is empty! " << G4endl; 
	  continue;
	}
        hFile = nameGen;
        if(nameGen == "binary")       { hFile = "bic"; }
        else if(nameGen == "bertini") { hFile = "bert"; }
        else if(nameGen == "lepar")   { nameGen = "lhep";  hFile = "lhep"; }
        else if(nameGen == "CHIPS")   { nameGen = "chips"; hFile = "chips"; }
	else if(nameGen == "elastic") { break; }

	char* c = getenv(nameGen);
        if(!c) {
	  G4cout << "Generator <" << nameGen << "> is not included in the "
		 << " list SetModels.csh, so is ignored!" 
		 << G4endl; 
	  continue;
	}
	G4String s(c);
	if(s=="1") { break; }	
      } else if(line == "#rad") {
	xssolang = false;
      } else if(line == "#xs_ghad") {
	xsbgg = false;
      } else if(line == "#verbose") {
        (*fin) >> verbose;
        G4cout << "### New verbose level " << verbose << G4endl;
	G4NistManager::Instance()->SetVerbose(verbose);
      } else if(line == "#saverand") {
        saverand = true;
      } else if(line == "#random") {
        G4String sss("");
        (*fin) >> sss;
	CLHEP::HepRandom::restoreEngineStatus(sss);
	if(verbose>0) G4cout << "Random Engine restored from file <"
			     << sss << ">" << G4endl;
	CLHEP::HepRandom::showEngineStatus();
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
      } else if(line == "#logx") {
        logx = true;
      } else if(line == "#exit") {
        end = false;
        break;
      } else if(line == "#HETCEmission") {
        gemis = true;
      } else if(line == "#GNASHTransition") {
        gtran = true;
      } else if(line == "#GEMEvaporation") {
        nevap = true;
      }
    } while(end);

    if(!end) { break; }

    cout << "###### Start new run # " << run << "     #####" << endl;

    material = mate->GetMaterial(nameMat);
    if(!material) {
      cout << "Material <" << nameMat
	     << "> is not found out"
	     << endl;
	     exit(1);
    }

    // -------- Projectile

    G4ParticleDefinition* part = 
      (G4ParticleTable::GetParticleTable())->FindParticle(namePart);
    G4double mass = part->GetPDGMass();
    if(m_pmax == 0.0) { m_pmax = emax; }
    else              { emax   = m_pmax; }
    energy = sqrt(m_p*m_p + mass*mass) - mass; 
    G4DynamicParticle dParticle(part,aDirection,energy);

    // ------- Select model
    G4VProcess* proc = phys->GetProcess(nameGen, namePart, material);
    G4QInelastic* chips = 0;
    if(nameGen == "chips") { chips = new G4QInelastic(); }

    if(!proc) { 
      cout << "For particle: " << part->GetParticleName()
	     << " generator " << nameGen << " is unavailable"
	     << endl;
      break;
    }

    // ------- Define target A
    const G4Element* elm = material->GetElement(0); 
    G4int A = (G4int)(elm->GetN()+0.5);
    G4int Z = (G4int)(elm->GetZ()+0.5);

    // ------- Binning 
    cout << "The particle:  " << part->GetParticleName() << endl;
    if(verbose > 0) {
      cout << "The material:  " << material->GetName() 
	   << "  Z= " << Z << "  A= " << A << endl;
      cout << "The step:      " << theStep/mm << " mm" << endl;
      cout << "The position:  " << aPosition/mm << " mm" << endl;
      cout << "The direction: " << aDirection << endl;
      cout << "The time:      " << aTime/ns << " ns" << endl;
    }

    G4double amass = phys->GetNucleusMass();
    G4double m_pmin = 0.0;
    G4double m_ptmax = 0.65;
    G4double m_pth = 0.2;
    G4int m_binp = 65;
    G4int m_bint = 20;
    G4double m_thetamax = 300.;
    G4double m_thetamin = 15.;
    G4double cosmin = cos(m_thetamax*0.001);
    G4double cosmax = cos(m_thetamin*0.001);

    cout << "Ekin = " << energy/GeV << " GeV" << endl;
    if(verbose > 0) {
      cout << "emax   = " << emax/GeV << " GeV" << endl;
      cout << "pmax   = " << m_pmax/GeV << " GeV" << endl;
    }
    if(usepaw && !isInitH) {

      isInitH = true;

      histo.add1D("10","Number of protons",10,-0.5,9.5);
      histo.add1D("11","Number of pions",10,-0.5,9.5);
      histo.add1D("12","Proton momentum",nbins,m_pmin/GeV,m_pmax/GeV);
      histo.add1D("13","Pion momentum",nbins,m_pmin/GeV,m_pmax/GeV);
      histo.add1D("14","Proton Pt",m_binp,0.0,m_ptmax/GeV);
      histo.add1D("15","Pion Pt",m_binp,0.0,m_ptmax/GeV);
      histo.add1D("16","Proton theta",m_bint,0.0,m_thetamax);
      histo.add1D("17","Pion theta",m_bint,0.0,m_thetamax);
      histo.add1D("18","Proton cos(theta)",m_bint,cosmin,1.0);
      histo.add1D("19","Proton cos(theta)",m_bint,cosmin,1.0);
    }
   
    if(usepaw) {
      G4String namef = hFile + ".hbook";
      histo.setFileName(namef);
      histo.book();
      G4cout << "Histograms are booked output file <" << hFile << "> "
	     << G4endl;
    }

    G4VCrossSectionDataSet* cs = 0;
    G4double cross_sec = 0.0;

    if(chips) {
      chips->SetParameters();
    } else if(nameGen == "LElastic" || nameGen == "BertiniElastic") {
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

    // -------- Track

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
    G4double safety = 10000.*cm;
    aPoint->SetSafety(safety);
    step->SetPreStepPoint(aPoint);

    bPoint = aPoint;
    G4ThreeVector bPosition = aDirection*theStep;
    bPosition += aPosition;
    bPoint->SetPosition(bPosition);
    step->SetPostStepPoint(bPoint);
    step->SetStepLength(theStep);

    if(chips) {
      G4ForceCondition condition = NotForced;
      cross_sec = 1.0/(material->GetTotNbOfAtomsPerVolume()*
		       chips->GetMeanFreePath(*gTrack, DBL_MAX, 
					      &condition));
    } else {
      cross_sec = (G4HadronCrossSections::Instance())->
        GetInelasticCrossSection(&dParticle, Z, A);
    }

    G4double factor = 
      cross_sec*MeV*1000.0*(G4double)nbinse/(energy*barn*(G4double)nevt);
    cout << "### factor  = " << factor
	 << "    cross(b)= " << cross_sec/barn << endl;

    double coeff = cross_sec*GeV*1000.0/(barn*(G4double)nevt);
    int    nmomtet [50][20];
    int    nmomtetm[50][20];
    int    nmomtet0[50][20];
    int    nmomtetp[50][20];
    double harpcs[50][20];
    double respip[50][20];
    double respin[50][20];
    double respi0[50][20];
    double resp[50][20];
    double errpip[50][20];
    double errpin[50][20];
    double errpi0[50][20];
    double errp[50][20];

    double mom0 = 0.0;
    int j, k;
    for(k=0; k<nmompi; k++) {
      double dp = mompi[k] - mom0;
      double ang0 = 0.0;
      for(j=0; j<nanglpi; j++) {
        nmomtet[k][j]  = 0;
        nmomtetm[k][j] = 0;
        nmomtet0[k][j] = 0;
        nmomtetp[k][j] = 0;
	if(xssolang)
	  harpcs[k][j]   = coeff/(twopi*dp*(cos(ang0) - cos(angpi[j])));
	else 
	  harpcs[k][j]   = coeff/(dp*(angpi[j] - ang0));
        ang0 = angpi[j];
      }
      mom0 = mompi[k];
    }    

    if(verbose > 0) {
      cout << " nbinTheta= " << nanglpi 
	   << " nbinP= " << nmompi << " ang0= " << angpi[0] << endl; 
    }

    if(!G4StateManager::GetStateManager()->SetNewState(G4State_Idle))
      cout << "G4StateManager PROBLEM! " << endl;

    G4Timer* timer = new G4Timer();
    timer->Start();
    const G4DynamicParticle* sec = 0;
    G4ParticleDefinition* pd;
    G4ThreeVector  mom;
    G4LorentzVector labv, fm;
    G4double e, p, m, px, py, pz, pt, theta, x;
    G4VParticleChange* aChange = 0;

    for (G4int iter=0; iter<nevt; ++iter) {

      if(verbose>=1 || iter == modu*(iter/modu)) { 
        G4cout << "### " << iter << "-th event start " << G4endl;
      }
      if(saverand) { CLHEP::HepRandom::saveEngineStatus("random.txt"); }
     
      G4double e0 = energy;
      if(sigmae > 0.0) {
        do {e0 = G4RandGauss::shoot(energy,sigmae);} while (e0 < 0.0);
      }
      dParticle.SetKineticEnergy(e0);

      gTrack->SetStep(step);
      gTrack->SetKineticEnergy(e0);

      labv = G4LorentzVector(0., 0., sqrt(e0*(e0 + 2.0*mass)), 
			     e0 + mass + amass);

      if(chips) { aChange = chips->PostStepDoIt(*gTrack,*step); }
      else      { aChange = proc->PostStepDoIt(*gTrack,*step); }

      G4int n = aChange->GetNumberOfSecondaries();

      G4int nbar = 0;
      G4int n_pr = 0;
      G4int n_pi = 0;
      
      for(G4int j=0; j<n; ++j) {

        sec = aChange->GetSecondary(j)->GetDynamicParticle();
        pd  = sec->GetDefinition();
        if(pd->GetPDGMass() > 100.*MeV) { ++nbar; }
      }

      for(G4int i=0; i<n; ++i) {

        sec = aChange->GetSecondary(i)->GetDynamicParticle();
        pd  = sec->GetDefinition();
        mom = sec->GetMomentumDirection();
        e   = sec->GetKineticEnergy();
	if (e < 0.0) { e = 0.0; }

        theta = mom.theta();
        G4double cost  = cos(theta);

        m = pd->GetPDGMass();
	p = sqrt(e*(e + 2.0*m));
	mom *= p;
	fm = G4LorentzVector(mom, e + m);
	labv -= fm;
        px = mom.x();
        py = mom.y();
        pz = mom.z();
        pt = sqrt(px*px +py*py);

        G4double thetamr = theta*1000.;

	if(usepaw && cost > cosmin && cost <= cosmax && p>m_pth) {
          if(pd == proton) {
            histo.fill(2,p/GeV,1.0);
            histo.fill(4,pt/GeV,1.0);
            if(p>m_pth) {
              histo.fill(6,thetamr,1.0);
              histo.fill(8,cost,1.0);
              ++n_pr;
            }
	  } else if(pd == pip || pd == pin ) {
            histo.fill(3,p/GeV,1.0);
            histo.fill(5,pt/GeV,1.0);
            if(p>m_pth) {
              histo.fill(7,thetamr,1.0);
              histo.fill(9,cost,1.0);
              ++n_pi;
            }
          }
          histo.fill(0,n_pr,1.0);
          histo.fill(1,n_pi,1.0);

	}
	if((verbose >2) ||
	   (verbose==2 && 
	    (pd == pip || pd == pin || pd == pi0 || pd == proton))) {
	  G4cout << i << "-th secondary: " 
		 << pd->GetParticleName() << "  p= "<<p
		 << "  theta= " << theta
		 << G4endl; 
	}
	if(p < mompi[nmompi-1] && theta < angpi[nanglpi-1] 
	   && (pd == pip || pd == pin || pd == pi0 || pd == proton)) {

	  G4int kang = -1;
	  G4int kp   = -1;
	  do {++kp;}   while (p > mompi[kp]);  
	  do {++kang;} while (theta > angpi[kang]);  
	  if(verbose>=2)
	    cout << " kp= " << kp << " kang= " << kang
		   <<endl; 

	  if(pd == proton)   { nmomtetp[kp][kang] += 1; }
	  else if(pd == pip) { nmomtet [kp][kang] += 1; }
	  else if(pd == pin) { nmomtetm[kp][kang] += 1; }
	  else               { nmomtet0[kp][kang] += 1; }

	}        
        delete aChange->GetSecondary(i);
      }
      aChange->Clear();
    }
  
    timer->Stop();
    G4cout << "  "  << *timer << G4endl;
    delete timer;

    // Committing the transaction with the tree
    if(usepaw) {
      std::cout << "###### Save histograms" << std::endl;
      histo.save();
    }

    std::cout.setf(std::ios::fixed);
    G4int prec = std::cout.precision(6);
    G4cout << "#################################################################" 
	   << G4endl;
    G4cout << "###### Cross section per bin for pi+" << std::setw(6) << G4endl;
    G4cout << G4endl;    

    for(k=0; k<nmompi; ++k) {
      for(j=0; j<nanglpi; ++j) {
        respip[k][j] = G4double(nmomtet[k][j]) *harpcs[k][j];
        respin[k][j] = G4double(nmomtetm[k][j])*harpcs[k][j];
        respi0[k][j] = G4double(nmomtet0[k][j])*harpcs[k][j];
        resp  [k][j] = G4double(nmomtetp[k][j])*harpcs[k][j];

        if(nmomtetp[k][j] > 0)
	  errp[k][j] = (resp[k][j])/sqrt(G4double(nmomtetp[k][j]));
        else  errp[k][j] = harpcs[k][j];

        if(nmomtet[k][j] > 0)
	  errpip[k][j] = (respip[k][j])/sqrt(G4double(nmomtet[k][j]));
        else  errpip[k][j] = harpcs[k][j];

        if(nmomtetm[k][j] > 0)
	  errpin[k][j] = (respin[k][j])/sqrt(G4double(nmomtetm[k][j]));
        else  errpin[k][j] = harpcs[k][j];

        if(nmomtet0[k][j] > 0)
	  errpi0[k][j] = (respi0[k][j])/sqrt(G4double(nmomtet0[k][j]));
        else  errpi0[k][j] = harpcs[k][j];
      }
    }
    std::ofstream* fout = new std::ofstream();
    std::string fname = hFile+"_pi+.dpdo";
    fout->open(fname.c_str(), std::ios::out|std::ios::trunc);
    ofstream* fout1 = new ofstream();
    string fname1 = hFile+"_pi-.dpdo";
    fout1->open(fname1.c_str(), std::ios::out|std::ios::trunc);
    ofstream* fout2 = new ofstream();
    string fname2 = hFile+"_pi0.dpdo";
    fout2->open(fname2.c_str(), std::ios::out|std::ios::trunc);
    ofstream* fout3 = new ofstream();
    string fname3 = hFile+"_p.dpdo";
    fout3->open(fname3.c_str(), std::ios::out|std::ios::trunc);

    mom0 = 0.0;
    G4double mom1;
    G4double dsdm[50]  = {0.0};    
    G4double dsda[20]  = {0.0};    
    G4double dsdmm[50] = {0.0};    
    G4double dsdam[20] = {0.0};    
    G4double dsdm0[50] = {0.0};    
    G4double dsda0[20] = {0.0};    
    G4double dsdmp[50] = {0.0};    
    G4double dsdap[20] = {0.0};    
    
    G4double cross;
    G4double crossm;
    G4double cross0;
    G4double crossp;

    for(k=0; k<nmompi; ++k) {
      mom1 = mompi[k];
      dsdm[k] = 0.0;
      G4cout << "## Next momentum bin " 
	     << mom0 << "  -  " << mom1 << "  MeV/c" 
	     << G4endl;

      G4double ang0 = 0.0;
      G4double ang1 = 0.0;

      for(j=0; j<nanglpi; ++j) {
        ang1  = angpi[j];
        cross = respip[k][j];
        crossm= respin[k][j];
        cross0= respi0[k][j];
        crossp= resp[k][j];
        G4cout << "  " << cross;
        (*fout) << mom0/GeV << " " << mom1/GeV << " " << ang0 << " " << ang1 
		<< " " << cross << " " << errpip[k][j] << std::endl; 
	(*fout1)<< mom0/GeV << " " << mom1/GeV << " " << ang0 << " " << ang1 
		<< " " << crossm << " " << errpin[k][j] << std::endl; 
	(*fout2)<< mom0/GeV << " " << mom1/GeV << " " << ang0 << " " << ang1 
		<< " " << cross0 << " " << errpi0[k][j] << std::endl; 
	(*fout3)<< mom0/GeV << " " << mom1/GeV << " " << ang0 << " " << ang1 
		<< " " << crossp << " " << errp[k][j] << std::endl; 

        if(k>0) {
          x = (mom1 - mom0)/GeV;
          dsda[j]  += cross*x;
	  dsdam[j] += crossm*x;
	  dsda0[j] += cross0*x;
	  dsdap[j] += crossp*x;
        } else {
          dsda[j]  = 0.0;
	  dsdam[j] = 0.0;
	  dsda0[j] = 0.0;
	  dsdap[j] = 0.0;
	}
        if(j>0) {
          x = twopi*(cos(ang0) - cos(ang1));
          dsdm[k]  += cross*x;
	  dsdmm[k] += crossm*x;
	  dsdm0[k] += cross0*x;
	  dsdmp[k] += crossp*x;
	}

        ang0 = ang1;
      }
      cout << endl;
      mom0 = mom1;
    }
    G4cout << G4endl;    
    G4cout << G4endl;    
    G4cout << "## ds/do(mb/strad) for pi+ with momentum cut: " 
           << mompi[0] << "  -  " << mompi[nmompi-1] 
           << "  MeV/c" << G4endl;
    for(j=0; j<nanglpi; ++j) {
      cout << "  " << dsda[j];
    }
    G4cout << G4endl;  
    G4cout << "## ds/dp(mb/GeV) for pi+ with theta cut " 
           << angpi[0] << "  -  " << angpi[nanglpi-1] << " radian" << G4endl;
    cross = 0.0;
    mom0  = 0.0;
    for(k=0; k<nmompi; ++k) {
      G4cout << "  " << dsdm[k];
      cross += (dsdm[k]*(mompi[k] - mom0)/GeV);
      mom0 = mompi[k];
    }
    G4cout << G4endl;
    G4cout << "## cross_tot(bn)= " << cross/1000. << G4endl;  
    G4cout << G4endl;    

    G4cout << "#################################################################" 
	   << G4endl;
    G4cout << "## ds/do(mb/strad) for pi0 with momentum cut: " 
           << mompi[0] << "  -  " << mompi[nmompi-1] 
           << "  MeV/c" << G4endl;
    cross = 0.0;
    for(j=0; j<nanglpi; ++j) {
      G4cout << "  " << dsda0[j];
    }
    G4cout << G4endl;    
    G4cout << G4endl;    
    G4cout << "## ds/dp(mb/GeV) for pi0 with theta cut " 
           << angpi[0] << "  -  " << angpi[nanglpi-1] << " radian" << G4endl;
    cross = 0.0;
    mom0  = 0.0;
    for(k=0; k<nmompi; ++k) {
      G4cout << "  " << dsdm0[k];
      cross += ((dsdm0[k])*(mompi[k] - mom0)/GeV);
      mom0 = mompi[k];
    }
    G4cout << G4endl;
    G4cout << "## cross_tot(bn)= " << cross/1000. << G4endl;  
    G4cout << G4endl;

    G4cout << "#################################################################" 
	   << G4endl;
    G4cout << "## ds/do(mb/strad) for protons with momentum cut: " 
           << mompi[0] << "  -  " << mompi[nmompi-1] 
           << "  MeV/c" << G4endl;
    cross = 0.0;
    for(j=0; j<nanglpi; ++j) {
      G4cout << "  " << dsdap[j];
    }
    G4cout << G4endl;    
    G4cout << G4endl;    
    G4cout << "## ds/dp(mb/GeV) for pi0 with theta cut " 
           << angpi[0] << "  -  " << angpi[nanglpi-1] << " radian" << G4endl;
    cross = 0.0;
    mom0  = 0.0;
    for(k=0; k<nmompi; ++k) {
      G4cout << "  " << dsdmp[k];
      cross += ((dsdmp[k])*(mompi[k] - mom0)/GeV);
      mom0 = mompi[k];
    }
    G4cout << G4endl;
    G4cout << "## cross_tot(bn)= " << cross/1000. << G4endl;  
    G4cout << G4endl;

    G4cout << "#################################################################" 
	   << G4endl;
    G4cout << "## ds/do(mb/strad) for pi- with momentum cut: " 
           << mompi[0] << "  -  " << mompi[nmompi-1] 
           << "  MeV/c" << G4endl;
    for(j=0; j<nanglpi; ++j) {
      cout << "  " << dsdam[j];
    }
    G4cout << G4endl;    
    G4cout << G4endl;    
    G4cout << "## ds/dp(mb/GeV) for pi- with theta cut " 
           << angpi[0] << "  -  " << angpi[nanglpi-1] << " radian" << G4endl;
    cross = 0.0;
    mom0  = 0.0;
    for(k=0; k<nmompi; ++k) {
      G4cout << "  " << dsdmm[k];
      cross += (dsdmm[k]*(mompi[k] - mom0)/GeV);
      mom0 = mompi[k];
    }
    G4cout << G4endl;
    G4cout << "## cross_tot(bn)= " << cross/1000. << G4endl;  
    G4cout << G4endl;

    mom0 = 0.0;
    if(mom0 == 0.0) {
      for(k=0; k<nmompi; ++k) {
	G4double mom1 = mompi[k];
	G4cout << "## Next momentum bin " 
	       << mom0 << "  -  " << mom1 << "  MeV/c" 
	       << G4endl;

	for(j=0; j<nanglpi; ++j) {
	  G4cout << "  " << respin[k][j];
	}
	G4cout << G4endl;
	mom0 = mom1;
      }
    }
    fout->close();
    fout1->close();
    fout2->close();
    fout3->close();
    if(verbose > 0) {
      cout << "###### End of run # " << run << "     ######" << endl;
    }
    cout.precision(prec);
  }  

  delete pFrame;
  delete lFrame;
  delete sFrame;
  delete mate;
  delete fin;
  delete phys;
  partTable->DeleteAllParticles();

  G4cout << "###### End of test #####" << G4endl;
}

