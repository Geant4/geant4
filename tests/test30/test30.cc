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
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//
//      File name:     Test30
//
//      Author:        V.Ivanchenko 
// 
//      Creation date: 12 March 2002
//
//      Modifications: 
//
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4ios.hh"
#include "g4std/fstream"
#include "g4std/iomanip"

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
#include "G4ForceCondition.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4Step.hh"
#include "G4GRSVolume.hh"

#include "G4UnitsTable.hh"
#include "CLHEP/Hist/TupleManager.h"
#include "CLHEP/Hist/HBookFile.h"
#include "CLHEP/Hist/Histogram.h"
#include "CLHEP/Hist/Tuple.h"


#include "G4Timer.hh"

int main(int argc, char** argv)
{
  HepTupleManager* hbookManager = 0;

  // -------------------------------------------------------------------
  // Setup

  G4String  namePart = "proton";
  G4String  nameMat  = "Si";
  G4String  nameGen  = "stringCHIPS";
  G4bool    usepaw   = false;
  G4int     verbose  = 0;
  G4double  energy   = 100.*MeV;
  G4double  elim     = 30.*MeV;
  G4int     nevt     = 1000;
  G4int     nbins    = 100;
  G4int     nbinsa   = 40;
  G4int     nbinse   = 80;
  G4int     nbinsd   = 20;
  G4String hFile     = "";
  G4double theStep   = 0.01*micrometer;
  G4double range     = 1.0*micrometer;
  G4double  emax     = 160.*MeV;
  G4Material* material = 0; 

  //  G4double ang[13] = {0.,11.,24.,35.,45.,56.,69.,82.,95.,106.,121.,134.,145.};
  G4double bng[14] = {0.,6.,18.,30.,40.,50.,62.,75.,88.,100.,114.,127.,140.,180.};
  G4double cng[13];


  // Track 
  G4ThreeVector aPosition = G4ThreeVector(0.,0.,0.);
  G4double      aTime     = 0. ;
  G4ThreeVector aDirection      = G4ThreeVector(0.0,0.0,1.0);
  G4double nx = 0.0, ny = 0.0, nz = 0.0;
 

  G4cout.setf( ios::scientific, ios::floatfield );


  // -------------------------------------------------------------------
  // Control on input

  if(argc < 2) {
    G4cout << "Input file is not specified! Exit" << G4endl;
    exit(1);
  }

  ifstream* fin = new ifstream();
  string fname = argv[1];
  fin->open(fname.c_str());
  if( !fin->is_open()) {
    G4cout << "Input file <" << fname << "> does not exist! Exit" << G4endl;
    exit(1);
  }

  // -------------------------------------------------------------------
  //--------- Materials definition ---------
 
  Test30Material*  mate = new Test30Material();
  Test30Physics*   phys = new Test30Physics();
	
  // Geometry 

  G4double dimX = 100.0*cm;
  G4double dimY = 100.0*cm;
  G4double dimZ = 100.0*cm;  

  G4Box* sFrame = new G4Box ("Box",dimX, dimY, dimZ);
  G4LogicalVolume* lFrame = new G4LogicalVolume(sFrame,material,"Box",0,0,0);
  G4PVPlacement* pFrame = new G4PVPlacement(0,G4ThreeVector(),"Box",
                                            lFrame,0,false,0);

  assert(pFrame);

  // -------------------------------------------------------------------
  // ---- Read input file
  G4cout << "Available commands are: " << G4endl;
  G4cout << "#events" << G4endl;
  G4cout << "#nbins" << G4endl;
  G4cout << "#nbinsa" << G4endl;
  G4cout << "#nbinse" << G4endl;
  G4cout << "#particle" << G4endl;
  G4cout << "#energy(MeV)" << G4endl;
  G4cout << "#emax(MeV)" << G4endl;
  G4cout << "#range(mm)" << G4endl;
  G4cout << "#step(mm)" << G4endl;
  G4cout << "#material" << G4endl;
  G4cout << "#generator" << G4endl;
  G4cout << "#paw" << G4endl;
  G4cout << "#verbose" << G4endl;
  G4cout << "#position(mm)" << G4endl;
  G4cout << "#direction" << G4endl;
  G4cout << "#time(ns)" << G4endl;
  G4cout << "#run" << G4endl;
  G4cout << "#exit" << G4endl;



  string line, line1;
  G4bool end = true;
  for(G4int run=0; run<100; run++) {
    do {
      (*fin) >> line;
      G4cout << "Next line " << line << G4endl;
      if(line == "#particle") {
        (*fin) >> namePart;
      } else if(line == "#energy(MeV)") {
        (*fin) >> energy;
        energy *= MeV;
      } else if(line == "#emax(MeV)") {
        (*fin) >> emax;
        emax *= MeV;
      } else if(line == "#events") {
        (*fin) >> nevt;
      } else if(line == "#nbins") {
        (*fin) >> nbins;
      } else if(line == "#nbinse") {
        (*fin) >> nbinse;
      } else if(line == "#nbinsa") {
        (*fin) >> nbinsa;
      } else if(line == "#range(mm)") {
        (*fin) >> range;
        range *= mm;
      } else if(line == "#step(mm)") {
        (*fin) >> theStep;
        theStep *= mm;
      } else if(line == "#material") {
        (*fin) >> nameMat;
      } else if(line == "#particle") {
        (*fin) >> namePart;				
      } else if(line == "#generator") {
        (*fin) >> nameGen;
      } else if(line == "#paw") {
        usepaw = true;
        (*fin) >> hFile;
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
    if(!material) {
      G4cout << "Material <" << nameMat  
	     << "> is not found out" 
	     << G4endl;			
	     exit(1);
    }

    G4ParticleDefinition* part = (G4ParticleTable::GetParticleTable())->FindParticle(namePart);

    G4VProcess* proc = phys->GetProcess(nameGen, namePart, material);
    G4double amass = phys->GetNucleusMass();
				
    const G4ParticleDefinition* proton = G4Proton::Proton();
    const G4ParticleDefinition* neutron = G4Neutron::Neutron();
    const G4ParticleDefinition* pin = G4PionMinus::PionMinus();
    const G4ParticleDefinition* pip = G4PionPlus::PionPlus();
    const G4ParticleDefinition* pi0 = G4PionZero::PionZero();
    const G4ParticleDefinition* deu = G4Deuteron::DeuteronDefinition();
    const G4ParticleDefinition* tri = G4Triton::TritonDefinition();	
    const G4ParticleDefinition* alp = G4Alpha::AlphaDefinition();
		
    if(!proc) {
      G4cout << "For particle: " << part->GetParticleName() 
	     << " generator " << nameGen << " is unavailable"
	     << G4endl;			
	     exit(1);
    }

    G4int maxn = (G4int)((*(material->GetElementVector()))[0]->GetN()) + 2; 
		
    G4cout << "The particle:  " << part->GetParticleName() << G4endl;
    G4cout << "The material:  " << material->GetName() << "  Amax= " << maxn << G4endl;
    G4cout << "The step:      " << theStep/mm << " mm" << G4endl;
    G4cout << "The position:  " << aPosition/mm << " mm" << G4endl;
    G4cout << "The direction: " << aDirection << G4endl;
    G4cout << "The time:      " << aTime/ns << " ns" << G4endl;


    // -------------------------------------------------------------------
    // ---- HBOOK initialization
    const G4int nhisto = 40; 
    HepHistogram* h[nhisto];
    G4double mass = part->GetPDGMass();
    G4double pmax = sqrt(energy*(energy + 2.0*mass));
		
    if(usepaw) {
      hbookManager = new HBookFile(hFile, 58);
      assert (hbookManager != 0);
  
      // ---- Book a histogram and ntuples
      G4cout << "Hbook file name: <" 
             << ((HBookFile*) hbookManager)->filename() << ">" << G4endl;

      h[0]=hbookManager->histogram("Number of Secondaries",50,-0.5,49.5);
      h[1]=hbookManager->histogram("Type of secondary",10,-0.5,9.5);
      h[2]=hbookManager->histogram("Phi(degrees) of Secondaries",90,-180.0,180.0);
      h[3]=hbookManager->histogram("Pz (MeV) for protons",100,-pmax,pmax);
      h[4]=hbookManager->histogram("Pz (MeV) for pi-",100,-pmax,pmax);
      h[5]=hbookManager->histogram("Pz (MeV) for pi+",100,-pmax,pmax);
      h[6]=hbookManager->histogram("Pz (MeV) for neutrons",100,-pmax,pmax);
      h[7]=hbookManager->histogram("Pt (MeV) for protons",100,0.,pmax);
      h[8]=hbookManager->histogram("Pt (MeV) for pi-",100,0.,pmax);
      h[9]=hbookManager->histogram("Pt (MeV) for pi+",100,0.,pmax);
      h[10]=hbookManager->histogram("Pt (MeV) for neutrons",100,0.,pmax);
      h[11]=hbookManager->histogram("E (MeV) for protons",100,0.,energy);
      h[12]=hbookManager->histogram("E (MeV) for pi-",100,0.,energy);
      h[13]=hbookManager->histogram("E (MeV) for pi+",100,0.,energy);
      h[14]=hbookManager->histogram("E (MeV) for neutrons",100,0.,energy);
      h[15]=hbookManager->histogram("delta E (MeV)",20,-1.,1.);
      h[16]=hbookManager->histogram("delta Pz (GeV)",20,-1.,1.);
      h[17]=hbookManager->histogram("delta Pt (GeV)",20,-1.,1.);
      
      h[18]=hbookManager->histogram("E (MeV) for pi0",100,0.,energy);
      h[19]=hbookManager->histogram("Pz (MeV) for pi0",100,-pmax,pmax);
      h[20]=hbookManager->histogram("Pt (MeV) for pi0",100,0.,pmax);
      
      h[21]=hbookManager->histogram("E(MeV) protons",nbinse,0.,emax);
      h[22]=hbookManager->histogram("E(MeV) neutrons",nbinse,0.,emax);

      h[23]=hbookManager->histogram("Phi(degrees) of neutrons",90,-180.0,180.0);

      h[24]=hbookManager->histogram("cos(theta) protons",nbinsa,-1.,1.);
      h[25]=hbookManager->histogram("cos(theta) neutrons",nbinsa,-1.,1.);

      h[26]=hbookManager->histogram("Baryon charge",maxn,-0.5,(G4double)maxn + 0.5);

      h[27]=hbookManager->histogram("ds/dE at theta = 0",nbinsd,0.,emax);
      h[28]=hbookManager->histogram("ds/dE at theta = 1",nbinsd,0.,emax);
      h[29]=hbookManager->histogram("ds/dE at theta = 2",nbinsd,0.,emax);
      h[30]=hbookManager->histogram("ds/dE at theta = 3",nbinsd,0.,emax);
      h[31]=hbookManager->histogram("ds/dE at theta = 4",nbinsd,0.,emax);
      h[32]=hbookManager->histogram("ds/dE at theta = 5",nbinsd,0.,emax);
      h[33]=hbookManager->histogram("ds/dE at theta = 6",nbinsd,0.,emax);
      h[34]=hbookManager->histogram("ds/dE at theta = 7",nbinsd,0.,emax);
      h[35]=hbookManager->histogram("ds/dE at theta = 8",nbinsd,0.,emax);
      h[36]=hbookManager->histogram("ds/dE at theta = 9",nbinsd,0.,emax);
      h[37]=hbookManager->histogram("ds/dE at theta = 10",nbinsd,0.,emax);
      h[38]=hbookManager->histogram("ds/dE at theta = 11",nbinsd,0.,emax);
      h[39]=hbookManager->histogram("ds/dE at theta = 12",nbinsd,0.,emax);
      	
      G4cout << "Histograms is initialised nbins=" << nbins
             << G4endl;
    }		
    // Create a DynamicParticle  
  
    G4DynamicParticle dParticle(part,aDirection,energy);
    G4VCrossSectionDataSet* cs = 0;
    G4double cross_sec = 0.0;

    if(part == proton) {
      cs = new G4ProtonInelasticCrossSection();
    } else if(part == neutron) {
      cs = new G4NeutronInelasticCrossSection();
    }
    if(cs) {
      cs->BuildPhysicsTable(*part);
      cross_sec = cs->GetCrossSection(&dParticle, material->GetElement(0));
    } else {
      cross_sec = (G4HadronCrossSections::Instance())->
        GetInelasticCrossSection(&dParticle, material->GetElement(0));
    }

    G4double factor = cross_sec*MeV*1000.0*(G4double)nbins/(energy*barn*(G4double)nevt);
    G4double factora= cross_sec*MeV*1000.0*(G4double)nbinsa/(twopi*pi*barn*(G4double)nevt);
    G4cout << "### factor  = " << factor
           << "### factora = " << factor 
           << "    cross(b)= " << cross_sec/barn << G4endl;
    G4double dtet = pi/(G4int)nbinsa;

    for(G4int k=0; k<13; k++) {
      cng[k] = cross_sec*MeV*1000.0*(G4double)nbinsd/
         (twopi*(cos(degree*bng[k]) - cos(degree*bng[k+1]))*barn*emax*(G4double)nevt);
    }

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

    G4RotationMatrix* rot  = new G4RotationMatrix();
    G4double phi0 = aDirection.phi();
    G4double theta0 = aDirection.theta();
    rot->rotateZ(-phi0);
    rot->rotateY(-theta0);

    G4cout << "Test rotation= " << (*rot)*(aDirection) << G4endl;

    G4Timer* timer = new G4Timer();
    timer->Start();
    const G4DynamicParticle* sec;
    G4ParticleDefinition* pd;
    G4ThreeVector  mom;
    G4LorentzVector labv, fm;
    G4double e, p, m, px, py, pz, pt, theta, sint;
    G4VParticleChange* aChange = 0;
			
    for (G4int iter=0; iter<nevt; iter++) {

      gTrack->SetStep(step); 
      gTrack->SetKineticEnergy(energy); 
      //G4double x = 0.0;

      labv = G4LorentzVector(0.0, 0.0, pmax, energy + mass + amass);
      aChange = proc->PostStepDoIt(*gTrack,*step);

      G4double de = aChange->GetLocalEnergyDeposit();
      G4int n = aChange->GetNumberOfSecondaries();
      G4int nn= n;
			
      if(verbose) {
        G4cout << "### " << iter << "-th event de(MeV) = " 
               << de/MeV << " n= " << n 
               << G4endl;
      }
			 					
      for(G4int i=0; i<n+1; i++) {

	if(i<n) {
	   sec = aChange->GetSecondary(i)->GetDynamicParticle();
	   pd  = sec->GetDefinition();					
	   mom = sec->GetMomentumDirection();
	   e   = sec->GetKineticEnergy();

	} else {
	   if(aChange->GetStatusChange() == fStopAndKill) break;
           nn++;
	   pd = part;
	   G4ParticleChange* bChange = dynamic_cast<G4ParticleChange*>(aChange);
	   mom = *(bChange->GetMomentumDirectionChange());
	   e   = bChange->GetEnergyChange();	
	}
	if (e < 0.0) {
           e = 0.0;
	}
        m = pd->GetPDGMass();
	p = sqrt(e*(e + 2.0*m));
	mom *= p;
        m  = pd->GetPDGMass();
        fm = G4LorentzVector(mom, e + m);
        labv -= fm;
        mom = (*rot)*mom;
        px = mom.x();
        py = mom.y();
        pz = mom.z();
        p  = sqrt(px*px +py*py + pz*pz);
        pt = sqrt(px*px +py*py);

        theta = mom.theta();
        G4int itet = (G4int)(theta/dtet);
        sint  = sin(dtet*(0.5 + (G4double)itet));
				
	if(usepaw && e > 0.0 && pt > 0.0) {
          h[2]->accumulate(mom.phi()/degree,1.0);
          if(pd == neutron) h[23]->accumulate(mom.phi()/degree,1.0);
	}				
	de += e;
        if(verbose>0 && abs(mom.phi()/degree - 90.) < 0.01) {
          G4cout << i << "-th secondary  " 
		 << pd->GetParticleName() << "   Ekin(MeV)= "
                 << e/MeV
		 << "   p(MeV)= " << mom/MeV
		 << "   m(MeV)= " << m/MeV
		 << "   Etot(MeV)= " << (e+m)/MeV
		 << "   pt(MeV)= " << pt/MeV
                 << G4endl;
        }
				
 

	if(usepaw) {

          if(pd) h[26]->accumulate((G4double)pd->GetBaryonNumber(), 1.0);

          if(pd == proton) { 
						
            h[1]->accumulate(1.0, 1.0);						
            h[3]->accumulate(pz/MeV, 1.0); 
            h[7]->accumulate(pt/MeV, 1.0);
            h[11]->accumulate(e/MeV, 1.0);
            h[11]->accumulate(e/MeV, 1.0);
	    //    h[18]->accumulate(e/MeV, 1.0);
	    h[21]->accumulate(e/MeV, factor);
	    h[24]->accumulate(cos(theta), factora/sint);
		
          } else if(pd == pin) {
    
	    h[1]->accumulate(4.0, 1.0);
            h[4]->accumulate(pz/MeV, 1.0); 
            h[8]->accumulate(pt/MeV, 1.0);
            h[12]->accumulate(e/MeV, 1.0);
						
          } else if(pd == pip) {
    
	    h[1]->accumulate(3.0, 1.0);		
            h[5]->accumulate(pz/MeV, 1.0); 
            h[9]->accumulate(pt/MeV, 1.0);
            h[13]->accumulate(e/MeV, 1.0);

	  } else if(pd == pi0) {
    
	    h[1]->accumulate(5.0, 1.0);	
	    h[18]->accumulate(e/MeV, 1.0);		
	    h[19]->accumulate(pz/MeV, 1.0); 
	    h[20]->accumulate(pt/MeV, 1.0);

	  } else if(pd == neutron) {
    
	    h[1]->accumulate(2.0, 1.0);	  
            h[6]->accumulate(pz/MeV, 1.0); 
            h[10]->accumulate(pt/MeV, 1.0);
            h[14]->accumulate(e/MeV, 1.0);
            //h[22]->accumulate(e/MeV, 1.0);
	    h[22]->accumulate(e/MeV, factor);
	    if(e >= elim) h[25]->accumulate(cos(theta), factora/sint);
            theta /= degree;
            for(k=0; k<13; k++) {
              if(theta <= bng[k+1]) break;
	    }
            h[27+k]->accumulate(e/MeV, cng[k]); 

	  } else if(pd == deu) {
	    h[1]->accumulate(6.0, 1.0);	
	  } else if(pd == tri) {
	    h[1]->accumulate(7.0, 1.0);	
	  } else if(pd == alp) {
	    h[1]->accumulate(8.0, 1.0);							
	  } else {
	    h[1]->accumulate(9.0, 1.0);	
	  }
	}
        if(i<n) delete aChange->GetSecondary(i);
      }
								
      if(verbose > 0) {
        G4cout << "Energy/Momentum balance= " << labv << G4endl;
      }	


      px = labv.px();	
      py = labv.py();
      pz = labv.pz();							
      p  = sqrt(px*px +py*py + pz*pz);
      pt = sqrt(px*px +py*py);

      if(usepaw) {
        h[0]->accumulate((float)nn,1.0);
	h[15]->accumulate(labv.e()/MeV, 1.0);
	h[16]->accumulate(pz/GeV, 1.0);
	h[17]->accumulate(pt/GeV, 1.0);
      }	
      //      delete aChange;
      aChange->Clear();
	
    }

    timer->Stop();
    G4cout << "  "  << *timer << G4endl;
    delete timer;

    if(usepaw) {
      hbookManager->write();
      G4cout << "# hbook is writed" << G4endl;
      for(G4int i=0; i<nhisto; i++) {
        if(h[i]) delete h[i];
      }
      delete hbookManager;    
      G4cout << "# hbook is deleted" << G4endl;
    }
    G4cout << "###### End of run # " << run << "     ######" << G4endl;

  } while(end);

  delete mate;
  delete fin;
  delete phys;

  G4cout << "###### End of test #####" << G4endl;
}
