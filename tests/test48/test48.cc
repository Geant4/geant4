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
#include "G4IsotopeVector.hh"
#include "G4ElementVector.hh"
#include "G4NistManager.hh"

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
#include "G4IonTable.hh"

#include "G4ParticleChange.hh"
#include "G4DynamicParticle.hh"

#include "G4AntiProton.hh"
#include "G4AntiNeutron.hh"
#include "G4PionMinus.hh"
#include "G4KaonMinus.hh"
#include "G4GenericIon.hh"
#include "G4ForceCondition.hh"

#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4Step.hh"
#include "G4GRSVolume.hh"
#include "G4GRSSolid.hh"

#include "G4TrackingManager.hh"
#include "G4Region.hh"

#include "G4StateManager.hh"
#include "G4Navigator.hh"
#include "G4Timer.hh"

// random engine/seed settings
#include "Randomize.hh"

#include "G4UImanager.hh"

#include "G4SystemOfUnits.hh"

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
  G4String  nameMatRaw  = "Cu";
  ostringstream osMat(ios_base::out|ios_base::app);// to enable appending in output operations
  G4String  nameGen  = "Bertini";
  G4double  energy   = 0.;
  G4double  m_p      = 125.*MeV;
  G4int     nevt     = 1000;
  G4double theStep   = 0.01*micrometer; // should it not be 0. ?
  G4Material* material = 0;
  G4int  verbose  = 0; 

//
//
  G4long myseed   = 135799753;
  G4int  jobid = -1;
  // Choose the Random engine
  //
  CLHEP::RanecuEngine* erndm = new CLHEP::RanecuEngine();
  CLHEP::HepRandom::setTheEngine( erndm );
//  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  CLHEP::HepRandom::setTheSeed(myseed);

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

  /*  G4UImanager* UI = G4UImanager::GetUIpointer();
      G4String command("/tracking/verbose 2");
      UI->ApplyCommand( command );
  */  

  // we need to be in the preinit state to create particles

  if(!G4StateManager::GetStateManager()->SetNewState(G4State_PreInit))
    G4cout << "G4StateManager PROBLEM! " << G4endl;


  // physics needs to be initialized before the 1st use of particle table,
  // because it constructs particles - otherwise the table is just empty
  //
  TestStoppingPhysics*   phys = new TestStoppingPhysics(verbose); // leak? - fixed
  // we need to get verbose from the inp file...

  G4GenericIon* gion = G4GenericIon::GenericIon();
  gion->SetProcessManager(new G4ProcessManager(gion));

  G4ParticleTable* partTable = G4ParticleTable::GetParticleTable();
  partTable->SetReadiness();

  G4IonTable* ions = partTable->GetIonTable();
  ions->CreateAllIon();
  ions->CreateAllIsomer();

  // Geometry
  //
  G4Box* sFrame = new G4Box ("Box", 100.0*cm, 100.0*cm, 100.0*cm);
  G4LogicalVolume* lFrame = new G4LogicalVolume(sFrame,material,"Box",0,0,0);
  G4PVPlacement* pFrame = new G4PVPlacement(0,G4ThreeVector(),"Box",
                                            lFrame,0,false,0);
  
  assert(pFrame);

  G4Region* rFrame = new G4Region("Region"); // needed by tracking manager
  lFrame->SetRegion(rFrame);
  rFrame->AddRootLogicalVolume(lFrame);

    G4Navigator* nav = new G4Navigator(); // FIXME leak - fixed: deleted at the end of the job
    nav->SetWorldVolume(pFrame);
    G4TouchableHandle touch(nav->CreateTouchableHistory());
    
  G4TrackingManager* trackManager = new G4TrackingManager; // deleted at the end of the job

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
//
// for parallel processing
//
  G4cout << "#randomSeed" << G4endl;
  G4cout << "#jobID" << G4endl;
//
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
        (*fin) >> line1;
        istringstream is(line1);
        is >> nevt;
        G4cout << "nevt : " << nevt << G4endl;
      } else if(line == "#step(mm)") {
        (*fin) >> theStep;
        theStep *= mm;
      } else if(line == "#material") {
        (*fin) >> nameMatRaw;
      } else if(line == "#generator") {
        (*fin) >> nameGen;
        G4cout << nameGen << G4endl;
      } else if(line == "#run") {
        break;
      } else if(line == "#verbose") {
        (*fin) >> verbose;
        if (verbose > 0) {
          G4cout <<"test48"
                 <<" : verbose set to "
                 << verbose
                 << G4endl;
        }
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
      } 
//
// needed for parallel processing
//
      else if ( line == "#randomSeed" )
      {
        (*fin) >> myseed;
	CLHEP::HepRandom::setTheSeed(myseed);
	G4cout << "###### Set Random Seed to " << myseed << "     #####" << G4endl;      
      }
      else if ( line == "#jobID" )
      {
        (*fin) >> jobid ;
      }
//
      else if(line == "#exit") 
      {
        end = false;
        break;
      }
    } while(end);

    if(!end) break;

    G4cout << "###### Start new run # " << run << "     #####" << G4endl;

    osMat.clear();
    osMat.str("G4_");
    osMat << nameMatRaw;
    G4String nameMat = osMat.str();

    G4cout << "###### Material: " << nameMat << " derived from " << nameMatRaw << G4endl;

    material = G4NistManager::Instance()->FindOrBuildMaterial(nameMat);

    if (!material) {
      G4cout << "Material <" << nameMat << "> is not found" << G4endl;
      exit(1);
    }

    G4ParticleDefinition* part = 
      (G4ParticleTable::GetParticleTable())->FindParticle(namePart);

    G4VProcess* proc = phys->GetProcess(nameGen, namePart);

    if (!proc) { 
      G4cout << "For particle: " << part->GetParticleName()
	     << " generator " << nameGen << " is unavailable"
	     << G4endl;
      exit(1);
    } 

    //    proc->SetVerboseLevel(2);

    size_t nEl = material->GetNumberOfElements();

    const G4Element* elm = material->GetElement(0); 


    G4int A = (G4int)(elm->GetN()+0.5);
    G4int Z = (G4int)(elm->GetZ()+0.5);
        
    G4int NIso = G4NistManager::Instance()->GetNumberOfNistIsotopes( Z );
    G4cout << " For materials " << material->GetName() << " Number of NISTIsotopes = " << NIso << G4endl;
    

    G4cout << "The particle:  " << part->GetParticleName() << G4endl;
    G4cout << "The material:  " << material->GetName() 
	   << "  Z= " << Z 
           << "  A of the first element= " << A 
           << "  A of the last  element= " << material->GetElement(nEl-1)->GetN() << G4endl;

    G4IsotopeVector*  iv = elm->GetIsotopeVector();
    G4cout << "Number of isotopes " << iv->size() << G4endl;

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
    // does not work for muons, commented out for now
    //    cross_sec = (G4HadronCrossSections::Instance())->
    //    GetCaptureCrossSection(&dParticle, Z ); // 2nd arg used to be G4Element*, but now it's G4int Z
    cross_sec /= millibarn;


/* move out of the loop over configs

    G4Navigator* nav = new G4Navigator(); // FIXME leak
    nav->SetWorldVolume(pFrame);
    G4TouchableHandle touch(nav->CreateTouchableHistory());
*/    
    if(!G4StateManager::GetStateManager()->SetNewState(G4State_Idle))
      G4cout << "G4StateManager PROBLEM! " << G4endl;

    G4Track* gTrack;
    gTrack = new G4Track(&dParticle,aTime,aPosition); // this may look like a leak but an attempt
                                                      // to delete gTrack messes something up
						      // at a higher level in the G4 infrastructure
    gTrack->SetTouchableHandle(touch); // the touch is defined before the loop over configs,
                                       // where the navigator is intantiated
				       // NOTE: I thought this darn stuff was needed by CHIPS only
				       //       but maybe it's needed also for something else...

    // Step
    G4Step *step; 
    step = new G4Step(); // FIXME leak =- fixed: deleted at the end of the loop
    step->InitializeStep(gTrack);
    step->GetPreStepPoint()->SetMaterial(material);// we somehow still need it
    step->SetStepLength(theStep);
    
    // OK, we're almost done with prep work
    // final step: book histo's
    //

    TestStoppingHisto histo(namePart, nameMatRaw, nameGen); 
    // ?? better object than a pointer because it's done per run;
    // so make the life cycle be run only, rather than delete ??
    if ( jobid > -1 ) histo.SetJobID(jobid);

// -->    G4LorentzVector labv; // not needed in this application ???
// amass would be needed to calculate labv )see above) 
// which is needed in some parts of test47 but not in this one
//    G4double amass = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(Z, A);    
//    std::cout << " Nucleus mass = " << amass << std::endl;

    G4ThreeVector   mom;
    G4VParticleChange* aChange = 0;

    if (verbose>=2) {
      G4cout << "Possible secondary status codes: " 
             << " fAlive " << fAlive 
             << " fStopButAlive " << fStopButAlive 
             << " fStopAndKill " << fStopAndKill 
             << " fKillTrackAndSecondaries " << fKillTrackAndSecondaries
             << " fSuspend " << fSuspend
             << " fPostponeToNextEvent " << fPostponeToNextEvent
             << G4endl;
    }

// the timer should start right before the loop
//
    G4Timer* timer = new G4Timer();
    timer->Start();
    
    for (G4int iter=0; iter<nevt; ++iter) {

      if(verbose>=0 and iter == 50000*(iter/50000))  
        G4cout << "### " << iter << "-th event start " << G4endl;

      G4double e0 = energy-mass;

      dParticle.SetKineticEnergy(e0);

      gTrack->SetStep(step);
      gTrack->SetKineticEnergy(e0);
      // gTrack->SetTouchableHandle(touch);// Set Box touchable history 
                                           // - this darn stuff used to be needed by CHIPS
      
// not needed in this application ???
//
//      labv = G4LorentzVector(0., 0., sqrt(e0*(e0 + 2.0*mass))/GeV,
//                             (e0 + mass + amass)/GeV);
      
      aChange = proc->AtRestDoIt(*gTrack,*step);
      
      if(verbose>=0 and iter == 50000*(iter/50000)) 
        G4cout << "##### " << iter << "-th event  #####" << G4endl;
      else if (verbose>=1) 
        G4cout << "##### " << iter << "-th event  #####" << G4endl;

      // loop over secondaries and cleanup
      G4int nsecnd = aChange->GetNumberOfSecondaries();

      if(verbose>=0 and iter == 50000*(iter/50000)) 
        G4cout << "##### " << nsecnd << " nsecondaries #####" << G4endl;
      else if (verbose>=1)
        G4cout << "##### " << nsecnd << " nsecondaries #####" << G4endl;

      for (G4int i=0; i<nsecnd; ++i) {
        G4Track*           sTrack = aChange->GetSecondary(i);
        G4TrackStatus sts =  	sTrack->GetTrackStatus ();
        const G4DynamicParticle* sdp = sTrack->GetDynamicParticle();
        const G4ParticleDefinition* sdpd = sdp->GetDefinition();
        G4String sName = sdpd->GetParticleName();

        if (verbose>=1) {
          G4cout << "Secondary: " << i 
                 << " is " << sName
                 << ", " << sdp->GetPDGcode()
                 << ", status: " << sts
                 << " with kinetic energy of " 
                 << sdp->GetKineticEnergy()/GeV 
                 << " and momentum of " 
                 << sdp->GetTotalMomentum()/GeV << " GeV"
                 << G4endl;
        }

        // HACK we need to "stop" muonic atoms
        // it looks like the mu- atom secondaries are always "born alive" i.e. not stopped

        //fTrackStatus = fAlive, fBelowThreshold = false, fGoodForTracking = false, 

        if (sdp->GetPDGcode() > 2000000000) {

          if (namePart != "mu-" || nsecnd !=1) {
          // check that there is only one secondary and that namePart is mu-
            G4cout << "!!!!! " << namePart << " nsecondaries " 
                   << "nsecnd " << " expected mu- and 1 secondary" 
                   << G4endl;
          }
          trackManager->ProcessOneTrack( sTrack );
          G4TrackVector * secondaries = trackManager->GimmeSecondaries();

          G4int nssecnd = secondaries->size();

          if(verbose>=0 and iter == 50000*(iter/50000)) 
            G4cout << "##### " << nssecnd 
                   << " nsecondaries ##### of " << sName << G4endl;
          else if (verbose>=1)
             G4cout << "##### " << nssecnd 
                   << " nsecondaries ##### of " << sName << G4endl;

          histo.FillEvtMuonMinusBeam(secondaries);

          for (G4int k=0; k<nssecnd; ++k) {
            const G4Track*           ssTrack = (*secondaries)[k];
            G4TrackStatus ssts =  	ssTrack->GetTrackStatus();
            const G4DynamicParticle* ssdp = ssTrack->GetDynamicParticle();
            const G4ParticleDefinition* ssdpd = ssdp->GetDefinition();
            G4String ssName = ssdpd->GetParticleName();

            if (verbose>=1) {
              G4cout << "Secondary: " << k
                     << " is " << ssName
                     << ", " << ssdp->GetPDGcode()
                     << ", status: " << ssts
                     << " with kinetic energy of " 
                     << ssdp->GetKineticEnergy()/GeV 
                     << " and momentum of " 
                     << ssdp->GetTotalMomentum()/GeV << " GeV"
                     << G4endl;
            }
            delete (*secondaries)[k];
          }
          secondaries->clear();
        }
      }

      // we do not treat the intermediate mu_atom as a secondary

      if ( namePart != "mu-" || 
           nsecnd!=1 ||
           aChange->GetSecondary(0)->GetDynamicParticle()->GetPDGcode() < 2000000000 ) {
        histo.FillEvt(aChange);
      } 
      for (G4int i=0; i<nsecnd; ++i) 
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
    
    delete step; // cleanup G4Step instantiated at the beginning of the loop over configs

  }

  //G4cout << " mate "  << mate << G4endl;
  //delete mate;
  
  delete fin;
  
  delete nav;
  
  delete phys;

  delete trackManager;
  
  delete erndm;

  G4cout << "###### End of test #####" << G4endl;
}
