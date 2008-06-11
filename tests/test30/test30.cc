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
//      14.11.03 Renamed to cascade
//      09.05.06 Return back to test30
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4ios.hh"
#include <fstream>
#include <iomanip>

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
#include "G4HadronElasticDataSet.hh"
#include "G4IonsShenCrossSection.hh"
#include "G4TripathiCrossSection.hh"
#include "G4HadronElasticDataSet.hh"
#include "G4PiNuclearCrossSection.hh"

#include "G4BGGNucleonElasticXS.hh"
#include "G4BGGPionElasticXS.hh"
#include "G4BGGNucleonInelasticXS.hh"
#include "G4BGGPionInelasticXS.hh"

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
#include "G4IonTable.hh"
#include "G4DecayPhysics.hh"

#include "G4ForceCondition.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4Step.hh"
#include "G4GRSVolume.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4ExcitationHandler.hh"
#include "G4PreCompoundModel.hh"
#include "G4Evaporation.hh"

#include "G4StateManager.hh"
#include "G4NistManager.hh"

#include "Histo.hh"
#include "G4Timer.hh"

int main(int argc, char** argv)
{
  G4cout << "========================================================" << G4endl;
  G4cout << "======        Cascade Test (test30) Start       ========" << G4endl;
  G4cout << "========================================================" << G4endl;

  //-----------------------------------------------------------------------------
  // ------- Initialisation 

  Histo     histo;
  G4NistManager::Instance()->SetVerbose(2);
  G4String  namePart = "proton";
  G4bool    ionParticle = false;
  G4bool    Shen     = false;
  G4int     ionZ(0), ionA(0);
  G4String  nameMat  = "G4_Al";
  G4String  nameGen  = "binary";
  G4bool    logx     = false;
  G4bool    usepaw   = false;
  G4bool    isInitH  = false;
  G4bool    inclusive= true;
  G4int     verbose  = 0;
  G4double  energy   = 100.*MeV;
  G4double  sigmae   = 0.0;
  G4double  elim     = 30.*MeV;
  G4double  energyg  = 40.*MeV;
  G4double  balance  = 5.0*MeV;
  G4double  dangl    = 5.0;
  G4int     nevt     = 1000;
  G4int     nbins    = 100;
  G4int     nbinsa   = 40;
  G4int     nbinse   = 80;
  G4int     nbinsd   = 20;
  G4int     nbinspi  = 20;
  G4int     nangl    = 0;
  G4int     nanglpr  = 0;
  G4int     nanglpi  = 0;
  G4int     modu     = 10000;
  G4int     targetA  = 0;
  G4String hFile     = "hbook.paw";
  G4double theStep   = 0.01*micrometer;
  G4double range     = 1.0*micrometer;
  G4double  emax     = 160.*MeV;
  G4double  emaxpi   = 200.*MeV;
  G4double ebinlog   = 2.0*MeV;
  G4double eBound    = 70.*MeV;
  G4double tetmax    = 180.;
  G4double costmax   = -1.0;
  G4double xxl       = 100.;
  G4double kBound    = 0.2;
  G4Material* material = 0;
  G4bool nevap = false;
  G4bool gtran = false;
  G4bool gemis = false;
  G4bool xsbgg = true;
  G4bool elastic = false;

  G4double ang[20] = {0.0};
  G4double bng1[20] = {0.0};
  G4double bng2[20] = {0.0};
  G4double cng[20] = {0.0};
  G4double angpr[20] = {0.0};
  G4double bng1pr[20] = {0.0};
  G4double bng2pr[20] = {0.0};
  G4double cngpr[20] = {0.0};
  G4double angpi[20] = {0.0};
  G4double bngpi1[20] = {0.0};
  G4double bngpi2[20] = {0.0};
  G4double cngpi[20] = {0.0};
  float bestZ[250] = {
    0.0, 1.0, 1.0, 2.0, 2.0, 0.0, 0.0, 4.0, 0.0, 0.0,   //0
    4.0, 0.0, 0.0, 0.0, 6.0, 0.0, 0.0, 0.0, 9.0, 0.0,   //10
    10.0, 10.0, 11.0, 0.0, 11.0, 0.0, 13.0, 0.0, 0.0, 0.0,   //20
    0.0, 0.0, 15.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   //30
    0.0, 0.0, 0.0, 0.0, 21.0, 0.0, 21.0, 21.0, 23.0, 0.0,   //40
    0.0, 24.0, 25.0, 0.0, 25.0, 0.0, 27.0, 27.0, 27.0, 26.0,   //50
    27.0, 0.0, 0.0, 0.0, 0.0, 30.0, 31.0, 31.0, 32.0, 32.0,   //60
    33.0, 33.0, 33.0, 34.0, 33.0, 34.0, 35.0, 35.0, 0.0, 36.0,   //70
    36.0, 37.0, 36.0, 37.0, 37.0, 38.0, 39.0, 39.0, 41.0, 40.0,   //80
    41.0, 39.0, 41.0, 38.0, 39.0, 40.0, 40.0, 39.0, 40.0, 0.0,   //90
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   //100
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   //110
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   //120
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   //130
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   //140
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   //150
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   //160
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,          //170
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,          //180
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,          //190
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,          //200
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,          //210
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,          //220
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,          //230
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };        //240

  // Track
  G4ThreeVector aPosition = G4ThreeVector(0.,0.,0.);
  G4double      aTime     = 0. ;
  G4ThreeVector aDirection      = G4ThreeVector(0.0,0.0,1.0);
  G4double nx = 0.0, ny = 0.0, nz = 0.0;

  G4cout.setf( std::ios::scientific, std::ios::floatfield );

  // -------------------------------------------------------------------
  // Control on input

  if(argc < 2) {
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

  // -------------------------------------------------------------------
  //--------- Materials definition ---------
  Test30Material*  mate = new Test30Material();

  //--------- Particles definition ---------
  Test30Physics*   phys = new Test30Physics();

  const G4ParticleDefinition* gamma = G4Gamma::Gamma();
  const G4ParticleDefinition* proton = G4Proton::Proton();
  const G4ParticleDefinition* neutron = G4Neutron::Neutron();
  const G4ParticleDefinition* pin = G4PionMinus::PionMinus();
  const G4ParticleDefinition* pip = G4PionPlus::PionPlus();
  const G4ParticleDefinition* pi0 = G4PionZero::PionZero();
  const G4ParticleDefinition* deu = G4Deuteron::DeuteronDefinition();
  const G4ParticleDefinition* tri = G4Triton::TritonDefinition();
  const G4ParticleDefinition* alp = G4Alpha::AlphaDefinition();
  const G4ParticleDefinition* ion = G4GenericIon::GenericIon();

  G4DecayPhysics decays;
  decays.ConstructParticle();  

  G4ParticleTable* partTable = G4ParticleTable::GetParticleTable();
  partTable->SetReadiness();
  assert(ion);

  //--------- Geometry definition

  G4double dimX = 100.0*cm;
  G4double dimY = 100.0*cm;
  G4double dimZ = 100.0*cm;

  G4Box* sFrame = new G4Box ("Box",dimX, dimY, dimZ);
  G4LogicalVolume* lFrame = new G4LogicalVolume(sFrame,material,"Box",0,0,0);
  G4PVPlacement* pFrame = new G4PVPlacement(0,G4ThreeVector(),"Box",
                                            lFrame,0,false,0);
  assert(pFrame);

  // -------------------------------------------------------------------
  // -------- Loop over run

  G4String line, line1;
  G4bool end = true;
  for(G4int run=0; run<100; run++) {

    // ---- Read input file
    do {
      (*fin) >> line;
      G4cout << "Next line " << line << G4endl;
      if(line == "#particle") {
        (*fin) >> namePart;
      } else if(line == "#ion") {
        ionParticle= true;
	namePart="GenericIon";
        (*fin) >> ionA >> ionZ;
      } else if(line == "#energy(MeV)") {
        (*fin) >> energy;
        energy *= MeV;
        emax    = energy;
      } else if(line == "#sigmae(MeV)") {
        (*fin) >> sigmae;
        sigmae *= MeV;
      } else if(line == "#emax(MeV)") {
        (*fin) >> emax;
        emax *= MeV;
      } else if(line == "#emaxg(MeV)") {
        (*fin) >> energyg;
        energyg *= MeV;
      } else if(line == "#emaxpi(MeV)") {
        (*fin) >> emaxpi;
        emaxpi *= MeV;
      } else if(line == "#elim(MeV)") {
        (*fin) >> elim;
        elim *= MeV;
      } else if(line == "#ebalance(MeV)") {
        (*fin) >> balance;
        balance *= MeV;
      } else if(line == "#ebinlog(MeV)") {
        (*fin) >> ebinlog;
	if (ebinlog < 1.1) ebinlog = 1.1;
        ebinlog *= MeV;
      } else if(line == "#events") {
        (*fin) >> nevt;
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
      } else if(line == "#nangle") {
        (*fin) >> nangl;
      } else if(line == "#nanglepr") {
        (*fin) >> nanglpr;
      } else if(line == "#nanglepi") {
        (*fin) >> nanglpi;
      } else if(line == "#anglespi") {
        for(int k=0; k<nanglpi; k++) {(*fin) >> angpi[k];}
      } else if(line == "#dangle") {
        (*fin) >> dangl;
      } else if(line == "#angles") {
        for(int k=0; k<nangl; k++) {(*fin) >> ang[k];}
      } else if(line == "#anglespr") {
        for(int k=0; k<nanglpr; k++) {(*fin) >> angpr[k];}
      } else if(line == "#thetamax") {
        (*fin) >> tetmax;
        if(tetmax <= 0.0 || tetmax >= 180.) tetmax = 180.;
        costmax = std::cos(tetmax*degree); 
      } else if(line == "#range(mm)") {
        (*fin) >> range;
        range *= mm;
      } else if(line == "#step(mm)") {
        (*fin) >> theStep;
        theStep *= mm;
      } else if(line == "#print") {
        (*fin) >> modu;
      } else if(line == "#eBound(MeV)") {
        (*fin) >> eBound;
        eBound *= MeV;
      } else if(line == "#kBound") {
        (*fin) >> kBound;
      } else if(line == "#material") {
        (*fin) >> nameMat;
      } else if(line == "#targetA") {
        (*fin) >> targetA;
      } else if(line == "#Shen") {
        Shen = true;
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
      } else if(line == "#logx") {
        logx = true;
      } else if(line == "#xs_ghad") {
	xsbgg = false;
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

    if(!end) break;

    // -------------------------------------------------------------------
    // -------- Start run processing

    G4StateManager* g4State=G4StateManager::GetStateManager();
    if (! g4State->SetNewState(G4State_Init)) 
      G4cout << "error changing G4state"<< G4endl;;   

    G4cout << "###### Start new run # " << run << "   for "
	   << nevt << " events  #####" << G4endl;
    
    // -------- Target 

    material = mate->GetMaterial(nameMat);
    if(!material) {
      G4cout << "Material <" << nameMat
	     << "> is not found out"
	     << G4endl;
	     exit(1);
    }
    const G4Element* elm = material->GetElement(0); 

    G4int A = (G4int)(elm->GetN()+0.5);
    G4int Z = (G4int)(elm->GetZ()+0.5);
    if(targetA > 0) A = targetA;

    // -------- Projectile

    G4ParticleDefinition* part(0);
    if (!ionParticle) {
      part = (G4ParticleTable::GetParticleTable())->FindParticle(namePart);
    } else {
      part = (G4ParticleTable::GetParticleTable())->GetIon(ionZ, ionA, 0.);
    }
    if (! part ) {
      G4cout << " Sorry, No definition for particle" <<namePart 
	     << " found" << G4endl;
      G4Exception(" ");  
    }
    G4DynamicParticle dParticle(part,aDirection,energy);

    // ------- Select model

    G4VProcess* proc = phys->GetProcess(nameGen, namePart, material);
    G4ExcitationHandler* theDeExcitation = phys->GetDeExcitation();
    G4PreCompoundModel* thePreCompound = phys->GetPreCompound();
    if (gtran && thePreCompound) thePreCompound->UseGNASHTransition();
    if (gemis && thePreCompound) thePreCompound->UseHETCEmission();
    if (nevap) {
      G4Evaporation* evp = new G4Evaporation();
      evp->SetGEMChannel();
      theDeExcitation->SetEvaporation(evp);
    }
    G4double amass = phys->GetNucleusMass();
    phys->SetA(A);

    if(!proc) {
      G4cout << "For particle: " << part->GetParticleName()
	     << " generator " << nameGen << " is unavailable"
	     << G4endl;
	     exit(1);
    }

    // ------- Binning 

    G4int maxn = A + 1;

    G4cout << "The particle:  " << part->GetParticleName() << G4endl;
    G4cout << "The material:  " << material->GetName() 
	   << " Z= " << Z << " A= " << A
	   << "  Amax= " << maxn << G4endl;
    G4cout << "The step:      " << theStep/mm << " mm" << G4endl;
    G4cout << "The position:  " << aPosition/mm << " mm" << G4endl;
    G4cout << "The direction: " << aDirection << G4endl;
    G4cout << "The time:      " << aTime/ns << " ns" << G4endl;

    // linear histograms
    G4double mass = part->GetPDGMass();
    G4double pmax = std::sqrt(energy*(energy + 2.0*mass));
    G4double bine = emax/(G4double)nbinse;
    G4double bind = emax/(G4double)nbinsd;

    G4cout << "energy = " << energy/MeV << " MeV" << G4endl;
    G4cout << "emax   = " << emax/MeV << " MeV" << G4endl;
    G4cout << "pmax   = " << pmax/MeV << " MeV" << G4endl;

    // logarithmic histograms
    G4double logmax = 2.0;
    G4double logmin = -2.0;
    G4int nbinlog   = 80;
    if(emax > 100.*MeV) {
      logmax  = 3.0;
      logmin  = -1.0;
      nbinlog = 80;
    }
    if(emax > 1000.*MeV) {
      logmax  = 4.0;
      logmin  = -1.0;
      nbinlog = 100;
    }
    G4double binlog = (logmax - logmin)/G4double(nbinlog);
    G4cout << "### Log10 scale from " << logmin << " to " << logmax 
	   << " in " << nbinlog << " bins" << G4endl; 

    if(nameGen == "LElastic" || 
       nameGen == "BertiniElastic" ||
       nameGen == "elastic" ||
       nameGen == "HElastic" ||
       nameGen == "DElastic") elastic = true;

    // ------- Histograms

    if(usepaw && !isInitH) {

      isInitH = true;
    
      histo.add1D("1","Number of Secondaries",100,-0.5,99.5);
      histo.add1D("2","Type of secondary",10,-0.5,9.5);
      histo.add1D("3","Phi(degrees) of Secondaries",90,-180.0,180.0);

      histo.add1D("4","ds/dE for protons at theta = 0",nbinsd,0.,emax);
      histo.add1D("5","ds/dE for protons at theta = 1",nbinsd,0.,emax);
      histo.add1D("6","ds/dE for protons at theta = 2",nbinsd,0.,emax);
      histo.add1D("7","ds/dE for protons at theta = 3",nbinsd,0.,emax);
      histo.add1D("8","ds/dE for protons at theta = 4",nbinsd,0.,emax);
      histo.add1D("9","ds/dE for protons at theta = 5",nbinsd,0.,emax);
      histo.add1D("10","ds/dE for protons at theta = 6",nbinsd,0.,emax);
      histo.add1D("11","ds/dE for protons at theta = 7",nbinsd,0.,emax);
      histo.add1D("12","ds/dE for protons at theta = 8",nbinsd,0.,emax);
      histo.add1D("13","ds/dE for protons at theta = 9",nbinsd,0.,emax);
      histo.add1D("14","ds/dE for protons at theta = 10",nbinsd,0.,emax);

      G4int i;
      // proton double differencial histograms are active by request
      for(i=nanglpr; i<11; i++) {histo.activate(3+i, false);}

      histo.add1D("15","E(MeV) for gamma",200,0.,energyg);
      histo.add1D("16","delta E (MeV) ",100,-balance,balance);
      histo.add1D("17","delta Pz (MeV/c)",100,-balance,balance);
      histo.add1D("18","delta Pt (MeV/c)",100,-balance,balance);

      histo.add1D("19","E(MeV) for pi0",nbinse,0.,emax);
      histo.add1D("20","E(MeV) for pi+",nbinse,0.,emax);
      histo.add1D("21","E(MeV) for pi-",nbinse,0.,emax);
      histo.add1D("22","E(MeV) protons",nbinse,0.,emax);
      histo.add1D("23","E(MeV) neutrons",nbinse,0.,emax);

      histo.add1D("24","Phi(degrees) of neutrons",90,-180.0,180.0);

      histo.add1D("25","cos(theta) protons",nbinsa,costmax,1.);
      histo.add1D("26","cos(theta) neutrons",nbinsa,costmax,1.);

      histo.add1D("27","Baryon number (mbn)",maxn,-0.5,(G4double)maxn + 0.5);

      histo.add1D("28","ds/dE for neutrons at theta = 0",nbinsd,0.,emax);
      histo.add1D("29","ds/dE for neutrons at theta = 1",nbinsd,0.,emax);
      histo.add1D("30","ds/dE for neutrons at theta = 2",nbinsd,0.,emax);
      histo.add1D("31","ds/dE for neutrons at theta = 3",nbinsd,0.,emax);
      histo.add1D("32","ds/dE for neutrons at theta = 4",nbinsd,0.,emax);
      histo.add1D("33","ds/dE for neutrons at theta = 5",nbinsd,0.,emax);
      histo.add1D("34","ds/dE for neutrons at theta = 6",nbinsd,0.,emax);
      histo.add1D("35","ds/dE for neutrons at theta = 7",nbinsd,0.,emax);
      histo.add1D("36","ds/dE for neutrons at theta = 8",nbinsd,0.,emax);
      histo.add1D("37","ds/dE for neutrons at theta = 9",nbinsd,0.,emax);
      histo.add1D("38","ds/dE for neutrons at theta = 10",nbinsd,0.,emax);
      histo.add1D("39","ds/dE for neutrons at theta = 11",nbinsd,0.,emax);
      histo.add1D("40","ds/dE for neutrons at theta = 12",nbinsd,0.,emax);

      // neutron double differencial histograms are active by request
      for(i=nangl; i<13; i++) {histo.activate(27+i, false);}
 
      histo.add1D("41","ds/dE for pi- at theta = 0",nbinspi,0.,emaxpi);
      histo.add1D("42","ds/dE for pi- at theta = 1",nbinspi,0.,emaxpi);
      histo.add1D("43","ds/dE for pi- at theta = 2",nbinspi,0.,emaxpi);
      histo.add1D("44","ds/dE for pi- at theta = 3",nbinspi,0.,emaxpi);
      histo.add1D("45","ds/dE for pi- at theta = 4",nbinspi,0.,emaxpi);
      histo.add1D("46","ds/dE for pi+ at theta = 0",nbinspi,0.,emaxpi);
      histo.add1D("47","ds/dE for pi+ at theta = 1",nbinspi,0.,emaxpi);
      histo.add1D("48","ds/dE for pi+ at theta = 2",nbinspi,0.,emaxpi);
      histo.add1D("49","ds/dE for pi+ at theta = 3",nbinspi,0.,emaxpi);
      histo.add1D("50","ds/dE for pi+ at theta = 4",nbinspi,0.,emaxpi);

      // pion double differencial histograms are active by request
      for(i=nanglpi; i<5; i++) {
	histo.activate(40+i, false);
	histo.activate(45+i, false);
      }
      histo.add1D("51","E(MeV) neutrons",nbinlog,logmin,logmax);
      histo.add1D("52","ds/dE for neutrons at theta = 0",nbinlog,logmin,logmax);
      histo.add1D("53","ds/dE for neutrons at theta = 1",nbinlog,logmin,logmax);
      histo.add1D("54","ds/dE for neutrons at theta = 2",nbinlog,logmin,logmax);
      histo.add1D("55","ds/dE for neutrons at theta = 3",nbinlog,logmin,logmax);
      histo.add1D("56","ds/dE for neutrons at theta = 4",nbinlog,logmin,logmax);
      histo.add1D("57","ds/dE for neutrons at theta = 5",nbinlog,logmin,logmax);

      // neutron double differencial histograms are active by request
      for(i=nangl; i<6; i++) {histo.activate(51+i, false);}

      if(elastic) {
	histo.add1D("58","Ekin (MeV) for primary particle",120,0.,energy*1.2/MeV);
	histo.add1D("59","cos(Theta) for recoil particle in Lab.Sys.",nbinsa,costmax,1.);
	histo.add1D("60","cos(Theta) for primary particle in Lab.Sys.",nbinsa,costmax,1.);
	histo.add1D("61","cos(Theta) for recoil particle in CM.Sys.",nbinsa,costmax,1.);
	histo.add1D("62","cos(Theta) for primary particle in CM.Sys.",nbinsa,costmax,1.);
	histo.add1D("63","Ekin (MeV) for recoil particle",120,0.,energy*1.2/MeV);
	histo.add1D("64","Theta (degree) for primary particle in Lab.Sys.",nbinsa,0.0,tetmax);
        G4double x2 = std::log10(tetmax);
        G4double x1 = x2 - std::log10(nbinsa);
        xxl = x2 - x1;
	histo.add1D("65","log10(theta (degree)) for primary particle in Lab.Sys.",nbinsa,x1,x2);
	// desactivate not needed hist for elastic
	histo.activate(50, false);
	for(i=0; i<13; i++) {histo.activate(14+i, false);}
      }
    }

    if(usepaw) {
      histo.setFileName(hFile);
      histo.book();
      G4cout << "Histograms are booked output file <" << hFile << "> "
	     << G4endl;
    }

    // -------- Normalisation

    G4VCrossSectionDataSet* cs = 0;
    G4double cross_sec = 0.0;

    if(nameGen == "LElastic" || nameGen == "BertiniElastic" ) {
      cs = new G4HadronElasticDataSet();
    } else if (nameGen == "elastic" || 
	       nameGen == "HElastic" || 
	       nameGen == "DElastic") {
      if(Z == 1) cs = new G4HadronElasticDataSet();
      else if(part == proton || part == neutron) {
	if(xsbgg) cs = new G4BGGNucleonElasticXS(part);
	else      cs = new G4HadronElasticDataSet();
      } else if(part == pip || part == pin) {
	if(xsbgg) cs = new G4BGGPionElasticXS(part);
	else      cs = new G4HadronElasticDataSet();
      } else {
	cs = new G4HadronElasticDataSet();
      }
    } else if(part == proton && Z > 1 && nameGen != "lepar") {
      if(xsbgg) cs = new G4BGGNucleonInelasticXS(part);
      else      cs = new G4ProtonInelasticCrossSection();
    } else if(part == neutron && Z > 1 && nameGen != "lepar") {
      if(xsbgg) cs = new G4BGGNucleonInelasticXS(part);
      else      cs = new G4NeutronInelasticCrossSection();
    } else if((part == pin || part == pip) && Z > 1 && nameGen != "lepar") {
      if(xsbgg) cs = new G4BGGPionInelasticXS(part);
      else cs = new G4PiNuclearCrossSection();
    } else if( ionParticle ) {
      if ( Shen ) {
        cs = new G4IonsShenCrossSection();
	G4cout << "Using Shen Cross section for Ions" << G4endl;
      } else {
	cs = new G4TripathiCrossSection();
	G4cout << "Using Tripathi Cross section for Ions" << G4endl;
      }
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

    G4double factor  = 
      cross_sec*MeV*1000.0*(G4double)nbinse/(energy*barn*(G4double)nevt);
    G4double factora = 
      cross_sec*1000.0*(G4double)nbinsa/(twopi*(1.0 - costmax)*barn*(G4double)nevt);
    G4double factoraa = 
      cross_sec*1000.0*(G4double)nbinsa/(twopi*tetmax*degree*barn*(G4double)nevt);
    G4double factoral = 
      cross_sec*1000.0*(G4double)nbinsa/(twopi*xxl*std::log(10.)*barn*(G4double)nevt);
    G4double factorb= cross_sec*1000.0/(barn*(G4double)nevt);
    G4cout << "### factor  = " << factor
           << "### factora = " << factora
           << "### factorb = " << factorb
           << "    cross(mb)= " << cross_sec*1000./barn << G4endl;

    // -------- Limit number of angles

    if(nangl   > 13) nangl   = 13;
    if(nanglpr > 11) nanglpr = 11;
    if(nanglpi >  5) nanglpi = 5;

    if(nangl > 0) {
      for(G4int k=0; k<nangl; k++) {

        if(nangl == 1) {
          bng1[0] = std::max(0.0,ang[0] - dangl);
          bng2[0] = std::min(180., ang[0] + dangl);
        } else if(k == 0) {
          bng1[0] = std::max(0.0,ang[0] - dangl);
          bng2[0] = std::min(0.5*(ang[0] + ang[1]), ang[0] + dangl);
        } else if(k < nangl-1) {
          bng1[k] = std::max(bng2[k-1], ang[k]-dangl);
          bng2[k] = std::min(0.5*(ang[k] + ang[k+1]), ang[k] + dangl);
        } else {
          bng1[k] = std::max(bng2[k-1], ang[k]-dangl);
          bng2[k] = std::min(180., ang[k] + dangl);
        }

        cng[k] = cross_sec*MeV*1000.0*(G4double)nbinsd/
         (twopi*(std::cos(degree*bng1[k]) - std::cos(degree*bng2[k]))*
                barn*emax*(G4double)nevt);
      }
    }

    if(nanglpr > 0) {
      for(G4int k=0; k<nanglpr; k++) {

        if(nanglpr == 1) {
          bng1pr[0] = std::max(0.0,angpr[0] - dangl);
          bng2pr[0] = std::min(180., angpr[0] + dangl);
        } else if(k == 0) {
          bng1pr[0] = std::max(0.0,angpr[0] - dangl);
          bng2pr[0] = std::min(0.5*(angpr[0] + angpr[1]), angpr[0] + dangl);
        } else if(k < nanglpr-1) {
          bng1pr[k] = std::max(bng2pr[k-1], angpr[k]-dangl);
          bng2pr[k] = std::min(0.5*(angpr[k] + angpr[k+1]), angpr[k] + dangl);
        } else {
          bng1pr[k] = std::max(bng2pr[k-1], angpr[k]-dangl);
          bng2pr[k] = std::min(180., angpr[k] + dangl);
        }

        cngpr[k] = cross_sec*MeV*1000.0*(G4double)nbinsd/
         (twopi*(std::cos(degree*bng1pr[k]) - std::cos(degree*bng2pr[k]))*
                barn*emax*(G4double)nevt);
      }
    }

    if(nanglpi > 0) {
      for(G4int k=0; k<nanglpi; k++) {

        if(nangl == 1) {
          bngpi1[0] = std::max(0.0,angpi[0] - dangl);
          bngpi2[0] = std::min(180., angpi[0] + dangl);
        } else if(k == 0) {
          bngpi1[0] = std::max(0.0,angpi[0] - dangl);
          bngpi2[0] = std::min(0.5*(angpi[0] + angpi[1]), angpi[0] + dangl);
        } else if(k < nanglpi-1) {
          bngpi1[k] = std::max(bngpi2[k-1], angpi[k]-dangl);
          bngpi2[k] = std::min(0.5*(angpi[k] + angpi[k+1]), angpi[k] + dangl);
        } else {
          bngpi1[k] = std::max(bngpi2[k-1], angpi[k]-dangl);
          bngpi2[k] = std::min(180., angpi[k] + dangl);
        }

        cngpi[k] = cross_sec*MeV*1000.0*(G4double)nbinspi/
         (twopi*(std::cos(degree*bngpi1[k]) - std::cos(degree*bngpi2[k]))*
                 barn*emaxpi*(G4double)nevt);
      }
    }

    // -------- Track

    G4Track* gTrack;
    gTrack = new G4Track(&dParticle,aTime,aPosition);

    // -------- Step

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
    G4RotationMatrix* rot = new G4RotationMatrix();
    G4double phi0 = aDirection.phi();
    G4double theta0 = aDirection.theta();
    rot->rotateZ(-phi0);
    rot->rotateY(-theta0);

    G4cout << "Test rotation= " << (*rot)*(aDirection) << G4endl;

    G4Timer* timer = new G4Timer();
    timer->Start();
    const G4DynamicParticle* sec = 0;
    G4ParticleDefinition* pd;
    G4ThreeVector  mom;
    G4LorentzVector labv, fm;
    G4double e, p, m, px, py, pz, pt, theta;
    G4VParticleChange* aChange = 0;
    G4int warn = 0;

    // -------- Event loop

    for (G4int iter=0; iter<nevt; iter++) {

      if(verbose>1) 
        G4cout << "### " << iter << "-th event start " << G4endl;

      G4double e0 = energy;
      do {
        if(sigmae > 0.0) e0 = G4RandGauss::shoot(energy,sigmae);
      } while (e0 < 0.0);

      dParticle.SetKineticEnergy(e0);

      gTrack->SetStep(step);
      gTrack->SetKineticEnergy(e0);

      labv = G4LorentzVector(0.0, 0.0, std::sqrt(e0*(e0 + 2.*mass)), 
			     e0 + mass + amass);
      G4ThreeVector bst = labv.boostVector();
      
      aChange = proc->PostStepDoIt(*gTrack,*step);

      // take into account local energy deposit
      G4double de = aChange->GetLocalEnergyDeposit();
      G4LorentzVector dee = G4LorentzVector(0.0, 0.0, 0.0, de); 
      labv -= dee;

      G4int n = aChange->GetNumberOfSecondaries();

      if(iter == modu*(iter/modu)) 
        G4cout << "##### " << iter << "-th event  #####" << G4endl;

      G4int nbar = 0;

      for(G4int j=0; j<n; j++) {

        sec = aChange->GetSecondary(j)->GetDynamicParticle();
        pd  = sec->GetDefinition();
        if(pd->GetPDGMass() > 100.*MeV) nbar++;
      }

      for(G4int i=0; i<n; i++) {

        sec = aChange->GetSecondary(i)->GetDynamicParticle();
        pd  = sec->GetDefinition();
        fm  = sec->Get4Momentum();
        mom = sec->GetMomentum();

	// for exclusive reaction 2 particles in final state
        if(!inclusive && nbar != 2) break;

        m = pd->GetPDGMass();
	p = mom.mag();
        labv -= fm;
        px = mom.x();
        py = mom.y();
        pz = mom.z();
        pt = std::sqrt(px*px +py*py);
        e  = fm.e() - m;

        theta = mom.theta();
        G4double cost  = std::cos(theta);
        G4double thetad = theta/degree;

        fm.boost(-bst);
        G4double costcm = std::cos(fm.theta());

	if(usepaw) {
          if(elastic) {
	    if(i==0)  {
	      histo.fill(57,e/MeV,1.0);
	      histo.fill(59,cost,factora);
	      histo.fill(61,costcm,factora);
	    } else if(i==1) {
	      histo.fill(58,e/MeV,1.0);
	      histo.fill(60,cost,factora);
	      histo.fill(62,costcm,factora);
	      histo.fill(63,thetad,factoraa/std::sin(theta));
	      histo.fill(64,std::log10(thetad),factoral*theta/std::sin(theta));
	    }
	  }

          histo.fill(2,mom.phi()/degree,1.0);
          if(pd == neutron) histo.fill(23,mom.phi()/degree,1.0);
	}

	if( (e == 0.0 || pt == 0.0) && warn < 50 ) {
          warn++;
          G4cout << "Warning! evt# " << iter 
	         << "  " << i << "-th sec  "
		 << pd->GetParticleName() << "   Ekin(MeV)= "
                 << e/MeV
                 << " Pt(MeV/c)= " << pt/MeV
		 << G4endl;
	}
	de += e;
        if((verbose>0 || std::fabs(mom.phi()/degree - 90.) < 0.001 ) && 
	   warn < 50 && verbose>1) {
          warn++;
          G4cout << "Warning! evt# " << iter 
                 << "  " << i << "-th sec  "
		 << pd->GetParticleName() << "  Ekin(MeV)= "
                 << e/MeV
		 << "  p(MeV)= " << mom/MeV
		 << "  m(MeV)= " << m/MeV
		 << "  Etot(MeV)= " << (e+m)/MeV
		 << "  pt(MeV)= " << pt/MeV
                 << "  std::sin(tet)= " << pt/p
                 << "  phi(deg)= " << mom.phi()/degree
                 << G4endl;
        }

	if(usepaw) {

          if(pd) {
            G4double N = pd->GetBaryonNumber();
            G4double Z1= pd->GetPDGCharge()/eplus;
            G4double Z0= bestZ[(G4int)N];
            if(std::fabs(Z0 - Z1) < 0.1 || Z0 == 0.0) 
	      histo.fill(26, N, factorb);
	  }

          if(pd == proton) {

            histo.fill(1,1.0, 1.0);
	    histo.fill(21,e/MeV, factor);
	    histo.fill(24,cost, factora);
            for(G4int kk=0; kk<nanglpr; kk++) {
              if(bng1pr[kk] <= thetad && thetad <= bng2pr[kk]) {
                histo.fill(3+kk,e/MeV, cngpr[kk]);
                break;
	      }
	    }

          } else if(pd == pin) {

	    histo.fill(1,4.0, 1.0);
            histo.fill(20,e/MeV, 1.0);
            for(G4int kk=0; kk<nanglpi; kk++) {
              if(bngpi1[kk] <= thetad && thetad <= bngpi2[kk]) {
                histo.fill(40+kk,e/MeV, cngpi[kk]);
                break;
	      }
	    }

          } else if(pd == pip) {

	    histo.fill(1,3.0, 1.0);
            histo.fill(19,e/MeV, 1.0);
            for(G4int kk=0; kk<nanglpi; kk++) {
              if(bngpi1[kk] <= thetad && thetad <= bngpi2[kk]) {
                histo.fill(45+kk,e/MeV, cngpi[kk]);
                break;
	      }
	    }

	  } else if(pd == pi0) {

	    histo.fill(1,5.0, 1.0);
	    histo.fill(18,e/MeV, 1.0);

	  } else if(pd == neutron) {

	    histo.fill(1,2.0, 1.0);
	    histo.fill(22,e/MeV, factor);
            G4double ee = std::log10(e/MeV);
	    G4double e2 = ee;
            G4bool islog= false;
            if(ee >= logmin && ee <= logmax) { 
	      islog = true;
	      G4int nbb = G4int(((ee - logmin)/binlog));
              G4double e1 = logmin + binlog*nbb;
	      e2 = std::pow(10.,e1 + binlog) - std::pow(10.,e1);
	      histo.fill(50, ee, factor*bine/e2);
	    } 
	    if(e >= elim) histo.fill(25, cost, factora);
            for(G4int kk=0; kk<nangl; kk++) {
              if(bng1[kk] <= thetad && thetad <= bng2[kk]) {
                histo.fill(27+kk,e/MeV, cng[kk]);
                if(islog && kk < 6) histo.fill(51+kk,ee,cng[kk]*bind/e2);
                break;
	      }
	    }

	  } else if(pd == gamma) {
	    histo.fill(14,e/MeV, 1.0);

	  } else if(pd == deu) {
	    histo.fill(1,6.0, 1.0);
	  } else if(pd == tri) {
	    histo.fill(1,7.0, 1.0);
	  } else if(pd == alp) {
	    histo.fill(1,8.0, 1.0);
	  } else {
	    histo.fill(1,9.0, 1.0);
	  }
	}
	//	delete sec;       	 
        delete aChange->GetSecondary(i);
      }

      if(verbose > 0) 
        G4cout << "Energy/Momentum balance= " << labv << G4endl;

      px = labv.px();
      py = labv.py();
      pz = labv.pz();
      p  = std::sqrt(px*px +py*py + pz*pz);
      pt = std::sqrt(px*px +py*py);

      if(usepaw) {
        histo.fill(0,(G4double)n,1.0);
	histo.fill(15,labv.e()/MeV, 1.0);
	histo.fill(16,pz/MeV, 1.0);
	histo.fill(17,pt/MeV, 1.0);
      }
      aChange->Clear();
    }

    timer->Stop();
    G4cout << "  "  << *timer << G4endl;
    delete timer;

    // -------- Committing the transaction with the tree

    if(usepaw) {
      G4cout << "###### Save histograms" << G4endl;
      histo.save();
    }
    G4cout << "###### End of run # " << run << "     ######" << G4endl;
  }

  delete mate;
  delete fin;
  delete phys;

  G4cout << "###### End of test #####" << G4endl;
}
