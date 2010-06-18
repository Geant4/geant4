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
//      GEANT 4 class file 
//      CERN Geneva Switzerland
//
//
//      File name:     Test43
//
//      Author:        Gunter Folger
//           created from test30 files originally by Vladimir Ivanchenko
//
//      Creation date: 2007
//
//      Modifications:
//
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4ios.hh"
#include <fstream>
#include <iomanip>

#include "G4Material.hh"
#include "G4ElementVector.hh"
#include "Test43Material.hh"
#include "Test43Physics.hh"
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
#include "G4PiNuclearCrossSection.hh"

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
#include "G4ForceCondition.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4Step.hh"
#include "G4GRSVolume.hh"

#include "G4UnitsTable.hh"
#include "G4ExcitationHandler.hh"
#include "G4PreCompoundModel.hh"
#include "G4Evaporation.hh"

#ifdef G4FPE_DEBUG
  #include "G4FPEDetection.hh"
#endif

#include "G4StateManager.hh"

#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TTree.h"

#include "G4Timer.hh"

G4double Rapidity(const G4double p, const G4double E)
{
	return 0.5*log((E+p)/(E-p));
} 

int main(int argc, char** argv)
{
  G4cout << "========================================================" << G4endl;
  G4cout << "======             Cascade Test Start           ========" << G4endl;
  G4cout << "========================================================" << G4endl;

#ifdef G4FPE_DEBUG
  InvalidOperationDetection();
#endif

  G4String  namePart = "proton";
  G4bool    ionParticle = false;
  G4bool    Shen = false;
  G4int     ionZ(0), ionA(0);
  G4String  nameMat  = "Si";
  G4String  nameGen  = "stringCHIPS";
  G4bool    logx     = false;
  G4bool    usepaw   = false;
  G4bool    inclusive= true;
  G4int     verbose  = 0;
  G4double  Kinetic_energy   = 100.*GeV;
  G4double  sigmae   = 0.0;
  G4double  elim     = 30.*GeV;
  G4double  dangl    = 5.0;
  G4int     nevt     = 1000;
  G4int     nbins    = 100;
  G4int     nbinsa   = 40;
  G4int     nbinse   = 80;
  G4int     nbinsd   = 20;
  G4int     nbinspi  = 20;
  G4int     nangl    = 0;
  G4int     nanglpi  = 0;
  G4String hFile     = "test43.root";
  G4double theStep   = 0.01*micrometer;
  G4double range     = 1.0*micrometer;
  G4double  emax     = 160.*GeV;
  G4double  emaxpi   = 200.*GeV;
  G4double ebinlog   = 2.0*GeV;
  G4double eBound    = 70.*GeV;
  G4double kBound    = 0.2;
  G4Material* material = 0;
  G4bool nevap = false;
  G4bool gtran = false;
  G4bool gemis = false;

  G4double ang[15] = {0.0};
  G4double bng1[15] = {0.0};
  G4double bng2[15] = {0.0};
  G4double cng[15] = {0.0};
  G4double angpi[10] = {0.0};
  G4double bngpi1[10] = {0.0};
  G4double bngpi2[10] = {0.0};
  G4double cngpi[10] = {0.0};



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

  Test43Material*  mate = new Test43Material();
  Test43Physics*   phys = new Test43Physics();

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
  G4cout << "#exclusive" << G4endl;
  G4cout << "#inclusive" << G4endl;
  G4cout << "#logx" << G4endl;
  G4cout << "#nbins" << G4endl;
  G4cout << "#nbinsa" << G4endl;
  G4cout << "#nbinse" << G4endl;
  G4cout << "#nbinsd" << G4endl;
  G4cout << "#nbinspi" << G4endl;
  G4cout << "#nanglepi" << G4endl;
  G4cout << "#anglespi" << G4endl;
  G4cout << "#nangle" << G4endl;
  G4cout << "#angles" << G4endl;
  G4cout << "#dangle" << G4endl;
  G4cout << "#particle" << G4endl;
  G4cout << "#ion Z A" << G4endl;
  G4cout << "#Shen (use Shen cross-section)" << G4endl;
  G4cout << "#Kinetic_energy(GeV)" << G4endl;
  G4cout << "#sigmae(GeV)" << G4endl;
  G4cout << "#emax(GeV)" << G4endl;
  G4cout << "#emaxpi(GeV)" << G4endl;
  G4cout << "#elim(GeV)" << G4endl;
  G4cout << "#ebinlog(GeV)" << G4endl;
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
  G4cout << "#HETCEmission" << G4endl;
  G4cout << "#GNASHTransition" << G4endl;
  G4cout << "#GEMEvaporation" << G4endl;
  G4cout << "#eBound" << G4endl;
  G4cout << "#kBound" << G4endl;

  const G4ParticleDefinition* gamma = G4Gamma::Gamma();
  const G4ParticleDefinition* proton = G4Proton::Proton();
  const G4ParticleDefinition* neutron = G4Neutron::Neutron();
  const G4ParticleDefinition* pim = G4PionMinus::PionMinus();
  const G4ParticleDefinition* pip = G4PionPlus::PionPlus();
  const G4ParticleDefinition* pi0 = G4PionZero::PionZero();
  const G4ParticleDefinition* deu = G4Deuteron::DeuteronDefinition();
  const G4ParticleDefinition* tri = G4Triton::TritonDefinition();
  const G4ParticleDefinition* alp = G4Alpha::AlphaDefinition();
  const G4ParticleDefinition* ion = G4GenericIon::GenericIon();
  G4ParticleTable* partTable = G4ParticleTable::GetParticleTable();
  partTable->SetReadiness();

  G4String line, line1;
  G4bool end = true;
  for(G4int run=0; run<100; run++) {
    do {
      (*fin) >> line;
      G4cout << "Next line " << line << G4endl;
      if(line == "#particle") {
        (*fin) >> namePart;
      } else if(line == "#ion") {
        ionParticle= true;
	namePart="GenericIon";
        (*fin) >> ionA >> ionZ;
      } else if(line == "#Kinetic_energy(GeV)") {
        (*fin) >> Kinetic_energy;
        Kinetic_energy *= GeV;
        emax    = Kinetic_energy;
      } else if(line == "#sigmae(GeV)") {
        (*fin) >> sigmae;
        sigmae *= GeV;
      } else if(line == "#emax(GeV)") {
        (*fin) >> emax;
        emax *= GeV;
      } else if(line == "#emaxpi(GeV)") {
        (*fin) >> emaxpi;
        emaxpi *= GeV;
      } else if(line == "#elim(GeV)") {
        (*fin) >> elim;
        elim *= GeV;
      } else if(line == "#ebinlog(GeV)") {
        (*fin) >> ebinlog;
	if (ebinlog < 1.1) ebinlog = 1.1;
        ebinlog *= GeV;
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
      } else if(line == "#nanglepi") {
        (*fin) >> nanglpi;
      } else if(line == "#anglespi") {
        for(int k=0; k<nanglpi; k++) {(*fin) >> angpi[k];}
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
      } else if(line == "#eBound(GeV)") {
        (*fin) >> eBound;
        eBound *= GeV;
      } else if(line == "#kBound") {
        (*fin) >> kBound;
      } else if(line == "#material") {
        (*fin) >> nameMat;
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

    G4cout << "###### Start new run # " << run << "     #####" << G4endl;

    if ( ionParticle ) {
       Kinetic_energy*=ionA;
    }
    material = mate->GetMaterial(nameMat);
    if(!material) {
      G4cout << "Material <" << nameMat
	     << "> is not found out"
	     << G4endl;
	     exit(1);
    }
    G4ParticleDefinition* part(0);
    if (namePart != "GenericIon") {
       part = (G4ParticleTable::GetParticleTable())->FindParticle(namePart);
    } else {
        G4StateManager* g4State=G4StateManager::GetStateManager();
	if (! g4State->SetNewState(G4State_Init)) G4cout << "error changing G4state"<< G4endl;;   
       G4IonTable ions;
       part = ions.GetIon(ionZ, ionA);
    }
   if (! part ) {
       G4cout << " Sorry, No definition for particle" <<namePart << " found" << G4endl;
       G4Exception(" "); 
       
    }
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
    G4double nucleus_mass = phys->GetNucleusMass();

    if(!proc) {
      G4cout << "For particle: " << part->GetParticleName()
	     << " generator " << nameGen << " is unavailable"
	     << G4endl;
	     exit(1);
    }

    G4int maxn = (G4int)((*(material->GetElementVector()))[0]->GetN()) + 1;
    // G4int maxz = (G4int)((*(material->GetElementVector()))[0]->GetZ()) + 1;

    G4cout << "The particle:  " << part->GetParticleName() << G4endl;
    G4cout << "The material:  " << material->GetName() << "  Amax= " << maxn << G4endl;
    G4cout << "The step:      " << theStep/mm << " mm" << G4endl;
    G4cout << "The position:  " << aPosition/mm << " mm" << G4endl;
    G4cout << "The direction: " << aDirection << G4endl;
    G4cout << "The time:      " << aTime/ns << " ns" << G4endl;


    // -------------------------------------------------------------------

    TFile * rootFile = new TFile(hFile,"CREATE");
    if ( ! rootFile ) 
    {
       std::cout << " Fail to create root file " << std::endl;
       exit(1);
    }   
    const G4int nhisto = 63;
    TH1 * h[nhisto];

    G4double mass = part->GetPDGMass();
    G4double pmax = std::sqrt(Kinetic_energy*(Kinetic_energy + 2.0*mass));
    G4double binlog = std::log10(ebinlog);
    G4int nbinlog = (G4int)(std::log10(2.0*emax)/binlog);
    G4double logmax = binlog*nbinlog;
    G4double bine = emax/(G4double)nbinse;
    G4double bind = emax/(G4double)nbinsd;
    G4double p_pmax=pmax*1.2;
    G4double p_Kinetic_energy=Kinetic_energy*1.2;


      // Creating an 1-dimensional histogram in the root directory of the tree

      // ---- Book a histogram and ntuples
      G4cout << "Hbook file name: <" << hFile << ">" << G4endl;
      G4cout << "Kinetic_energy = " << Kinetic_energy/GeV << " GeV" << G4endl;
      G4cout << "emax   = " << emax/GeV << " GeV" << G4endl;
      G4cout << "pmax   = " << pmax/GeV << " GeV" << G4endl;

      h[0]=new TH1F("1","Number of Secondaries",90,-0.5,89.5);
      h[1]=new TH1F("2","Type of secondary",10,-0.5,9.5);
      h[2]=new TH1F("3","Phi(degrees) of Secondaries",90,-180.0,180.0);
      h[3]=new TH1F("4","Pz (GeV) for protons",100,-p_pmax/GeV,p_pmax/GeV);
      h[4]=new TH1F("5","Pz (GeV) for pi-",100,-p_pmax/GeV,p_pmax/GeV);
      h[5]=new TH1F("6","Pz (GeV) for pi+",100,-p_pmax/GeV,p_pmax/GeV);
      h[6]=new TH1F("7","Pz (GeV) for neutrons",100,-p_pmax/GeV,p_pmax/GeV);
      h[7]=new TH1F("8","Pt (GeV) for protons",100,0.,p_pmax/GeV);
      h[8]=new TH1F("9","Pt (GeV) for pi-",100,0.,p_pmax/GeV);
      h[9]=new TH1F("10","Pt (GeV) for pi+",100,0.,p_pmax/GeV);
      h[10]=new TH1F("11","Pt (GeV) for neutrons",100,0.,p_pmax/GeV);
      h[11]=new TH1F("12","E (GeV) for protons",100,0.,p_Kinetic_energy/GeV);
      h[12]=new TH1F("13","E (GeV) for pi-",100,0.,p_Kinetic_energy/GeV);
      h[13]=new TH1F("14","E (GeV) for pi+",100,0.,p_Kinetic_energy/GeV);
      h[14]=new TH1F("15","E (GeV) for neutrons",100,0.,p_Kinetic_energy/GeV);
      h[15]=new TH1F("16","delta E (GeV)",20,-1.,1.);
      h[16]=new TH1F("17","delta Pz (GeV)",20,-1.,1.);
      h[17]=new TH1F("18","delta Pt (GeV)",20,-1.,1.);

      h[18]=new TH1F("19","E (GeV) for pi0",100,0.,p_Kinetic_energy/GeV);
      h[19]=new TH1F("20","Pz (GeV) for pi0",100,-p_pmax/GeV,p_pmax/GeV);
      h[20]=new TH1F("21","Pt (GeV) for pi0",100,0.,p_pmax/GeV);

      h[21]=new TH1F("22","E(GeV) protons",nbinse,0.,emax/GeV);
      h[22]=new TH1F("23","E(GeV) neutrons",nbinse,0.,emax/GeV);

      h[23]=new TH1F("24","Phi(degrees) of neutrons",90,-180.0,180.0);

      h[24]=new TH1F("25","cos(theta) protons",nbinsa,-1.,1.);
      h[25]=new TH1F("26","cos(theta) neutrons",nbinsa,-1.,1.);

      h[26]=new TH1F("27","Baryon number (mbn)",maxn,-0.5,(G4double)maxn + 0.5);

      if(nangl>0)
       h[27]=new TH1F("28","ds/dE for neutrons at theta = 0",nbinsd,0.,emax/GeV);
      if(nangl>1)
       h[28]=new TH1F("29","ds/dE for neutrons at theta = 1",nbinsd,0.,emax/GeV);
      if(nangl>2)
       h[29]=new TH1F("30","ds/dE for neutrons at theta = 2",nbinsd,0.,emax/GeV);
      if(nangl>3)
       h[30]=new TH1F("31","ds/dE for neutrons at theta = 3",nbinsd,0.,emax/GeV);
      if(nangl>4)
       h[31]=new TH1F("32","ds/dE for neutrons at theta = 4",nbinsd,0.,emax/GeV);
      if(nangl>5)
       h[32]=new TH1F("33","ds/dE for neutrons at theta = 5",nbinsd,0.,emax/GeV);
      if(nangl>6)
       h[33]=new TH1F("34","ds/dE for neutrons at theta = 6",nbinsd,0.,emax/GeV);
      if(nangl>7)
       h[34]=new TH1F("35","ds/dE for neutrons at theta = 7",nbinsd,0.,emax/GeV);
      if(nangl>8)
       h[35]=new TH1F("36","ds/dE for neutrons at theta = 8",nbinsd,0.,emax/GeV);
      if(nangl>9)
       h[36]=new TH1F("37","ds/dE for neutrons at theta = 9",nbinsd,0.,emax/GeV);
      if(nangl>10)
       h[37]=new TH1F("38","ds/dE for neutrons at theta = 10",nbinsd,0.,emax/GeV);
      if(nangl>11)
       h[38]=new TH1F("39","ds/dE for neutrons at theta = 11",nbinsd,0.,emax/GeV);
      if(nangl>12)
       h[39]=new TH1F("40","ds/dE for neutrons at theta = 12",nbinsd,0.,emax/GeV);

      if(nanglpi>0)
       h[40]=new TH1F("41","ds/dE for pi- at theta = 0",nbinspi,0.,emaxpi/GeV);
      if(nanglpi>1)
       h[41]=new TH1F("42","ds/dE for pi- at theta = 1",nbinspi,0.,emaxpi/GeV);
      if(nanglpi>2)
       h[42]=new TH1F("43","ds/dE for pi- at theta = 2",nbinspi,0.,emaxpi/GeV);
      if(nanglpi>3)
       h[43]=new TH1F("44","ds/dE for pi- at theta = 3",nbinspi,0.,emaxpi/GeV);
      if(nanglpi>4)
       h[44]=new TH1F("45","ds/dE for pi- at theta = 4",nbinspi,0.,emaxpi/GeV);
      if(nanglpi>0)
       h[45]=new TH1F("46","ds/dE for pi+ at theta = 0",nbinspi,0.,emaxpi/GeV);
      if(nanglpi>1)
       h[46]=new TH1F("47","ds/dE for pi+ at theta = 1",nbinspi,0.,emaxpi/GeV);
      if(nanglpi>2)
       h[47]=new TH1F("48","ds/dE for pi+ at theta = 2",nbinspi,0.,emaxpi/GeV);
      if(nanglpi>3)
       h[48]=new TH1F("49","ds/dE for pi+ at theta = 3",nbinspi,0.,emaxpi/GeV);
      if(nanglpi>4)
       h[49]=new TH1F("50","ds/dE for pi+ at theta = 4",nbinspi,0.,emaxpi/GeV);

      h[50]=new TH1F("51","E(GeV) neutrons",nbinlog,0.,logmax);
      if(nangl>0)
	h[51]=new TH1F("52","ds/dE for neutrons at theta = 0",nbinlog,0.,logmax);
      if(nangl>1)
	h[52]=new TH1F("53","ds/dE for neutrons at theta = 1",nbinlog,0.,logmax);
      if(nangl>2)
	h[53]=new TH1F("54","ds/dE for neutrons at theta = 2",nbinlog,0.,logmax);
      if(nangl>3)
	h[54]=new TH1F("55","ds/dE for neutrons at theta = 3",nbinlog,0.,logmax);
      if(nangl>4)
	h[55]=new TH1F("56","ds/dE for neutrons at theta = 4",nbinlog,0.,logmax);

      h[56]=new TH1F("57","Ekin (GeV) for 1st particle",120,0.,p_Kinetic_energy/GeV);
      h[57]=new TH1F("58","Ekin (GeV) for 2nd particle",120,0.,p_Kinetic_energy/GeV);
      h[58]=new TH1F("59","cos(Theta) for 1st particle in Lab.Sys.",nbinsa,-1.,1.);
      h[59]=new TH1F("60","cos(Theta) for 2nd particle in Lab.Sys.",nbinsa,-1.,1.);
      h[60]=new TH1F("61","cos(Theta) for 1st particle in CM.Sys.",nbinsa,-1.,1.);
      h[61]=new TH1F("62","cos(Theta) for 2nd particle in CM.Sys.",nbinsa,-1.,1.);
      h[62]=new TH1F("63","cos(Theta) for 1 & 2 particle in CM.Sys.",nbinsa,-1.,1.);

      G4cout << "Histograms is initialised nbins=" << nbins
             << G4endl;
	TH1* Mpim = new TH1F("100","pi-",30,-0.25,14.75);
	TH1* Mpip = new TH1F("101","pi+",30,-0.25,14.75);
	TH1* Mpi0 = new TH1F("102","pi0",30,-0.25,14.75);
	
	TH1* MKm  = new TH1F("110","K-",30,-0.25,14.75);
	TH1* MKp  = new TH1F("112","K+",30,-0.25,14.75);
	TH1* MK0s = new TH1F("113","K0s",30,-0.25,14.75);
	TH1* MK0l = new TH1F("114","K0l",30,-0.25,14.75);
	
	TH1* MLambda    = new TH1F("120","Lambda",30,-0.25,14.75);
	TH1* MLambdaBar = new TH1F("121","LambdaBar",30,-0.25,14.75);
	
	TH1* MSigma0    = new TH1F("130","Sigma0",30,-0.25,14.75);
	TH1* MSigma0Bar = new TH1F("131","Sigma0Bar",30,-0.25,14.75);
	
	TH1* Mproton    = new TH1F("140","proton",30,-0.25,14.75);
	TH1* MprotonBar = new TH1F("141","antiproton",30,-0.25,14.75);
	
	TH1* Mneutron    = new TH1F("142","neutron",30,-0.25,14.75);
	TH1* MneutronBar = new TH1F("143","antineutron",30,-0.25,14.75);	

// 	std::string description = "int PDG; float pt; float fx; float rapidity; float E; float Ekin; float cmsE; float p; float theta";
        TTree * Tinfo = new TTree("SEC","Secondary info");



        struct secondary_info {
	     int PDG;
	     float pt;
	     float fx;
	     float rapidity;
	     float E;
	     float Ekin;
	     float cmsE;
	     float p;
	     float theta;
	};
	secondary_info secs;     
	Tinfo->Branch("secondaries",&secs.PDG,"PDG/I:pt/F:fx/F:rapidity/F:E/F:Ekin/F:cmsE/F:p/F:theta/F");



	TH1* Info = new TH1D("Info","information",10,.5,10.5);

	G4double beta = pmax / ( Kinetic_energy + mass);
	G4double delta_eta = atanh ( beta );        
	G4cout << " offset for rapidity = " << delta_eta << G4endl;
        Info->Fill(1.,double(nevt));
//    //Info->Fill(2.,    will get cross section in use.
	Info->Fill(3.,delta_eta);
	Info->Fill(4., material->GetElement(0)->GetA());
//    //Info->Fill(5.,    will get default (LHEP) cross section.
//            



	G4int Spim = 0;
	G4int Spip = 0;
	G4int Spi0 = 0;

	G4int SKm = 0;
	G4int SKp = 0;
	G4int SK0s = 0;
	G4int SK0l = 0;

	G4int SLambda = 0;
	G4int SLambdaBar = 0;

	G4int SSigma0 = 0;
	G4int SSigma0Bar = 0;

	G4int Sproton = 0;
	G4int SprotonBar = 0;

	G4int Sneutron = 0;
	G4int SneutronBar = 0;
	
	G4int Seta = 0;
	G4int Seta_prime = 0;

    // Create a DynamicParticle
    G4DynamicParticle dParticle(part,aDirection,Kinetic_energy);
    
    G4double proj_momentum= sqrt(Kinetic_energy*(Kinetic_energy + 2*part->GetPDGMass()));



    G4VCrossSectionDataSet* cs = 0;
    G4double cross_sec = 0.0;

    if(nameGen == "LElastic" || nameGen == "LElasticB" || 
       nameGen == "HElastic" || nameGen == "BertiniElastic") {
      cs = new G4HadronElasticDataSet();
    } else if ( nameGen.substr(0,3) == "QGS" || 
	        nameGen.substr(0,3) == "FTF" ) {
      G4cout << G4endl << " Using cross sections from QGS physics list: ";
      if(part == proton && material->GetElement(0)->GetZ() > 1.5)
      {
          cs = new G4ProtonInelasticCrossSection();
	  G4cout << "G4ProtonInelasticCrossSection" << G4endl;
      }
      if(part == neutron && material->GetElement(0)->GetZ() > 1.5)
      {
          cs = new G4NeutronInelasticCrossSection();
	  G4cout << "G4NeutronInelasticCrossSection" << G4endl;
      }  
      if( (part == pip || part == pi0 || part== pim) && material->GetElement(0)->GetZ() > 1.5)
      {
          cs = new G4PiNuclearCrossSection();
	  G4cout << "G4PiNuclearCrossSection" << G4endl;
      }
	   
    } else if(part == proton && material->GetElement(0)->GetZ() > 1.5) {
      cs = new G4ProtonInelasticCrossSection();
    } else if(part == neutron && material->GetElement(0)->GetZ() > 1.5) {
      cs = new G4NeutronInelasticCrossSection();
    } else if( ionParticle ) {
      if ( Shen ) {
        cs = new G4IonsShenCrossSection();
	G4cout << "Using Shen Cross section for Ions" << G4endl;
      }
      if ( ! cs ) {
      cs = new G4TripathiCrossSection();
      G4cout << "Using Tripathi Cross section for Ions" << G4endl;
      }
    } else {
      G4cout << G4endl<< "Using default (LHEP) cross sections" << G4endl;
      cs = new G4HadronInelasticDataSet();
    }

    if(cs) {
      cs->BuildPhysicsTable(*part);
      cross_sec = cs->GetCrossSection(&dParticle, material->GetElement(0));
    } else {
      G4cout << G4endl << " No cross section set up -->> using default G4HadronCrossSections" << G4endl;
      cross_sec = (G4HadronCrossSections::Instance())->
        GetInelasticCrossSection(&dParticle, material->GetElement(0));
    }
    G4cout << "    cross(b)= " << cross_sec/barn << G4endl;
    Info->Fill(2.,cross_sec/barn);
    Info->Fill(5.,(G4HadronCrossSections::Instance())->GetInelasticCrossSection(&dParticle, material->GetElement(0))/barn);
    G4double factor = cross_sec*MeV*1000.0*(G4double)nbinse/(Kinetic_energy*barn*(G4double)nevt);
    G4double factora= cross_sec*MeV*1000.0*(G4double)nbinsa/(twopi*2.0*barn*(G4double)nevt);
    G4double factorb= cross_sec*1000.0/(barn*(G4double)nevt);
    G4cout << "### factor  = " << factor
           << "### factora = " << factora
           << "### factorb = " << factorb
           << "    cross(b)= " << cross_sec/barn << G4endl;

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
                 barn*emax*(G4double)nevt);
      }
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

    if(!G4StateManager::GetStateManager()->SetNewState(G4State_Idle))
      G4cout << "G4StateManager PROBLEM! " << G4endl;
    /*
    G4cout << "### proton:" << G4endl;
    proton->DumpTable();
    G4cout << "### piminus:" << G4endl;
    pim->DumpTable();
    */
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

    for (G4int iter=0; iter<nevt; iter++) {

      if(verbose>1) {
        G4cout << "### " << iter << "-th event start " << G4endl;
      }

      G4double e0 = Kinetic_energy;
      do {
        if(sigmae > 0.0) e0 = G4RandGauss::shoot(Kinetic_energy,sigmae);
      } while (e0 < 0.0);

      dParticle.SetKineticEnergy(e0);

      gTrack->SetStep(step);
      gTrack->SetKineticEnergy(Kinetic_energy);

      labv = G4LorentzVector(0.0, 0.0, pmax, Kinetic_energy + mass + nucleus_mass);
      G4double Ecms=labv.mag2();
      G4ThreeVector bst = labv.boostVector();

      G4LorentzVector labNN(0,0,pmax,Kinetic_energy + mass + G4Proton::Proton()->GetPDGMass());
      G4double NucleonNucleoncms=labNN.mag2();
      G4ThreeVector boostNN = labNN.boostVector();
      
      aChange = proc->PostStepDoIt(*gTrack,*step);

      G4double de = aChange->GetLocalEnergyDeposit();
      G4int n = aChange->GetNumberOfSecondaries();
      
      if(iter == 1000*(iter/1000)) {
        G4cerr << "##### " << iter << "-th event  #####" << G4endl;
      }

      G4int nbar = 0;
      G4int Npim = 0;
      G4int Npip = 0;
      G4int Npi0 = 0;

      G4int NKm = 0;
      G4int NKp = 0;
      G4int NK0s = 0;
      G4int NK0l = 0;

      G4int NLambda = 0;
      G4int NLambdaBar = 0;

      G4int NSigma0 = 0;
      G4int NSigma0Bar = 0;

      G4int Nproton = 0;
      G4int NprotonBar = 0;

      G4int Nneutron = 0;
      G4int NneutronBar = 0;
      
      G4int Neta = 0;
      G4int Neta_prime = 0;

      for(G4int j=0; j<n; j++) {

        sec = aChange->GetSecondary(j)->GetDynamicParticle();
        pd  = sec->GetDefinition();
        if(pd->GetPDGMass() > 100.*MeV) nbar++;
	G4String pname=pd->GetParticleName();
	if ( pname == "pi-" ) Npim++;
	else if ( pname == "pi+" ) Npip++;
	else if ( pname == "pi0" ) Npi0++;

	else if ( pname == "kaon-" ) NKm++;
	else if ( pname == "kaon+" ) NKp++;
	else if ( pname == "kaon0S" ) NK0s++;
	else if ( pname == "kaon0L" ) NK0l++;

	else if ( pname == "lambda" ) NLambda++;
	else if ( pname == "anti_lambda" ) NLambdaBar++;

	else if ( pname == "sigma0" ) NSigma0++;
	else if ( pname == "anti_sigma0" ) NSigma0Bar++;
	else if ( pname == "sigma+" ) ;
	else if ( pname == "sigma-" ) ;
	else if ( pname == "anti_sigma+" ) ;
	else if ( pname == "anti_sigma-" ) ;
	
	else if ( pname == "xi0" ) ;
	else if ( pname == "xi-" ) ;
	else if ( pname == "anti_xi0" ) ;
	else if ( pname == "anti_xi-" ) ;
	
	else if ( pname == "omega-" ) ;
	else if ( pname == "anti_omega-" ) ;

	else if ( pname == "proton" ) Nproton++;
	else if ( pname == "anti_proton" ) NprotonBar++;

	else if ( pname == "neutron" ) Nneutron++;
	else if ( pname == "anti_neutron" ) NneutronBar++;
	else if ( pname == "eta" ) Neta++;
	else if ( pname == "eta_prime" ) Neta_prime++;
	else if ( pname == "gamma" ) ;
	else if (pd->GetParticleType() == "nucleus" ) ;
	else {
	   G4cout << "****Found " << pname ;
	   if ( pd->IsShortLived() ) G4cout << "  is Shortlived" ;
	   G4cout << G4endl<< " .... width, Shortlived... " << pd->GetPDGWidth() 
	          << " " << pd->IsShortLived() << G4endl;
		pd->DumpTable();  
	}

      }
	   Mpim->Fill(Npim);
	   Mpip->Fill(Npip);
	   Mpi0->Fill(Npi0);

	   MKm->Fill(NKm);
	   MKp->Fill(NKp);
	   MK0s->Fill(NK0s);
	   MK0l->Fill(NK0l);

	   MLambda->Fill(NLambda);
	   MLambdaBar->Fill(NLambdaBar);

	   MSigma0->Fill(NSigma0);
	   MSigma0Bar->Fill(NSigma0Bar);

	   Mproton->Fill(Nproton);
	   MprotonBar->Fill(NprotonBar);

	   Mneutron->Fill(Nneutron);
	   MneutronBar->Fill(NneutronBar);

	   Spim += Npim;
	   Spip += Npip;
	   Spi0 += Npi0;

	   SKm += NKm;
	   SKp += NKp;
	   SK0s += NK0s;
	   SK0l += NK0l;

	   SLambda += NLambda;
	   SLambdaBar += NLambdaBar;

	   SSigma0 += NSigma0;
	   SSigma0Bar += NSigma0Bar;

	   Sproton += Nproton;
	   SprotonBar += NprotonBar;

	   Sneutron += Nneutron;
	   SneutronBar += NneutronBar;
	   
	   Seta += Neta;
	   Seta_prime += Neta_prime;

      for(G4int i=0; i<n; i++) {

        sec = aChange->GetSecondary(i)->GetDynamicParticle();
        pd  = sec->GetDefinition();
        mom = sec->GetMomentumDirection();
        e   = sec->GetKineticEnergy();

	if (e < 0.0) {
           e = 0.0;
	}

	// for exclusive reaction 2 particles in final state
        if(!inclusive && nbar != 2) break;

        m = pd->GetPDGMass();
	p = std::sqrt(e*(e + 2.0*m));   // e is T - kinetic energy
	mom *= p;
        fm = G4LorentzVector(mom, e + m);
        labv -= fm;
        mom = (*rot)*mom;
        px = mom.x();
        py = mom.y();
        pz = mom.z();
        pt = std::sqrt(px*px +py*py);

        theta = mom.theta();
        G4double cost  = std::cos(theta);
        G4double thetad = theta/degree;

        fm.boost(-bst);
        G4double costcm = std::cos(fm.theta());

	G4LorentzVector fmNN(mom, e + m);
	fmNN.boost(-boostNN);
	G4double charge=pd->GetPDGCharge();
	G4double feynmanX=2*fmNN.z()/std::sqrt(NucleonNucleoncms);
	G4double rapidity=Rapidity(pz,e+m);


	G4double ptSquare= sqr(pt)/sqr(GeV);
	G4double theta=mom.theta();

	secs.PDG = pd->GetPDGEncoding();
	secs.pt = float(pt/GeV);
	secs.fx = float(feynmanX);
	secs.rapidity = float(rapidity);
	secs.E = float((e + m )/GeV);
	secs.Ekin = float(e/GeV);
	secs.cmsE = float(fmNN.e()/GeV);
	secs.p    = float(p/GeV);
	secs.theta = float(theta);
	Tinfo->Fill();




	if(usepaw) {
	  if(i==0)  {
	    h[56]->Fill(e/GeV,1.0);
	    h[58]->Fill(cost,factora);
	    h[60]->Fill(costcm,factora);
	    h[62]->Fill(costcm,factora);
	  } else if(i==1) {
	    h[57]->Fill(e/GeV,1.0);
	    h[59]->Fill(cost,factora);
	    h[61]->Fill(costcm,factora);
	    h[62]->Fill(costcm,factora);
	  }
          h[2]->Fill(mom.phi()/degree,1.0);
          if(pd == neutron) h[23]->Fill(mom.phi()/degree,1.0);
	}

	if( e == 0.0 || pt == 0.0) {
          G4cout << "Warning! in event # " << iter 
	         << i << "-th secondary  "
		 << pd->GetParticleName() << "   Ekin(GeV)= "
                 << e/GeV
                 << " Pt(GeV/c)= " << pt/GeV
		 << G4endl;
	}
	de += e;
        if(verbose>0 || std::fabs(mom.phi()/degree - 90.) < 0.001) {
          G4cout << i << "-th secondary  "
		 << pd->GetParticleName() << "   Ekin(GeV)= "
                 << e/GeV
		 << "   p(GeV)= " << mom/GeV
		 << "   m(GeV)= " << m/GeV
		 << "   Etot(GeV)= " << (e+m)/GeV
		 << "   pt(GeV)= " << pt/GeV
                 << " has deg = " << mom.phi()/degree
                 << G4endl;
        }

	if(usepaw) {

          if(pd) {
            float N = pd->GetBaryonNumber();
            float Z = pd->GetPDGCharge()/eplus;
	     h[26]->Fill(N, factorb);
	  }

          if(pd == proton) {

            h[1]->Fill(1.0, 1.0);
            h[3]->Fill(pz/GeV, 1.0);
            h[7]->Fill(pt/GeV, 1.0);
            h[11]->Fill(e/GeV, 1.0);
	    h[21]->Fill(e/GeV, factor);
	    h[24]->Fill(cost, factora);

          } else if(pd == pim) {

	    h[1]->Fill(4.0, 1.0);
            h[4]->Fill(pz/GeV, 1.0);
            h[8]->Fill(pt/GeV, 1.0);
            h[12]->Fill(e/GeV, 1.0);
            for(G4int kk=0; kk<nanglpi; kk++) {
              if(bngpi1[kk] <= thetad && thetad <= bngpi2[kk]) {
                h[40+kk]->Fill(e/GeV, cngpi[kk]);
                break;
	      }
	    }

          } else if(pd == pip) {

	    h[1]->Fill(3.0, 1.0);
            h[5]->Fill(pz/GeV, 1.0);
            h[9]->Fill(pt/GeV, 1.0);
            h[13]->Fill(e/GeV, 1.0);
            for(G4int kk=0; kk<nanglpi; kk++) {
              if(bngpi1[kk] <= thetad && thetad <= bngpi2[kk]) {
                h[45+kk]->Fill(e/GeV, cngpi[kk]);
                break;
	      }
	    }

	  } else if(pd == pi0) {

	    h[1]->Fill(5.0, 1.0);
	    h[18]->Fill(e/GeV, 1.0);
	    h[19]->Fill(pz/GeV, 1.0);
	    h[20]->Fill(pt/GeV, 1.0);

	  } else if(pd == neutron) {

	    h[1]->Fill(2.0, 1.0);
            h[6]->Fill(pz/GeV, 1.0);
            h[10]->Fill(pt/GeV, 1.0);
            h[14]->Fill(e/GeV, 1.0);
	    h[22]->Fill(e/GeV, factor);
            G4double ee = std::log10(e/GeV);
            G4int    nb = (G4int)(ee/binlog);
            G4double e1 = binlog*nb;
            G4double e2 = e1 + binlog;
            e1 = std::pow(10., e1);
            e2 = std::pow(10., e2) - e1;
            G4double f  = factor*bine/e2;
	    h[50]->Fill(ee, f);
	    if(e >= elim) h[25]->Fill(cost, factora);
            for(G4int kk=0; kk<nangl; kk++) {
              if(bng1[kk] <= thetad && thetad <= bng2[kk]) {
                h[27+kk]->Fill(e/GeV, cng[kk]);
                if(kk < 5) h[51+kk]->Fill(ee, cng[kk]*bind/e2);
                break;
	      }
	    }

	  } else if(pd == deu) {
	    h[1]->Fill(6.0, 1.0);
	  } else if(pd == tri) {
	    h[1]->Fill(7.0, 1.0);
	  } else if(pd == alp) {
	    h[1]->Fill(8.0, 1.0);
	  } else {
	    h[1]->Fill(9.0, 1.0);
	  }
	}
	//	delete sec;       	 
        delete aChange->GetSecondary(i);
      }

      if(verbose > 0) {
        G4cout << "Energy/Momentum balance= " << labv << G4endl;
      }

      px = labv.px();
      py = labv.py();
      pz = labv.pz();
      p  = std::sqrt(px*px +py*py + pz*pz);
      pt = std::sqrt(px*px +py*py);

      if(usepaw) {
        h[0]->Fill((float)n,1.0);
	h[15]->Fill(labv.e()/GeV, 1.0);
	h[16]->Fill(pz/GeV, 1.0);
	h[17]->Fill(pt/GeV, 1.0);
      }
      aChange->Clear();

    }

    timer->Stop();
    G4cout << "  "  << *timer << G4endl;
    delete timer;

    // Committing the transaction with the tree
    if(usepaw) {
      std::cout << "Committing..." << std::endl;
      rootFile->Write();

    }
     
	std::cout << "  pim : " <<        G4double(Spim) /nevt 
			 << " (" << sqrt(G4double(Spim))/nevt << ")" << G4endl;
	std::cout << "  pip : " <<        G4double(Spip) /nevt 		
			 << " (" << sqrt(G4double(Spip))/nevt << ")" << G4endl;
	std::cout << "  pi0 : " <<        G4double(Spi0) /nevt 		
			 << " (" << sqrt(G4double(Spi0))/nevt << ")" << G4endl;

	std::cout << "  Kp : " <<  	     G4double(SKp) /nevt 		
			 << " (" << sqrt(G4double(SKp))/nevt << ")" << G4endl;
	std::cout << "  Km : " <<  	     G4double(SKm) /nevt 		
			 << " (" << sqrt(G4double(SKm))/nevt << ")" << G4endl;
	std::cout << "  K0s : " <<        G4double(SK0s) /nevt 		
			 << " (" << sqrt(G4double(SK0s))/nevt << ")" << G4endl;
	std::cout << "  K0l : " << 	     G4double(SK0l) /nevt 		
			 << " (" << sqrt(G4double(SK0l))/nevt << ")" << G4endl;
	
	std::cout << "  Lambda : " <<     G4double(SLambda) /nevt 	
			 << " (" << sqrt(G4double(SLambda))/nevt << ")" << G4endl;
	std::cout << "  LambdaBar : " <<  G4double(SLambdaBar) /nevt 	
			 << " (" << sqrt(G4double(SLambdaBar))/nevt << ")" << G4endl;
	
	std::cout << "  Sigma0 : " <<     G4double(SSigma0) /nevt 	
			 << " (" << sqrt(G4double(SSigma0))/nevt << ")" << G4endl;
	std::cout << "  Sigma0Bar : " <<  G4double(SSigma0Bar) /nevt 	
			 << " (" << sqrt(G4double(SSigma0Bar))/nevt << ")" << G4endl;

	std::cout << "  proton : " <<     G4double(Sproton) /nevt 			
			 << " (" << sqrt(G4double(Sproton))/nevt << ")" << G4endl;
	std::cout << "  protonBar : " <<  G4double(SprotonBar) /nevt 			
			 << " (" << sqrt(G4double(SprotonBar))/nevt << ")" << G4endl;

	std::cout << "  neutron : " <<    G4double(Sneutron) /nevt 			
			 << " (" << sqrt(G4double(Sneutron))/nevt << ")" << G4endl;
	std::cout << "  neutronBar : " << G4double(SneutronBar) /nevt 			
			 << " (" << sqrt(G4double(SneutronBar))/nevt << ")" << G4endl;

	std::cout << " "  << G4double(Spim) /nevt
	     << " "  << G4double(Spip) /nevt
	     << " "  << G4double(Spi0) /nevt

	     << " "  << G4double(SKp) /nevt
	     << " "  << G4double(SKm) /nevt
	     << " "  << G4double(SK0s) /nevt
	
	     << " "  << G4double(SLambda + SSigma0) /nevt
	     << " "  << G4double(SLambdaBar + SSigma0Bar) /nevt
	
	     << " "  <<  G4double(Sproton) /nevt
	     << " "  <<  G4double(SprotonBar) /nevt
	     << "  | "
	     << proj_momentum/GeV 
	     << G4endl;
    G4cerr << "###### End of run # " << run << "     ######" << G4endl;

  } while(end);

  delete mate;
  delete fin;
  delete phys;

  G4cout << "###### End of test #####" << G4endl;
}
