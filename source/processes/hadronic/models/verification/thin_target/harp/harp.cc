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
//      File name:     Harp
//
//      Author:        V.Ivanchenko
//
//      Creation date: 12 March 2002
//
//      Modifications:
//      14.11.03 Renamed to cascade
//
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
#include "G4ExcitationHandler.hh"
#include "G4PreCompoundModel.hh"
#include "G4Evaporation.hh"

#include <memory> // for the auto_ptr(T>
#include "AIDA/AIDA.h"

#include "G4Timer.hh"

int main(int argc, char** argv)
{
  G4cout << "========================================================" << G4endl;
  G4cout << "======             Cascade Test Start           ========" << G4endl;
  G4cout << "========================================================" << G4endl;

  G4String  namePart = "proton";
  G4String  nameMat  = "Si";
  G4String  nameGen  = "stringCHIPS";
  G4bool    logx     = false;
  G4bool    usepaw   = false;
  G4bool    inclusive= true;
  G4int     verbose  = 0;
  G4double  energy   = 100.*MeV;
  G4double  sigmae   = 0.0;
  G4double  elim     = 30.*MeV;
  G4double  dangl    = 5.0;
  G4int     nevt     = 1000;
  G4int     nbins    = 100;
  G4int     nbinsa   = 40;
  G4int     nbinse   = 80;
  G4int     nbinsd   = 20;
  G4int     nbinspi  = 20;
  G4int     nangl    = 0;
  G4int     nanglpi  = 0;
  G4String hFile     = "hbook.paw";
  G4double theStep   = 0.01*micrometer;
  G4double range     = 1.0*micrometer;
  G4double  emax     = 160.*MeV;
  G4double  emaxpi   = 200.*MeV;
  G4double ebinlog   = 2.0*MeV;
  G4double eBound    = 70.*MeV;
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
  G4cout << "#energy(MeV)" << G4endl;
  G4cout << "#sigmae(MeV)" << G4endl;
  G4cout << "#emax(MeV)" << G4endl;
  G4cout << "#emaxpi(MeV)" << G4endl;
  G4cout << "#elim(MeV)" << G4endl;
  G4cout << "#ebinlog(MeV)" << G4endl;
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



  G4String line, line1;
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
        emax    = energy;
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

    material = mate->GetMaterial(nameMat);
    if(!material) {
      G4cout << "Material <" << nameMat
	     << "> is not found out"
	     << G4endl;
	     exit(1);
    }

    G4ParticleDefinition* part = (G4ParticleTable::GetParticleTable())->FindParticle(namePart);

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
    /*
    if(nameGen == "binary") {
      phys->setCutOnP(eBound);
      phys->setCutOnPPP(kBound);
    }
    */
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

    G4int maxn = (G4int)((*(material->GetElementVector()))[0]->GetN()) + 1;
    // G4int maxz = (G4int)((*(material->GetElementVector()))[0]->GetZ()) + 1;

    G4cout << "The particle:  " << part->GetParticleName() << G4endl;
    G4cout << "The material:  " << material->GetName() << "  Amax= " << maxn << G4endl;
    G4cout << "The step:      " << theStep/mm << " mm" << G4endl;
    G4cout << "The position:  " << aPosition/mm << " mm" << G4endl;
    G4cout << "The direction: " << aDirection << G4endl;
    G4cout << "The time:      " << aTime/ns << " ns" << G4endl;


    // -------------------------------------------------------------------

    // Creating the analysis factory
    std::auto_ptr< AIDA::IAnalysisFactory > af( AIDA_createAnalysisFactory() );

    // Creating the tree factory
    std::auto_ptr< AIDA::ITreeFactory > tf( af->createTreeFactory() );

    // Creating a tree mapped to a new hbook file.
    std::auto_ptr< AIDA::ITree > tree( tf->create( hFile,  "hbook", false, true) );
    std::cout << "Tree store : " << tree->storeName() << std::endl;

    // Creating a tuple factory, whose tuples will be handled by the tree
    //   std::auto_ptr< AIDA::ITupleFactory > tpf( af->createTupleFactory( *tree ) );

    const G4int nhisto = 10;
    AIDA::IHistogram1D* h[nhisto];
    //    AIDA::IHistogram2D* h2;
    //AIDA::ITuple* ntuple1 = 0;


    G4double mass = part->GetPDGMass();
    G4double pmax = sqrt(energy*(energy + 2.0*mass));
    G4double binlog = log10(ebinlog);
    G4int nbinlog = (G4int)(log10(2.0*emax)/binlog);
    G4double logmax = binlog*nbinlog;
    G4double bine = emax/(G4double)nbinse;
    G4double bind = emax/(G4double)nbinsd;

    double m_pmin = 0.0;
    double m_pmax = emax/GeV;
    double m_ptmax = 0.65;
    double m_pth = 0.2;
    int m_binp = 65;
    int m_bint = 20;
    double m_thetamax = 300.;
    double m_thetamin = 15.;
    double cosmin = cos(m_thetamax*0.001);
    double cosmax = cos(m_thetamin*0.001);

    if(usepaw) {

      // Creating a histogram factory, whose histograms will be handled by the tree
      std::auto_ptr< AIDA::IHistogramFactory > hf( af->createHistogramFactory( *tree ) );

      // Creating an 1-dimensional histogram in the root directory of the tree

      // ---- Book a histogram and ntuples
      G4cout << "Hbook file name: <" << hFile << ">" << G4endl;
      G4cout << "energy = " << energy/GeV << " GeV" << G4endl;
      G4cout << "emax   = " << emax/GeV << " GeV" << G4endl;
      G4cout << "pmax   = " << m_pmax << " GeV" << G4endl;

      h[0]=hf->createHistogram1D("10","Number of protons",10,-0.5,9.5);
      h[1]=hf->createHistogram1D("11","Number of pions",10,-0.5,9.5);
      h[2]=hf->createHistogram1D("12","Proton momentum",m_binp,m_pmin,m_pmax);
      h[3]=hf->createHistogram1D("13","Pion momentum",m_binp,m_pmin,m_pmax);
      h[4]=hf->createHistogram1D("14","Proton Pt",m_binp,0.0,m_ptmax);
      h[5]=hf->createHistogram1D("15","Pion Pt",m_binp,0.0,m_ptmax);
      h[6]=hf->createHistogram1D("16","Proton theta",m_bint,0.0,m_thetamax);
      h[7]=hf->createHistogram1D("17","Pion theta",m_bint,0.0,m_thetamax);
      h[8]=hf->createHistogram1D("18","Proton cos(theta)",m_bint,cosmin,1.0);
      h[9]=hf->createHistogram1D("19","Proton cos(theta)",m_bint,cosmin,1.0);

      G4cout << "Histograms is initialised nbins=" << nbins
             << G4endl;
    }
    // Create a DynamicParticle

    G4DynamicParticle dParticle(part,aDirection,energy);
    G4VCrossSectionDataSet* cs = 0;
    G4double cross_sec = 0.0;

    if(part == proton && material->GetElement(0)->GetZ() > 1.5) {
      cs = new G4ProtonInelasticCrossSection();
    } else if(part == neutron && material->GetElement(0)->GetZ() > 1.5) {
      cs = new G4NeutronInelasticCrossSection();
    } else {
      cs = new G4HadronInelasticDataSet();
    }

    if(cs) {
      cs->BuildPhysicsTable(*part);
      cross_sec = cs->GetCrossSection(&dParticle, material->GetElement(0));
    } else {
      cross_sec = (G4HadronCrossSections::Instance())->
        GetInelasticCrossSection(&dParticle, material->GetElement(0));
    }

    G4double factor = cross_sec*MeV*1000.0*(G4double)nbinse/(energy*barn*(G4double)nevt);
    G4double factora= cross_sec*MeV*1000.0*(G4double)nbinsa/(twopi*2.0*barn*(G4double)nevt);
    G4double factorb= cross_sec*1000.0/(barn*(G4double)nevt);
    G4cout << "### factor  = " << factor
           << "### factora = " << factor
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
         (twopi*(cos(degree*bng1[k]) - cos(degree*bng2[k]))*
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
         (twopi*(cos(degree*bngpi1[k]) - cos(degree*bngpi2[k]))*
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

    G4RotationMatrix* rot  = new G4RotationMatrix();
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

      G4double e0 = energy;
      do {
        if(sigmae > 0.0) e0 = G4RandGauss::shoot(energy,sigmae);
      } while (e0 < 0.0);

      dParticle.SetKineticEnergy(e0);

      gTrack->SetStep(step);
      gTrack->SetKineticEnergy(energy);

      labv = G4LorentzVector(0.0, 0.0, pmax, energy + mass + amass);
      aChange = proc->PostStepDoIt(*gTrack,*step);

      G4double de = aChange->GetLocalEnergyDeposit();
      G4int n = aChange->GetNumberOfSecondaries();

      if(iter == 1000*(iter/1000)) {
        G4cerr << "##### " << iter << "-th event  #####" << G4endl;
      }

      G4int nbar = 0;
      int n_pr = 0;
      int n_pi = 0;
      

      for(G4int j=0; j<n; j++) {

        sec = aChange->GetSecondary(j)->GetDynamicParticle();
        pd  = sec->GetDefinition();
        if(pd->GetPDGMass() > 100.*MeV) nbar++;
      }

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
        double cost  = cos(theta);
        double sint  = sin(theta);
        G4double thetad = theta/degree;
        G4double thetamr = theta*1000.;

        if(thetamr>= m_thetamin && thetamr <= m_thetamax && p>m_pth) {
          std::string nam = pd->GetParticleName();
          if(nam == "proton") {
            h[2]->fill(float(p/GeV),1.0);
            h[4]->fill(float(pt/GeV),1.0);
//            h2[0]->fill(p,cost);
            if(p>m_pth) {
              h[6]->fill(float(thetamr),1.0);
              h[8]->fill(float(cost),1.0);
              n_pr++;
            }
          } else if(nam == "pi+" || nam == "kaon+" || nam == "pi-" || nam == "kaon-") {
            h[3]->fill(float(p/GeV),1.0);
            h[5]->fill(float(pt/GeV),1.0);
//            h2[1]->fill(p,cost);
            if(p>m_pth) {
              h[7]->fill(float(thetamr),1.0);
              h[9]->fill(float(cost),1.0);
              n_pi++;
            }
          }
          h[0]->fill(float(n_pr),1.0);
          h[1]->fill(float(n_pi),1.0);
	}
      
	//	delete sec;       	 
        delete aChange->GetSecondary(i);
      }
    
      aChange->Clear();
    }
  
    timer->Stop();
    G4cout << "  "  << *timer << G4endl;
    delete timer;

    // Committing the transaction with the tree
    if(usepaw) {
      std::cout << "Committing..." << std::endl;
      tree->commit();
      std::cout << "Closing the tree..." << std::endl;
      tree->close();
    }

    G4cerr << "###### End of run # " << run << "     ######" << G4endl;
  }  
    //  } while(end);

  delete mate;
  delete fin;
  delete phys;

  G4cout << "###### End of test #####" << G4endl;
}
