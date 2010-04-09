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
//      File name:     Test19 (Comparison of G4QString Results with Data)
//
//      Author:        M.Kossov (made of test30 created by V.Ivanchenko)
// 
//      Creation date: 27 Jan 2004
//
//      Modifications: G4electroweak CHIPS interface (G4QCaptureAtRest)
//                     which is independent of the "hadronic" package
//                     classes is added to this test for the debugging,
//                     physics tests and development of the CHIPS model.
//
// -------------------------------------------------------------------
//       1         2         3         4         5         6         7         8         9
//34567890123456789012345678901234567890123456789012345678901234567890123456789012345678901

//#define pdebug
#define nout
//#define inter
//#define pscan
//#define csdebug
//#define masstest
//#define smear
//#define escan
//#define idebug
//#define tdebug
//#define pidebug
//#define debug
//#define ppdebug
//#define hdebug
//#define lhepdbg
//#define meandbg
//#define ekindbg
//#define histdbg
//--- Random seed
//#define ranseed
//--- Flags of models (only one must be chosen), CHIPS is a default for System Testing ---
#define chips
//#define synch
//#define lep
//#define hep
//#define preco
//#define berti
//#define binar
//#define qgsc
// ------------------------------------- FLAGS ------------------
#include "G4UIterminal.hh"
#include "globals.hh"
#include "G4ios.hh"
#include <iostream>
#include <fstream>
#include <iomanip>

// Fake classes for this test
#include "G4RunManager.hh"
#include "G4VUserPhysicsList.hh"
//#include "G4VUserDetectorConstruction.hh"
#include "G4Material.hh"
#include "G4ElementVector.hh"
#include "G4VContinuousDiscreteProcess.hh"
#include "G4ProcessManager.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleChange.hh"
#include "G4GenericIon.hh"

#include "G4QInelastic.hh"
#include "G4QElastic.hh"
#include "G4QSynchRad.hh"
#include "G4QDiffraction.hh"
#include "G4QLowEnergy.hh"
#include "G4TheoFSGenerator.hh"
#include "G4StringChipsParticleLevelInterface.hh"
#include "G4QGSModel.hh"
#include "G4QGSParticipants.hh"
#include "G4QGSMFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4QuasiElasticChannel.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4AntiProtonInelasticProcess.hh"
#include "G4ProtonInelasticCrossSection.hh"
#include "G4HadronInelasticDataSet.hh"
#include "G4HadronCaptureDataSet.hh"
#include "G4HadronFissionDataSet.hh"
#include "G4PionPlusInelasticProcess.hh"
#include "G4BGGNucleonElasticXS.hh"
#include "G4BGGPionElasticXS.hh"
#include "G4NucleonNuclearCrossSection.hh"
#include "G4GlauberGribovCrossSection.hh"
#include "G4PiNuclearCrossSection.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4KaonPlusInelasticProcess.hh"
#include "G4KaonMinusInelasticProcess.hh"
// CHIPS Inelastic Cross-Sections
#include "G4QMuonNuclearCrossSection.hh"
#include "G4QProtonNuclearCrossSection.hh"
#include "G4QNeutronNuclearCrossSection.hh"
#include "G4QNeutronCaptureRatio.hh"
#include "G4QPionMinusNuclearCrossSection.hh"
#include "G4QPionPlusNuclearCrossSection.hh"
#include "G4QKaonPlusNuclearCrossSection.hh"
#include "G4QKaonMinusNuclearCrossSection.hh"
#include "G4QKaonZeroNuclearCrossSection.hh"
#include "G4QHyperonNuclearCrossSection.hh"
#include "G4QHyperonPlusNuclearCrossSection.hh"
#include "G4QAntiBaryonPlusNuclearCrossSection.hh"
#include "G4QAntiBaryonNuclearCrossSection.hh"
// CHIPS ElasticCross-Sections
#include "G4QProtonElasticCrossSection.hh"
#include "G4QNeutronElasticCrossSection.hh"
#include "G4QPionMinusElasticCrossSection.hh"
#include "G4QPionPlusElasticCrossSection.hh"
#include "G4QKaonPlusElasticCrossSection.hh"
#include "G4QKaonMinusElasticCrossSection.hh"
//#include "G4QKaonZeroElasticCrossSection.hh"
#include "G4QHyperonElasticCrossSection.hh"
#include "G4QHyperonPlusElasticCrossSection.hh"
//#include "G4QAntiBaryonPlusElasticCrossSection.hh"
#include "G4QAntiBaryonElasticCrossSection.hh"
// GHAD Cross-Sections
#include "G4PiNuclearCrossSection.hh"
#include "G4CrossSectionDataStore.hh"
#include "G4HadronElasticDataSet.hh"
#include "G4HadronInelasticDataSet.hh"
#include "G4NeutronHPElasticData.hh"
#include "G4NeutronHPCaptureData.hh"
#include "G4NeutronHPFissionData.hh"

#include "G4ApplicationState.hh"
#include "G4StateManager.hh"

#include "G4UnitsTable.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleChange.hh"
#include "G4DynamicParticle.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4AntiProton.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4MuonMinus.hh"
#include "G4Gamma.hh"
#include "G4PionZero.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4Lambda.hh"
#include "G4SigmaMinus.hh"
#include "G4SigmaZero.hh"
#include "G4SigmaPlus.hh"
#include "G4Alpha.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4ForceCondition.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4Step.hh"
#include "G4GRSVolume.hh"

#include "G4UImanager.hh" 

//#include "Test19Physics.hh"
#include "Test19PhysicsList.hh"
#include "Test19MagneticField.hh"

//#include "G4HadronCrossSections.hh"
//#include "G4VCrossSectionDataSet.hh"
//#include "G4ProtonInelasticCrossSection.hh"
//#include "G4NeutronInelasticCrossSection.hh"
//#include "G4HadronInelasticDataSet.hh"

#include "G4SynchrotronRadiation.hh"
#include "G4SynchrotronRadiationInMat.hh"

// *** M o d e l s ****
#include "G4LEProtonInelastic.hh"
#include "G4HEProtonInelastic.hh"
#include "G4LENeutronInelastic.hh"
#include "G4HENeutronInelastic.hh"
#include "G4LEPionMinusInelastic.hh"
#include "G4HEPionMinusInelastic.hh"
#include "G4LEPionPlusInelastic.hh"
#include "G4HEPionPlusInelastic.hh"
#include "G4LEKaonPlusInelastic.hh"
#include "G4HEKaonPlusInelastic.hh"
#include "G4LEKaonMinusInelastic.hh"
#include "G4HEKaonMinusInelastic.hh"
#include "G4LEAntiProtonInelastic.hh"
#include "G4HEAntiProtonInelastic.hh"
#include "G4BinaryCascade.hh"
#include "G4CascadeInterface.hh"
#include "G4ExcitationHandler.hh"
#include "G4PreCompoundModel.hh"
#include "G4Evaporation.hh"

// New Histogramming (from AIDA and Anaphe):
//#include <memory> // for the auto_ptr(T>

// Temporary AIDA can be not used in test19
//#include "AIDA/AIDA.h"

// ASCII Pseudo NTUPLE 
#ifndef nout
#include "G4QHBook.hh"
#endif

#include "G4Timer.hh"
#include "time.h"

//int main(int argc, char** argv)
int main()
{
  static const G4double mNuc=939*MeV;  // Nucleon mass for the visible energy estimate
#ifdef csdebug
  //GHAD temperature//
  static const G4double T0=273*kelvin; // Room temperature for the XS tests
#endif
#ifdef ranseed
  //G4UImanager* UI = G4UImanager::GetUIpointer();
  //UI->ApplyCommand( "" );
  CLHEP::RanluxEngine defaultEngine( 1234567, 4 ); 
  CLHEP::HepRandom::setTheEngine( &defaultEngine ); 
  G4int seed = time( NULL ); 
  CLHEP::HepRandom::setTheSeed( seed ); 
#endif
#ifdef masstest
  G4int tstPDG=91000000;
  G4cout<<"Test19: testPDG="<<tstPDG<<", M = "<<G4QPDGCode(tstPDG).GetMass()<<G4endl;
  return EXIT_SUCCESS;
#endif
#ifdef meandbg
  G4double pE=0.;
  G4double sE=0.;
  G4double sE2=0.;
  G4double nE=0.;
  G4double nD=0.;
#endif
#ifdef ekindbg
  G4double pE=0.;
  G4double sE=0.;
  G4double sE2=0.;
  G4double nE=0.;
  G4double nD=0.;
#endif
#ifdef histdbg
  const G4int nHst=40;
  G4int nHst1=nHst-1;
  G4int dQnM=nHst/2;
  G4int dQnM1=dQnM-1;
  G4int dQnE=nHst/2;
  //G4int dQnE1=dQnE-1;
  G4int rEnH[nHst];
  G4int rPzH[nHst];
  G4int rPrH[nHst];
  G4int dEnH[nHst];
  G4int dPzH[nHst];
  G4int dPrH[nHst];
  G4int dChH[nHst];
  G4int dBnH[nHst];
  G4int nGaH[nHst];
  G4int GaSE[nHst];
  for(G4int ih=0; ih<nHst; ih++)
  {
    rEnH[ih]=0; rPzH[ih]=0; rPrH[ih]=0; dEnH[ih]=0; dPzH[ih]=0; dPrH[ih]=0;
    dChH[ih]=0; dBnH[ih]=0; nGaH[ih]=0; GaSE[ih]=0;
  }
#endif
#ifdef hdebug
  G4StringChipsParticleLevelInterface::SetMaxB(20.); // Impact parameter limit
  G4StringChipsParticleLevelInterface::SetMaxE(20.); // Energy deposition limit
  G4StringChipsParticleLevelInterface::Reset();   // Initialize historgamming in the class
#endif
  //const G4int nTg=4;   // Length of the target list for the Performance test
  //G4int tli[nTg]={90001000,90002002,90007007,90027032}; // PDG Codes of targets
  //G4String tnm[nTg]={"Hydrogen","Helium","Nitrogen","Cobalt"}; // Target names
  //G4String tsy[nTg]={"1H","He","14N","59Co"}; // Target symbols for the Target Loop
  //G4Material* mat[nTg]={0,0,0,0}; // Material pointers for the Target Loop
  const G4int nTg=5;   // Length of the target list for the Performance test
  G4int tli[nTg]={90001000,90002002,90007007,90027032,90092146}; // PDG Codes of targets
  G4String tnm[nTg]={"Hydrogen","Helium","Nitrogen","Cobalt","Uranium"}; // Target names
  G4String tsy[nTg]={"1H","He","14N","59Co","238U"};// Target symbols for the Target Loop
  G4Material* mat[nTg]={0,0,0,0,0}; // Material pointers for the Target Loop
  const G4int nPr=12;          // Length of the projectile list for the Performance test
  G4int pli[nPr] = {2212, 2112, 211, -211, 3122, 321, -321, -2212, -2112, 11, 13, 22};//Pro
  const G4int nEn=4;           // Length of the kin. energy list for the Performance test
  G4double eli[nEn] = {27., 227., 999., 9999.}; // Kinetic energy values
  // ^^^ End of the Performance On Flight test definition for targets/projectiles/energies
#ifdef tdebug
  const G4int nT=20;           // Dimension of the t-distribution vectors
  //const G4int nT1=nT-1;        // The last bin of the t-distribution vector
  const G4double maxT=200000.; // -t_max for the T-vectors
  const G4double dT=maxT/nT;   // Step of the t-hystogram
  const G4double fT=dT/2;      // Center of the first bin of the t-histogram
  G4double tVal[nT];           // -t values for centers of bins of the t-hystogram
  G4int tSig[nT];              // t-histogram (filled and reset many times
#endif
  // Run manager
  G4RunManager* runManager = new G4RunManager;
  runManager->SetUserInitialization(new Test19PhysicsList);
  Test19MagneticField* fpMagField = new Test19MagneticField();
  G4double fieldValue=10.*tesla/3.;  // For all tests the B-field=(3.333*tesla, 0., 0.) 
  fpMagField->SetFieldValue(fieldValue);
  G4StateManager::GetStateManager()->SetNewState(G4State_Init); // To let create ions
#ifdef debug
  G4cout<<"Test19: Prepare G4QHBook files or ntuples"<<G4endl;
#endif
#ifndef nout
  G4QHBook* ntp = new G4QHBook;
#endif
  //-------- set standard output format-------
  G4cout.setf( std::ios::scientific, std::ios::floatfield );
  //-------- take parameters from a file -------
  std::ifstream inFile("chipstest.in", std::ios::in);
  G4double temperature;
  G4double ssin2g;
  G4double eteps;
  G4double momb;
  G4double enb;
  G4double cP;
  G4double fN;
  G4double fD;
  G4double rM;
  G4double sA;
  G4int    nop;
  G4int    pPDG;
  G4int    tPDG;
  G4int    nEvt;
  //G4int    nofdecays;
  //G4int    decmask=0;
  inFile>>temperature>>ssin2g>>eteps>>nop>>momb>>enb>>pPDG>>tPDG>>nEvt>>fN>>fD>>cP>>rM>>sA;
  G4cout<<"Test19:Par's: T="<<temperature<<",S="<<ssin2g<<",Eta="<<eteps<<",nop="
        <<nop<<",p="<<momb<<",e="<<enb<<",pPDG="<<pPDG<<",tPDG="<<tPDG<<",nEv="
        <<nEvt<<",fN="<<fN<<",fD="<<fD<<",cP="<<cP<<",rM="<<rM<<",sA="<<sA<<G4endl;
  //-------- Initialize CHIPS
  G4QCHIPSWorld* theW=G4QCHIPSWorld::Get();
  theW->GetParticles(nop);           // Create CHIPS World of nop particles
  //G4Exception("***CHIPStest: TMP");
  G4QInelastic::SetParameters(temperature,ssin2g,eteps,fN,fD,cP,rM,nop,sA);
  G4QInelastic::SetManual();
  // ********** Now momb is a momentum of the incident particle, if =0 => LOOP ************
//#ifdef chips
//#else
  // Not necessary for CHIPS, only for GHAD: fake force condition
  G4ForceCondition* cond = new G4ForceCondition;
  *cond=NotForced;
//#endif
  G4double mp=G4QPDGCode(pPDG).GetMass();
  G4double ep=mp;
  G4int cnE=1;
  if(momb==0.) cnE=nEn;
  else
  {
    ep=std::sqrt(mp*mp+momb*momb);
    if(enb>0.) ep=enb;
  }
  //G4int tPDG=90000000+tgZ*1000+tgN;
  G4int    tgZ=(tPDG-90000000)/1000;
  G4int    tgN=tPDG-90000000-tgZ*1000;
#ifdef debug
  G4cout<<"Test19: projPDG="<<pPDG<<",tgZ ="<<tgZ<<",tgN = "<<tgN<<",tgPDG="<<tPDG<<G4endl;
#endif
  // ---------- Define material for the simulation ------------------
  G4int tgA        = tgZ+tgN; // Mass number - fake
  // The material can be copied from the commented cMaterial Factory above
  G4Isotope* isotope=0;
  G4Element* element=0;
  G4Material* material=0;
  // LOOPs for the wide performance test
  G4double      aTime      = 0. ;
  //G4ThreeVector aDirection = G4ThreeVector(0.,0.,1.);
  G4ThreeVector aDirection = G4ThreeVector(0.,1.,0.); // @@ Temporary "not along Z"
  G4int tgm=1;                                        // By default only one target
  if(!tPDG) // Make max for the LOOP over all targets and define materials
  {
#ifdef debug
    G4cout<<"Test19: targetPDG="<<tPDG<<", tgm = "<<tgm<<", nTg = "<<nTg<<G4endl;
#endif
    tgm=nTg;
    for(G4int tgi=0; tgi<tgm; tgi++)
    {
      tPDG=tli[tgi];
      tgZ = (tPDG-90000000)/1000;
      tgN = tPDG-90000000-tgZ*1000;
      tgA = tgZ+tgN; // Baryon number
      // The material can be copied from the commented cMaterial Factory above
      isotope = new G4Isotope(tsy[tgi], tgZ, tgA);
      element = new G4Element(tnm[tgi], tsy[tgi], 1);
      element->AddIsotope(isotope, 100.*perCent);
      material = new G4Material(tnm[tgi], 1.*g/cm3, 1);
      material->AddElement(element, 1);
      G4cout<<"Test19:-->Material("<<tgZ<<","<<tgN<<"):"<<tnm[tgi]<<","<<tsy[tgi]<<G4endl;
      mat[tgi]=material;
    }
  }
  else
  {
    material = new G4Material("ZA_Isomer", 1.*g/cm3, 1);
    if(tgZ < 93 && tgN < 151)
    {
#ifdef debug
    G4cout<<"Test19: Before isotope tgZ="<<tgZ<<", tgN = "<<tgN<<G4endl;
#endif
      element = new G4Element("ZA_Isotop", "ZA", 1);
      isotope = new G4Isotope("Isotop", tgZ, tgA);
      element->AddIsotope(isotope, 100.*perCent);
      material->AddElement(element, 1);
    }
  }
#ifdef debug
  G4cout<<"Test19:--- Material "<<material->GetName()<<" is defined ---" << G4endl;
#endif
  if(!material)
  {
    G4cout<<"Test19: Last Material "<<material->GetName()<<" is not defined."<<G4endl;
    exit(1);
  }
  // ---------------------- PSEUDO-TRACK DEFINISIONS -------------------------
  G4double nx = 0.;
  G4double ny = 0.;
  G4double nz = 0.;
  G4int npart=1;                                       // By default only one particle
  if(!pPDG) npart=nPr;                                 // Make a LOOP ove all particles
#ifdef chips
  // *************** CHIPS process definition starts here *******************************
  G4QInelastic* proc = new G4QInelastic;               // A general CHIPS process
  ///G4QDiffraction* proc = new G4QDiffraction;           // A diffraction CHIPS process
  ///G4QLowEnergy* proc = new G4QLowEnergy;               // fragment-nucleus universal
#endif
  // **************** GHAD process definition starts here *******************************
#ifdef preco
  G4PreCompoundModel* aModel = new G4PreCompoundModel(new G4ExcitationHandler); // Preco
#endif
#ifdef berti
  G4CascadeInterface* aModel = new G4CascadeInterface; // Bertini
#endif
#ifdef binar
  G4BinaryCascade* aModel = new G4BinaryCascade;       // Binary
#endif
#ifdef lep
  G4HadronicInteraction* aModel;                       // LEP generators
  if     (pPDG==2212) aModel = new G4LEProtonInelastic;
  else if(pPDG==2112) aModel = new G4LENeutronInelastic;
  else if(pPDG==-211) aModel = new G4LEPionMinusInelastic;
  else if(pPDG== 211) aModel = new G4LEPionPlusInelastic;
  else if(pPDG==-321) aModel = new G4LEKaonMinusInelastic;
  else if(pPDG== 321) aModel = new G4LEKaonPlusInelastic;
  else if(pPDG==-2112)aModel = new G4LEAntiProtonInelastic;
  else
  {
    G4cout<<"-Error-Test19: Process is not defined in LEP for PDG="<<pPDG<<G4endl;
    aModel = new G4LEProtonInelastic;
  }
#endif
#ifdef hep
  G4HadronicInteraction* aModel;                       // HEP generators
  if     (pPDG==2212) aModel = new G4HEProtonInelastic;
  else if(pPDG==2112) aModel = new G4HEProtonInelastic;
  else if(pPDG==-211) aModel = new G4HEPionMinusInelastic;
  else if(pPDG== 211) aModel = new G4HEPionPlusInelastic;
  else if(pPDG==-321) aModel = new G4HEKaonPlusInelastic;
  else if(pPDG== 321) aModel = new G4HEKaonPlusInelastic;
  else if(pPDG==-2112)aModel = new G4HEAntiProtonInelastic;
  else
  {
    aModel = new G4HEProtonInelastic;
    G4cout<<"-Error-Test19: Process is not defined in HEP for PDG="<<pPDG<<G4endl;
  }
#endif
#ifdef qgsc
  G4TheoFSGenerator* aModel = new G4TheoFSGenerator;           // The same for QGS & FTF
  G4StringChipsParticleLevelInterface* theCHIPS=new G4StringChipsParticleLevelInterface;
  //G4cout<<"Tst19:*> Nuclear fragmentation model is defined"<<G4endl;
  //// ------------- Defines a Kind of nuclear fragmentation model--------
  aModel->SetTransport(theCHIPS);
  G4QGSModel<G4QGSParticipants>* aStringModel = new G4QGSModel<G4QGSParticipants>;
  //G4cout<<"Tst19:*> Intranuclear transport model is defined"<<G4endl;
  //// ----------- Defines a Kind of the QGS model -------------
  G4QGSMFragmentation aFragmentation;       // @@ Can be a general solution (move up)
  G4ExcitedStringDecay* aStringDecay = new G4ExcitedStringDecay(&aFragmentation);
  aStringModel->SetFragmentationModel(aStringDecay);
  aModel->SetHighEnergyGenerator(aStringModel);
  G4QuasiElasticChannel* theQuasiElastic = new G4QuasiElasticChannel;
  aModel->SetQuasiElasticChannel(theQuasiElastic);
  //G4cout<<"Tst19:*> String model is defined"<<G4endl;
  //// ----------- Defines energy limits of the model ----------
  ///aModel->SetMinEnergy(8*GeV);                // Do we need this ?
  ///aModel->SetMaxEnergy(100*TeV);              // Do we need that ?
#endif
#ifdef chips
  if(!proc)
  {
    G4cout<<"Tst19: there is no G4QInelastic process"<<G4endl;
    exit(1);
  }
  // Comment for G4QLowEnergies or G4QDiffraction (not G4QInelastic)
  //proc->SetParameters(temperature, ssin2g, eteps, fN, fD, cP, rM, nop, sA);
#else
#ifdef synch
  G4double s_n=0.;
  G4double s_s=0.;
  G4double s_d=0.;
  const G4int nBin=50;
  G4double minom=.000000001;
  G4double logmin=std::log(minom);
  G4double maxom=10.;
  G4double logmax=std::log(maxom);
  G4double dlog=(logmax-logmin)/nBin;
  G4double eBin[nBin];
  G4double hBin[nBin];
  for(G4int i=0; i < nBin; ++i)
  {
    eBin[i]=0.;
    hBin[i]=0.;
  }
  //
  //G4QSynchRad* proc = new G4QSynchRad;           // CHIPS Synchrotron Radiation process
  G4SynchrotronRadiation* proc = new G4SynchrotronRadiation; // GHadSynchratRad process 1
  //G4SynchrotronRadiationInMat* proc = new G4SynchrotronRadiationInMat; // GHadSR process2
#else
  G4HadronInelasticProcess* proc = 0;
  if     (pPDG==2212) proc = new G4ProtonInelasticProcess;
  else if(pPDG==2112) proc = new G4ProtonInelasticProcess;
  else if(pPDG==-211) proc = new G4PionMinusInelasticProcess;
  else if(pPDG== 211) proc = new G4PionPlusInelasticProcess;
  else if(pPDG==-321) proc = new G4KaonMinusInelasticProcess;
  else if(pPDG== 321) proc = new G4KaonPlusInelasticProcess;
  else if(pPDG==-2112) proc = new G4AntiProtonInelasticProcess;
  else G4cout<<"-Error-Test19: Process is not defined for PDG="<<pPDG<<G4endl;
  proc->RegisterMe(aModel); // from G4HadronicProcess
  G4VCrossSectionDataSet* theCS;
  if(pPDG==2212 || pPDG==2112) theCS = new G4ProtonInelasticCrossSection;
  else if(pPDG==211||pPDG==-211) theCS = new G4PiNuclearCrossSection;
  else theCS = new G4HadronInelasticDataSet; // For kaons anti-protons and others...
  proc->AddDataSet(theCS);   // Can not be skipped for the event generator
#endif
#endif
#ifdef debug
  G4cout<<"Test19:--***-- process is created --***--" << G4endl; // only one run
#endif
#ifdef hdebug
  G4Timer* timer = new G4Timer();
  timer->Start();
#endif
  G4int nTot=npart*tgm*cnE;
  G4int nCur=0;
  for(G4int pnb=0; pnb<npart; pnb++)           // LOOP over particles
  {
   if (npart>1) pPDG=pli[pnb];
   G4QContent pQC=G4QPDGCode(pPDG).GetQuarkContent();
   G4int        cp = pQC.GetCharge();          // Charge of the projectile
   G4int        sp = pQC.GetStrangeness();     // Strangeness of the projectile
   if(pPDG==22 || pPDG==12 || pPDG==14 || pPDG== 16 || pPDG==-11 || pPDG==-13 || pPDG==-15)
   {                                           // QCont of gamma&neutrinos is not supported
     cp=0;
     sp=0;
   }
   if(pPDG==11 || pPDG==13 || pPDG== 15)       // QuarkContent of leptons is not supported
   {
     cp=-1;
     sp=0;
   }
   if(pPDG==-11|| pPDG==-13|| pPDG==-15)       // QuarkContent of antilept isn't supported
   {
     cp=1;
     sp=0;
   }
   G4ParticleDefinition* part=G4AntiSigmaPlus::AntiSigmaPlus(); // DefaultProj definition
   if(pPDG==2212) part=G4Proton::Proton();           // Definition of Proton projectile
   else if(pPDG==22) part=G4Gamma::Gamma();          // Definition of gamma projectile
   else if(pPDG==11) part=G4Electron::Electron();    // Definition of electron projectile
   else if(pPDG==-11) part=G4Positron::Positron();   // Definition of positron projectile
   else if(pPDG==13) part=G4MuonMinus::MuonMinus();  // Definition of mu- projectile
   else if(pPDG==-13) part=G4MuonPlus::MuonPlus();   // Definition of the mu+ projectile
   else if(pPDG==2112) part=G4Neutron::Neutron();    // Definition of Neutron projectile
   else if(pPDG==-211) part=G4PionMinus::PionMinus();// Definition of Pi- projectile
   else if(pPDG==211)  part=G4PionPlus::PionPlus();  // Definition of Pi+ projectile
   else if(pPDG==321)  part=G4KaonPlus::KaonPlus();  // Definition of K+ projectile
   else if(pPDG==-321) part=G4KaonMinus::KaonMinus();// Definition of K- projectile
   else if(pPDG==130)part=G4KaonZeroLong::KaonZeroLong();//Definition of K0_L projectile
   else if(pPDG==310)part=G4KaonZeroShort::KaonZeroShort();//Definition of K0S projectile
   else if(pPDG==-2212)part=G4AntiProton::AntiProton();//Define antiProton projectile
   else if(pPDG==-2112)part=G4AntiNeutron::AntiNeutron();// Define AntiNeutron projectile
   else if(pPDG>999999)
   {
     G4int ZN=pPDG%1000000;
     G4int Z=ZN/1000;
     G4int A=Z+ZN%1000;
#ifdef debug
     G4cout<<"Test19: Nucleus with Z="<<Z<<", A="<<A<<G4endl; // only one run
#endif
     part = G4ParticleTable::GetParticleTable()->FindIon(Z,A,0,Z);// Define the Generic Ion
   }
   else if(pPDG==3122) part=G4Lambda::Lambda();      // Definition of Lambda projectile
   else if(pPDG==3112) part=G4SigmaMinus::SigmaMinus();// Definition of Sigma- projectile
   else if(pPDG==3122) part=G4SigmaZero::SigmaZero();// Definition of Sigma0 projectile
   else if(pPDG==3222) part=G4SigmaPlus::SigmaPlus();// Definition of Sigma+ projectile
   else if(pPDG==3322) part=G4XiMinus::XiMinus();    // Definition of the Xi- projectile
   else if(pPDG==3312) part=G4XiZero::XiZero();      // Definition of the Xi- projectile
   else if(pPDG==3334) part=G4OmegaMinus::OmegaMinus(); // Definition of Omega- projectile
   else if(pPDG==-3122) part=G4AntiLambda::AntiLambda(); // Definition of AntiLambda
   else if(pPDG==-3112) part=G4AntiSigmaMinus::AntiSigmaMinus();// Definition of AntiSigma-
   else if(pPDG==-3122) part=G4AntiSigmaZero::AntiSigmaZero();// Definition of AntiSigma0
   else if(pPDG==-3322) part=G4AntiXiMinus::AntiXiMinus(); // Definition of the AntiXi-
   else if(pPDG==-3312) part=G4AntiXiZero::AntiXiZero(); // Definition of the Xi0
   else if(pPDG==-3334) part=G4AntiOmegaMinus::AntiOmegaMinus();// Definition of AntiOmega-
   else if(pPDG==15) part=G4TauMinus::TauMinus();    // Definition of tau- projectile
   else if(pPDG==-15) part=G4TauPlus::TauPlus();     // Definition of tau+ projectile
   else if(pPDG==12) part=G4NeutrinoE::NeutrinoE();  // Definition of e_neutrino projectile
   else if(pPDG==-12) part=G4AntiNeutrinoE::AntiNeutrinoE(); // anti-e_neutrino projectile
   else if(pPDG==14) part=G4NeutrinoMu::NeutrinoMu();  // Definition of mu-neutrino proj.
   else if(pPDG==-14) part=G4AntiNeutrinoMu::AntiNeutrinoMu(); // anti-mu_neutrino proj.
   //else if(pPDG==16) part=G4NeutrinoTau::NeutrinoTau(); // Definition of tau-neutrino
   //else if(pPDG==-16) part=G4AntiNeutrinoTau::AntiNeutrinoTau(); // anti-tau_neutrino
   else if(pPDG!=-3222) // Leave defaulf definition for Anti Sigma+ projectile
   {
     G4cerr<<"***Test19: "<<pPDG<<" is a PDG code of not supported particle"<<G4endl;
     G4Exception("***Test19: OnFlight Process is called for not supported particle");
   }
   G4double pMass = part->GetPDGMass();                 // Mass of the projectile
   //
   G4ThreeVector aPosition(nx*mm, ny*mm, nz*mm);
   G4DynamicParticle* dParticle = new G4DynamicParticle(part,aDirection,0.);// Dummy Energy

   //G4cout<<"Test19: before the new Track definition"<<G4endl;
   // General Track definition
   G4Track* gTrack = new G4Track(dParticle,aTime,aPosition); // Track definition (NoTarg)

   G4int    bnp=pQC.GetBaryonNumber();
   G4int    bns=pQC.GetStrangeness();
   // ---------- Define material for the simulation ------------------
   //G4double tgA     = 26.98; // @@ Important? Can it be just tgZ+tgN?
   G4int tgA        = tgZ+tgN; // @@ Temporary, not good
   G4String nameMat = "Thin Target";  // @@ Not important can be an arbitrary name
   // The material can be copied from the commented cMaterial Factory above
   G4Material* material = new G4Material("ZA_Isomer", 1.*g/cm3, 1);
   if(tgZ<93 && tgN<151)
   {
     G4Isotope* isotope = new G4Isotope("Isotop", tgZ, tgA);
     G4Element* element = new G4Element("ZA_Isotop", "ZA", 1);
     element->AddIsotope(isotope, 100.*perCent);
     material->AddElement(element, 1);
   }
#ifdef debug
   G4cout<<"Test19:--- Material is defined ---" << G4endl; // only one run
#endif
   // For the Material Factory case
   //G4String  nameMat  = "Si";
   //G4Material* material = G4Material::GetMaterial(nameMat); // Get definition of Material
   if(!material)
   {
     G4cout<<"Test19:Material "<<nameMat<<" is not defined in the Test19Material"<<G4endl;
     exit(1);
   }
   // ----------- Geometry definition (simple BOX) ----------------
   G4double dimX = 100.*cm;
   G4double dimY = 100.*cm;
   G4double dimZ = 100.*cm;
   G4Box* sFrame = new G4Box ("Box",dimX, dimY, dimZ);
   // Zero,Zero,Zero position
   G4LogicalVolume* lFrame = new G4LogicalVolume(sFrame,material,"Box",0,0,0);
   G4PVPlacement* pFrame = new G4PVPlacement(0,G4ThreeVector(),"Box", lFrame,0,false,0);
   assert(pFrame);
#ifdef pverb
   G4cout<<"Test19:###### Start new run #####" << G4endl; // only one run
#endif
#ifdef debug
   G4double totCN  = tgZ+part->GetPDGCharge();
   G4double totSN  = tgZ+part->GetPDGCharge();
   G4int    totBNN = tgZ+tgN+part->GetBaryonNumber();
#endif
   G4double theStep   = 0.01*micrometer;
   // G4int maxz = (G4int)((*(material->GetElementVector()))[0]->GetZ()) + 1;
#ifdef pdebug
   G4double sumKE=0.;
   G4double su2KE=0.;
   G4double sEnAB=0.;
   G4double sEnH=0.;
   G4int    sumAB=0;
   G4int targPDG=90000000+1000*tgZ+tgN;                 // PDG Code of the target
   G4cout<<"Test19: The particle: "<<part->GetParticleName()<<G4endl;
   G4cout<<"Test19: The material: "<<material->GetName()<<G4endl;
   G4cout<<"Test19: The target:   "<<targPDG<<G4endl;
   G4cout<<"Test19: The step:     "<<theStep/mm<<" mm"<<G4endl;
   G4cout<<"Test19: The position: "<<aPosition/mm<<" mm"<<G4endl;
   G4cout<<"Test19: The direction:"<<aDirection<<G4endl;
   G4cout<<"Test19: The time:     "<<aTime/ns<<" ns"<<G4endl;
#endif
#ifdef tdebug
   for(G4int ip=0; ip<nT; ip++) tVal[ip]=(fT+dT*ip)/1000000.;// Fill the t-histogram
#endif
   // Create a DynamicParticle
   for(G4int nen=0; nen<cnE; nen++)                          // LOOP over projectile energy
   {
    G4double  energy = (ep-mp)*MeV;                          // projectile kinetic energy
    if(cnE>1) energy = eli[nen]*MeV;                         // Kin energy from the vector
    G4double  pvE = energy+mp;                               // projectile "visible" energy
    if(bnp>0)
    {
      if(bns>0) pvE-= bnp*mNuc;
      else      pvE-= mp;
    }
#ifdef debug
    G4cout<<"Test19: M="<<mp<<", T="<<energy<<", MeV"<<G4endl;
#endif
    dParticle->SetKineticEnergy(energy); // Fill the Kinetic Energy of the projectile

    for(G4int tgi=0; tgi<tgm; tgi++) // Loop over materials
    {
     nCur++;
     if (tgm>1)
     {
      tPDG=tli[tgi];
      tgZ = (tPDG-90000000)/1000;
      tgN = tPDG-90000000-tgZ*1000;
      material=mat[tgi];
      G4Element* curEl=(*(material->GetElementVector()))[0];
      G4cout<<"Test19: Material="<<material->GetName()<<", Element[0]="<<curEl->GetName()
            <<",A[0]="<<(*(curEl->GetIsotopeVector()))[0]->GetN()<<" is selected."<<G4endl;
     }
     G4cout<<"Test19:NewRun: Targ="<<tPDG<<",M="<<G4QPDGCode(tPDG).GetMass()<<", Proj="
           <<pPDG<<", E="<<energy<<" MeV, Run #"<<nCur<<" of "<<nTot<<G4endl;
     //G4double mt=G4QPDGCode(tPDG).GetMass();             // @@ just for check
     G4QContent tQC=G4QPDGCode(tPDG).GetQuarkContent();
     G4int    ct=tQC.GetCharge();
     G4int    st=tQC.GetStrangeness();
     G4int    bnt=tQC.GetBaryonNumber();
#ifdef debug
     G4cout<<"Test19: pQC"<<pQC<<", pch="<<cp<<", tQC"<<tQC<<", tch="<<ct<<G4endl;
#endif
     G4int    totC=cp+ct;
     G4int    totS=sp+st;
     G4int    totBN=bnp+bnt;
#ifdef debug
     G4cout<<"Test19:tC="<<totC<<"?="<<totCN<<",tB="<<totBN<<"?="<<totBNN<<",tS="<<totS
           <<"?="<<totSN<<G4endl;
#endif
     // Step Definition
     G4Step* step = new G4Step();
     step->SetTrack(gTrack);          // Step is initialized by the Track (?)

     G4StepPoint* aPoint = new G4StepPoint(); // It cant be initialized right away (!?)
     aPoint->SetPosition(aPosition);
     aPoint->SetMaterial(material);
     G4double safety = 10000.*cm;
     aPoint->SetSafety(safety);
     step->SetPreStepPoint(aPoint);   // Begin of the step

     //G4StepPoint* bPoint = aPoint;
     G4StepPoint* bPoint = new G4StepPoint(); // It cant be initialized right away (!?)
     bPoint->SetMaterial(material);
     bPoint->SetSafety(safety);
     G4ThreeVector bPosition = aDirection*theStep+aPosition;// aDirection is defined byCard
     bPoint->SetPosition(bPosition);
     step->SetPostStepPoint(bPoint);  // End of the step
     step->SetStepLength(theStep);    // Step is set byCard above
#ifdef pverb
     G4cout<<"Test19: The end point is defined and filled in the step "<<G4endl;
#endif
     G4Navigator* nav = new G4Navigator;
#ifdef pverb
     G4cout<<"Test19: The Navigator is defined "<<G4endl;
#endif
     nav->SetWorldVolume(pFrame);
#ifdef pverb
     G4cout<<"Test19: The Box frame is set "<<G4endl;
#endif
     G4TouchableHandle touch(nav->CreateTouchableHistory());
#ifdef pverb
     G4cout<<"Test19: The TouchableHandle is defined "<<G4endl;
#endif
     G4Timer* timer = new G4Timer();
     timer->Start();
#ifdef debug
     G4cout<<"Test19: Run is started, timer is started, kinEnergy="<<energy<<G4endl;
#endif
     const G4DynamicParticle* sec = 0;
     G4ParticleDefinition* pd = 0;
     G4ThreeVector  mom(0.,0.,0.);
     G4LorentzVector totSum, lorV;
     G4double e, p, m, totKE;
     // @@ G4double px, py, pz, pt;
     //G4VParticleChange* aChange = 0;
     G4ParticleChange* aChange = 0;
     G4double e0 = energy+pMass;
     G4ThreeVector pmax=std::sqrt(e0*e0-pMass*pMass)*aDirection;
     //G4double et=e0+mt;
     //G4int nEvt=100;
     // Randomization loop: cycle random generator, using 2 lower digits in nEvt
     G4int    iRandCount = nEvt%100;
     G4double vRandCount = 0.;
     while (iRandCount>0)                // Shift of the RNDN values 
     {
      vRandCount = G4UniformRand();     // Fake calls
      iRandCount--;
     }
#ifdef tdebug
     for(G4int it=0; it<nT; it++) tSig[it]=0; // Reset the output value
#endif
#ifdef csdebug
     //G4cout<<"Test19: before new Mu-Nuclear process"<<G4endl;
     // -----> For GHAD muons
     ///G4MuNuclearInteraction* HadrPR = new G4MuNuclearInteraction;// GHAD MuNuclear Proc.
     ///HadrPR->SetPhysicsTableBining(.01*GeV, 7.e8*GeV, 1000); // For the table
     ///HadrPR->BuildPhysicsTable(*part);      //NotNecessary for CHIPS G4QInelastic
     // -----> For GHAD hadrons
     //G4BGGPionElasticXS    barGGPiAElXS(part);//Barashenkov for <100GeV, GG for >100GeV
     //barGGPiAElXS.BuildPhysicsTable(*part);    //NotNecessary for CHIPS G4QInelastic
     //G4BGGNucleonElasticXS    barGGNAElXS(part);//Barashenkov for <100GeV, GG for >100GeV
     //barGGNAElXS.BuildPhysicsTable(*part);    //NotNecessary for CHIPS G4QInelastic
     //G4NucleonNuclearCrossSection BarashNAElXS;// NA elastic Barashenkov approximation
     //G4GlauberGribovCrossSection  GGNAElXS;    // NA GlaGrib approximation: high energies
     //G4PiNuclearCrossSection  barashPiXS;      // Barashenkov parameterization of pi-A XS
     G4CrossSectionDataStore  theElasticXS;     // GEISHA Elastic hadron-A XS
     G4CrossSectionDataStore  theCaptureXS;     // GEISHA neutron capture XS
     G4CrossSectionDataStore  theFissionXS;     // GEISHA fission XS
     G4HadronElasticDataSet   theElasticData;   // GEISHA Elastic XS Tables
     G4HadronCaptureDataSet   theCaptureData;   // GEISHA Capture XS Tables
     G4HadronFissionDataSet   theFissionData;   // GEISHA Fission XS Tables
     G4NeutronHPElasticData   theHPElasticData; // NeutronHP Elastic XS Tables
     G4NeutronHPCaptureData   theHPCaptureData; // NeutronHP Capture XS Tables
     G4NeutronHPFissionData   theHPFissionData; // NeutronHP Fission XS Tables
     //theElasticXS.AddDataSet(&theElasticData);  // Put LHEP parameterization in XS class
     theCaptureXS.AddDataSet(&theCaptureData);  // Put LHEP parameterization in XS class
     //theFissionXS.AddDataSet(&theFissionData);  // Put LHEP parameterization in XS class
     theElasticXS.AddDataSet(&theHPElasticData);// Put HP parameterization in XS class
     //theCaptureXS.AddDataSet(&theHPElasticData);// Put HP parameterization in XS class
     theFissionXS.AddDataSet(&theHPFissionData);// Put HP parameterization in XS class
     //G4CrossSectionDataStore  theInelasticXS;   // GEISHA Inelastic hadron-A XS
     //G4HadronInelasticDataSet theInelasticData; // GEISHA Inelastic XS Tables
     //theInelasticXS.AddDataSet(&theInelasticData);// Put parameterization in the XS class
     // -----> For CHIPS
     // ..... CHIPS on the Process level
     ///G4QInelastic* HadrPR = new G4QInelastic(); // CHIPS Inelastic
     ///G4QInelastic* HadrPR = new G4QElastic(); // CHIPS Inelastic
     // ..... CHIPS on the Cross-Section level
     G4VQCrossSection* HadrCS = 0;             // ProtoPointer to CHIPS CrossSections
     
     // CHIPS Inelastic
     G4QNeutronCaptureRatio* theNCaptureRatio=G4QNeutronCaptureRatio::GetPointer();
     if(pPDG==13 || pPDG==-13) HadrCS = G4QMuonNuclearCrossSection::GetPointer();// MuNuc
     else if(pPDG==2212)    HadrCS = G4QProtonNuclearCrossSection::GetPointer(); // pA
     else if(pPDG==2112)    HadrCS = G4QNeutronNuclearCrossSection::GetPointer(); // nA
     else if(pPDG==-211)    HadrCS = G4QPionMinusNuclearCrossSection::GetPointer();//pi-A
     else if(pPDG== 211)    HadrCS = G4QPionPlusNuclearCrossSection::GetPointer(); //pi+A
     else if(pPDG==-321)    HadrCS = G4QKaonMinusNuclearCrossSection::GetPointer();// K-A
     else if(pPDG== 321)    HadrCS = G4QKaonPlusNuclearCrossSection::GetPointer(); // K+A
     else if(pPDG== 311 || pPDG==-311 || pPDG== 310 || pPDG== 130)
                            HadrCS = G4QKaonZeroNuclearCrossSection::GetPointer(); // K0A
     else if(pPDG==3222)    HadrCS = G4QHyperonPlusNuclearCrossSection::GetPointer();
     else if(pPDG>3000 && pPDG<4000 && pPDG!=-3222)                                 // @@
                            HadrCS = G4QHyperonNuclearCrossSection::GetPointer();
     else if(pPDG==-3112 || pPDG==-3312 || pPDG==-3334)
                            HadrCS = G4QAntiBaryonPlusNuclearCrossSection::GetPointer();
     else if(pPDG>-4000 && pPDG<-2000 && pPDG!=-3112 && pPDG!=-3312 && pPDG!=-3334) // @@
                            HadrCS = G4QAntiBaryonNuclearCrossSection::GetPointer();
     // CHIPS Elastic
     //if     (pPDG==2212)    HadrCS = G4QProtonElasticCrossSection::GetPointer(); // pA
     //else if(pPDG==2112)    HadrCS = G4QNeutronElasticCrossSection::GetPointer(); // nA
     //else if(pPDG==-211)    HadrCS = G4QPionMinusElasticCrossSection::GetPointer();//pi-A
     //else if(pPDG== 211)    HadrCS = G4QPionPlusElasticCrossSection::GetPointer(); //pi+A
     //else if(pPDG==-321)    HadrCS = G4QKaonMinusElasticCrossSection::GetPointer();// K-A
     //else if(pPDG== 321)    HadrCS = G4QKaonPlusElasticCrossSection::GetPointer(); // K+A
     ////else if(pPDG== 311 || pPDG==-311 || pPDG== 310 || pPDG== 130)
     ////                       HadrCS = G4QKaonZeroElasticCrossSection::GetPointer();//K0A
     //else if(pPDG==3222)    HadrCS = G4QHyperonPlusElasticCrossSection::GetPointer();//SiP
     //else if(pPDG>3000 && pPDG<4000 && pPDG!=-3222)                                 // @@
     //                       HadrCS = G4QHyperonElasticCrossSection::GetPointer();
     ////else if(pPDG==-3112 || pPDG==-3312 || pPDG==-3334)
     ////                      HadrCS = G4QAntiBaryonPlusElasticCrossSection::GetPointer();
     //else if(pPDG>-4000 && pPDG<-2000)
     //                       HadrCS = G4QAntiBaryonElasticCrossSection::GetPointer();
     // ______ End of CHIPS/GHAD ___________
     //G4cout<<"Test19: Cross-section process is defined pPDG="<<pPDG<<G4endl;
     // --- A temporary LOOP for calculation of total cross section ------------
     //G4double pMin=.02;                       // in GeV --> for protons
     //G4double pMax=370.;                      // in GeV ==> for HE (CHIPS/LHEP)
     //G4double pMax=.22;                       // in GeV ==> for HP
     //G4double pMax=10000000.;                 // in GeV --> np
     //G4double pMax=1000.;                     // in GeV --> np->inelastic
     G4double pMax=1.;                        // in GeV --> LHEP/CHIPS (Capture/Fission)
     //G4double pMin=.03;                       // in GeV --> np->dg
     //G4double pMin=.002;                      // in GeV --> for HE/HP
     G4double pMin=.000000162;                // in GeV --> for neutrons (CHIPS,LHEP)
     //G4double pMin=.000000177;                 // in GeV --> for neutrons (HP)
     G4int nic=50;                            // Number of points
     G4double lpMin=std::log(pMin);
     G4double lpMax=std::log(pMax);
     G4double dlp=(lpMax-lpMin)/nic;
     G4double lmic=lpMin-dlp/2;
     G4double hMa=0.;                         // Mass of a hadron in GeV
     if(pPDG==2212) hMa=.938272;                  // Mass of a proton in GeV
     else if(pPDG==2112) hMa=.93957;              // Mass of a neutron in GeV
     else if(pPDG==211 || pPDG==-211) hMa=.13957; // Mass of a charged pion in GeV
     else G4cout<<"-Warning-Test19: *** M=0. *** Add mass for PDG="<<pPDG<<G4endl;
     G4double hMa2=hMa*hMa;
     G4cout<<"Test19: mi="<<lmic+dlp<<",ma="<<lmic+dlp*nic<<",d="<<dlp<<",n="<<nic<<", Z="
           <<tgZ<<", A="<<tgZ+tgN<<G4endl;
     for(G4int ic=0; ic<nic; ic++)
     {
       lmic+=dlp;
       G4double mic=std::exp(lmic);            // current Momentum in GeV
       G4double p2=mic*mic;
       G4double ken=std::sqrt(p2+hMa2)-hMa;
       // CHIPS calculation by G4QElasticCrossSection___ Only for CHIPS CS level
       //G4double CS = HadrCS->GetCrossSection(false, mic*GeV, tgZ, tgN, pPDG);
       G4double CS = HadrCS->GetCrossSection(false, mic*GeV, tgZ, tgN, pPDG)*
                     theNCaptureRatio->GetRatio(mic*GeV, tgZ, tgN);
       //
       // === From here CHECK of energy-momentum conservation starts ===
       //G4double den=0.;
       //gTrack->SetStep(step);             // Now step is included in theTrack (see above)
       //gTrack->SetKineticEnergy(ken*GeV); // Duplication of KinEnergy for the Track (?!)
       //gTrack->SetTouchableHandle(touch); // Set Box touchable history
       ////
       ////G4cout<<"Test19: before MeanFreePath out of Loop"<<G4endl;
       //G4double MFP = 2.e99;
       //MFP = HadrPR->GetMeanFreePath(*gTrack,.1,cond);
       ////G4cout<<"Test19: after MeanFreePath out of Loop MFP="<<MFP<<G4endl;
       //G4int nen=10000;
       //G4int cen=0;
       /////if(CS>0.) for(G4int ie=0; ie<nen; ie++) //@@ for CHIPS only
       //if(MFP<1.e99)for(G4int ie=0; ie<nen; ie++)
       //{
       //  //
       //  gTrack->SetStep(step);            // Now step is included inTheTrack (see above)
       //  gTrack->SetKineticEnergy(ken*GeV);// Duplication of KinEnergy for the Track (?!)
       //  gTrack->SetTouchableHandle(touch);// Set Box touchable history
       //  //
       //  //G4cout<<"Test19: before MeanFreePath"<<G4endl;
       //  HadrPR->GetMeanFreePath(*gTrack,.1,cond);
       //  // GHAD/CHIPS Energy Deposition
       //  //G4cout<<"Test19: before PostStepDoIt ie="<<ie<<G4endl;
       //  G4ParticleChange* aChange =
       //              static_cast<G4ParticleChange*>(HadrPR->PostStepDoIt(*gTrack,*step));
       //  //G4cout<<"Test19:afterPostStepDoIt, Alive="<<aChange->GetTrackStatus()<<G4endl;
       //  //if(aChange->GetTrackStatus()==fAlive) den += ken*GeV-aChange->GetEnergy();//dE
       //  if(aChange->GetTrackStatus()==fAlive && 1.-aChange->GetEnergy()/ken/GeV>.000001)
       //  {
       //    G4double dz=aChange->GetMomentumDirection()->z();
       //    den += std::sqrt(1.-dz*dz);
       //    cen++;
       //  }
       //  G4int nS = aChange->GetNumberOfSecondaries();
       //  //G4cout<<"Test19: before delete secondaries, nS="<<nS<<G4endl;
       //  if(nS) for(G4int ides=0; ides<nS; ides++) delete aChange->GetSecondary(ides);
       //  //G4cout<<"Test19: before aCange->Clear, nS="<<nS<<G4endl;
       //  aChange->Clear();
       //  //
       dParticle->SetKineticEnergy(ken*GeV);   // Fill Kinetic Energy of thep rojectile
       //  // GHAD Cross-section on process level
       //  //G4double MFP = HadrPR->GetMeanFreePath(*gTrack,.1,cond);
       // ..... GHAD on the Cross-Section level .... mu-nuclear
       //G4double aZ=tgZ;
       //G4double aA=(tgZ+tgN)*g/mole;
       //G4double CS = HadrPR->ComputeMicroscopicCrossSection(part,ken*GeV,aZ,aA);
       //G4double CS = HadrPR->GetElasticCrossSection(dParticle,element);
       // ..... GHAD on the Cross-Section level .... pions
       // === Elastic ===
       ///G4double CS = theElasticXS.GetCrossSection(dParticle,element,T0);//GHAD/HPElastic
       ///G4double CS = theCaptureXS.GetCrossSection(dParticle,element,T0);//GHAD/HPCapture
       ///G4double CS = theFissionXS.GetCrossSection(dParticle,element,T0);//GHAD/HPFission
       //barashPiXS.GetCrossSection(dParticle,element,T0);// BarashenkovPiAin (initializes)
       //G4double CS= barashPiXS.GetElasticXsc();// BarashenkovPiAel (call after GetCrSec)
       //G4double CS = barGGNAElXS.GetCrossSection(dParticle,element,T0);//BarashenkGG_NAel
       //G4double CS = barGGPiAElXS.GetCrossSection(dParticle,element,T0);//BarashenGGPiAel
       //BarashNAElXS.GetCrossSection(dParticle,element,T0);//BarashenkovNAel
       //G4double CS = BarashNAElXS.GetElasticXsc();//BarashenkovNAel (call after GetCrSec)
       //GGNAElXS.GetCrossSection(dParticle,element,T0);// GlaGribHadronAin
       //G4double CS = GGNAElXS.GetElasticGlauberGribovXsc();//GlaGribHAEl(callAfterGetCS)
       // === Inelastic ===
       //G4double CS = theElasticXS.GetCrossSection(dParticle,element,T0) +
       //              theInelasticXS.GetCrossSection(dParticle,element,T0);// GEISHA total
       //G4double CS = theInelasticXS.GetCrossSection(dParticle,element,T0);// GEISHA total
       //G4double ICS= barashPiXS.GetCrossSection(dParticle,element,T0);// BarashenkovPiAin
       //G4double ECS= barashPiXS.GetElasticXsc();// BarashenkovPiAin (call after GetCrSec)
       //G4double CS = barashPiXS.GetTotalXsc();// BarashenkovPiAtot (call after GetCrosSe)
       //  // ==================== End of GHAD ==============================
       //  // CHIPS calculation by G4QElasticCrossSection
       //  // ..... CHIPS CS on the CS level
       //  //G4double CS = HadrCS->GetCrossSection(true, mic*GeV, tgZ, tgN, pPDG); 
       //  // ..... CHIPS dE on the CS level
       //  //den += HadrCS->GetExchangeEnergy(); 
       //  // ..... CHIPS CS on the Process level
       //  //G4double MFP = HadrPR->GetMeanFreePath(*gTrack,.1,cond);
       //} // End of energy loop
       // === End of CHECK of energy and momentum conservation
       //G4cout<<"Test19: P="<<mic" (GeV/c), MeanFreePath="<<MFP<<G4endl;
       //////////G4cout<<"Test19: P="<<mic<<" "<<CS/millibarn<<G4endl;
       G4cout<<mic<<" "<<CS/millibarn<<G4endl;
       //G4cout<<"Test19: P="<<mic<<" (GeV/c), XS="<<CS/millibarn/tgA<<G4endl; // Muons
       //G4cout<<"Test19: P="<<mic<<" (GeV/c), <dE>/E="<<den/nen/ken/GeV<<G4endl;
       //G4cout<<"Test19: P="<<mic<<" (GeV/c), <std::sin(t)>="<<den/cen<<G4endl;
     }
     // --- End of the temporary LOOP for calculation of total cross section ------------
     return EXIT_SUCCESS;
#endif
#ifdef idebug
     G4cout<<"Test19: Before the event loop, nEvents= "<<nEvt<<G4endl;
#endif
#ifdef pidebug
     // Tmp "Event counting for pions pi-N-N"
     G4int nP0=0;
     G4int nPP=0;
     G4int nPN=0;
     // End of Tmp
#endif
     G4double dTot=0.;
     G4double dEl=0.;
#ifdef pdebug
     sumKE=0.;
     su2KE=0.;
     sEnAB=0.;
     sEnH=0.;
     sumAB=0;
#endif
     G4int    nDoNothing=0;
     for (G4int iter=0; iter<nEvt; iter++)
     {
#ifdef debug
      G4cout<<"Test19: ### "<<iter<< "-th event starts.### energy="<<energy<<G4endl;
#endif

      if(!(iter%1000)&&iter)G4cout<<">TEST19: "<<iter<<" events are simulated"<<G4endl;

      gTrack->SetStep(step);            // Now step is included in the Track (see above)
      gTrack->SetKineticEnergy(energy); // Duplication of Kin. Energy for the Track (?!)
      gTrack->SetTouchableHandle(touch);// Set Box touchable history
#ifdef debug
      G4cout<<"Test19: Before the fake proc->GetMeanFreePath call"<<G4endl;
#endif
      proc->GetMeanFreePath(*gTrack,0.1,cond);
#ifdef pdebug
      G4double mfp=proc->GetMeanFreePath(*gTrack,0.1,cond);
      G4cout<<"Test19: Before PostStepDoIt MeanFreePath="<<mfp/meter<<" [m]"<<G4endl;
#endif
      aChange = static_cast<G4ParticleChange*>(proc->PostStepDoIt(*gTrack,*step)); // Up
      G4int nSec = aChange->GetNumberOfSecondaries();
      G4TrackStatus lead = aChange->GetTrackStatus();
#ifdef debug
      G4cout<<"T19:AfterPostStepDoIt,nS="<<nSec<<",l="<<lead<<" =? fA="<<fAlive<<G4endl;
#endif
      if(!nSec && lead==fAlive)
      {
        G4double eOut=aChange->GetEnergy();
        if(std::fabs(eOut-energy) < .01) nDoNothing++; // Calculate # of DoNothing events
        else G4cout<<"**Tst19: DoNothing E="<<eOut<<" # "<<energy<<G4endl;
      }
      else
      {
#ifdef debug
       G4cout<<"T19: "<<nSec<<",L="<<lead<<","<<fAlive<<",tC="<<totC<<",tS="<<totS<<G4endl;
#endif
       G4LorentzVector totSumM(0.,0.,0.,0.);
       G4int begi=0;
       if(nSec || lead==fAlive)
       {
        G4int totCharge = totC; 
        G4int totStran  = totS; 
#ifdef chips
        G4int    curN=proc->GetNumberOfNeutronsInTarget();
#else
        // for GHAD
        G4int    curN = tgN; 
#endif
        G4int    dBN = curN-tgN;
        G4int    totBaryN = totBN+dBN;
        G4int    curPDG=tPDG+dBN;
        G4double curM=G4QPDGCode(curPDG).GetMass(); // Update mass of the TargetNucleus
        totSum = G4LorentzVector(pmax, e0+curM);
        totKE=0.;
        if(lead==fAlive)
        {
          begi=-1;                                  // Means that leading particle is alive
          totCharge-=static_cast<G4int>(aChange->GetCharge());
          totStran-=sp;
          totBaryN-=bnp;
          G4double sen=aChange->GetEnergy();
          G4double sma=aChange->GetMass();
          const G4ThreeVector* smo=aChange->GetMomentumDirection();
          G4double ten=sen+sma;
          G4LorentzVector s4m(aChange->CalcMomentum(sen,*smo,sma),ten);
#ifdef debug
          G4cout<<"Test19: E="<<sen+sma<<",p="<<(*smo)*std::sqrt(ten*ten-sma*sma)<<G4endl;
          G4cout<<"Test19: L4m="<<s4m<<",T4m="<<totSum<<", N="<<curN<<",M="<<curM<<G4endl;
#endif
          totKE+=s4m.e();
          if(bnp>0)
          {
            if(bns>0) totKE-= bnp*mNuc;
            else      totKE-= mp;
          }
          else if (bnp<0) totKE-=bnp*mNuc; 
          totSum-=s4m;
          totSumM=totSum;
        }
#ifdef debug
        G4cout<<"Test19: tChg="<<totC<<", T4m="<<totSum<<",tBar="<<totBaryN<<G4endl;
#endif
#ifdef chips
        G4LorentzVector Residual=proc->GetEnegryMomentumConservation();
#else
        // for GHAD
        G4LorentzVector Residual(0.,0.,0.,0.);
#endif
#ifdef debug
        G4double de = aChange->GetLocalEnergyDeposit();// Init TotalEnergy by EnergyDeposit
        G4cout<<"Test19: "<<nSec<<" secondary particles are generated, dE="<<de<<G4endl;
#endif
        // @@ ----------------------- Begin
        G4double weight=aChange->GetSecondary(0)->GetDynamicParticle()->GetKineticEnergy();
#ifdef debug
        G4cout<<"Test19:----------------: Weigh="<<weight<<G4endl;
#endif
        dTot+=weight;
        if(nSec>1) dEl+=weight;
        // @@ ----------------------- End
        //G4int nbar = 0;

        // @@ G4int npt=0;
        G4int    c=0;    // Prototype of the PDG Code of the particle
        G4int nGamma=0;
        G4double EGamma=0;
#ifdef pidebug
        // fake line
#else
        G4int nP0=0;
        G4int nPP=0;
        G4int nPN=0;
#endif
        G4int nKaons=0;
        G4int nEta=0;
        // @@ G4int nAlphas=0;
        G4int nPhotons=0;
        G4int nProtons=0;
        G4int nNeutrons=0;
        G4int nSpNeut=0;
        // @@ G4int nSpAlph=0;
        G4int nOmega=0;
        // @@ G4int nDec=0;
        // @@ G4int dirN=0;
#ifdef pdebug
        G4cout<<"Test19:----DONE^^^^^^^*******^^^^^^^^:ir="<<iter<<": #ofH="<<nSec<<G4endl;
        if(!(iter%1000)) G4cerr<<"#"<<iter<<G4endl;
#endif
        G4bool alarm=false;
        // @@ G4bool rad=false;
        // @@ G4bool hyp=false;
        // @@ G4bool badPDG=false;
        // ------- LOOP over secondary particles -------
        G4int cST=0;
        G4int cCG=0;
        G4int cBN=0;
        for(G4int i=begi; i<nSec; i++)
        {
          if(i<0) // Leading particle
          {
            e  = aChange->GetEnergy();
            m  = aChange->GetMass();
            mom= *(aChange->GetMomentumDirection());
            c  = pPDG;
#ifdef debug
            G4cout<<"Test19:LeadingParticle is found, E="<<e<<",m="<<m<<",PDG="<<c<<G4endl;
#endif
            cST= sp;
            cCG= cp;
            cBN= bnp;

          }
          else    // Secondaries
          {
#ifdef debug
            G4cout<<"Test19:Secondary i="<<i<<",*Ses="<<aChange->GetSecondary(i)<<G4endl;
#endif
            sec = aChange->GetSecondary(i)->GetDynamicParticle();
#ifdef debug
            G4cout<<"Test19:Secondary i="<<i<<",*SPD="<<sec<<G4endl;
#endif
            pd  = sec->GetDefinition();
#ifdef debug
            G4cout<<"Test19:Secondary i="<<i<<",*PD="<<pd<<G4endl;
#endif
            c   = pd->GetPDGEncoding();
            cST = static_cast<G4int>(pd->GetQuarkContent(3)-pd->GetAntiQuarkContent(3));
            cCG = static_cast<G4int>(pd->GetPDGCharge());
            cBN = static_cast<G4int>(pd->GetBaryonNumber());
            if(!c) c=90000000+cST*999999+cCG*999+cBN;
            m   = pd->GetPDGMass();
            mom = sec->GetMomentumDirection();
            e   = sec->GetKineticEnergy();
#ifdef synch
            s_n+=1.;
            s_s+=e;
            s_d+=e*e;
            G4int ebin=(G4int)((std::log(e)-logmin)/dlog);
            if(ebin >= 0 && ebin < nBin)
            {
              hBin[ebin]+=1.;
              eBin[ebin]+=e;
            }
#endif
#ifdef debug
            G4cout<<"Test19:SecPt#"<<i<<",T="<<e<<",m="<<m<<",PDG="<<c<<",C="<<cCG<<G4endl;
#endif
          }
          if (e < 0.0)
          {
            G4cout<<"**Test19:Event#"<<iter<<",Hadron#"<<i<<",T="<<e<<"<0 (Set 0)"<<G4endl;
            e = 0.0;
          }
          // for exclusive reaction 2 particles in final state
          p = std::sqrt(e*(e + m + m));
          mom *= p;
#ifdef debug
          G4cout<<"Test19:Secondary #"<<i<<" E="<<e+m<<", P="<<mom<<G4endl;
#endif
          //if(i==1) // Means the target secondary for the ellastic scattering
          //{
          //  G4double t=e*m;
          //  t=t+t;
          //  G4int tb=static_cast<G4int>(t/dT);
          //  if(tb>nT1) tb=nT1;
          //  tSig[tb]++;
          //}
          lorV = G4LorentzVector(mom, e+m);    // "e" is a Kinetic energy!
          G4double visE=0.;
          if(i>=0)                             // Do not count LV & charges of leptons
          {
            totStran-=cST;
            totCharge-=cCG;
            totBaryN-=cBN;
            totSum -= lorV;
            visE=lorV.e();
            if(cBN>0)
            {
              if(cST>0) visE-=cBN*mNuc;
              else      visE-=m;
            }
            else if(cBN<0) visE-=cBN*mNuc;
            totKE +=  visE;
#ifdef debug
            G4cout<<"Test19:i="<<i<<",tS="<<totStran<<",tC="<<totC<<",t4="<<totSum<<",tB="
                  <<totBaryN<<G4endl;
#endif
          }
#ifdef debug
          G4cout<<"Test19: *** M="<<lorV.m()<<G4endl;
#endif
          //if(fabs(m-lorV.m())>.005&&1>2) // @@ Temporary closed
          if(std::fabs(m-lorV.m())>.005) // @@ Temporary closed
          {
            G4cout<<"***Test19: m="<<lorV.m()<<" # "<<m<<", d="<<lorV.m()-m<<G4endl;
            alarm=true;
          }
          if(!(lorV.e()>=0||lorV.e()<0)   || !(lorV.px()>=0||lorV.px()<0) ||
             !(lorV.py()>=0||lorV.py()<0) || !(lorV.pz()>=0||lorV.pz()<0))
          {
            G4cerr<<"***Test19: NAN in LorentzVector="<<lorV<<G4endl;
            alarm=true;
          }
          if(c==90000002||c==90002000||c==92000000)
          {
            G4cout<<"***Test19:***Dibaryon *** i="<<i<<", PDG="<<c<<G4endl;
            alarm=true;
          }
          if(c==90000003||c==90003000||c==93000000)
          {
            G4cout<<"***Test19:***Tribaryon *** i="<<i<<", PDG="<<c<<G4endl;
            alarm=true;
          }
#ifdef debug
          G4cout<<"Test19: ===> alarm="<<alarm<<G4endl;
#endif
          if(c==223) nOmega++;
          if(c==22) nPhotons++;
          if(c==311||c==321||c==-311||c==-321) nKaons++; // kaons
          if(c==221) nEta++;                             // etas
          //if(c==90002002) nAlphas++;                     // Alphas
          if(c==2212) nProtons++;                        // Protons
          if(c==2112) nNeutrons++;                       // Neutrons
          if(c==2112 && std::fabs(e-1005.)<3.) nSpNeut++;// Dibar-Neutrons
          //if(c==90002002 && e-m<7.) nSpAlph++;           // Special Alphas
#ifdef pidebug
          if(nSec==3 && c== 111) nP0++;                  // Neutral  pions
          if(nSec==3 && c==-211) nPN++;                  // Negative pions
          if(nSec==3 && c== 211) nPP++;                  // Positive pions
#else
          if(c== 111) nP0++;                             // Neutral  pions
          if(c==-211) nPN++;                             // Negative pions
          if(c== 211) nPP++;                             // Positive pions
#endif
          if(c==22)                                      // Gammas
          {
            nGamma++;                                    // # of gammas
            EGamma+=e;                                   // Total energy in gammas
          }
#ifdef ppdebug
          if(lead==fAlive && (i==-1 || c==2112 || c==2212 || c==90000001 || c==90001000))
          //if(lead==fAlive)
            G4cout<<"Test19:#"<<i<<",PDG="<<c<<",S="<<cST<<",C="<<cCG<<",B="<<cBN<<",4M="
                  <<lorV<<m<<",T="<<lorV.e()-m<<",vE="<<visE<<G4endl;
#endif
          //delete aChange->GetSecondary(i);
#ifdef pdebug
          if(cBN>0 && cST>0)
          {
            G4cout<<"---Test19:--- Hyperon with PDG="<<c<<",S="<<cST<<G4endl;
            sEnH+=visE;
          }
          if(cBN<0)
          {
            //if(cST<0) G4cout<<"***Test19:*** Anti-Hyperon with PDG="<<c<<G4endl;
            G4cout<<"***Test19:*** Anti-Baryon with PDG="<<c<<",S="<<cST<<G4endl;
            sEnAB+=visE;
            ++sumAB;
          }
          //if(lead==fAlive && (i==-1 || c==2112 || c==2212 || c==90000001 || c==90001000))
          G4cout<<"Test19:#"<<i<<",PDG="<<c<<",C="<<cCG<<",B="<<cBN<<",4M="<<lorV<<m<<",T="
                <<lorV.e()-m<<",vE="<<visE<<G4endl;
#endif
        } // End of the LOOP over secondaries
        // delete secondaries in the end of the event         
#ifdef chips
        G4double ss=std::fabs(totSum.t())+std::fabs(totSum.x())+
                    std::fabs(totSum.y())+std::fabs(totSum.z());
#ifdef inter
	  G4double misr=.27;
#else
	  G4double misr=-DBL_MAX;
#endif
        G4double sr=std::fabs(Residual.t())+std::fabs(Residual.x())+
                    std::fabs(Residual.y())+std::fabs(Residual.z());    
#endif
#ifdef lhepdbg
        // --- Start detailed correction
        G4double resGSM =
                   G4QNucleus(90000000+totStran*999999+totCharge*999+totBaryN).GetGSMass();
        G4double resMom = totSum.rho();
        G4double resE= std::sqrt(resGSM*resGSM+resMom*resMom);
        G4double redE= totSum.e()-resE;
        // *** Debug print
        //G4double resT= resE-resGSM;
        //G4cout<<"Test19: LHEP_COR_E="<<redE<<", T="<<resT<<", P="<<resMom<<G4endl;
        // *** End of debug print
        totSum=G4LorentzVector(0.,0.,0.,redE);
        // --- Stop detailed correction
        totStran=0;
        totCharge=0;
        totBaryN=0;
        // --- Print flag for
        G4bool prtf=false;
        //if(std::fabs(totSum.e())>.1) 
        //{
        //  prtf=true;
        //  G4cout<<"Test19:----->LHEP_COR_E="<<redE<<", T="<<resT<<", P="<<resMom<<G4endl;
        //}
#endif
#ifdef meandbg
        //G4double eng=energy+mp;
        G4double eng=energy;
        G4double cE=eng-totSum.t();
        if(cE<0.) G4cout<<"Test19: **Negative**, E="<<cE<<G4endl;
        if(cE>eng+eng) nD++;
        nE++;
        sE+=cE;
        sE2+=cE*cE;
        pE+=eng;
#endif
#ifdef ekindbg
        if(totKE<0.) G4cout<<"Test19: **Negative**, E="<<totKE<<G4endl;
        if(totKE>energy+energy) nD++;
        nE++;
        sE+=totKE;
        sE2+=totKE*totKE;
        pE+=energy+mp;
#endif
#ifdef histdbg
       
        if(nGamma>=nHst) nGamma=nHst-1;
        nGaH[nGamma]++;
        G4int nESG=static_cast<G4int>(EGamma+.00001);
        if(nESG>=nHst) nESG=nHst-1;
        GaSE[nESG]++;
        if(totCharge>=dQnM) totCharge=dQnM1;
        if(totCharge<-dQnM) totCharge=-dQnM;
        dChH[totCharge+dQnM]++;
        if(totBaryN>=dQnM) totBaryN=dQnM1;
        if(totBaryN<-dQnM) totBaryN=-dQnM;
        dBnH[totBaryN+dQnM]++;
        G4double dE=totSum.t();
        G4int rdE=static_cast<G4int>(dQnM*20.*dE/e0-17+dQnM);
        //G4int dEn=static_cast<G4int>(dE*10.+dQnM);
        G4int dEn=static_cast<G4int>(dE/10.+dQnE);
        if(rdE>=nHst) rdE=nHst1;
        if(rdE<0) rdE=0;
        if(dEn>=nHst)
        {
          //G4cout<<"Test19: Too big energy nonconservation dE="<<dE<<G4endl;
          dEn=nHst1;
        }
        if(dEn<0) dEn=0;
        rEnH[rdE]++;
        dEnH[dEn]++;
        //if(nPN || nP0 || nPP) dEnH[dEn]++;
        G4double dZ=totSum.z();
        G4int rdZ=static_cast<G4int>(dQnM*20.*dZ/pmax-17+dQnM);
        G4int dZn=static_cast<G4int>(dZ*10.+dQnM);
        if(rdZ>=nHst) rdZ=nHst1;
        if(rdZ<0) rdZ=0;
        if(dZn>=nHst) dZn=nHst1;
        if(dZn<0) dZn=0;
        rPzH[rdZ]++;
        dPzH[dZn]++;
        G4double dR=std::sqrt(totSum.x()*totSum.x()+totSum.y()*totSum.y());
        G4int rdR=static_cast<G4int>(dQnM*20.*dR/pmax);
        G4int dRn=static_cast<G4int>(dR*100.);
        if(rdR>=nHst) rdR=nHst1;
        if(rdR<0) rdR=0;
        if(dRn>=nHst) dRn=nHst1;
        if(dRn<0) dRn=0;
        rPrH[rdR]++;
        dPrH[dRn]++;
        //if(!(iter%100)&&iter) G4cout<<"BN="<<totBaryN<<",dE="<<dE<<",E="<<e0<<",dZ="<<dZ
        //                            <<",P="<<pmax<<",dR="<<dR<<G4endl;
#endif
#ifdef pdebug
        sumKE+=totKE;
        su2KE+=totKE*totKE;
        G4cout<<">TEST19: r4M="<<totSum<<", rCh="<<totCharge<<", rBN="<<totBaryN<<", vE="
              <<totKE<<", vR="<<totKE/pvE<<G4endl;
#endif
#ifdef lhepdbg 
        if(prtf)
#endif
#ifdef chips
        if (totCharge || totBaryN || ss>.27 || alarm || (nGamma && !EGamma))// OnlyForCHIPS
        {
          G4cout<<"*Warning*Test19:#"<<iter<<":n="<<nSec<<",4M="<<totSum<<",Ch="<<totCharge
                <<",BaryN="<<totBaryN<<", R="<<Residual<<",D2="<<ss<<",nN="<<curN<<G4endl;
          totSum = totSumM;
          if(nGamma && !EGamma) G4cout<<"***Test19: Egamma=0"<<G4endl;
          for (G4int indx=begi; indx<nSec; indx++)
          {
            if(indx<0) // Leading particle
            {
              e  = aChange->GetEnergy();
              m  = aChange->GetMass();
              mom= *(aChange->GetMomentumDirection());
              c  = pPDG;
              cST= sp;
              cCG= cp;
              cBN= bnp;
            }
            else    // Secondaries
            {
              sec = aChange->GetSecondary(indx)->GetDynamicParticle();
              pd  = sec->GetDefinition();
              c   = pd->GetPDGEncoding();
              cST = static_cast<G4int>(pd->GetQuarkContent(3)-pd->GetAntiQuarkContent(3));
              cCG = static_cast<G4int>(pd->GetPDGCharge());
              cBN = static_cast<G4int>(pd->GetBaryonNumber());
              if(!c) c=90000000+cST*999999+cCG*999+cBN;
              m   = pd->GetPDGMass();
              //G4cout<<"***Tmp***Test19:#"<<indx<<",S="<<cST<<",C="<<cCG<<",B="<<cBN
              //      <<",M="<<m<<G4endl;
              mom = sec->GetMomentumDirection();
              e   = sec->GetKineticEnergy();
            }
            p = std::sqrt(e*(e + m + m));
            mom *= p;
            lorV = G4LorentzVector(mom, e + m);    // "e" is a Kinetic energy!
            totSum -= lorV;
            G4cout<<"Test19:#"<<indx<<",PDG="<<c<<",m="<<m<<",4M="<<lorV<<",T="<<e
                  <<", d4M="<<totSum<<", S="<<cST<<", C="<<cCG<<", B="<<cBN<<G4endl;
          }
          if(sr>misr)G4Exception("***Test19:ALARM/baryn/chrg/energy/mom isn't conserved");
          //G4Exception("***Test19: ALARM/baryn/chrg/energy/mom is not conserved");
        }
#endif
#ifndef nout
        ntp->FillEvt(aChange,dParticle); // Fill the simulated event in the ASCII "ntuple"
#endif
        // =============== May be print it here if it is not zero...
        for(G4int ides=0; ides<nSec; ides++) delete aChange->GetSecondary(ides);
        aChange->Clear();
#ifdef debug
        G4cout<<"Test19:--->>> Ev # "<<iter<<", After ntp.FillEvt"<<G4endl;
#endif
       } // End of nSec>0 IF
      } // End of Check for DoNothing
     } // End of the LOOP over events
#ifdef pdebug
     G4double meanVE=sumKE/nEvt;
     G4double rmsVE=su2KE/nEvt;
     rmsVE=std::sqrt(rmsVE-meanVE*meanVE);
     G4cout<<">Mean>Test19:RVE="<<meanVE/pvE<<",RMSVE="<<rmsVE/pvE<<",NAB="
           <<sumAB<<",RAB="<<sEnAB/nEvt/pvE<<",HAB="<<sEnH/nEvt/pvE<<G4endl;
#endif
     // Stop the timer to estimate the speed of the generation
#ifdef pidebug
     G4cout<<"Test19: nPi0="<<nP0<<", nPiP="<<nPP<<", nPiM="<<nPN<<G4endl;
#endif
     G4double dn=100.* nDoNothing / nEvt;
     timer->Stop();
     G4cout<<"Test19:Time/Ev="<<timer->GetUserElapsed()/nEvt<<" s, DN="<<dn<<" %"<<G4endl;
     delete timer;
#ifdef tdebug
     G4double pGeV=pmax/1000.;
     G4double alp=std::log(pGeV/(1.+1./pGeV/pGeV/pGeV));
     G4double exr=1.-.822/(1.+std::exp(-(alp-.2)*1.15));
     G4cout<<"Test19:EndOfRun,p="<<pmax<<",dT="<<dTot<<",r="<<dEl/dTot<<",e="<<exr
           <<",d="<<exr-dEl/dTot<<",ra="<<(exr-.178)/.822<<G4endl;
     for(G4int ir=0; ir<nT; ir++) G4cout<<tVal[ir]<<" "<<tSig[ir]<<G4endl;// Print t-vectors
#endif
#ifdef pverb
     G4cout<<"Test19: ########## End of run ##########"<<G4endl;
#endif
     delete step;   // The G4Step delets aPoint and bPoint (can be necessary in Loop)
    } // End of the target LOOP
   } // End of the Energy LOOP
   delete gTrack; // The G4Track delets the G4DynamicParticle
  } // End of the projectile LOOP
  delete proc;
  //delete man;      // == Should not be uncommented unless definition above is commented!
  //delete pFrame;   // This
  //delete lFrame;   // is
  //delete sFrame;   // not
  //delete material; // necessary
  //delete element;  // -> (automaticaly
  //delete isotope;  // in the end of main)
  //G4cout << "Test19: After delete process etc." << G4endl;
  //delete phys;                        // Absolete physics class
  //delete cmaterial;                   // Temporary material definition pointer (?)
#ifndef nout
  //G4cout << "Test19: Before ntp" << G4endl;
  delete ntp; // Delete the class to fill the#of events
  //G4cout << "Test19: After ntp" << G4endl;
#endif
#ifdef hdebug
  timer->Stop();
  G4int    nbnh=G4StringChipsParticleLevelInterface::GetNbn();
  G4double hbdb=G4StringChipsParticleLevelInterface::GetDB();
  //G4double hede=G4StringChipsParticleLevelInterface::GetDE(); // The same now
  G4double toth=G4StringChipsParticleLevelInterface::GetTot();
  G4double ovrb=G4StringChipsParticleLevelInterface::GetBov();
  G4double ovre=G4StringChipsParticleLevelInterface::GetEov();
  G4cout<<"Test19:TimePerEvent="<<timer->GetUserElapsed()/toth<<", N="<<toth
        <<", overB="<<ovrb/toth<<", overE="<<ovre/toth<<G4endl;
  delete timer;
  G4double bc=-hbdb/2;
  for(G4int ih=0; ih<nbnh; ih++)
  {
    bc+=hbdb;
    G4double nE=G4StringChipsParticleLevelInterface::GetE(ih);
    G4double nB=G4StringChipsParticleLevelInterface::GetB(ih);
    G4cout<<"Test19:be="<<bc<<" "<<nE/toth<<" "<<std::sqrt(nE)/toth
          <<" "<<nB/toth/bc<<" "<<std::sqrt(nB)/toth/bc<<G4endl;
  }
#endif
  delete runManager;
  delete cond;
#ifdef meandbg
  G4cout<<"E="<<pE/nE<<":sE/E="<<sE/pE<<",dE/E="<<std::sqrt(sE2*nE-sE*sE)/pE
        <<",D%="<<100*nD/nE<<",sE="<<sE<<",sE2="<<sE2<<",nE="<<nE<<G4endl;
#endif
#ifdef ekindbg
  G4cout<<"Ei="<<pE/nE<<": Ef="<<sE/nE<<", dE="<<std::sqrt(sE2*nE-sE*sE)/nE<<", D%="
        <<100*nD/nE<<",sE="<<sE<<",sE2="<<sE2<<",nE="<<nE<<G4endl;
#endif
#ifdef histdbg
  G4cout<<std::setw(2)<<"#"<<" : "<<std::setw(9)<<"#gamma"<<std::setw(9)<<"SumEGam"
        <<std::setw(9)<<"dBrN+20"<<std::setw(9)<<"dChg+20" <<std::setw(9)<<"rE(5%)"
    //<<std::setw(9)<<"rPz(5%)"<<std::setw(9)<<"rPr(5%)" <<std::setw(9)<<"dE(.1)"
        <<std::setw(9)<<"rPz(5%)"<<std::setw(9)<<"rPr(5%)" <<std::setw(9)<<"dE(10.)"
        <<std::setw(9)<<"dPz(.1)"<<std::setw(9)<<"dPr(.01)"<<G4endl;
  for(G4int jh=0; jh<nHst; jh++)
  {
    G4cout<<std::setw(2)<<jh<<" : "<<std::setw(9)<<nGaH[jh]<<std::setw(9)<<GaSE[jh]
          <<std::setw(9)<<dBnH[jh]<<std::setw(9)<<dChH[jh]<<std::setw(9)<<rEnH[jh]
          <<std::setw(9)<<rPzH[jh]<<std::setw(9)<<rPrH[jh]<<std::setw(9)<<dEnH[jh]
          <<std::setw(9)<<dPzH[jh]<<std::setw(9)<<dPrH[jh]<<G4endl;
  }
#endif
#ifdef synch
  s_s/=s_n;
  s_d=std::sqrt(s_d/s_n-s_s*s_s);
  G4cout<<"SynchroRadiation: Mean="<<s_s<<", RMS="<<s_d<<", nEv="<<s_n<<G4endl;
  G4double low=minom;
  G4double lgx=logmin;
  for(G4int i=0; i<nBin; ++i)
  {
    lgx+=dlog;
    G4double hx=std::exp(lgx);
    G4double dx=hx-low;
    G4double z=hBin[i];
    G4double x=(hx+low)/2.;
    G4double y=0.;
    if(z>0.)
    {
      x=eBin[i]/z;
      y=z*x/dx;
    }
    G4double d=std::sqrt(z+1.)*x/dx;
    G4cout<<x<<" "<<y<<" "<<d<<G4endl;
    low=hx;
  }
#endif
#ifdef pverb
  G4cout<<"###### End of Test19 #####"<<G4endl;
#endif
  //exit(1); // Never do this !
  //return no_of_errors;
  //return 0;
  //abort();
  //return EXIT_SUCCESS;
}
