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
//      File name:     Test49 (Comparison of G4QString Results with Data)
//
//      Author:        M.Kossov (made of test19)
// 
//      Creation date: 26 Aug 2009
//
//      Modifications: G4electroweak CHIPS interface (G4QCaptureAtRest)
//                     which is independent of the "hadronic" package
//                     classes is added to this test for the debugging,
//                     physics tests and development of the CHIPS model.
//
// -------------------------------------------------------------------
//       1         2         3         4         5         6         7         8         9
//34567890123456789012345678901234567890123456789012345678901234567890123456789012345678901

#define devel
//#define pdebug
//#define debug
//#define inter
//#define pscan
//#define csdebug
//#define smear
//#define escan
//#define idebug
//#define tdebug
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

#include "G4QInelastic.hh"
#include "G4QDiffraction.hh"
#include "G4QLowEnergy.hh"
#include "G4TheoFSGenerator.hh"
#include "G4StringChipsParticleLevelInterface.hh"
#include "G4ExcitationHandler.hh"
#include "G4PreCompoundModel.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4QGSModel.hh"
#include "G4FTFModel.hh"
#include "G4QGSParticipants.hh"
#include "G4QGSMFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4QuasiElasticChannel.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4AntiProtonInelasticProcess.hh"
#include "G4ProtonInelasticCrossSection.hh"
#include "G4NeutronInelasticCrossSection.hh"
#include "G4HadronInelasticDataSet.hh"
#include "G4PionPlusInelasticProcess.hh"
#include "G4PiNuclearCrossSection.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4KaonPlusInelasticProcess.hh"
#include "G4KaonMinusInelasticProcess.hh"
//Does not work with shared libraries//#include "G4MuNuclearInteraction.hh"
#include "G4QMuonNuclearCrossSection.hh"
#include "G4PiNuclearCrossSection.hh"
#include "G4CrossSectionDataStore.hh"
#include "G4HadronElasticDataSet.hh"
#include "G4HadronInelasticDataSet.hh"

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

#include "Test49Physics.hh"
#include "Test49PhysicsList.hh"

#include "G4HadronCrossSections.hh"
#include "G4VCrossSectionDataSet.hh"
#include "G4ProtonInelasticCrossSection.hh"
#include "G4NeutronInelasticCrossSection.hh"
#include "G4HadronInelasticDataSet.hh"

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

// Temporary AIDA can be not used in test49
//#include "AIDA/AIDA.h"

// ASCII Pseudo NTUPLE 
#include "G4QHBook.hh"
#include "G4Timer.hh"
#include "time.h"

//int main(int argc, char** argv)
int main()
{
#ifdef csdebug
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
  runManager->SetUserInitialization(new Test49PhysicsList);
  G4StateManager::GetStateManager()->SetNewState(G4State_Init); // To let create ions
#ifdef debug
  G4cout<<"Test49: Prepare G4QHBook files or ntuples"<<G4endl;
#endif
  G4QHBook* ntp = new G4QHBook;
  //-------- set standard output format-------
  G4cout.setf( std::ios::scientific, std::ios::floatfield );
  //-------- take parameters from a file -------
  std::ifstream inFile("chipstest.in", std::ios::in);
  G4double temp;
  G4double ssin2g;
  G4double eteps;
  G4double momb;
  G4double enb;
  G4double cP;
  G4double fN;
  G4double fD;
  G4double rM;
  G4double sA;
  G4String mName;
  G4int    nop;
  G4int    pPDG;
  G4int    tPDG;
  G4int    nEvt;
  //G4int    nofdecays;
  //G4int    decmask=0;
  inFile>>temp>>ssin2g>>eteps>>nop>>momb>>enb>>pPDG>>tPDG>>nEvt>>fN>>fD>>cP>>rM>>sA>>mName;
  G4cout<<"Test49:Par's: T="<<temp<<",S="<<ssin2g<<",Eta="<<eteps<<",nop="<<nop<<",p="
        <<momb<<",e="<<enb<<",pPDG="<<pPDG<<",tPDG="<<tPDG<<",nEv="<<nEvt<<",fN="<<fN
        <<",fD="<<fD<<",cP="<<cP<<",rM="<<rM<<",sA="<<sA<<",modName="<<mName<<G4endl;
  //-------- Initialize CHIPS
  G4QCHIPSWorld* theW=G4QCHIPSWorld::Get();
  theW->GetParticles(nop);           // Create CHIPS World of nop particles
  //G4Exception("***CHIPStest: TMP");
  G4QInelastic::SetParameters(temp,ssin2g,eteps,fN,fD,cP,rM,nop,sA);
  G4QInelastic::SetManual();
  // ********** Now momb is a momentum of the incident particle, if =0 => LOOP ************
//#ifdef chips
//#else
  // Not necessary for CHIPS, only for GHAD: fake force condition
  G4ForceCondition* cond = new G4ForceCondition;
  *cond=NotForced;
//#endif
  if     (pPDG == 1000010010) pPDG = 2212;
  else if(pPDG == 1000000010) pPDG = 2112;
  else if(pPDG  > 1000000000)
  {
    pPDG = pPDG/10 - 10000000;
    pPDG-= (pPDG-90000000)/1000;
  }
#ifdef debug
  G4cout<<"Test49: proj_CHIPS_PDG = "<<pPDG<<G4endl;
#endif
  G4double mp=G4QPDGCode(pPDG).GetMass();
  G4double ep=mp;
  G4int cnE=1;
  if(momb==0.) cnE=nEn;
  else
  {
    ep=std::sqrt(mp*mp+momb*momb);
    if(enb>0.) ep=enb;
  }
  tPDG = tPDG/10 - 10000000;
  G4int tgZ = (tPDG-90000000)/1000;
  G4int tgN = tPDG-90000000-tgZ*1001;
  tPDG-= tgZ;
#ifdef debug
  G4cout<<"Test49: targ_CHIPS_PDG = "<<tPDG<<", tgZ = "<<tgZ<<", tgN = "<<tgN<<G4endl;
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
  G4ThreeVector aDirection = G4ThreeVector(0.,0.,1.); // @@ Temporary "not along Z"
  G4int tgm=1;                                        // By default only one target
  if(!tPDG) // Make max for the LOOP over all targets and define materials
  {
#ifdef debug
    G4cout<<"Test49: targetPDG="<<tPDG<<", tgm = "<<tgm<<", nTg = "<<nTg<<G4endl;
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
      G4cout<<"Test49:-->Material("<<tgZ<<","<<tgN<<"):"<<tnm[tgi]<<","<<tsy[tgi]<<G4endl;
      mat[tgi]=material;
    }
  }
  else
  {
    isotope = new G4Isotope("Isotop", tgZ, tgA);
    element = new G4Element("ZA_Isotop", "ZA", 1);
    element->AddIsotope(isotope, 100.*perCent);
    material = new G4Material("ZA_Isomer", 1.*g/cm3, 1);
    material->AddElement(element, 1);
  }
#ifdef debug
  G4cout<<"Test49:--- Material "<<material->GetName()<<" is defined ---" << G4endl;
#endif
  if(!material)
  {
    G4cout<<"Test49: Last Material "<<material->GetName()<<" is not defined."<<G4endl;
    exit(1);
  }
  // ---------------------- PSEUDO-TRACK DEFINISIONS -------------------------
  G4double nx = 0.;
  G4double ny = 0.;
  G4double nz = 0.;
  G4int npart=1;                                       // By default only one particle
  if(!pPDG) npart=nPr;                                 // Make a LOOP ove all particles
  G4VDiscreteProcess* proc = 0;                        // Process to be simulated
  G4QInelastic* cproc = 0;                             // CHIPS
  G4HadronInelasticProcess*     hproc   = 0;           // GHAD
  G4ProtonInelasticProcess*     prhproc = 0;
  G4NeutronInelasticProcess*    nehproc = 0;
  G4PionMinusInelasticProcess*  pmhproc = 0;
  G4PionPlusInelasticProcess*   pphproc = 0;
  G4KaonMinusInelasticProcess*  kmhproc = 0;
  G4KaonPlusInelasticProcess*   kphproc = 0;
  G4AntiProtonInelasticProcess* aphproc = 0;
  if(mName == "chips")  // ************** CHIPS process definition starts here ************
  {
    cproc = new G4QInelastic;
#ifdef devel
    cproc->SetParameters(temp, ssin2g, eteps, fN, fD, cP, rM, nop, sA); // Development
#endif
    ///G4QDiffraction* proc = new G4QDiffraction;           // A diffraction CHIPS process
    ///G4QLowEnergy* proc = new G4QLowEnergy;               // fragment-nucleus universal
    proc=cproc;
  }
  else
  { // **************** GHAD process definition starts here *******************************
    G4HadronicInteraction*      aModel = 0;
    G4PreCompoundModel*      precModel = 0;
    G4CascadeInterface*      bertModel = 0;
    G4BinaryCascade*         binaModel = 0;
    G4LEProtonInelastic*     lephModel = 0;
    G4LENeutronInelastic*    lenhModel = 0;
    G4LEPionMinusInelastic*  lepmModel = 0;
    G4LEPionPlusInelastic*   leppModel = 0;
    G4LEKaonMinusInelastic*  lekmModel = 0;
    G4LEKaonPlusInelastic*   lekpModel = 0;
    G4LEAntiProtonInelastic* leapModel = 0;
    G4HEProtonInelastic*     hephModel = 0;
    G4HENeutronInelastic*    henhModel = 0;
    G4HEPionMinusInelastic*  hepmModel = 0;
    G4HEPionPlusInelastic*   heppModel = 0;
    G4HEKaonMinusInelastic*  hekmModel = 0;
    G4HEKaonPlusInelastic*   hekpModel = 0;
    G4HEAntiProtonInelastic* heapModel = 0;
    G4TheoFSGenerator*       hModel    = 0;
    G4VCrossSectionDataSet*  theCS     = 0;
    if     (mName == "preco")
    {
      precModel = new G4PreCompoundModel(new G4ExcitationHandler);
#ifdef debug
      G4cout<<"Test49: Preco Model is defined ="<<precModel<<G4endl;
#endif
      aModel = precModel;
    }
    else if(mName == "bertini")
    {
      bertModel = new G4CascadeInterface;
#ifdef debug
      G4cout<<"Test49: Preco Model is defined ="<<bertModel<<G4endl;
#endif
      aModel = bertModel;
    }
    else if(mName == "binary")
    {
      binaModel = new G4BinaryCascade;
#ifdef debug
      G4cout<<"Test49: Preco Model is defined ="<<binaModel<<G4endl;
#endif
      aModel = binaModel;
    }
    else if(mName == "lep")
    {
      if     (pPDG==2212)
      {
        lephModel = new G4LEProtonInelastic;
#ifdef debug
        G4cout<<"Test49: LEP p Model is defined ="<<lephModel<<G4endl;
#endif
        aModel = lephModel;
      }
      else if(pPDG==2112)
      {
        lenhModel = new G4LENeutronInelastic;
#ifdef debug
        G4cout<<"Test49: LEP n Model is defined ="<<lenhModel<<G4endl;
#endif
        aModel = lenhModel;
      }
      else if(pPDG==-211)
      {
        lepmModel = new G4LEPionMinusInelastic;
#ifdef debug
        G4cout<<"Test49: LEP pim Model is defined ="<<lepmModel<<G4endl;
#endif
        aModel = lepmModel;
      }
      else if(pPDG== 211)
      {
        leppModel = new G4LEPionPlusInelastic;
#ifdef debug
        G4cout<<"Test49: LEP pip Model is defined ="<<leppModel<<G4endl;
#endif
        aModel = leppModel;
      }
      else if(pPDG==-321)
      {
        lekmModel = new G4LEKaonMinusInelastic;
#ifdef debug
        G4cout<<"Test49: LEP km Model is defined ="<<lekmModel<<G4endl;
#endif
        aModel = lekmModel;
      }
      else if(pPDG== 321)
      {
        lekpModel = new G4LEKaonPlusInelastic;
#ifdef debug
        G4cout<<"Test49: LEP kp Model is defined ="<<lekpModel<<G4endl;
#endif
        aModel = lekpModel;
      }
      else if(pPDG==-2112)
      {
        leapModel = new G4LEAntiProtonInelastic;
#ifdef debug
        G4cout<<"Test49: LEP ap Model is defined ="<<leapModel<<G4endl;
#endif
        aModel = leapModel;
      }
      else
      {
        G4cout<<"-Error-Test49: Process is not defined in LEP for PDG="<<pPDG<<G4endl;
        aModel = new G4LEProtonInelastic;
      }
    }
    else if(mName == "hep")
    {
      if     (pPDG==2212)
      {
        hephModel = new G4HEProtonInelastic;
#ifdef debug
        G4cout<<"Test49: HEP p Model is defined ="<<hephModel<<G4endl;
#endif
        aModel = hephModel;
      }
      else if(pPDG==2112)
      {
        henhModel = new G4HENeutronInelastic;
#ifdef debug
        G4cout<<"Test49: HEP n Model is defined ="<<henhModel<<G4endl;
#endif
        aModel = henhModel;
      }
      else if(pPDG==-211)
      {
        hepmModel = new G4HEPionMinusInelastic;
#ifdef debug
        G4cout<<"Test49: HEP pim Model is defined ="<<hepmModel<<G4endl;
#endif
        aModel = hepmModel;
      }
      else if(pPDG== 211)
      {
        heppModel = new G4HEPionPlusInelastic;
#ifdef debug
        G4cout<<"Test49: HEP pip Model is defined ="<<heppModel<<G4endl;
#endif
        aModel = heppModel;
      }
      else if(pPDG==-321)
      {
        hekmModel = new G4HEKaonMinusInelastic;
#ifdef debug
        G4cout<<"Test49: HEP km Model is defined ="<<hekmModel<<G4endl;
#endif
        aModel = hekmModel;
      }
      else if(pPDG== 321)
      {
        hekpModel = new G4HEKaonPlusInelastic;
#ifdef debug
        G4cout<<"Test49: HEP kp Model is defined ="<<hekpModel<<G4endl;
#endif
        aModel = hekpModel;
      }
      else if(pPDG==-2112)
      {
        heapModel = new G4HEAntiProtonInelastic;
#ifdef debug
        G4cout<<"Test49: HEP ap Model is defined ="<<heapModel<<G4endl;
#endif
        aModel = heapModel;
      }
      else
      {
        aModel = new G4HEProtonInelastic;
        G4cout<<"-Error-Test49: Process is not defined in HEP for PDG="<<pPDG<<G4endl;
      }
    }
    else if(mName == "qgsc")
    {
      hModel = new G4TheoFSGenerator("QGSC"); // The same for QGS & FTF
#ifdef debug
      G4cout<<"Test49: QGSC Model is defined ="<<hModel<<G4endl;
#endif
      G4QGSModel<G4QGSParticipants>* aStringModel = new G4QGSModel<G4QGSParticipants>;
      G4QGSMFragmentation* aFragmentation = new G4QGSMFragmentation;
      G4ExcitedStringDecay* aStringDecay = new G4ExcitedStringDecay(aFragmentation);
      aStringModel->SetFragmentationModel(aStringDecay);

      G4StringChipsParticleLevelInterface*theCHIPS=new G4StringChipsParticleLevelInterface;

      hModel->SetTransport(theCHIPS);
      hModel->SetHighEnergyGenerator(aStringModel);

      G4QuasiElasticChannel* theQuasiElastic = new G4QuasiElasticChannel;
      hModel->SetQuasiElasticChannel(theQuasiElastic);

      G4ProjectileDiffractiveChannel* theProjectileDiffraction =
                                                       new G4ProjectileDiffractiveChannel;
      hModel->SetProjectileDiffraction(theProjectileDiffraction);

      hModel->SetMinEnergy(0.);
      hModel->SetMaxEnergy(100*TeV);
      aModel = hModel;
#ifdef debug
      G4cout<<"Test49: model="<<mName<<"(QGSC) is defined"<<G4endl;
#endif
    }
    else if(mName == "qgsp")
    {
      hModel = new G4TheoFSGenerator("QGSP"); // The same for QGS & FTF
#ifdef debug
      G4cout<<"Test49: QGSP Model is defined ="<<hModel<<G4endl;
#endif
      G4QGSModel<G4QGSParticipants>* aStringModel = new G4QGSModel<G4QGSParticipants>;
      G4QGSMFragmentation* aFragmentation = new G4QGSMFragmentation;
      G4ExcitedStringDecay* aStringDecay = new G4ExcitedStringDecay(aFragmentation);
      aStringModel->SetFragmentationModel(aStringDecay);

      G4GeneratorPrecompoundInterface* theCascade = new G4GeneratorPrecompoundInterface;
      G4PreCompoundModel* thePreEquilib = new G4PreCompoundModel(new G4ExcitationHandler);
      theCascade->SetDeExcitation(thePreEquilib);  

      hModel->SetTransport(theCascade);
      hModel->SetHighEnergyGenerator(aStringModel);

      G4QuasiElasticChannel* theQuasiElastic = new G4QuasiElasticChannel;
      hModel->SetQuasiElasticChannel(theQuasiElastic);

      G4ProjectileDiffractiveChannel* theProjectileDiffraction =
                                                       new G4ProjectileDiffractiveChannel;
      hModel->SetProjectileDiffraction(theProjectileDiffraction);

      hModel->SetMinEnergy(0.);
      hModel->SetMaxEnergy(100*TeV);
      aModel = hModel;
#ifdef debug
      G4cout<<"Test49: model="<<mName<<"(QGSP) is defined"<<G4endl;
#endif
    }
    else if(mName == "ftfc")
    {
      hModel = new G4TheoFSGenerator("FTFC"); // The same for QGS & FTF
#ifdef debug
      G4cout<<"Test49: FTFC Model is defined ="<<hModel<<G4endl;
#endif
      G4FTFModel* aStringModel = new G4FTFModel;
      G4LundStringFragmentation* aFragmentation = new G4LundStringFragmentation;
      G4ExcitedStringDecay* aStringDecay = new G4ExcitedStringDecay(aFragmentation);
      aStringModel->SetFragmentationModel(aStringDecay);

      G4StringChipsParticleLevelInterface*theCHIPS=new G4StringChipsParticleLevelInterface;

      hModel->SetTransport(theCHIPS);
      hModel->SetHighEnergyGenerator(aStringModel);

      G4QuasiElasticChannel* theQuasiElastic = new G4QuasiElasticChannel;
      hModel->SetQuasiElasticChannel(theQuasiElastic);

      G4ProjectileDiffractiveChannel* theProjectileDiffraction =
                                                       new G4ProjectileDiffractiveChannel;
      hModel->SetProjectileDiffraction(theProjectileDiffraction);

      hModel->SetMinEnergy(0.);
      hModel->SetMaxEnergy(100*TeV);
      aModel = hModel;
#ifdef debug
      G4cout<<"Test49: model="<<mName<<"(FTFC) is defined"<<G4endl;
#endif
    }
    else if(mName == "ftfp")
    {
      hModel = new G4TheoFSGenerator("FTFP"); // The same for QGS & FTF
#ifdef debug
      G4cout<<"Test49: FTFP Model is defined ="<<hModel<<G4endl;
#endif
      G4FTFModel* aStringModel = new G4FTFModel;
      G4LundStringFragmentation* aFragmentation = new G4LundStringFragmentation;
      G4ExcitedStringDecay* aStringDecay = new G4ExcitedStringDecay(aFragmentation);
      aStringModel->SetFragmentationModel(aStringDecay);

      G4GeneratorPrecompoundInterface* theCascade = new G4GeneratorPrecompoundInterface;
      G4PreCompoundModel* thePreEquilib = new G4PreCompoundModel(new G4ExcitationHandler);
      theCascade->SetDeExcitation(thePreEquilib);  

      hModel->SetTransport(theCascade);
      hModel->SetHighEnergyGenerator(aStringModel);

      G4QuasiElasticChannel* theQuasiElastic = new G4QuasiElasticChannel;
      hModel->SetQuasiElasticChannel(theQuasiElastic);

      G4ProjectileDiffractiveChannel* theProjectileDiffraction =
                                                       new G4ProjectileDiffractiveChannel;
      hModel->SetProjectileDiffraction(theProjectileDiffraction);

      hModel->SetMinEnergy(0.);
      hModel->SetMaxEnergy(100*TeV);
      aModel = hModel;
#ifdef debug
      G4cout<<"Test49: model="<<mName<<"(FTFP) is defined"<<G4endl;
#endif
    }
    else G4cout<<"-Error-Test49: The model name = "<<mName<<" isn't defined"<<G4endl;
    if     (pPDG == 2212)
    {
      prhproc = new G4ProtonInelasticProcess;
      hproc   = prhproc;
      proc   = hproc;
      if     (mName == "preco")   prhproc->RegisterMe(precModel);
      else if(mName == "bertini") prhproc->RegisterMe(bertModel);
      else if(mName == "binary")  prhproc->RegisterMe(binaModel);
      else if(mName == "lep")     prhproc->RegisterMe(lephModel);
      else if(mName == "hep")     prhproc->RegisterMe(hephModel);
      else if(mName=="qgsc" || mName=="qgsp" || mName=="ftfc" || mName=="ftfp")
                                  prhproc->RegisterMe(hModel);
      else G4cout<<"Test49: No Proton Model for the Name ="<<mName<<G4endl;
      theCS = new G4ProtonInelasticCrossSection;
      hproc->AddDataSet(theCS);   // Can not be skipped for the GHAD event generator
#ifdef debug
      G4cout<<"Test49: Process is defined for proton Model="<<mName<<G4endl;
#endif
    }
    else if(pPDG== 2112)
    {
      nehproc = new G4NeutronInelasticProcess;
      hproc   = nehproc;
      proc   = hproc;
      if     (mName == "preco")   nehproc->RegisterMe(precModel);
      else if(mName == "bertini") nehproc->RegisterMe(bertModel);
      else if(mName == "binary")  nehproc->RegisterMe(binaModel);
      else if(mName == "lep")     nehproc->RegisterMe(lenhModel);
      else if(mName == "hep")     nehproc->RegisterMe(henhModel);
      else if(mName=="qgsc" || mName=="qgsp" || mName=="ftfc" || mName=="ftfp")
                                  nehproc->RegisterMe(hModel);
      else G4cout<<"Test49: No PionPlus Model for the Name ="<<mName<<G4endl;
      theCS = new G4NeutronInelasticCrossSection;
      nehproc->AddDataSet(theCS);   // Can not be skipped for the GHAD event generator
#ifdef debug
      G4cout<<"Test49: Process is defined for neutron PDG="<<pPDG<<G4endl;
#endif
    }
    else if(pPDG== -211)
    {
      pmhproc = new G4PionMinusInelasticProcess;
      hproc   = pmhproc;
      proc   = hproc;
      if     (mName == "preco")   pmhproc->RegisterMe(precModel);
      else if(mName == "bertini") pmhproc->RegisterMe(bertModel);
      else if(mName == "binary")  pmhproc->RegisterMe(binaModel);
      else if(mName == "lep")     pmhproc->RegisterMe(lepmModel);
      else if(mName == "hep")     pmhproc->RegisterMe(hepmModel);
      else if(mName=="qgsc" || mName=="qgsp" || mName=="ftfc" || mName=="ftfp")
                                  pmhproc->RegisterMe(hModel);
      else G4cout<<"Test49: No PionMinus Model for the Name ="<<mName<<G4endl;
      theCS = new G4PiNuclearCrossSection; // The same for both charges of pions
      pmhproc->AddDataSet(theCS);   // Can not be skipped for the GHAD event generator
#ifdef debug
      G4cout<<"Test49: Process is defined for piminus PDG="<<pPDG<<G4endl;
#endif
    }
    else if(pPDG==  211)
    {
      pphproc = new G4PionPlusInelasticProcess;
      proc   = pphproc;
      if     (mName == "preco")   pphproc->RegisterMe(precModel);
      else if(mName == "bertini") pphproc->RegisterMe(bertModel);
      else if(mName == "binary")  pphproc->RegisterMe(binaModel);
      else if(mName == "lep")     pphproc->RegisterMe(leppModel);
      else if(mName == "hep")     pphproc->RegisterMe(heppModel);
      else if(mName=="qgsc" || mName=="qgsp" || mName=="ftfc" || mName=="ftfp")
                                  pphproc->RegisterMe(hModel);
      else G4cout<<"Test49: No PionPlus Model for the Name ="<<mName<<G4endl;
      theCS = new G4PiNuclearCrossSection; // The same for both charges of pions
      pmhproc->AddDataSet(theCS);   // Can not be skipped for the GHAD event generator
#ifdef debug
      G4cout<<"Test49: Process is defined for piplus PDG="<<pPDG<<G4endl;
#endif
    }
    else if(pPDG== -321)
    {
      kmhproc = new G4KaonMinusInelasticProcess;
      hproc   = kmhproc;
      proc   = hproc;
      if     (mName == "preco")   kmhproc->RegisterMe(precModel);
      else if(mName == "bertini") kmhproc->RegisterMe(bertModel);
      else if(mName == "binary")  kmhproc->RegisterMe(binaModel);
      else if(mName == "lep")     kmhproc->RegisterMe(lekmModel);
      else if(mName == "hep")     kmhproc->RegisterMe(hekmModel);
      else if(mName=="qgsc" || mName=="qgsp" || mName=="ftfc" || mName=="ftfp")
                                  kmhproc->RegisterMe(hModel);
      else G4cout<<"Test49: No KaonMinus Model for the Name ="<<mName<<G4endl;
      theCS = new G4HadronInelasticDataSet; // For kaons, anti-protons, and other hadrons
      kmhproc->AddDataSet(theCS);   // Can not be skipped for the GHAD event generator
#ifdef debug
      G4cout<<"Test49: Process is defined for kaminus PDG="<<pPDG<<G4endl;
#endif
    }
    else if(pPDG==  321)
    {
      kphproc = new G4KaonPlusInelasticProcess;
      proc   = kphproc;
      if     (mName == "preco")   kphproc->RegisterMe(precModel);
      else if(mName == "bertini") kphproc->RegisterMe(bertModel);
      else if(mName == "binary")  kphproc->RegisterMe(binaModel);
      else if(mName == "lep")     kphproc->RegisterMe(lekpModel);
      else if(mName == "hep")     kphproc->RegisterMe(hekpModel);
      else if(mName=="qgsc" || mName=="qgsp" || mName=="ftfc" || mName=="ftfp")
                                  kphproc->RegisterMe(hModel);
      else G4cout<<"Test49: No KaonPlus Model for the Name ="<<mName<<G4endl;
      theCS = new G4HadronInelasticDataSet; // For kaons, anti-protons, and other hadrons
      kphproc->AddDataSet(theCS);   // Can not be skipped for the GHAD event generator
#ifdef debug
      G4cout<<"Test49: Process is defined for kaplus PDG="<<pPDG<<G4endl;
#endif
    }
    else if(pPDG==-2112)
    {
      aphproc = new G4AntiProtonInelasticProcess;
      hproc   = prhproc;
      proc   = hproc;
      if     (mName == "preco")   aphproc->RegisterMe(precModel);
      else if(mName == "bertini") aphproc->RegisterMe(bertModel);
      else if(mName == "binary")  aphproc->RegisterMe(binaModel);
      else if(mName == "lep")     aphproc->RegisterMe(leapModel);
      else if(mName == "hep")     aphproc->RegisterMe(heapModel);
      else if(mName=="qgsc" || mName=="qgsp" || mName=="ftfc" || mName=="ftfp")
                                  aphproc->RegisterMe(hModel);
      else G4cout<<"Test49: No AntiProton Model for the Name ="<<mName<<G4endl;
      theCS = new G4HadronInelasticDataSet; // For kaons, anti-protons, and other hadrons
      aphproc->AddDataSet(theCS);   // Can not be skipped for the GHAD event generator
#ifdef debug
      G4cout<<"Test49: Process is defined for anti-proton PDG="<<pPDG<<G4endl;
#endif
    }
    else G4cout<<"-Error-Test49: Process is not defined for PDG="<<pPDG<<G4endl;
#ifdef debug
    G4cout<<"Test49: GHAD proc="<<proc<<" (model="<<mName<<") is defined"<<G4endl;
#endif
  }
#ifdef debug
  G4cout<<"Test49:--***-- process is created --***--, proc="<<proc<<G4endl; // only one run
#endif
#ifdef hdebug
  G4Timer* timer = new G4Timer(); // ??
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
   else if(pPDG==15) part=G4TauMinus::TauMinus();    // Definition of tau- projectile
   else if(pPDG==-15) part=G4TauPlus::TauPlus();     // Definition of tau+ projectile
   else if(pPDG==12) part=G4NeutrinoE::NeutrinoE();  // Definition of e_neutrino projectile
   else if(pPDG==-12) part=G4AntiNeutrinoE::AntiNeutrinoE(); // anti-e_neutrino projectile
   else if(pPDG==14) part=G4NeutrinoMu::NeutrinoMu();  // Definition of mu-neutrino proj.
   else if(pPDG==-14) part=G4AntiNeutrinoMu::AntiNeutrinoMu(); // anti-mu_neutrino proj.
   //else if(pPDG==16) part=G4NeutrinoTau::NeutrinoTau(); // Definition of tau-neutrino
   //else if(pPDG==-16) part=G4AntiNeutrinoTau::AntiNeutrinoTau(); // anti-tau_neutrino
   else if(pPDG==2112) part=G4Neutron::Neutron();    // Definition of Neutron projectile
   else if(pPDG==-211) part=G4PionMinus::PionMinus();// Definition of Pi- projectile
   else if(pPDG==211)  part=G4PionPlus::PionPlus();  // Definition of Pi+ projectile
   else if(pPDG==321)  part=G4KaonPlus::KaonPlus();  // Definition of K+ projectile
   else if(pPDG==-321) part=G4KaonMinus::KaonMinus();// Definition of K- projectile
   else if(pPDG==130)part=G4KaonZeroLong::KaonZeroLong();//Definition of K0_L projectile
   else if(pPDG==310)part=G4KaonZeroShort::KaonZeroShort();//Definition of K0S projectile
   else if(pPDG==-2212)part=G4AntiProton::AntiProton();//Define antiProton projectile
   else if(pPDG==-2112)part=G4AntiNeutron::AntiNeutron();// Define AntiNeutron projectile
   else if(pPDG==3122) part=G4Lambda::Lambda();      // Definition of Lambda projectile
   else if(pPDG==3112) part=G4SigmaMinus::SigmaMinus();// Definition of Sigma- projectile
   else if(pPDG==3222) part=G4SigmaPlus::SigmaPlus();// Definition of Sigma+ projectile
   else if(pPDG==3322) part=G4XiMinus::XiMinus();    // Definition of the Xi- projectile
   else if(pPDG==3312) part=G4XiZero::XiZero();      // Definition of the Xi- projectile
   else if(pPDG==3334) part=G4OmegaMinus::OmegaMinus(); // Definition of Omega- projectile
   else if(pPDG==2112) part=G4Neutron::Neutron();    // Definition of Neutron projectile
   else if(pPDG!=-3222) // Leave defaulf definition for Anti Sigma+ projectile
   {
     G4cerr<<"***Test49: "<<pPDG<<" is a PDG code of not supported particle"<<G4endl;
     G4Exception("***Test49: OnFlight Process is called for not supported particle");
   }
   G4double pMass = part->GetPDGMass();                 // Mass of the projectile
   //
   G4ThreeVector aPosition(nx*mm, ny*mm, nz*mm);
   G4DynamicParticle* dParticle = new G4DynamicParticle(part,aDirection,0.);// Dummy Energy

   //G4cout<<"Test49: before the new Track definition"<<G4endl;
   // General Track definition
   G4Track* gTrack = new G4Track(dParticle,aTime,aPosition); // Track definition (NoTarg)

   G4int    bnp=pQC.GetBaryonNumber();
   // ---------- Define material for the simulation ------------------
   //G4double tgA     = 26.98; // @@ Important? Can it be just tgZ+tgN?
   G4int tgA        = tgZ+tgN; // @@ Temporary, not good
   G4String nameMat = "Thin Target";  // @@ Not important can be an arbitrary name
   // The material can be copied from the commented cMaterial Factory above
   G4Isotope* isotope = new G4Isotope("Isotop", tgZ, tgA);
   G4Element* element = new G4Element("ZA_Isotop", "ZA", 1);
   element->AddIsotope(isotope, 100.*perCent);
   G4Material* material = new G4Material("ZA_Isomer", 1.*g/cm3, 1);
   material->AddElement(element, 1);
#ifdef debug
   G4cout<<"Test49:--- Material is defined ---" << G4endl; // only one run
#endif
   // For the Material Factory case
   //G4String  nameMat  = "Si";
   //G4Material* material = G4Material::GetMaterial(nameMat); // Get definition of Material
   if(!material)
   {
     G4cout<<"Test49:Material "<<nameMat<<" is not defined in the Test49Material"<<G4endl;
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
   G4cout<<"Test49:###### Start new run #####" << G4endl; // only one run
#endif
#ifdef debug
   G4double totCN  = tgZ+part->GetPDGCharge();
   G4double totSN  = tgZ+part->GetPDGCharge();
   G4int    totBNN = tgZ+tgN+part->GetBaryonNumber();
#endif
   G4double theStep   = 0.01*micrometer;
   // G4int maxz = (G4int)((*(material->GetElementVector()))[0]->GetZ()) + 1;
#ifdef pdebug
   G4int targPDG=90000000+1000*tgZ+tgN;                 // PDG Code of the target
   G4cout<<"Test49: The particle: "<<part->GetParticleName()<<G4endl;
   G4cout<<"Test49: The material: "<<material->GetName()<<G4endl;
   G4cout<<"Test49: The target:   "<<targPDG<<G4endl;
   G4cout<<"Test49: The step:     "<<theStep/mm<<" mm"<<G4endl;
   G4cout<<"Test49: The position: "<<aPosition/mm<<" mm"<<G4endl;
   G4cout<<"Test49: The direction:"<<aDirection<<G4endl;
   G4cout<<"Test49: The time:     "<<aTime/ns<<" ns"<<G4endl;
   G4cout<<"Test49: Process:      "<<mName<<" (proc="<<proc<<")"<<G4endl;
#endif
#ifdef tdebug
   for(G4int ip=0; ip<nT; ip++) tVal[ip]=(fT+dT*ip)/1000000.;// Fill the t-histogram
#endif
   // Create a DynamicParticle
   for(G4int nen=0; nen<cnE; nen++)                          // LOOP over projectile energy
   {
    G4double  energy = (ep-mp)*MeV;                          // projectile kinetic energy
    if(cnE>1) energy = eli[nen]*MeV;
#ifdef debug
    G4cout<<"Test49: M="<<mp<<", T="<<energy<<" MeV, proc="<<proc<<G4endl;
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
      G4cout<<"Test49: Material="<<material->GetName()<<", Element[0]="<<curEl->GetName()
            <<",A[0]="<<(*(curEl->GetIsotopeVector()))[0]->GetN()<<" is selected."<<G4endl;
     }
     G4cout<<"Test49:NewRun:Targ="<<tPDG<<",M="<<G4QPDGCode(tPDG).GetMass()<<",Proj="<<pPDG
           <<",E="<<energy<<" MeV, Run #"<<nCur<<" of "<<nTot<<", Model="<<mName<<G4endl;
     //G4double mt=G4QPDGCode(tPDG).GetMass();             // @@ just for check
     G4QContent tQC=G4QPDGCode(tPDG).GetQuarkContent();
     G4int    ct=tQC.GetCharge();
     G4int    st=tQC.GetStrangeness();
     G4int    bnt=tQC.GetBaryonNumber();
#ifdef debug
     G4cout<<"Test49: pQC"<<pQC<<", pch="<<cp<<", tQC"<<tQC<<", tch="<<ct<<G4endl;
#endif
     G4int    totC=cp+ct;
     G4int    totS=sp+st;
     G4int    totBN=bnp+bnt;
#ifdef debug
     G4cout<<"Test49:tC="<<totC<<"?="<<totCN<<",tB="<<totBN<<"?="<<totBNN<<",tS="<<totS
           <<"?="<<totSN<<", proc="<<proc<<G4endl;
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
     G4cout<<"Test49: The end point is defined and filled in the step "<<G4endl;
#endif
     G4Navigator* nav = new G4Navigator;
#ifdef pverb
     G4cout<<"Test49: The Navigator is defined "<<G4endl;
#endif
     nav->SetWorldVolume(pFrame);
#ifdef pverb
     G4cout<<"Test49: The Box frame is set "<<G4endl;
#endif
     G4TouchableHandle touch(nav->CreateTouchableHistory());
#ifdef pverb
    G4cout<<"Test49: The TouchableHandle is defined "<<G4endl;
#endif
     G4Timer* timer = new G4Timer();
     timer->Start();
#ifdef idebug
     G4cout<<"Test49: Run is started, timer is started, kinEnergy="<<energy<<G4endl;
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
#ifdef idebug
     G4cout<<"Test49: Before the event loop, nEvents= "<<nEvt<<G4endl;
#endif
     G4double dTot=0.;
     G4double dEl=0.;
     G4int    nDoNothing=0;
     for (G4int iter=0; iter<nEvt; iter++)
     {
#ifdef debug
	 G4cout<<"Test49: ## "<<iter<<"-th event ##,energy="<<energy<<",proc="<<proc<<G4endl;
#endif

      if(!(iter%1000)&&iter)G4cout<<"***=>TEST19: "<<iter<<" events are simulated"<<G4endl;

      gTrack->SetStep(step);            // Now step is included in the Track (see above)
      gTrack->SetKineticEnergy(energy); // Duplication of Kin. Energy for the Track (?!)
      gTrack->SetTouchableHandle(touch);// Set Box touchable history
#ifdef debug
      G4cout<<"Test49: Before the fake proc->GetMeanFreePath call, proc="<<proc<<G4endl;
#endif
      if(mName != "chips")hproc->GetMeanFreePath(*gTrack,0.1,cond);//ToAvoid GHAD complains
#ifdef debug
      G4cout<<"Test49:Before PostStepDoIt, p="<<proc<<",hp="<<hproc<<",cp="<<cproc
            <<", model = "<<mName<<G4endl;
#endif
      if(mName == "chips")
      {
#ifdef debug
        G4cout<<"Test49: Before PostStepDoIt, CHIPS proc="<<hproc<<G4endl;
#endif
         aChange = static_cast<G4ParticleChange*>(cproc->PostStepDoIt(*gTrack,*step));
      }
      else
      {
#ifdef debug
        G4cout<<"Test49: Before PostStepDoIt, GHAD proc="<<hproc<<", pPDG="<<pPDG<<G4endl;
#endif
        if     (pPDG == 2212)
        {
          aChange = static_cast<G4ParticleChange*>(prhproc->PostStepDoIt(*gTrack,*step));
#ifdef debug
          G4cout<<"Test49: Simulation is done for proton PDG="<<pPDG<<G4endl;
#endif
        }
        else if(pPDG== 2112)
        {
          aChange = static_cast<G4ParticleChange*>(nehproc->PostStepDoIt(*gTrack,*step));
#ifdef debug
          G4cout<<"Test49: Simulation is done for neutron PDG="<<pPDG<<G4endl;
#endif
        }
        else if(pPDG== -211)
        {
          aChange = static_cast<G4ParticleChange*>(pmhproc->PostStepDoIt(*gTrack,*step));
#ifdef debug
          G4cout<<"Test49: Simulation is done for piminus PDG="<<pPDG<<G4endl;
#endif
        }
        else if(pPDG==  211)
        {
          aChange = static_cast<G4ParticleChange*>(pphproc->PostStepDoIt(*gTrack,*step));
#ifdef debug
          G4cout<<"Test49: Simulation is done for piplus PDG="<<pPDG<<G4endl;
#endif
        }
        else if(pPDG== -321)
        {
          aChange = static_cast<G4ParticleChange*>(kmhproc->PostStepDoIt(*gTrack,*step));
#ifdef debug
          G4cout<<"Test49: Simulation is done for kaminus PDG="<<pPDG<<G4endl;
#endif
        }
        else if(pPDG==  321)
        {
          aChange = static_cast<G4ParticleChange*>(kphproc->PostStepDoIt(*gTrack,*step));
#ifdef debug
          G4cout<<"Test49: Simulation is done for kaplus PDG="<<pPDG<<G4endl;
#endif
        }
        else if(pPDG==-2112)
        {
          aChange = static_cast<G4ParticleChange*>(aphproc->PostStepDoIt(*gTrack,*step));
#ifdef debug
          G4cout<<"Test49: Simulation is done for anti-proton PDG="<<pPDG<<G4endl;
#endif
        }
        else G4cout<<"-Error-Test49: Process is not defined for PDG="<<pPDG<<G4endl;
      }
#ifdef debug
      G4cout<<"Test49: Before PostStepDoIt, proc="<<proc<<G4endl;
#endif
      G4int nSec = aChange->GetNumberOfSecondaries();
      G4TrackStatus lead = aChange->GetTrackStatus();
#ifdef debug
      G4cout<<"T19:AfterPostStepDoIt,nS="<<nSec<<",l="<<lead<<" =? fA="<<fAlive<<G4endl;
#endif
      if(!nSec && lead==fAlive)
      {
        G4double eOut=aChange->GetEnergy();
        if(std::fabs(eOut-energy) < .01) nDoNothing++; // Calculate # of DoNothing events
        else G4cerr<<"**Tst19: DoNothing E="<<eOut<<" # "<<energy<<G4endl;
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
        G4int    curN = tgN;
        if(mName == "chips") curN = cproc->GetNumberOfNeutronsInTarget();
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
          G4cout<<"Test49: E="<<sen+sma<<",p="<<(*smo)*std::sqrt(ten*ten-sma*sma)<<G4endl;
          G4cout<<"Test49: L4m="<<s4m<<",T4m="<<totSum<<", N="<<curN<<",M="<<curM<<G4endl;
#endif
          totKE+=s4m.e()-sma;
          totSum-=s4m;
          totSumM=totSum;
        }
#ifdef debug
        G4cout<<"Test49: tChg="<<totC<<", T4m="<<totSum<<",tBar="<<totBaryN<<G4endl;
#endif
        G4LorentzVector Residual(0.,0.,0.,0.);
        if(mName == "chips") Residual = cproc->GetEnegryMomentumConservation();
#ifdef debug
        G4double de = aChange->GetLocalEnergyDeposit();// Init TotalEnergy by EnergyDeposit
        G4cout<<"Test49: "<<nSec<<" secondary particles are generated, dE="<<de<<G4endl;
#endif
        // @@ ----------------------- Begin
        G4double weight=aChange->GetSecondary(0)->GetDynamicParticle()->GetKineticEnergy();
#ifdef debug
        G4cout<<"Test49:----------------: Weigh="<<weight<<G4endl;
#endif
        dTot+=weight;
        if(nSec>1) dEl+=weight;
        // @@ ----------------------- End
        //G4int nbar = 0;

        // @@ G4int npt=0;
        G4int    c=0;    // Prototype of the PDG Code of the particle
        G4int nGamma=0;
        G4double EGamma=0;
        G4int nP0=0;
        G4int nPP=0;
        G4int nPN=0;
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
        G4cout<<"Test49:----DONE^^^^^^^*******^^^^^^^^:ir="<<iter<<": #ofH="<<nSec<<G4endl;
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
            G4cout<<"Test49:LeadingParticle is found, E="<<e<<",m="<<m<<",PDG="<<c<<G4endl;
#endif
            cST= sp;
            cCG= cp;
            cBN= bnp;

          }
          else    // Secondaries
          {
#ifdef debug
            G4cout<<"Test49:Secondary i="<<i<<",*Ses="<<aChange->GetSecondary(i)<<G4endl;
#endif
            sec = aChange->GetSecondary(i)->GetDynamicParticle();
#ifdef debug
            G4cout<<"Test49:Secondary i="<<i<<",*SPD="<<sec<<G4endl;
#endif
            pd  = sec->GetDefinition();
#ifdef debug
            G4cout<<"Test49:Secondary i="<<i<<",*PD="<<pd<<G4endl;
#endif
            c   = pd->GetPDGEncoding();
            cST = static_cast<G4int>(pd->GetQuarkContent(3)-pd->GetAntiQuarkContent(3));
            cCG = static_cast<G4int>(pd->GetPDGCharge());
            cBN = static_cast<G4int>(pd->GetBaryonNumber());
            if(!c) c=90000000+cST*999999+cCG*999+cBN;
            m   = pd->GetPDGMass();
            mom = sec->GetMomentumDirection();
            e   = sec->GetKineticEnergy();
#ifdef debug
            G4cout<<"Test49:SecPt#"<<i<<",T="<<e<<",m="<<m<<",PDG="<<c<<",C="<<cCG<<G4endl;
#endif
          }
          if (e < 0.0)
          {
            G4cerr<<"**Test49:Event#"<<iter<<",Hadron#"<<i<<",T="<<e<<"<0 (Set 0)"<<G4endl;
            e = 0.0;
          }
          // for exclusive reaction 2 particles in final state
          p = std::sqrt(e*(e + m + m));
          mom *= p;
#ifdef debug
          G4cout<<"Test49:Secondary #"<<i<<" E="<<e+m<<", P="<<mom<<G4endl;
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
          if(i>=0)                             // Do not count LV & charges of leptons
          {
            totStran-=cST;
            totCharge-=cCG;
            totBaryN-=cBN;
            totSum -= lorV;
            G4double visE=lorV.e();
            if(cBN>0) visE-=m;
            if(cBN<0) visE+=m;
            totKE +=  visE;
#ifdef debug
            G4cout<<"Test49:i="<<i<<",tS="<<totStran<<",tC="<<totC<<",t4="<<totSum<<",tB="
                  <<totBaryN<<G4endl;
#endif
          }
#ifdef debug
          G4cout<<"Test49: *** M="<<lorV.m()<<G4endl;
#endif
          if(fabs(m-lorV.m())>.005&&1>2)
          //if(std::fabs(m-lorV.m())>.005) // @@ Temporary check
          {
            G4cerr<<"***Test49: m="<<lorV.m()<<" # "<<m<<", d="<<lorV.m()-m<<G4endl;
            alarm=true;
          }
          if(!(lorV.e()>=0||lorV.e()<0)   || !(lorV.px()>=0||lorV.px()<0) ||
             !(lorV.py()>=0||lorV.py()<0) || !(lorV.pz()>=0||lorV.pz()<0))
          {
            G4cerr<<"***Test49: NAN in LorentzVector="<<lorV<<G4endl;
            alarm=true;
          }
          if(c==90000002||c==90002000||c==92000000)
          {
            G4cout<<"***Test49:***Dibaryon *** i="<<i<<", PDG="<<c<<G4endl;
            alarm=true;
          }
          if(c==90000003||c==90003000||c==93000000)
          {
            G4cout<<"***Test49:***Tribaryon *** i="<<i<<", PDG="<<c<<G4endl;
            alarm=true;
          }
#ifdef debug
          G4cout<<"Test49: ===> alarm="<<alarm<<G4endl;
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
          if(c==111) nP0++;                              // Neutral  pions
          if(c==-211) nPN++;                             // Negative pions
          if(c==211) nPP++;                              // Positive pions
          if(c==22)                                      // Gammas
          {
            nGamma++;                                    // # of gammas
            EGamma+=e;                                   // Total energy in gammas
          }
#ifdef ppdebug
          if(lead==fAlive && (i==-1 || c==2112 || c==2212 || c==90000001 || c==90001000))
          //if(lead==fAlive)
            G4cout<<"Test49:#"<<i<<",PDG="<<c<<",S="<<cST<<",C="<<cCG<<",B="<<cBN<<",4M="
                  <<lorV<<m<<",T="<<lorV.e()-m<<G4endl;
#endif
          //delete aChange->GetSecondary(i);
#ifdef pdebug
          //if(lead==fAlive && (i==-1 || c==2112 || c==2212 || c==90000001 || c==90001000))
          G4cout<<"Test49:#"<<i<<",PDG="<<c<<",C="<<cCG<<",B="<<cBN<<",4M="<<lorV<<m<<",T="
                <<lorV.e()-m<<G4endl;
#endif
        } // End of the LOOP over secondaries
        // delete secondaries in the end of the event         
	  G4double misr=-DBL_MAX;
#ifdef inter
	  misr=.27;
#endif
#ifdef lhepdbg
        if(mName=="lep" || mName=="hep")
        {
          // --- Start detailed correction
          G4double resGSM =
                   G4QNucleus(90000000+totStran*999999+totCharge*999+totBaryN).GetGSMass();
          G4double resMom = totSum.rho();
          G4double resE= std::sqrt(resGSM*resGSM+resMom*resMom);
          G4double redE= totSum.e()-resE;
          // *** Debug print
          //G4double resT= resE-resGSM;
          //G4cout<<"Test49: LHEP_COR_E="<<redE<<", T="<<resT<<", P="<<resMom<<G4endl;
          // *** End of debug print
          totSum=G4LorentzVector(0.,0.,0.,redE);
          // --- Stop detailed correction
          totStran=0;
          totCharge=0;
          totBaryN=0;
          //if(std::fabs(totSum.e())>.1) 
          //  G4cout<<"Test49:--->LHEP_COR_E="<<redE<<", T="<<resT<<", P="<<resMom<<G4endl;
        }
#endif
#ifdef meandbg
        //G4double eng=energy+mp;
        G4double eng=energy;
        G4double cE=eng-totSum.t();
        if(cE<0.) G4cout<<"Test49: **Negative**, E="<<cE<<G4endl;
        if(cE>eng+eng) nD++;
        nE++;
        sE+=cE;
        sE2+=cE*cE;
        pE+=eng;
#endif
#ifdef ekindbg
        if(totKE<0.) G4cout<<"Test49: **Negative**, E="<<totKE<<G4endl;
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
          //G4cout<<"Test49: Too big energy nonconservation dE="<<dE<<G4endl;
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
        G4cout<<">TEST19: r4M="<<totSum<<",rCh="<<totCharge<<",rBN="<<totBaryN<<G4endl;
#endif
        if(mName == "chips")
        {
          G4double ss=std::fabs(totSum.t())+std::fabs(totSum.x())+
                      std::fabs(totSum.y())+std::fabs(totSum.z());
          G4double sr=std::fabs(Residual.t())+std::fabs(Residual.x())+
                      std::fabs(Residual.y())+std::fabs(Residual.z());    

          if (totCharge || totBaryN || ss>.27 || alarm || (nGamma && !EGamma))
          {
            G4cerr<<"**Test49:#"<<iter<<":n="<<nSec<<",4M="<<totSum<<",Charge="<<totCharge
                  <<",BaryN="<<totBaryN<<",R="<<Residual<<",D2="<<ss<<",nN="<<curN<<G4endl;
            totSum = totSumM;
            if(nGamma&&!EGamma)G4cerr<<"***Test49: Egamma=0"<<G4endl;
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
                cST =static_cast<G4int>(pd->GetQuarkContent(3)-pd->GetAntiQuarkContent(3));
                cCG = static_cast<G4int>(pd->GetPDGCharge());
                cBN = static_cast<G4int>(pd->GetBaryonNumber());
                if(!c) c=90000000+cST*999999+cCG*999+cBN;
                m   = pd->GetPDGMass();
                //G4cerr<<"***Tmp***Test49:#"<<indx<<",S="<<cST<<",C="<<cCG<<",B="<<cBN
                //      <<",M="<<m<<G4endl;
                mom = sec->GetMomentumDirection();
                e   = sec->GetKineticEnergy();
              }
              p = std::sqrt(e*(e + m + m));
              mom *= p;
              lorV = G4LorentzVector(mom, e + m);    // "e" is a Kinetic energy!
              totSum -= lorV;
              G4cerr<<"Test49:#"<<indx<<",PDG="<<c<<",m="<<m<<",4M="<<lorV<<",T="<<e
                    <<", d4M="<<totSum<<", S="<<cST<<", C="<<cCG<<", B="<<cBN<<G4endl;
            }
            if(sr>misr)G4Exception("**Test49:ALARM/baryn/chrg/energy/mom isn't conserved");
          }
        }
        ntp->FillEvt(aChange,dParticle); // Fill the simulated event in the ASCII "ntuple"
        // =============== May be print it here if it is not zero...
        for(G4int ides=0; ides<nSec; ides++) delete aChange->GetSecondary(ides);
        aChange->Clear();
#ifdef debug
        G4cout<<"Test49:--->>> Ev # "<<iter<<", After ntp.FillEvt"<<G4endl;
#endif
       } // End of nSec>0 IF
      } // End of Check for DoNothing
     } // End of the LOOP over events
     // Stop the timer to estimate the speed of the generation
     G4double dn=100.* nDoNothing / nEvt;
     timer->Stop();
     G4cout<<"Test49:Time/Ev="<<timer->GetUserElapsed()/nEvt<<" s, DN="<<dn<<" %"<<G4endl;
     delete timer;
#ifdef tdebug
     G4double pGeV=pmax/1000.;
     G4double alp=std::log(pGeV/(1.+1./pGeV/pGeV/pGeV));
     G4double exr=1.-.822/(1.+std::exp(-(alp-.2)*1.15));
     G4cout<<"Test49:EndOfRun,p="<<pmax<<",dT="<<dTot<<",r="<<dEl/dTot<<",e="<<exr
           <<",d="<<exr-dEl/dTot<<",ra="<<(exr-.178)/.822<<G4endl;
     for(G4int ir=0; ir<nT; ir++) G4cout<<tVal[ir]<<" "<<tSig[ir]<<G4endl;// Print t-vectors
#endif
#ifdef pverb
     G4cerr<<"Test49: ########## End of run ##########"<<G4endl;
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
  //G4cout << "Test49: After delete process etc." << G4endl;
  //delete phys;                        // Absolete physics class
  //delete cmaterial;                   // Temporary material definition pointer (?)
  //G4cout << "Test49: Before ntp" << G4endl;
  delete ntp; // Delete the class to fill the#of events
  //G4cout << "Test49: After ntp" << G4endl;
#ifdef hdebug
  timer->Stop();
  G4int    nbnh=G4StringChipsParticleLevelInterface::GetNbn();
  G4double hbdb=G4StringChipsParticleLevelInterface::GetDB();
  //G4double hede=G4StringChipsParticleLevelInterface::GetDE(); // The same now
  G4double toth=G4StringChipsParticleLevelInterface::GetTot();
  G4double ovrb=G4StringChipsParticleLevelInterface::GetBov();
  G4double ovre=G4StringChipsParticleLevelInterface::GetEov();
  G4cout<<"Test49:TimePerEvent="<<timer->GetUserElapsed()/toth<<", N="<<toth
        <<", overB="<<ovrb/toth<<", overE="<<ovre/toth<<G4endl;
  delete timer;
  G4double bc=-hbdb/2;
  for(G4int ih=0; ih<nbnh; ih++)
  {
    bc+=hbdb;
    G4double nE=G4StringChipsParticleLevelInterface::GetE(ih);
    G4double nB=G4StringChipsParticleLevelInterface::GetB(ih);
    G4cout<<"Test49:be="<<bc<<" "<<nE/toth<<" "<<std::sqrt(nE)/toth
          <<" "<<nB/toth/bc<<" "<<std::sqrt(nB)/toth/bc<<G4endl;
  }
#endif
  delete runManager;
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
#ifdef pverb
  G4cout << "###### End of Test49 #####" << G4endl;
#endif
}
