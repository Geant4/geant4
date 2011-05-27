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
//      File name:     Test39 (Comparison of G4QElastic Results with Data)
//
//      Author:        M.Kossov (made of test19)
// 
//      Creation date: 23 Jan 2006
//
//      Modifications:
//
// -------------------------------------------------------------------
//       1         2         3         4         5         6         7         8         9
//34567890123456789012345678901234567890123456789012345678901234567890123456789012345678901

#define nout
//#define pscan
//#define smear
//#define escan
//#define idebug
//#define tdebug
//#define csdebug
#define tdisthist
//#define pverb
//#define debug
//#define pdebug
//#define fdebug
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
#include "G4VParticleChange.hh"        // @@ should we keep both? (this is not necessary)
#include "G4ParticleChange.hh"
#include "G4QCaptureAtRest.hh"
#include "G4QElastic.hh"               // CHIPS
#include "G4LElastic.hh"
#include "G4WHadronElasticProcess.hh"
#include "G4HadronElastic.hh"
//#include "G4LElasticB.hh"
#include "G4LEnp.hh"
#include "G4LEpp.hh"
#include "G4NeutronHPChannel.hh"       // HP cross-section
#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPElasticData.hh"
#include "G4NeutronHPElasticFS.hh"
#include "G4NeutronHPThermalBoost.hh"
#include "G4ElasticHadrNucleusHE.hh"
#include "G4HadronElasticProcess.hh"
#include "G4HadronCrossSections.hh"
#include "G4QProtonElasticCrossSection.hh"  // CHIPS pA elastic Cross-Sections
#include "G4QNeutronElasticCrossSection.hh" // CHIPS nA elasticCross-Sections

#include "G4ApplicationState.hh"
#include "G4StateManager.hh"

#include "G4UnitsTable.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleChange.hh"
#include "G4DynamicParticle.hh"
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

#include "Test39Physics.hh"
#include "Test39PhysicsList.hh"

//#include "G4HadronCrossSections.hh"
//#include "G4VCrossSectionDataSet.hh"
//#include "G4ProtonInelasticCrossSection.hh"
//#include "G4NeutronInelasticCrossSection.hh"
//#include "G4HadronInelasticDataSet.hh"

//#include "G4ExcitationHandler.hh"
//#include "G4PreCompoundModel.hh"
//#include "G4Evaporation.hh"

// New Histogramming (from AIDA and Anaphe):
//#include <memory> // for the auto_ptr(T>

// Temporary AIDA can be not used in test39
//#include "AIDA/AIDA.h"

// ASCII Pseudo NTUPLE 
#ifndef nout
#include "G4QHBook.hh"
#endif

#include "G4Timer.hh"
#include "time.h"

//#define debug

//int main(int argc, char** argv)
int main()
{
  const G4int nTg=3;   // Length of the target list for the Performance test
  G4int tli[nTg]={90001000,90013014,90082126};     // Targets CHIPS codes
  G4String tnm[nTg]={"Hydrogen","Aluminum","Lead"};// Target names
  G4String tsy[nTg]={"H","Al","Pb"};               // Target symbols
  G4Material* mat[nTg]={0,0,0,};                   // Material pointers for the Target Loop
  const G4int nPr=10;  // Length of the projectile list for the Performance test
  G4int pli[nPr] = {211, 2112, 2212, -211, 321,-321, 310, 3122, 3222, -2212}; // projPDGs
  const G4int nMom=3;  // Length of the projectile momentum list for the Performance test
  G4double mom[nMom] = {100., 1000., 10000.}; // Set of momenta in MeV/c !
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
  const G4int nAZ=270;         // Dimension of the table
  const G4int mAZ=266;         // Maximum filled A (at present). Must be mAZ<nAZ
  // Best Z for the given A - changed by MK (@@Not one-to-one correspondance! Make alt "-")
  const G4int bestZ[nAZ] = {
     0,  1,  1,  2,  2,  0,  3,  3,  4,  4,   //0
     5,  5,  6,  6,  7,  7,  8,  8,  8,  9,   //10
    10, 10, 10, 11, 12, 12, 12, 13, 14, 14,   //20
    14, 15, 16, 16, 16, 17, 18, 17, 18, 19,   //30
    18, 19, 20, 20, 20, 21, 22, 22, 23, 23,   //40
    22, 23, 24, 24, 26, 25, 26, 26, 28, 27,   //50
    28, 28, 28, 29, 30, 29, 30, 30, 30, 31,   //60
    32, 31, 32, 32, 32, 33, 34, 34, 34, 35,   //70
    34, 35, 36, 36, 36, 37, 39, 36, 38, 39,   //80
    40, 40, 41, 40, 40, 42, 42, 42, 42, 44,   //90
    44, 44, 44, 45, 44, 46, 46, 47, 46, 47,   //100
    48, 48, 48, 48, 48, 49, 50, 50, 50, 50,   //110
    50, 51, 50, 51, 50, 52, 52, 53, 52, 54,   //120
    52, 54, 54, 55, 54, 56, 54, 56, 56, 57,   //130
    58, 59, 60, 60, 60, 60, 60, 62, 62, 62,   //140
    62, 63, 62, 63, 62, 64, 64, 64, 64, 65,   //150
    64, 66, 66, 66, 66, 67, 68, 68, 68, 69,   //160
    68, 70, 70, 70, 70, 71, 70, 72, 72, 72,   //170
    72, 73, 74, 74, 74, 75, 74, 75, 76, 76,   //180
    76, 77, 76, 77, 78, 78, 78, 79, 80, 80,   //190
    80, 80, 80, 81, 80, 81, 82, 82, 82, 83,   //200
    82,  0, 82,  0, 82,  0, 84,  0,  0,  0,   //210
    86,  0, 86, 87, 88,  0, 88, 89, 88, 89,   //220
    89, 91, 90,  0, 92, 92,  0, 93, 92, 94,   //230
     0,  0,  0, 95, 94,  0,  0, 96,  0,  0,   //240
     0, 98, 99,  0,  0,  0,  0,100,101,102,   //250
   103,104,105,106,  0,108,109,  0,  0,  0};  //260
  // 0   1   2   3   4   5   6   7   8   9
  // Second candidate
  const G4int secoZ[nAZ] = {
     0,  0,  0,  1,  0,  0,  0,  4,  0,  0,   //0
     4,  6,  5,  7,  0,  8,  0,  0,  0,  8,   //10
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,   //20
     0,  0, 15,  0,  0, 16,  0,  0,  0,  0,   //30
    20, 20,  0,  0,  0,  0, 20,  0,  0,  0,   //40
    24,  0,  0,  0, 24,  0,  0,  0, 26,  0,   //50
    27,  0,  0,  0, 28,  0,  0,  0,  0,  0,   //60
    30,  0,  0,  0, 34,  0, 32,  0,  0,  0,   //70
    36,  0, 34,  0, 38,  0, 38, 38,  0,  0,   //80
     0,  0, 42,  0, 42,  0, 44,  0, 44,  0,   //90
    42,  0, 46,  0, 46,  0, 48,  0, 48,  0,   //100
    46,  0, 50, 49, 50, 50, 48,  0,  0,  0,   //110
    52,  0, 52,  0, 52,  0, 54,  0, 54,  0,   //120
    54,  0, 56,  0, 56,  0, 56,  0, 58,  0,   //130
    54,  0, 58,  0, 62, 61,  0,  0, 60,  0,   //140
    60,  0, 64,  0, 64,  0, 66,  0, 66,  0,   //150
    66,  0, 68,  0, 68,  0,  0,  0, 70,  0,   //160
    70,  0,  0,  0, 72,  0, 72,  0,  0,  0,   //170
    74,  0,  0,  0, 76,  0, 76, 76,  0,  0,   //180
    78,  0, 78,  0,  0,  0, 80,  0, 78,  0,   //190
     0,  0,  0,  0, 82,  0,  0,  0,  0, 84,   //200
    84,  0, 83,  0, 83,  0,  0,  0,  0,  0,   //210
     0,  0,  0,  0,  0,  0,  0,  0, 89,  0,   //220
     0,  0,  0,  0, 93,  0,  0,  0, 93,  0,   //230
     0,  0,  0,  0,  0,  0,  0, 97,  0,  0,   //240
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,   //250
     0,  0,107,  0,  0,  0,  0,  0,  0,  0};  //260
  // 0   1   2   3   4   5   6   7   8   9
  // Second candidate
  const G4int thrdZ[nAZ] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //0
    0, 0, 7, 0, 0, 0, 0, 0, 0, 0,   //10
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //20
    0, 0, 0, 0, 0,20, 0, 0, 0, 0,   //30
   19, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //40
   23, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //50
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //60
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //70
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //80
    0, 0,36, 0,38, 0,40, 0, 0, 0,   //90
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //100
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //110
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //120
   56, 0,50, 0, 0, 0,58, 0,57, 0,   //130
    0, 0, 0, 0, 0, 0, 0, 0,65, 0,   //140
    0, 0,66, 0, 0, 0, 0, 0, 0, 0,   //150
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //160
    0, 0, 0, 0, 0, 0,71, 0, 0, 0,   //170
   73, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //180
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //190
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //200
   83, 0,84, 0,84, 0, 0, 0, 0, 0,   //210
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //220
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //230
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //240
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //250
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0};  //260
  //0  1  2  3  4  5  6  7  8  9
  // Fourth candidate (only two isotopes)
  const G4int quadZ[nAZ] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //0
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //10
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //20
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //30
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //40
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //50
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //60
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //70
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //80
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //90
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //100
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //110
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //120
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //130
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //140
    0, 0,67, 0, 0, 0, 0, 0, 0, 0,   //150
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //160
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //170
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //180
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //190
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //200
   85, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //210
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //220
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //230
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //240
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   //250
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0};  //260
  //0  1  2  3  4  5  6  7  8  9

  // Test of not overlaping four tables above (can be commented)
  if(mAZ>=nAZ) G4cout<<"***Test39: Too big mAZ="<<mAZ<<", nAZ="<<nAZ<<G4endl;

  // Run manager
  G4RunManager* runManager = new G4RunManager;
  runManager->SetUserInitialization(new Test39PhysicsList);
  G4StateManager::GetStateManager()->SetNewState(G4State_Init); // To let create ions
  G4ParticleDefinition* ionDefinition=0;
  ionDefinition=G4ParticleTable::GetParticleTable()->FindIon(6,12,0,6);
  if(!ionDefinition)
  {
    G4cerr<<"*** Error! *** Test39:(6,6) ion can not be defined"<<G4endl;
    return 0;
  }
  else G4cout<<"Test39: (6,6) ion is OK, Run State="<<G4StateManager::GetStateManager()->
              GetStateString(G4StateManager::GetStateManager()->GetCurrentState())<<G4endl;

  for(G4int a=1; a<nAZ; a++)
  {
    G4int z1=bestZ[a]; // First table
    G4int z2=secoZ[a]; // Second table
    G4int z3=thrdZ[a]; // Third table
    G4int z4=quadZ[a]; // Fourth table
    if(z2)
    {
      if(z1==z2)G4cout<<"#"<<a<<": z1="<<z1<<" = z2="<<z2<<", z3="<<z3<<",z4="<<z4<<G4endl;
      if(z3)
      {
        if(z1==z3)G4cout<<"#"<<a<<",z1="<<z1<<" = z3="<<z3<<",z2="<<z2<<",z4="<<z4<<G4endl;
        if(z2==z3)G4cout<<"#"<<a<<",z1="<<z1<<",z2="<<z2<<" = z3="<<z3<<",z4="<<z4<<G4endl;
        if(z4)
        {
          if(z1==z4)G4cout<<"#"<<a<<",z1="<<z1<<"=z4="<<z4<<",z2="<<z2<<",z3="<<z3<<G4endl;
          if(z2==z4)G4cout<<"#"<<a<<",z1="<<z1<<",z2="<<z2<<"=z4="<<z4<<",z3="<<z3<<G4endl;
          if(z3==z4)G4cout<<"#"<<a<<",z1="<<z1<<",z2="<<z2<<",z3="<<z3<<"=z4="<<z4<<G4endl;
        }
      }
      else if(z4) G4cout<<"#"<<a<<",z1="<<z1<<",z2="<<z2<<",z3=0.(!) & z4="<<z4<<G4endl;
    }
    else if(z3||z4)G4cout<<"#"<<a<<",z1="<<z1<<",z2=0.(!) & z3="<<z3<<" & z4="<<z4<<G4endl;
  }
#ifdef debug
  G4cout<<"Test39: Prepare G4QHBook files or ntuples"<<G4endl;
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
  G4cout<<"Test39:Par's: T="<<temperature<<",S="<<ssin2g<<",Eta="<<eteps<<",nop="
        <<nop<<",p="<<momb<<",e="<<enb<<",pPDG="<<pPDG<<",tPDG="<<tPDG<<",nEv="
        <<nEvt<<",fN="<<fN<<",fD="<<fD<<",cP="<<cP<<",rM="<<rM<<",sA="<<sA<<G4endl;
  //-------- Initialize CHIPS interaction (not necessary)
  //G4QCHIPSWorld* theW=G4QCHIPSWorld::Get();
  //theW->GetParticles(nop);           // Create CHIPS World of nop particles
  // ********** Now momb is a momentum of the incident particle, if =0 => LOOP ************
  G4double mp=G4QPDGCode(pPDG).GetMass();
  G4double ep=mp;
  G4int cnE=1;
  if(momb==0.) cnE=nMom;
  else
  {
    ep=std::sqrt(mp*mp+momb*momb);
    if(enb>0.) ep=enb;
  }
  //G4int tPDG=90000000+tgZ*1000+tgN;
  G4int    tgZ=(tPDG-90000000)/1000;
  G4int    tgN=tPDG-90000000-tgZ*1000;
  // ---------- Define material for the simulation ------------------
  G4int tgA        = tgZ+tgN; // Mass number - fake
  // The material can be copied from the commented cMaterial Factory above
  G4Isotope* isotope=0;
  G4Element* element=0;
  G4Material* material=0;
  // LOOPs for the wide performance test
  G4double      aTime      = 0. ;
  G4ThreeVector aDirection = G4ThreeVector(0.,0.,1.);
  G4int tgm=1;                                        // By default only one target
  if(!tPDG) // Make max for the LOOP over all targets and define materials
  {
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
      G4cout<<"Test39:-->Material("<<tgZ<<","<<tgN<<"):"<<tnm[tgi]<<","<<tsy[tgi]
            <<"; index="<<material->GetIndex()<<G4endl;
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
  G4cout<<"Test39:--- Material "<<material->GetName()<<" is defined ---" << G4endl;
#endif
  if(!material)
  {
    G4cout<<"Test39: Last Material "<<material->GetName()<<" is not defined."<<G4endl;
    exit(1);
  }
  // ---------------------- PSEUDO-TRACK DEFINISIONS -------------------------
  G4double nx = 0.;
  G4double ny = 0.;
  G4double nz = 0.;
  G4int npart=1;                                      // By default only one particle
  if(!pPDG) npart=nPr;                                // Make a LOOP ove all particles
  G4QElastic* proc = new G4QElastic;
#ifdef debug
  G4cout<<"G4GH_ElasticPhysics::ConstructProcess: the elastic model is defined"<<G4endl;
#endif
  if(!proc)
  {
    G4cout<<"Tst39: there is no G4...Elastic process"<<G4endl;
    exit(1);
  }
#ifdef debug
  G4cout<<"Test39:--***-- process is created --***--" << G4endl; // only one run
#endif
  G4int nTot=npart*tgm*cnE;
  G4int nCur=0;
  for(G4int pnb=0; pnb<npart; pnb++) // LOOP over particles
  {
   if (npart>1) pPDG=pli[pnb];
   G4QContent pQC=G4QPDGCode(pPDG).GetQuarkContent();
   G4int        cp = pQC.GetCharge();          // @@ must be the same as chp
   if(pPDG==22) cp = 0;                        // QuarkContent of photon is not supported
   if(pPDG==11 || pPDG==13 || pPDG==15) cp=-1; // QuarkContent of leptons is not supported
   if(pPDG==-11|| pPDG==-13|| pPDG==-15) cp=1; // QuarkContent of antilept isn't supported
   G4ParticleDefinition* part=0;               // DefaultProj is particle is empty
   if(pPDG==2212) part=G4Proton::Proton();           // Definition of Proton projectile
   else if(pPDG==2112) part=G4Neutron::Neutron();    // Definition of Neutron projectile
   else if(pPDG==22) part=G4Gamma::Gamma();          // Definition of gamma projectile
   else if(pPDG==11) part=G4Electron::Electron();    // Definition of electron projectile
   else if(pPDG==-11) part=G4Positron::Positron();   // Definition of positron projectile
   else if(pPDG==13) part=G4MuonMinus::MuonMinus();  // Definition of mu- projectile
   else if(pPDG==-13) part=G4MuonPlus::MuonPlus();   // Definition of the mu+ projectile
   else if(pPDG==15) part=G4TauMinus::TauMinus();    // Definition of tau- projectile
   else if(pPDG==-15) part=G4TauPlus::TauPlus();     // Definition of tau+ projectile
   else if(pPDG==-211) part=G4PionMinus::PionMinus();// Definition of Pi- projectile
   else if(pPDG==211)  part=G4PionPlus::PionPlus();// Definition of Pi+ projectile
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
     G4cerr<<"***Test39: "<<pPDG<<" is a PDG code of not supported particle"<<G4endl;
     G4Exception("***Test39: OnFlight Process is called for not supported particle");
   }
   // Not for CHIPS & not for GHAD/G4LElasic: only for G4UHadronElastic of V. Ivanchenko
   ///proc->BuildPhysicsTable(*part);//NotNecessary forG4LElastic&CHIPS. Only forNewEl V.I.
   // ----------------------------
   G4double pMass = part->GetPDGMass();                 // Mass of the projectile in IU
   //
   G4ThreeVector aPosition(nx*mm, ny*mm, nz*mm);
   G4DynamicParticle* dParticle = new G4DynamicParticle(part,aDirection,0.);// Dummy Energy
   // General Track definition
   G4Track* gTrack = new G4Track(dParticle,aTime,aPosition); // Track definition (NoTarg)
   // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   //G4DynamicParticle* pPart = new G4DynamicParticle(G4Proton::Proton(),aDirection,0.);
   //G4DynamicParticle* nPart = new G4DynamicParticle(G4Neutron::Neutron(),aDirection,0.);
   //G4Track* pTrack = new G4Track(pPart,aTime,aPosition); // Track definition (NoTarg)
   //G4Track* nTrack = new G4Track(nPart,aTime,aPosition); // Track definition (NoTarg)
   // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   G4int    bnp=pQC.GetBaryonNumber();                  // Baryon number of the projectile
   G4LorentzVector proj4Mom(0.,0.,momb,ep);
   // ---------- Define material for the simulation ------------------
   //G4double tgA     = 26.98; // @@ Important? Can it be just tgZ+tgN?
   G4int tgA        = tgZ+tgN;
   G4String nameMat = "Thin Target";  // @@ Not important can be an arbitrary name
   // The material can be copied from the commented cMaterial Factory above
   G4Isotope* isotope = new G4Isotope("Isotop", tgZ, tgA);
   G4Element* element = new G4Element("ZA_Isotop", "ZA", 1);
   element->AddIsotope(isotope, 100.*perCent);
   G4Material* material = new G4Material("ZA_Isomer", 2.*g/cm3, 1);
   material->AddElement(element, 1);
#ifdef debug
   G4cout<<"Test39:Material=Element: Z="<<element->GetZ()<<"("<<tgZ<<"), A="
         <<element->GetA()*mole/g<<"("<<tgA<<")"<<element->GetN()<<G4endl;
   G4cout << *material << G4endl;
#endif
   // For the Material Factory case
   //G4String  nameMat  = "Si";
   //G4Material* material = G4Material::GetMaterial(nameMat); // Get definition of Material
   if(!material)
   {
     G4cout<<"Test39:Material "<<nameMat<<" is not defined in the Test39Material"<<G4endl;
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
   G4cout<<"Test39:###### Start new run #####" << G4endl; // only one run
#endif
   G4int    chp=static_cast<G4int>(part->GetPDGCharge()); // Charge of the projectile
#ifdef debug
   G4double totCN  = tgZ+chp;
   G4int    totBNN = tgA+bnp;
#endif
   G4double theStep   = 0.01*micrometer;
   // G4int maxz = (G4int)((*(material->GetElementVector()))[0]->GetZ()) + 1;
#ifdef pdebug
   G4int targPDG=90000000+1000*tgZ+tgN;                 // PDG Code of the target
   G4cout<<"Test39: The particle: "<<part->GetParticleName()<<G4endl;
   G4cout<<"Test39: The material: "<<material->GetName()<<G4endl;
   G4cout<<"Test39: The target:   "<<targPDG<<G4endl;
   G4cout<<"Test39: The step:     "<<theStep/mm<<" mm"<<G4endl;
   G4cout<<"Test39: The position: "<<aPosition/mm<<" mm"<<G4endl;
   G4cout<<"Test39: The direction:"<<aDirection<<G4endl;
   G4cout<<"Test39: The time:     "<<aTime/ns<<" ns"<<G4endl;
#endif
#ifdef tdebug
   for(G4int ip=0; ip<nT; ip++) tVal[ip]=(fT+dT*ip)/1000000.;// Fill the t-histogram
#endif
   // Create a DynamicParticle
   for(G4int nen=0; nen<cnE; nen++)                         // LOOP over projectile energy
   {
    if(cnE>1)
    {
      ep = std::sqrt(mom[nen]*mom[nen]+pMass*pMass);
      mp = pMass;
    }
    G4double  energy = (ep-mp)*MeV;                         // kinetic particle energy IU
#ifdef pdebug
    G4cout<<"Test39: M1="<<mp<<", MM="<<pMass<<", T="<<energy<<", p="<<mom[nen]<<G4endl;
#endif
    dParticle->SetKineticEnergy(energy);// Fill the Kinetic Energy of the projectile

    for(G4int tgi=0; tgi<tgm; tgi++) // Loop over materials
    {
     nCur++;
     if (tgm>1)
     {
      tPDG=tli[tgi];
      tgZ = (tPDG-90000000)/1000;
      tgN = tPDG-90000000-tgZ*1000;
      tgA = tgN+tgZ;
      material=mat[tgi];
      G4Element* curEl=(*(material->GetElementVector()))[0];
      G4cout<<"Test39: Material="<<material->GetName()<<", Element[0]="<<curEl->GetName()
            <<",A[0]="<<(*(curEl->GetIsotopeVector()))[0]->GetN()<<" is selected."<<G4endl;
     }
     G4cout<<"Test39:*NewRun* Targ="<<tPDG<<",Proj="<<pPDG<<",P="<<mom[nen]<<" MeV/c, =>> "
           <<nCur<<" of "<<nTot<<G4endl;
     G4double mt=G4QPDGCode(tPDG).GetMass()*MeV;             // Target mass in IU (MeV?)
     G4QContent tQC=G4QPDGCode(tPDG).GetQuarkContent();
     G4int    ct=tQC.GetCharge();
     G4int    bnt=tQC.GetBaryonNumber();
#ifdef debug
     G4cout<<"Test39: pQC"<<pQC<<", pch="<<cp<<", tQC"<<tQC<<", tch="<<ct<<G4endl;
#endif

     G4int    totC=cp+ct;
     G4int    totBN=bnp+bnt;
#ifdef debug
     G4cout<<"Test39:tC="<<totC<<"?="<<totCN<<",tB="<<totBN<<"?="<<totBNN<<G4endl;
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
     G4cout<<"Test39: The end point is defined and filled in the step "<<G4endl;
#endif
     G4Navigator* nav = new G4Navigator;
#ifdef pverb
     G4cout<<"Test39: The Navigator is defined "<<G4endl;
#endif
     nav->SetWorldVolume(pFrame);
#ifdef pverb
     G4cout<<"Test39: The Box frame is set "<<G4endl;
#endif
     //G4VTouchable* vtouch = nav->CreateTouchableHistory();
     G4TouchableHandle touch(nav->CreateTouchableHistory());
#ifdef pverb
     G4cout<<"Test39: The TouchableHandle is defined "<<G4endl;
#endif
     G4Timer* timer = new G4Timer();
     timer->Start();
#ifdef idebug
     G4cout<<"Test39: Run's started, timer's started, kinEnergy="<<energy<<"(MeV)"<<G4endl;
#endif
     const G4DynamicParticle* sec = 0;
     G4ParticleDefinition* pd = 0;
     G4ThreeVector  mom;
     G4LorentzVector totSum, lorV;
     G4double e, p, m, sm;
     // @@ G4double px, py, pz, pt;
     G4ParticleChange* aChange = 0;
     G4double e0 = energy+pMass;                  // Total energy of the projectile in IU
     G4double pMass2=pMass*pMass;
     G4double pmax=std::sqrt(e0*e0-pMass2);       // Momentum of the projectile in IU
     G4LorentzVector incident4M(0.,0.,pmax,e0);
     G4double et=e0+mt;                           // Total energy in the reaction in IU
     //G4int nEvt=100;
     // Randomization loop: cycle random generator, using 2 lower digits in nEvt
     G4int    iRandCount = nEvt%100;
     G4double vRandCount = 0.;
     while (iRandCount>0)                         // Shift of the RNDN values 
     {
      vRandCount = G4UniformRand();               // Fake calls
      iRandCount--;
     }
#ifdef tdebug
     for(G4int it=0; it<nT; it++) tSig[it]=0;     // Reset the output value
#endif
#ifdef idebug
     G4cout<<"Test39: Before EventLoop, nEvs="<<nEvt<<",M="<<pMass<<",T="<<energy<<G4endl;
#endif
     G4double dTot=0.;
     G4double dEl=0.;
     G4int goodE=0;                              // a#of good events
     G4int misTG=0;                               // a#of events with missing recoil
     G4int badTG=0;                               // a#of events with wrong recoil
     G4int zeroO=0;                               // a#of events with noLead & 0 in OUT
     G4int alonO=0;                               // a#of events with noLead & 1 in OUT
     G4int badOT=0;                               // a#of events with wrong secondaries
     const G4int ntpt=3;
     G4double thist[ntpt];                        // Values of t in the histogram
#ifdef tdisthist
     G4double shist[ntpt];                        // Collected values in the histogram
#endif
     //G4double tMin=.0001;                         // in GeV^2  (pPb)
     //G4double tMin=.00001;                        // in GeV^2  (pd)
     //G4double tMin=.000005;                       // in GeV^2 (pp_le)
     //G4double tMin=.0005;                         // in GeV^2 (pp_he4)
     G4double tMin=.0005;                         // in GeV^2 (pp_he3)
     ///G4double tMax=energy*pMass/GeV/GeV;          // 1/2 max -t value in GeV^2 (pp)
     ///if(tMax>9.)tMax=9.;
     G4double gev2=GeV*GeV;
     G4double tmt=(e0-mp)*mt/gev2;                // T_proj*M_targ (GeV^2)
     //G4double emt=e0*mt/gev2;                     // E_proj*M_targ (GeV^2) - Wrong Calc
     //G4double mp2=mp*mp/gev2;                     // M_proj^2 (GeV^2)
     //G4double mt2=mt*mt/gev2;                     // M_trag^2 (GeV^2)
     //G4double pl2=pmax*pmax/gev2;                 // P_proj^2 (GeV^2)
     //G4double tMax=4.8*pl2*mt2/(mt2+mp2+emt+emt); // max -t value in GeV^2 (pd?-Wrong)
     G4double tMax=1.*tmt;                        // max -t value in GeV^2 (pA)
     G4double tMaM=1.;                            // max_max -t value in GeV^2 (pA)
     G4cout<<"Test39: Mt="<<mt<<", Mp="<<mp<<", tmt="<<tmt<<", E="<<e0<<G4endl;
     if(tMax>tMaM) tMax=tMaM;
     G4double ltMin=std::log(tMin);
     G4double ltMax=std::log(tMax);
     G4double dlt=(ltMax-ltMin)/ntpt;
#ifdef debug
     G4cout<<"Test39: n="<<ntpt<<", ti="<<tMin<<", ta="<<tMax<<", dl="<<dlt<<G4endl;
#endif
     G4double hdlt=dlt/2;
     G4double beglt=ltMin+hdlt;
     for(G4int ti=0; ti<ntpt; ti++) thist[ti]=std::exp(beglt+ti*dlt);
#ifdef debug
     G4cout<<"Test39: beg="<<beglt<<", dl="<<dlt<<", t0="<<thist[0]<<G4endl;
#endif
     // Fake call to avoib GHAD complains
     G4ForceCondition* cond = new G4ForceCondition;
     *cond=NotForced;
     // End of the fake call
     // ************************** CROSS-SECTION *******************************
     // -----> For GHAD
     ///G4HadronCrossSections* HadrCS = new G4HadronCrossSections;// GHAD/DB of CrossSec's
     // -----> For HP
     //G4NeutronHPElasticFS* theFS = new G4NeutronHPElasticFS;
     //G4String dirName=getenv("NeutronHPCrossSections");
     //G4String tString = "/Elastic/";
     //dirName = dirName + tString;
     //G4NeutronHPChannel HadrCS;
     //HadrCS.Init((*(G4Element::GetElementTable()))[0], dirName);
     //while(!HadrCS.Register(theFS));
     //delete theFS;
     // ------> for CHIPS
     G4VQCrossSection* HadrCS = 0; // Prototype of CHIPS Elastic manager
     if     (pPDG==2212) HadrCS=G4QProtonElasticCrossSection::GetPointer();
     else if(pPDG==2112) HadrCS=G4QNeutronElasticCrossSection::GetPointer();  
     else if(pPDG== 211) HadrCS=G4QPionPlusElasticCrossSection::GetPointer();  
     else if(pPDG==-211) HadrCS=G4QPionMinusElasticCrossSection::GetPointer();  
     else if(pPDG== 321) HadrCS=G4QKaonPlusElasticCrossSection::GetPointer();  
     else if(pPDG==-321) HadrCS=G4QKaonMinusElasticCrossSection::GetPointer();  
     else if(pPDG== 130 || pPDG==310) HadrCS=G4QKaonMinusElasticCrossSection::GetPointer();
     else if(pPDG==3222) HadrCS=G4QHyperonPlusElasticCrossSection::GetPointer();  
     else if(pPDG>3110 && pPDG<3334) HadrCS=G4QHyperonElasticCrossSection::GetPointer();
     else if(pPDG>-3334&&pPDG<-1110) HadrCS=G4QAntiBaryonElasticCrossSection::GetPointer();
     else G4cout<<"*Warning*G4QElastic::Test39: wrong PDG="<<pPDG<<G4endl;
#ifdef csdebug
     // --- A temporary LOOP for calculation of total cross section ------------
     //G4double pMin=.02;                            // in GeV --> for protons
     G4double pMax=6000.;                            // in GeV
     //G4double pMax=600000000.;                     // in GeV
     //G4double pMax=10000000.;                      // in GeV --> np
     //G4double pMax=1000.;                           // in GeV --> np->inelastic
     ////G4double pMin=.03;                            // in GeV --> np->dg
     G4double pMin=.000004;                        // in GeV --> for neutrons
     //G4double pMax=700.;                           // in GeV
     G4int nic=50;
     G4double lpMin=std::log(pMin);
     G4double lpMax=std::log(pMax);
     G4double dlp=(lpMax-lpMin)/nic;
     G4double lmic=lpMin-dlp/2;
     G4double hMa=pMass/GeV;                       //Mass of a projectile in GeV
     G4double hMa2=hMa*hMa;
     G4cout<<"Test39: mi="<<lmic+dlp<<",ma="<<lmic+dlp*nic<<",d="<<dlp<<",n="<<nic<<", Z="
           <<tgZ<<", N="<<tgN<<", Element="<<*element<<G4endl;
     for(G4int ic=0; ic<nic; ic++)
     {
       lmic+=dlp;
       G4double mic=std::exp(lmic);                // current Momentum in GeV
       G4double p2=mic*mic;
       G4double ken=(std::sqrt(p2+hMa2)-hMa)*GeV;
       dParticle->SetKineticEnergy(ken);           // Fill Kinetic Energy of the projectile
       // GHAD Cross-section
       ///G4double CS = HadrCS->GetElasticCrossSection(dParticle,element);
       //HadrCS->SetCorrectInelasticNearZero(true);
       //G4double CS = HadrCS->GetInelasticCrossSection(dParticle,element);
       // HP_GHAD Cross-section
       //G4double CS = HadrCS.GetXsec(ken*GeV);
       //G4NeutronHPThermalBoost aThermalE;
       //G4double CS = HadrCS.GetXsec(aThermalE.GetThermalEnergy(aTrack,element,
       //                                                     material->GetTemperature()));
       // ==================== End of GHAD ==============================
       // CHIPS calculation by G4Q...ElasticCrossSection
       G4double CS = HadrCS->GetCrossSection(true, mic*GeV, tgZ, tgN, pPDG); 
       // ------ direct CHIPS approximation of elastic cross section --------
       //G4double sp=std::sqrt(mic);
       //G4double dl=lmic-3.;
       //G4double CS=2.648/p2/sp+(18.73+.6351*dl*dl+9./mic)/(1+.4186*lmic)/(1+.3953/p2/p2);
       //CS*=millibarn;
       // ------ end of direct CHIPS approximation (temporary, not necessary)--------
       gTrack->SetStep(step);            // Now step is included in the Track (see above)
       gTrack->SetKineticEnergy(ken);    // Duplication of Kin. Energy for the Track
       gTrack->SetTouchableHandle(touch);// Set Box touchable history
       G4double mfp=proc->GetMeanFreePath(*gTrack,0.1,cond); // Calculate the meanFreePath
       G4cout<<"Test39: P="<<mic<<" (GeV/c), CrosSec="<<CS/millibarn<<",MFP="<<mfp<<G4endl;
     }
     // --- End of the temporary LOOP for calculation of total cross section ------------
#endif
     if(HadrCS) dParticle->SetKineticEnergy(energy);// Fill the KineticEnergy of projectile
#ifdef tdisthist
     // --> For GHAD
     ///G4double reactCS = HadrCS->GetElasticCrossSection(dParticle,element);
     // --> For HP_GHAD
     //G4double reactCS = HadrCS.GetXsec(energy);
     // ==================== End of GHAD ==============================
     // --> For CHIPS ("false" of onlyCS to prepare parameters for differential cross-sect)
     // P=pmax is in GeV, but this is not necessary
     G4double reactCS = HadrCS->GetCrossSection(false, pmax, tgZ, tgN, pPDG);
     // --- End
     for(G4int is=0; is<ntpt; is++) shist[is]=0.; // Reset Collected values
#endif
     // ********************* END OF CROSS-SECTION ***************************************
#ifdef debug
     G4cout<<"Test39:---->, E="<<energy<<"(IU=MeV), CrossSect="<<reactCS/millibarn<<G4endl;
#endif
     gTrack->SetStep(step);            // Now step is included in the Track (see above)
     gTrack->SetKineticEnergy(energy); // Duplication of Kin. Energy for the Track (?!)
     gTrack->SetTouchableHandle(touch);// Set Box touchable history
     G4double mfp=proc->GetMeanFreePath(*gTrack,0.1,cond);
#ifdef debug
     G4cout<<"Test39: Mean free path "<<mfp/meter<<"(m)"<<G4endl;
#endif
     if(mfp<1.e20) for (G4int iter=0; iter<nEvt; iter++)
     {
#ifdef debug
      G4cout<<"Test39: ### "<<iter<< "-th event starts.### energy(IU)="<<energy<<G4endl;
#endif
      if(!(iter%100000)&&iter)G4cout<<"*=>TEST39: "<<iter<<" events are simulated"<<G4endl;

      // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ RANDOM ENERGY AND RANDOM n/p SEQUENCE
      //energy=G4UniformRand()*.3*GeV;
      //gTrack=pTrack;
      //if(G4UniformRand()>.5) gTrack=nTrack;
      // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      gTrack->SetStep(step);            // Now step is included in the Track (see above)
      gTrack->SetKineticEnergy(energy); // Duplication of Kin. Energy for the Track (?!)
      gTrack->SetTouchableHandle(touch);// Set Box touchable history
#ifdef debug
      G4cout<<"Test39: Before PostStepDoIt, T="<<energy<<G4endl;
#endif
      G4bool flead=false;   // existance of the leading particle
      G4bool fakeel=true;   // bad event (must be changed to good inside while)
      G4int fcn=0;
      while(fakeel)
      {
        fcn++;
#ifdef debug
        G4cout<<"Test39: Before the fake proc->GetMeanFreePath call"<<G4endl;
#endif
        // CHIPS does not need it: keep only for GHAD
        proc->GetMeanFreePath(*gTrack,0.1,cond); // Fake call to avoid complains of GHAD
#ifdef debug
        G4cout<<"Test39: After the fake proc->GetMeanFreePath call"<<G4endl;
#endif
        //..............................................................................
        aChange = static_cast<G4ParticleChange*>(proc->PostStepDoIt(*gTrack,*step));
        G4int nS = aChange->GetNumberOfSecondaries();
        G4TrackStatus lead = aChange->GetTrackStatus();
#ifdef fdebug
        G4cout<<"Test39:AfterPostStepDoIt #"<<fcn<<",NSec="<<nS<<",Lead=0?="<<lead<<G4endl;
#endif
        flead=false;
        G4int bn1=-2727;
        G4int ch1=-2727;
        if(nS)
        {
          const G4DynamicParticle* sec0=aChange->GetSecondary(0)->GetDynamicParticle();
          bn1=sec0->GetDefinition()->GetBaryonNumber();
          ch1=static_cast<G4int>(sec0->GetDefinition()->GetPDGCharge()); 
#ifdef fdebug
          G4double Ekin0=sec0->GetKineticEnergy();
          G4int PDG0=sec0->GetDefinition()->GetPDGEncoding();
          G4cout<<"Test39:Rec, CH="<<ch1<<",BN="<<bn1<<",PDG="<<PDG0<<",T="<<Ekin0<<G4endl;
#endif
        }
        G4int ch2=-2727;
        G4int bn2=-2727;
        if(lead==fAlive)
        {
          flead=true;
#ifdef fdebug
          G4double len=aChange->GetEnergy();  // Energy of the leading particle in MeV
          G4double lma=aChange->GetMass();
          G4double lch=aChange->GetCharge();
          G4cout<<"Test39: Leading particle Ekin="<<len<<",M="<<lma<<",Chrg="<<lch<<G4endl;
          G4cout<<"Test39: b1="<<bn1<<" = tgA="<<tgA<<", c1="<<ch1<<" = tgZ="<<tgZ<<G4endl;
#endif
          if(bn1==tgA&&ch1==tgZ) fakeel=false;
          else if(bn1==-2727) misTG++;
          else badTG++;
        }
        else
        {
          if(nS>1)
          {
            const G4DynamicParticle* sec1=aChange->GetSecondary(1)->GetDynamicParticle();
            ch2=static_cast<G4int>(sec1->GetDefinition()->GetPDGCharge());
            bn2=sec1->GetDefinition()->GetBaryonNumber();
#ifdef fdebug
            G4double Ekin1=sec1->GetKineticEnergy();
            G4int PDG1=sec1->GetDefinition()->GetPDGEncoding();
            G4cout<<"Test39: 2nd secondary BN="<<bn2<<",PDG="<<PDG1<<",Ek="<<Ekin1<<G4endl;
#endif
          }
          if( (bn1==bnp&&ch1==chp && bn2==tgA&&ch2==tgZ) ||
              (bn2==bnp&&ch2==chp && bn1==tgA&&ch1==tgZ) )     fakeel=false;
          else if(bn1==-2727) zeroO++;
          else if(bn2==-2727) alonO++;
          else badOT++;
        }
#ifdef debug
        G4cout<<"Test39:b1="<<bn1<<",b2="<<bn2<<",p="<<bnp<<",t="<<tgA<<",c="<<fcn
              <<", fake="<<fakeel<<G4endl;
#endif
        if(fakeel) // Cleaning of the bad secondaries
        {
          G4int nS = aChange->GetNumberOfSecondaries();
          for(G4int ides=0; ides<nS; ides++) delete aChange->GetSecondary(ides);
          aChange->Clear();
          if(alonO+zeroO+badOT+misTG+badTG>100000+nEvt) break;
        }
      } // End of while
      if(fakeel) break;
      goodE++;
#ifdef debug
      G4cout<<"Test39: ***OUT OF LOOP*** c="<<fcn<<" ****************************"<<G4endl;
#endif
      G4int nSec = aChange->GetNumberOfSecondaries();
#ifdef debug
      G4cout<<"Test39: "<<nSec<<" secondaries in the event"<<G4endl;
#endif
      G4double totCharge = totC;
      //............................ONLY for CHIPS................................
      G4int    curN=proc->GetNumberOfNeutronsInTarget(); // Works only for CHIPS
      G4int    dBN = curN-tgN;
      G4int    totBaryN = totBN+dBN;
      G4int    curPDG=tPDG+dBN;
      G4LorentzVector Residual=proc->GetEnegryMomentumConservation(); // Only for CHIPS
      //............................ONLY FOR GHAD...............................
      ///G4int    totBaryN = totBN;                   // Substitute for not CHIPS
      ///G4int    curPDG=tPDG;                        // Substitute for not CHIPS
      //........................................................................
      G4double curM=G4QPDGCode(curPDG).GetMass();  // Update mass of the TargetNucleus
      totSum = G4LorentzVector(0., 0., pmax, et+curM-mt);
      //
#ifdef debug
      G4double de = aChange->GetLocalEnergyDeposit();// Init TotalEnergy by EnergyDeposit
      G4cout<<"Test39: "<<nSec<<" secondary particles are generated, dE="<<de<<G4endl;
#endif
      // @@ ----------------------- Begin
      G4double weight = aChange->GetSecondary(0)->GetDynamicParticle()->GetKineticEnergy();
#ifdef debug
      G4cout<<"Test39:----------------: Weight="<<weight<<G4endl;
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
      G4cout<<"Test39:----DONE^^^^^^^*******^^^^^^^^:ir="<<iter<<": #ofH="<<nSec<<G4endl;
      if(!(iter%100)) G4cerr<<"#"<<iter<<G4endl;
#endif
      G4bool alarm=false;
      // @@ G4bool rad=false;
      // @@ G4bool hyp=false;
      // @@ G4bool badPDG=false;
      // ------- LOOP over secondary particles -------
      G4int ii=0;
      if(flead) ii=-1;
      for(G4int i=ii; i<nSec; i++)
      {
        if(i>-1)
        {
          sec = aChange->GetSecondary(i)->GetDynamicParticle();
          pd  = sec->GetDefinition();
          c   = pd->GetPDGEncoding();
          if(!c)
          {
            G4int chrg=static_cast<G4int>(pd->GetPDGCharge());
            G4int bary=static_cast<G4int>(pd->GetBaryonNumber());
            c=90000000+chrg*999+bary;
          }
          m   = pd->GetPDGMass();
          mom = sec->GetMomentumDirection();
          e   = sec->GetKineticEnergy();
        }
        else                                // leading particle
        {
          pd  = part;                       // @@ if there is a mistake...?
          c   = pd->GetPDGEncoding();
          if(!c)
          {
            G4int chrg=static_cast<G4int>(pd->GetPDGCharge());
            G4int bary=static_cast<G4int>(pd->GetBaryonNumber());
            c=90000000+chrg*999+bary;
          }
          m   = aChange->GetMass();
          mom = *(aChange->GetMomentumDirection());
          e   = aChange->GetEnergy();
        }
        if (e < 0.0)
        {
          G4cerr<<"**Test39:Event#"<<iter<<",Hadron#"<<i<<", E="<<e<<" <0 (Set 0)"<<G4endl;
          e = 0.0;
        }
        // for exclusive reaction 2 particles in final state
        p = std::sqrt(e*(e + m + m));
        mom *= p;
        //if(i==1) // Means the target secondary for the ellastic scattering
        //{
        //  G4double t=e*m;
        //  t=t+t;
        //  G4int tb=static_cast<G4int>(t/dT);
        //  if(tb>nT1) tb=nT1;
        //  tSig[tb]++;
        //}
        lorV = G4LorentzVector(mom, e+m);    // "e" is a Kinetic energy!
        totSum -= lorV;
        //if(fabs(m-lorV.m())>.005&&1>2) // @@ Temporary closed
        if(std::fabs(m-lorV.m())>.005) // @@ Temporary closed
        {
          G4cerr<<"***Test39: m="<<lorV.m()<<" # "<<m<<", d="<<lorV.m()-m<<G4endl;
          alarm=true;
        }
        if(!(lorV.e()>=0||lorV.e()<0)   || !(lorV.px()>=0||lorV.px()<0) ||
           !(lorV.py()>=0||lorV.py()<0) || !(lorV.pz()>=0||lorV.pz()<0))
        {
          G4cerr<<"***Test39: NAN in LorentzVector="<<lorV<<G4endl;
          alarm=true;
        }
        if(c==90000002||c==90002000||c==92000000)
        {
          G4cout<<"***Test39:***Dibaryon *** i="<<i<<", PDG="<<c<<G4endl;
          alarm=true;
        }
        if(c==90000003||c==90003000||c==93000000)
        {
          G4cout<<"***Test39:***Tribaryon *** i="<<i<<", PDG="<<c<<G4endl;
          alarm=true;
        }
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
        if(c==22) nGamma++;                            // Gammas
        if(c==22) EGamma+=e;                           // Energy of gammas
        G4int cCG=0;
        G4int cBN=0;
        if(std::abs(c)>99)                               // Do not count charge of leptons
        {
          cCG=static_cast<G4int>(pd->GetPDGCharge());
          cBN=static_cast<G4int>(pd->GetBaryonNumber());
          totCharge-=cCG;
          totBaryN-=cBN;
        }
#ifdef pdebug
        G4cout<<"Test39:#"<<i<<",PDG="<<c<<",C="<<cCG<<",B="<<cBN<<",4M="<<lorV<<m<<",T="
              <<lorV.e()-m<<G4endl;
#endif
        // ********************* calculation of t and filling histogram *****************
#ifdef tdisthist
        G4LorentzVector t4M=incident4M-lorV;
        G4double tval=-t4M.m2()/GeV/GeV;
#ifdef pdebug
        G4cout<<"Ts39:c="<<c<<",t="<<tval<<",i="<<tMin<<",a="<<tMax<<",s="<<lorV<<G4endl;
#endif
        if(c==pPDG && tval>tMin && tval<tMax) // fill the histogram
        {
          G4int nchan=static_cast<G4int>((std::log(tval)-ltMin)/dlt);
          if(nchan<0) nchan=0;
          if(nchan>=ntpt) nchan=ntpt-1;
          shist[nchan]+=1./tval; // 1/t in GeV^2
#ifdef pdebug
          G4cout<<"Test39: shist["<<nchan<<"]="<<shist[nchan]<<G4endl;
#endif
        }
#endif
        // *****************End of calculation of t and filling histogram ***************
        //delete aChange->GetSecondary(i);
      } // End of the LOOP over secondaries
      // delete secondaries in the end of the event         
      G4double ss=std::fabs(totSum.t())+std::fabs(totSum.x())+std::fabs(totSum.y())+
                  std::fabs(totSum.z());
      //G4double sr=std::fabs(Residual.t())+std::fabs(Residual.x())+
      //            std::fabs(Residual.y())+std::fabs(Residual.z());    
#ifdef pdebug
      G4cout<<">TEST39:r4M="<<totSum<<ss<<",rCh="<<totCharge<<",rBaryN="<<totBaryN<<G4endl;
#endif
      ///if (1>2) // @@ The check is temporary closed LHEP
      //if (ss>.27) // Not CHIPS, but conservs energy
      if (totCharge ||totBaryN || ss>.27 || alarm || (nGamma && !EGamma)) // for CHIPS
      {
        totSum = G4LorentzVector(0., 0., pmax, et);
        G4cerr<<"**Test39:#"<<iter<<": n="<<nSec<<", 4M="<<totSum<<", Charge="<<totCharge
        ///      <<", BaryN="<<totBaryN<<",D2="<<ss<<G4endl; // Not CHIPS
              <<",BaryN="<<totBaryN<<", R="<<Residual<<",D2="<<ss<<",nN="<<curN<<G4endl;
        if(nGamma&&!EGamma)G4cerr<<"***Test39: Egamma=0"<<G4endl;
        G4int indi=0;
        if(flead) indi=-1;
        for (G4int indx=indi; indx<nSec; indx++)
        {
          if(indx>-1)
          {
            sec = aChange->GetSecondary(indx)->GetDynamicParticle();
            pd  = sec->GetDefinition();
            c   = pd->GetPDGEncoding();
            if(!c)
            {
              G4int chrg=static_cast<G4int>(pd->GetPDGCharge());
              G4int bary=static_cast<G4int>(pd->GetBaryonNumber());
              c=90000000+chrg*999+bary;
            }
            m   = pd->GetPDGMass();
            mom = sec->GetMomentumDirection();
            G4QPDGCode cQPDG(c);
            sm   = cQPDG.GetMass();
            e   = sec->GetKineticEnergy();
          }
          else                                // leading particle
          {
            pd  = part;
            c   = pd->GetPDGEncoding();
            if(!c)
            {
              G4int chrg=static_cast<G4int>(pd->GetPDGCharge());
              G4int bary=static_cast<G4int>(pd->GetBaryonNumber());
              c=90000000+chrg*999+bary;
            }
            m   = aChange->GetMass();
            mom = *(aChange->GetMomentumDirection());
            G4QPDGCode cQPDG(c);
            sm  = cQPDG.GetMass();
            e   = aChange->GetEnergy();
          }
          p = std::sqrt(e*(e + m + m));
          mom *= p;
          lorV = G4LorentzVector(mom, e + m);    // "e" is a Kinetic energy!
          totSum -= lorV;
          G4cerr<<"Test39:#"<<indx<<",PDG="<<c<<",m="<<m<<",4M="<<lorV<<",T="<<e
                <<", d4M="<<totSum<<G4endl;
        }
        if(sr>.27)
          G4Exception("***Test39: ALARM or baryn/charge/energy/momentum is not conserved");
      }
#ifndef nout
      ntp->FillEvt(aChange); // Fill the simulated event in the ASCII "ntuple"
#endif
      // =============== May be print it here if it is not zero...
      for(G4int ides=0; ides<nSec; ides++) delete aChange->GetSecondary(ides);
      aChange->Clear();
#ifdef debug
      G4cout<<"Test39:--->>> After ntp.FillEvt, ntpt="<<ntpt<<G4endl;
#endif
     } // End of the LOOP over events
     else
     {
       G4cout<<"Test39: Simulation is skipped because of big MeanFreePath="<<mfp<<G4endl;
       goodE=1;
     }
     // Stop the timer to estimate the speed of the generation
     timer->Stop();
     G4cout<<"Test39:CalculationTimePerEvent="<<timer->GetUserElapsed()/nEvt<<" s"<<G4endl;
     delete timer;
     G4double tS=goodE+misTG+badTG+zeroO+alonO+badOT;
     G4cout<<"Test39: g="<<goodE/tS<<" mTG="<<misTG/tS<<" bTG="<<badTG/tS
           <<" 0="<<zeroO/tS<<" 1="<<alonO/tS<<" b="<<badOT/tS<<G4endl;
     // *************** Resulting print for the collected t-histograms ***************
#ifdef tdisthist
     G4double csmb=reactCS/millibarn/goodE/dlt;
     for(G4int ti=0; ti<ntpt; ti++)
     {
       shist[ti]*=csmb;
       G4cout<<"Test39: t="<<thist[ti]<<" ,dsig/dt="<<shist[ti]<<G4endl;
     }
#endif
     // ************ End of resulting print for the collected t-histograms ***************
#ifdef ppdebug
     G4double pGeV=pmax/1000.;
     G4double alp=std::log(pGeV/(1.+1./pGeV/pGeV/pGeV));
     G4double exr=1.-.822/(1.+std::exp(-(alp-.2)*1.15));
     G4cout<<"Test39:EndOfRun,p="<<pmax<<",dT="<<dTot<<",r="<<dEl/dTot<<",e="<<exr
           <<",d="<<exr-dEl/dTot<<",ra="<<(exr-.178)/.822<<G4endl;
     for(G4int ir=0; ir<nT; ir++) G4cout<<tVal[ir]<<" "<<tSig[ir]<<G4endl;// Print t-vectors
#endif
#ifdef pverb
     G4cerr<<"Test39: ########## End of run ##########"<<G4endl;
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
  //G4cout << "Test39: After delete process etc." << G4endl;
  //delete phys;                        // Absolete physics class
  //delete cmaterial;                   // Temporary material definition pointer (?)
#ifndef nout
  //G4cout << "Test39: Before ntp" << G4endl;
  delete ntp; // Delete the class to fill the#of events
  //G4cout << "Test39: After ntp" << G4endl;
#endif
  delete runManager;
#ifdef pverb
  G4cout << "###### End of Test39 #####" << G4endl;
#endif
  //exit(1); // Never do this !
  //return no_of_errors;
  //return 0;
  //abort();
  //return EXIT_SUCCESS;
}
