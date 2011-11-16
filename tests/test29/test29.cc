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
//      File name:     Test29 (Comparison of G4QString Results with Data)
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
#define inter
//#define pverb
//#define pscan
//#define smear
//#define escan
//#define debug
//#define rdebug
//#define mtst
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
#include "G4QCaptureAtRest.hh"
#include "G4MuonMinusCaptureAtRest.hh"
//#include "G4QuasmonString.hh"

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
#include "G4Ellipsoid.hh"
#include "G4PVPlacement.hh"
#include "G4Step.hh"
#include "G4GRSVolume.hh"
#include "G4GRSSolid.hh"

#include "Test29Physics.hh"
#include "Test29PhysicsList.hh"

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

// Temporary AIDA can be not used in test29
//#include "AIDA/AIDA.h"

// ASCII Pseudo NTUPLE 
#ifndef nout
#include "G4QHBook.hh"
#endif

#include "G4Timer.hh"
#include "time.h"
#include <cmath>

//int main(int argc, char** argv)
extern "C" double drand();
int main()
{
  const G4int nTg=8;   // Length of the target list for the Performance test
  G4int tli[nTg]={90001000,90002002,90003004,90007007,90013014,90027032,90047060,90092146};
  G4String tnm[nTg]={"Hydrogen","Helium","Lithium","Nitrogen","Aluminum","Cobalt","Silver"
                     ,"Uranium"};
  G4String tsy[nTg]={"1H","He","7Li","14N","27Al","59Co","107Ag","238U"};
  G4Material* mat[nTg]={0,0,0,0,0,0,0,0};
  const G4int nPr=11;  // Length of the projectile list for the Performance test
  G4int pli[nPr] = {-2212, -211, -321, 13, 15, 3112, 3312, 3334, 2112, -2112, -3222};
  const G4int nAZ=270;  // Dimension of the table
  const G4int mAZ=266;  // Maximum filled A (at present). Must be mAZ<nAZ
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
  if(mAZ>=nAZ) G4cout<<"***Test29: Too big mAZ="<<mAZ<<", nAZ="<<nAZ<<G4endl;

  // Run manager
  G4RunManager* runManager = new G4RunManager;
  runManager->SetUserInitialization(new Test29PhysicsList);
  G4StateManager::GetStateManager()->SetNewState(G4State_Init); // To let create ions
  G4ParticleDefinition* ionDefinition=0;
  ionDefinition=G4ParticleTable::GetParticleTable()->FindIon(6,12,0,6);
  if(!ionDefinition)
  {
    G4cerr<<"*** Error! *** Test29:(6,6) ion can not be defined"<<G4endl;
    return 0;
  }
  else G4cout<<"Test29: (6,6) ion is OK, Run State="<<G4StateManager::GetStateManager()->
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
  G4cout<<"Test29: Prepare G4QHBook files or ntuples"<<G4endl;
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
  inFile>>temperature>>ssin2g>>eteps>>nop>>momb>>enb>>pPDG>>tPDG>>nEvt>>fN>>fD>>cP>>rM>>sA;
  G4cout<<"Test29:Par's: T="<<temperature<<",S="<<ssin2g<<",Eta="<<eteps<<",nop="
        <<nop<<",p="<<momb<<",e="<<enb<<",pPDG="<<pPDG<<",tPDG="<<tPDG<<",nEv="
        <<nEvt<<",fN="<<fN<<",fD="<<fD<<",cP="<<cP<<",rM="<<rM<<",sA="<<sA<<G4endl;
  //-------- Initialize CHIPS
  G4QCHIPSWorld* theW=G4QCHIPSWorld::Get(); // (?)
  theW->GetParticles(nop);           // Create CHIPS World of nop particles (?)
  G4QCaptureAtRest::SetParameters(temperature,ssin2g,eteps,fN,fD,cP,rM,nop,sA);
  G4QCaptureAtRest::SetManual();
#ifdef mtst
  // @@ Temporary mass test -- Begin
  G4int vZ=1;
  G4int vN=6;
  G4int vA=vZ+vN;
  G4cout<<"Test29:MTst,Z="<<vZ<<",N="<<vN<<",T="<<G4NucleiPropertiesTable::IsInTable(vZ,vA)
        <<",tV="<<G4NucleiProperties::GetNuclearMass(vA,vZ)<<",calcV="
        <<G4ParticleTable::GetParticleTable()->FindIon(vZ,vA,0,vZ)->GetPDGMass()<<",chipM="
        <<G4QPDGCode(90000000+vZ*1000+vN).GetMass()<<G4endl;
  // @@ Temporary mass test -- End
#endif
  // Not necessary for CHIPS, only for GHAD: fake force condition
  G4ForceCondition* cond = new G4ForceCondition;
  *cond=NotForced;
//#endif
  // Construct all particles of G4 instead of sraight forward one by one
  ///G4ParticlePhysics* allParticles = new G4ParticlePhysics(); // Short cut from IntPhysL
  ///allParticles->ConstructParticle();
  //
  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();
  
  // gamma
  G4Gamma::GammaDefinition();
  
  // optical photon
  G4OpticalPhoton::OpticalPhotonDefinition();

  // leptons


  // --------------------------------------- End of the Particle definition ---***---
  // *********** Find tgZ and tgN from tPDG *************
  G4int tgZ=(tPDG-90000000)/1000;
  G4int tgN=tPDG-90000000-tgZ*1000;
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
      G4cout<<"Test29:-->Material("<<tgZ<<","<<tgN<<"):"<<tnm[tgi]<<","<<tsy[tgi]<<G4endl;
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
  G4cout<<"Test29:--- Material "<<material->GetName()<<" is defined ---" << G4endl;
#endif
  if(!material)
  {
    G4cout<<"Test29: Last Material "<<material->GetName()<<" is not defined."<<G4endl;
    exit(1);
  }
  G4double nx = 0.;
  G4double ny = 0.;
  G4double nz = 0.;
  G4int npart=1;                                      // By default only one particle
  if(!pPDG) npart=nPr;                                // Make a LOOP ove all particles
  // Different Process Managers are used for the atRest and onFlight processes
  //G4ParticleDefinition* part=G4GenericIon::GenericIonDefinition();
  //G4ProcessManager* man = new G4ProcessManager(part);
  //part->SetProcessManager(man);
  //
  //G4VRestProcess* proc = new G4QCaptureAtRest;
  G4QCaptureAtRest* proc = new G4QCaptureAtRest;   // CHIPS
  ///G4MuonMinusCaptureAtRest* proc = new G4MuonMinusCaptureAtRest; // GHAD
  if(!proc)
  {
    G4cout<<"Tst29: there is no G4QCaptureAtRest process"<<G4endl;
    exit(1);
  }
#ifdef debug
  G4cout<<"Test29:--***-- process is created --***--" << G4endl; // only one run
#endif
  // Only for CHIPS
  proc->SetParameters(temperature, ssin2g, eteps, fN, fD, cP, rM, nop, sA);
  //
  // Geometry

  G4double dimX = 100.0*cm;
  G4double dimY = 100.0*cm;
  G4double dimZ = 100.0*cm;

  G4Box* sFrame = new G4Box ("Box",dimX, dimY, dimZ);
  G4LogicalVolume* lFrame = new G4LogicalVolume(sFrame,material,"Box",0,0,0);
  G4PVPlacement* pFrame = new G4PVPlacement(0,G4ThreeVector(),"Box",lFrame,0,false,0);

  assert(pFrame);
  //#ifdef pverb
  G4cout<<"Test29: Geometry is defined "<<G4endl;
  //#endif
  //const G4RotationMatrix* rotation=pFrame->GetFrameRotation();
  G4Navigator* nav = new G4Navigator;
  //#ifdef pverb
  G4cout<<"Test29: The Navigator is defined "<<G4endl;
  //#endif
  nav->SetWorldVolume(pFrame);
  //#ifdef pverb
  G4cout<<"Test29: The Box frame is set "<<G4endl;
  //#endif
  //G4VTouchable* vtouch = nav->CreateTouchableHistory();
  G4TouchableHandle touch(nav->CreateTouchableHistory());
  //#ifdef pverb
  G4cout<<"Test29: The TouchableHandle is defined "<<G4endl;
  //#endif
  G4int nTot=npart*tgm;
  G4int nCur=0;
  for(G4int pnb=0; pnb<npart; pnb++) // LOOP over particles
  {
   if (npart>1) pPDG=pli[pnb];
   G4QContent pQC=G4QPDGCode(pPDG).GetQuarkContent();
   G4int      cp = pQC.GetCharge();
   if (pPDG==11||pPDG==13||pPDG==15)cp=-1;// QC/Charge of leptons isn't supported by CHIPS
   if (pPDG==-11||pPDG==-13||pPDG==-15) cp=1;// QC/Charge isn't supported by CHIPS
   G4ParticleDefinition* part=G4AntiSigmaPlus::AntiSigmaPlus();//DefaultProj is AntiSigma+
   if(pPDG==-2212) part=G4AntiProton::AntiProton();  //Definition of antiProton projectile
   else if(pPDG==-211) part=G4PionMinus::PionMinus();// Definition of the Pi- projectile
   else if(pPDG==-321) part=G4KaonMinus::KaonMinus();// Definition of the K- projectile
   else if(pPDG==13) part=G4MuonMinus::MuonMinus();  // Definition of the mu- projectile
   else if(pPDG==15) part=G4TauMinus::TauMinus();    // Definition of the tau- projectile
   else if(pPDG==3112) part=G4SigmaMinus::SigmaMinus(); // Definition of Sigma- projectile
   else if(pPDG==3312) part=G4XiMinus::XiMinus();    // Definition of the Xi- projectile
   else if(pPDG==3334) part=G4OmegaMinus::OmegaMinus(); // Definition of Omega- projectile
   else if(pPDG==2112) part=G4Neutron::Neutron();    // Definition of Neutron projectile
   else if(pPDG==-2112) part=G4AntiNeutron::AntiNeutron();// Define AntiNeutron projectile
   else if(pPDG!=-3222) // Leave defaulf definition for Anti Sigma+ projectile
   {
     G4ExceptionDescription desc;
     desc<<"***Test29: At Rest Process is called for not negative particle"<<G4endl;
     desc<<"***Test29: pPDG="<<pPDG<<" is a PDG code of not negative particle"<<G4endl;
     G4Exception("***Test29", "Test29-01",FatalException, desc);
   }
   G4double pMass = part->GetPDGMass();                 // Mass of the projectile
   //
   G4ThreeVector aPosition(nx*mm, ny*mm, nz*mm);
   // Create a DynamicParticle
   G4double  energy   = 0.*MeV;                              // 0 GeV particle energy(Cap)
   G4DynamicParticle* dParticle = new G4DynamicParticle(part,aDirection,energy);
   //dParticle->SetKineticEnergy(energy);// Fill the Kinetic Energy of the projectile
   // General Track definition
   G4Track* gTrack = new G4Track(dParticle,aTime,aPosition);// Track definition (NoTarg)

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
      G4cout<<"Test29: Material="<<material->GetName()<<", Element[0]="<<curEl->GetName()
            <<",A[0]="<<(*(curEl->GetIsotopeVector()))[0]->GetN()<<" is selected."<<G4endl;
    }
    G4cout<<"Test29:NewRun: Targ="<<tPDG<<",Proj="<<pPDG<<", "<<nCur<<" of "<<nTot<<G4endl;
    G4int    bnp=pQC.GetBaryonNumber();
    G4QContent tQC=G4QPDGCode(tPDG).GetQuarkContent();
    G4int    ct=tQC.GetCharge();
    G4int    bnt=tQC.GetBaryonNumber();
#ifdef debug
    G4cout<<"Test29: pQC"<<pQC<<", pch="<<cp<<", tQC"<<tQC<<", tch="<<ct<<G4endl;
#endif
    G4int    totC=cp+ct;
    G4int    totBN=bnp+bnt;
    // ---------------------- PSEUDO-STEP DEFINISIONS -------------------------
#ifdef debug
    G4double totCN  = tgZ+part->GetPDGCharge();
    G4int    totBNN = tgZ+tgN+part->GetBaryonNumber();
    G4cout<<"Test29:tC="<<totC<<"?="<<totCN<<",tB="<<totBN<<"?="<<totBNN<<G4endl;
#endif
    G4double theStep   = 0.01*micrometer;
#ifdef pverb
    G4cout<<"Test29: The particle: "<<part->GetParticleName()<<G4endl;
    G4cout<<"Test29: The material: "<<material->GetName()<<G4endl;
    G4cout<<"Test29: The step:     "<<theStep/mm<<" mm"<<G4endl;
    G4cout<<"Test29: The position: "<<aPosition/mm<<" mm"<<G4endl;
    G4cout<<"Test29: The direction:"<<aDirection<<G4endl;
    G4cout<<"Test29: The time:     "<<aTime/ns<<" ns"<<G4endl;
#endif
    // Step Definition
    G4Step* step = new G4Step();
    step->SetTrack(gTrack);          // Step is initialized by the Track (?)

    G4StepPoint* aPoint = new G4StepPoint(); // It cant be initialized right away (!?)
    /////G4TouchableHandle touch = aPoint->GetTouchableHandle();
    aPoint->SetPosition(aPosition);
    aPoint->SetMaterial(material);
    G4double safety = 10000.*cm;
    aPoint->SetSafety(safety);
    step->SetPreStepPoint(aPoint);   // Begin of the step
#ifdef pverb
    G4cout<<"Test29: Starting point is defined and filled in the step "<<G4endl;
#endif
    //G4StepPoint* bPoint = new G4StepPoint(*aPoint); // @@ This does not work
#ifdef pverb
    //G4cout<<"Test29: A copy of the start point to the end point is made "<<G4endl;
#endif
    G4StepPoint* bPoint = new G4StepPoint(); // It cant be initialized right away (!?)
    bPoint->SetMaterial(material);
    bPoint->SetSafety(safety);
    G4ThreeVector bPosition = aDirection*theStep+aPosition; // aDirection is defined byCard
    bPoint->SetPosition(bPosition);
    step->SetPostStepPoint(bPoint);  // End of the step
    step->SetStepLength(theStep);    // Step is set byCard above
#ifdef pverb
    G4cout<<"Test29: The end point is defined and filled in the step "<<G4endl;
#endif
    G4Timer* timer = new G4Timer();
    timer->Start();
#ifdef pverb
    G4cout<<"Test29: Rub is started, timer is started, kinEnergy="<<energy<<G4endl;
#endif
    const G4DynamicParticle* sec = 0;
    G4ParticleDefinition* pd = 0;
    G4ThreeVector  mom;
    G4LorentzVector totSum, lorV;
    G4double e, p, m, sm;
    // @@ G4double px, py, pz, pt;
    G4VParticleChange* aChange = 0;
    G4double e0 = energy+pMass;
    G4double pmax=std::sqrt(e0*e0-pMass*pMass);
    // Randomization loop: cycle random generator, using 2 lower digits in nEvt
    G4int    iRandCount = nEvt%100;
    while (iRandCount>0)                // Shift of the RNDN values 
    {
      G4int nr = static_cast<G4int>(10*G4UniformRand())+1;        // Fake cicle number
      for(G4int j=0; j<nr; j++)  G4UniformRand();                 // Fake calls
      iRandCount--;
    }
    iRandCount = nEvt%100;
#ifdef pverb
    G4cout<<"Test29: Before the event loop, nEvents= "<<nEvt<<G4endl;
#endif
    for (G4int iter=0; iter<nEvt; iter++)
    {
      G4int nr = static_cast<G4int>(iRandCount*G4UniformRand())+1;// Fake cicle number
      for(G4int j=0; j<nr; j++) G4UniformRand();                  // Fake calls
#ifdef debug
      G4cout<<"Test29: ### "<<iter<< "-th event starts.###"<<G4endl;
#endif
      if(!(iter%1000) && iter)
         G4cout<<"****************>>TEST29: "<<iter<<" events are simulated"<<G4endl;

      gTrack->SetStep(step);            // Now step is included in the Track (see above)
      gTrack->SetKineticEnergy(energy); // Duplication of Kin. Energy for the Track (?!)
      gTrack->SetTouchableHandle(touch);// Set Box touchable history

      //aChange = proc->PostStepDoIt(*gTrack,*step); // For onFlight
      aChange = proc->AtRestDoIt(*gTrack,*step);  // For At Rest (step is defined twice?)
#ifdef debug
      G4cout<<"Test29:--***@@@***--After AtRestDoIt--***@@@***--" << G4endl;// only one run
#endif
      G4double totCharge = totC;
      // For CHIPS
      G4int    curN=proc->GetNumberOfNeutronsInTarget();
      // --- for GHAD ---
      ///G4int    curN=tgN;
      // --- End
      if(curN!=tgN) G4cerr<<"******Test29: tgN="<<tgN<<" # curN="<<curN<<G4endl;
      G4int    dBN = curN-tgN;
      G4int    totBaryN = totBN+dBN;
      G4int    curPDG=tPDG+dBN;
      G4double curM=G4QPDGCode(curPDG).GetMass(); // Update #ofNeutrons in theTargetNucleus
#ifdef debug
      G4cout<<"Test29:tN="<<tgN<<",cN="<<curN<<",tPDG="<<tPDG<<",cPDG="<<curPDG<<",cM="
            <<curM<<", E0="<<e0<<G4endl;
#endif
      // Only for CHIPS
      G4LorentzVector Residual=proc->GetEnegryMomentumConservation();
      // for GHAD
      ///G4LorentzVector Residual(0.,0.,0.,0.);
      //
      //G4double de = aChange->GetLocalEnergyDeposit();// Init TotalEnergy by EnergyDeposit
      G4int nSec = aChange->GetNumberOfSecondaries();
      G4double EnergyDep = aChange->GetLocalEnergyDeposit();
      totSum = G4LorentzVector(0., 0., pmax, e0+curM-EnergyDep); // For the En/Mom check
#ifdef debug
      G4cout<<"Test29:"<<nSec<<" secondaries, dE="<<EnergyDep<<",4M="<<totSum<<G4endl;
#endif
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
      G4int    c=0;    // PDG Code of the particle
#ifdef pdebug
      G4cout<<"Test29:----DONE^^^^^^^*********^^^^^^^: ir="<<iter<<", #ofH="<<nSec<<G4endl;
      if(!(iter%1000)) G4cerr<<"#"<<iter<<G4endl;
#endif
      G4bool alarm=false;
      // @@ G4bool rad=false;
      // @@ G4bool hyp=false;
      // @@ G4bool badPDG=false;
      //-- LOOP over secondary particles (test of En/Mom & Quantum Numbers conservation) --
      for(G4int i=0; i<nSec; i++)
      {
        sec = aChange->GetSecondary(i)->GetDynamicParticle();
        pd  = sec->GetDefinition();
        c   = pd->GetPDGEncoding();
        if(!c)
        {
          G4int chrg=static_cast<G4int>(pd->GetPDGCharge());
          G4int bary=static_cast<G4int>(pd->GetBaryonNumber());
          c=1000000000+chrg*10000+bary*10; // New PDG2006 code
        }
        m   = pd->GetPDGMass();
#ifdef rdebug
        G4int ac=std::abs(c);
        if(ac<10000 && ac>100)
        {
          G4int dc=0;
          if(ac<1000)                  // Mesons
          {
            dc=ac%100;
            dc=ac-dc*100;
          }
          else                         // Baryons
          {
            dc=ac%1000;
            dc=ac-dc*1000;
          }
          if(dc>2) G4cout<<"Test29: Resonance PDG="<<c<<", m="<<m<<G4endl;
        }
#endif
        mom = sec->GetMomentumDirection();
        e   = sec->GetKineticEnergy();
        if (e < -0.0)
        {
          G4cerr<<"**Test29:Event#"<<iter<<",Hadron#"<<i<<", E="<<e<<" <0 (Set 0)"<<G4endl;
          e = 0.0;
        }
        // for exclusive reaction 2 particles in final state
        p = std::sqrt(e*(e + m + m));
        mom *= p;
        lorV = G4LorentzVector(mom, e + m);    // "e" is a Kinetic energy!
        totSum -= lorV;
#ifdef debug
        G4cout<<"Test29: PDG="<<c<<", 4M="<<lorV<<", resid4M="<<totSum<<", m="<<m<<G4endl;
#endif
        if(std::fabs(m-lorV.m())>.005)
        {
          G4cerr<<"***Test29: m="<<lorV.m()<<" # "<<m<<", d="<<lorV.m()-m<<G4endl;
          alarm=true;
        }
        if(!(lorV.e()>=0||lorV.e()<0)   || !(lorV.px()>=0||lorV.px()<0) ||
           !(lorV.py()>=0||lorV.py()<0) || !(lorV.pz()>=0||lorV.pz()<0))
        {
          G4cerr<<"***Test29: NAN in LorentzVector="<<lorV<<G4endl;
          alarm=true;
        }
        if(c==90000002||c==90002000||c==92000000)
        {
          G4cout<<"***Test29:***Dibaryon *** i="<<i<<", PDG="<<c<<G4endl;
          alarm=true;
        }
        if(c==90000003||c==90003000||c==93000000)
        {
          G4cout<<"***Test29:***Tribaryon *** i="<<i<<", PDG="<<c<<G4endl;
          alarm=true;
        }
        //if((c==90000001||c==2112)&&.5*(e+p)>800.) alarm=true; // High energy neutron
        //if((c==90001000||c==2212)&&.5*(e+p)>800.) alarm=true; // High energy proton
        if(c==223) nOmega++;
        if(c==22) nPhotons++;
        if(c==311||c==321||c==-311||c==-321) nKaons++; // kaons
        if(c==221) nEta++;                             // etas
        //if(c==90002002) nAlphas++;                     // Alphas
        if(c==2212) nProtons++;                        // Protons
        if(c==2112) nNeutrons++;                       // Neutrons
        if(c==2112 && std::fabs(e-1005.)<3.) nSpNeut++;     // Dibar-Neutrons
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
        G4cout<<"Test29:#"<<i<<",PDG="<<c<<",C="<<cCG<<",B="<<cBN<<",4M="<<lorV<<m<<",T="
              <<lorV.e()-m<<G4endl;
#endif
        //delete aChange->GetSecondary(i);
      } // End of the LOOP over secondaries
      // delete secondaries in the end of the event         
      if(std::abs(pPDG)<99&&totCharge==totC&&totBaryN==totBN)//In lepton decay targetMass=0
      {
        totCharge=0;
        totBaryN=0;
        totSum-=G4LorentzVector(0., 0., 0., curM);
      }
      G4double ss=std::fabs(totSum.t())+std::fabs(totSum.x())+std::fabs(totSum.y())+
                  std::fabs(totSum.z());    
#ifdef inter
      G4double sr=std::fabs(Residual.t())+std::fabs(Residual.x())+std::fabs(Residual.y())+
                  std::fabs(Residual.z());    
#endif
#ifdef pdebug
      G4cout<<"TEST29:r4M="<<totSum<<ss<<",rChg="<<totCharge<<",rBaryN="<<totBaryN<<G4endl;
#endif
      // Only for CHIPS
      if (totBaryN || ss>.5 || alarm || (nGamma && !EGamma)) // muCap: chargeNoncons onHlev
      //if (totCharge ||totBaryN || ss>.5 || alarm || (nGamma && !EGamma))
      // for others no conservation checks
      //if (1>2)
      //
      {
        G4cout<<"*Warning*Test29:#"<<iter<<":n="<<nSec<<",4M="<<totSum<<",Chrg="<<totCharge
              <<",BaryN="<<totBaryN<<", R="<<Residual<<",D2="<<ss<<",nN="<<curN<<G4endl;
        totSum = G4LorentzVector(0., 0., pmax, e0+curM);
        if(nGamma&&!EGamma)G4cerr<<"***Test29: Egamma=0"<<G4endl;
        for (G4int indx=0; indx<nSec; indx++)
        {
          sec  = aChange->GetSecondary(indx)->GetDynamicParticle();
          pd   = sec->GetDefinition();
          c    = pd->GetPDGEncoding();
          if(!c)
          {
            G4int chrg=static_cast<G4int>(pd->GetPDGCharge());
            G4int bary=static_cast<G4int>(pd->GetBaryonNumber());
            c=90000000+chrg*999+bary;
          }
          m    = pd->GetPDGMass();
          mom  = sec->GetMomentumDirection();
          G4QPDGCode cQPDG(c);
          sm   = cQPDG.GetMass();
          e    = sec->GetKineticEnergy();
          if (e < -0.0)
          {
            G4cerr<<"**Test29:Ev#"<<iter<<",Hadr#"<<indx<<": E="<<e<<" <0 (Set 0)"<<G4endl;
            e = 0.0;
          }
          p    = std::sqrt(e*(e + m + m));
          mom *= p;
          lorV = G4LorentzVector(mom, e + m);    // "e" is a Kinetic energy!
          totSum -= lorV;
          G4cout<<"Test29:#"<<indx<<",PDG="<<c<<",m="<<m<<"("<<sm<<"),4M="<<lorV<<",T="<<e
                <<",r4M="<<totSum<<G4endl;
        }
#ifdef inter
        if(sr>.5)
#endif
          G4Exception("***Test29", "Test29-02",FatalException, "ALARM or baryn/charge/energy/momentum is not conserved");
      }
#ifndef nout
      if(npart==1) ntp->FillEvt(aChange); // Fill the simulated event in the ASCII "ntuple"
#endif
      for(G4int ides=0; ides<nSec; ides++) delete aChange->GetSecondary(ides);
      aChange->Clear();
#ifdef debug
      G4cout<<"Test29:--->>> After ntp.FillEvt"<<G4endl;
#endif
    } // End of the LOOP over events

    // Stop the timer to estimate the speed of the generation
    timer->Stop();
    G4cout<<"Test29:CalculationTimePerEvent= "<<timer->GetUserElapsed()/nEvt<<" s"<<G4endl;
    //delete timer;
    //#ifdef pverb
    G4cerr<<"Test29: ########## End of run ##########"<<G4endl;
    //#endif
    //delete aPoint;
    //delete bPoint;
    //delete step;  // The G4Step delets aPoint and bPoint
   }
   //delete gTrack; // The G4Track delets the G4DynamicParticle
  } // End of the projectile LOOP
  //delete proc;
#ifndef nout
  G4cout << "###### NTUPLE is written: Before delete ntp #####" << G4endl;
  //delete ntp; // Delete the class to fill the#of events
#endif
  //delete nav;
  //delete runManager;
  G4cout << "###### End of Test29 #####" << G4endl;
  //exit(1); // Never do this !
  return EXIT_SUCCESS;
}
