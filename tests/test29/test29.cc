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

#define nout
//#define pscan
//#define smear
//#define escan
//#define debug
//#define pdebug
// ------------------------------------- FLAGS ------------------
#include "G4UIterminal.hh"
#include "globals.hh"
#include "G4ios.hh"
#include <iostream>
#include <fstream>
#include <iomanip>

#include "G4Material.hh"
#include "G4ElementVector.hh"
#include "G4VContinuousDiscreteProcess.hh"
#include "G4ProcessManager.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleChange.hh"
#include "G4QCaptureAtRest.hh"
#include "G4QString.hh"

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

//#define debug

//int main(int argc, char** argv)
int main()
{
  const G4int nAZ=270;  // Dimension of the table
  const G4int mAZ=266;  // Mafimum filled A (at present). Must be mAZ<nAZ
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
  G4cout<<"Prepare G4QHBook files or ntuples"<<G4endl;
#endif
#ifndef nout
  G4QHBook ntp;
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
  G4int    tgZ;
  G4int    tgN;
  G4int    nEvt;
  //G4int    nofdecays;
  //G4int    decmask=0;
  inFile>>temperature>>ssin2g>>eteps>>nop>>momb>>enb>>pPDG>>tgZ>>tgN>>nEvt>>fN>>fD>>cP>>rM>>sA;
  G4cout<<"CHIPStest Par's: Temp="<<temperature<<",SSs="<<ssin2g<<",P/S="<<eteps<<",nop="
        <<nop<<",p="<<momb<<",e="<<enb<<",pPDG="<<pPDG<<",tZ="<<tgZ<<",tN="<<tgN<<",nEv="
        <<nEvt<<",fN="<<fN<<",fD="<<fD<<",cP="<<cP<<",rM="<<rM<<",sA="<<sA<<G4endl;
  //-------- Initialize CHIPS
  G4QCHIPSWorld* theW=G4QCHIPSWorld::Get();
  theW->GetParticles(nop);           // Create CHIPS World of nop particles
  //G4Exception("***CHIPStest: TMP");
  G4QNucleus::SetParameters(fN,fD,cP,rM);
  G4Quasmon::SetParameters(temperature,ssin2g,eteps);
  G4QEnvironment::SetParameters(sA);                 // SolAngle (pbar-A secondary capture)
  // *********** Now momb is a momentum of the incident particle *************
  G4double mp=G4QPDGCode(pPDG).GetMass();             // @@ just for the check
  momb=momb;
  G4double ep=sqrt(mp*mp+momb*momb);                  // @@ just for the check
  if(enb>0.) ep=enb;
  G4int tPDG=90000000+tgZ*1000+tgN;
  G4double mt=G4QPDGCode(tPDG).GetMass();             // @@ just for the check
  G4QContent pQC=G4QPDGCode(pPDG).GetQuarkContent();
  G4int    cp=pQC.GetCharge();
  G4int    bnp=pQC.GetBaryonNumber();
  G4QContent tQC=G4QPDGCode(tPDG).GetQuarkContent();
  G4int    ct=tQC.GetCharge();
  G4int    bnt=tQC.GetBaryonNumber();
#ifdef debug
  G4cout<<"CHIPStest: pQC"<<pQC<<", pch="<<cp<<", tQC"<<tQC<<", tch="<<ct<<G4endl;
#endif
  G4int    totC=cp+ct;
  G4int    totBN=bnp+bnt;
  G4LorentzVector proj4Mom(0.,0.,momb,ep);
  // @@ G4double fEvt=nEvt;
  // @@ G4double sumE=0.;
  // @@ G4double sumK=0.;
  // @@ G4double sumG=0.;
  // @@ G4double sumT=0.;
  // @@ G4double sumN=0.;
  // @@ G4double sum0=0.;
  // @@ G4double sumP=0.;
  // @@ G4double sum1N=0.;
  // @@ G4double sumNN=0.;
  // @@ G4double sumPP=0.;
  // @@ G4double sumAL=0.;
  // @@ G4double time=clock()/CLOCKS_PER_SEC;
  // ---------------------- PSEUDO-TRACK DEFINISIONS -------------------------
  G4double nx = 0.;
  G4double ny = 0.;
  G4double nz = 0.;
  G4ThreeVector aPosition(nx*mm, ny*mm, nz*mm);
  G4double      aTime      = 0. ;
  G4ThreeVector aDirection = G4ThreeVector(0.,0.,1.);

  // ====================== CREATE G4 ENVIRONMENT =========================
  // ----------- Material Factory: Definition of Elements & Materials ----------

  //G4Element*   H  = new G4Element ("Hydrogen", "H",   1. ,  1.01*g/mole);
  //G4Element*   N  = new G4Element ("Nitrigen", "N",   7. , 14.00*g/mole);
  //G4Element*   O  = new G4Element ("Oxygen"  , "O",   8. , 16.00*g/mole);
  //G4Element*   C  = new G4Element ("Carbon"  , "C",   6. , 12.00*g/mole);
  //G4Element*  Cs  = new G4Element ("Cesium"  , "Cs", 55. , 132.905*g/mole);
  //G4Element*   I  = new G4Element ("Iodide"  , "I",  53. , 126.9044*g/mole);

  //G4Material* cmaterial = 0;
  //cmaterial  = new G4Material("H",       1.,   1.0000*g/mole,  1.000*g/cm3);
  //cmaterial  = new G4Material("D",       1.,   2.0000*g/mole,  1.000*g/cm3);
  //cmaterial  = new G4Material("Li",      3.,   6.9410*g/mole,  1.000*g/cm3);
  //cmaterial  = new G4Material("Be",      4.,   9.0100*g/mole,  1.848*g/cm3);
  //cmaterial  = new G4Material("C",       6.,  12.0000*g/mole,  2.000*g/cm3);
  //cmaterial  = new G4Material("Al",     13.,  26.9800*g/mole,  2.700*g/cm3);
  //cmaterial  = new G4Material("Si",     14.,  28.0550*g/mole,  2.330*g/cm3);
  //cmaterial  = new G4Material("LAr",    18.,  39.9500*g/mole,  1.393*g/cm3);
  //cmaterial  = new G4Material("Zr",     40.,  91.2240*g/mole,  4.000*g/cm3);
  //cmaterial  = new G4Material("LXe",    54., 131.2900*g/mole,  3.020*g/cm3);
  //cmaterial  = new G4Material("Fe",     26.,  55.8500*g/mole,  7.870*g/cm3);
  //cmaterial  = new G4Material("Ni",     28.,  58.6934*g/mole,  8.000*g/cm3);
  //cmaterial  = new G4Material("Cu",     29.,  63.5500*g/mole,  8.960*g/cm3);
  //cmaterial  = new G4Material("Au",     79., 196.9700*g/mole, 19.320*g/cm3);
  //cmaterial  = new G4Material("Ta",     73., 180.9479*g/mole, 16.670*g/cm3);
  //cmaterial  = new G4Material("W",      74., 183.8500*g/mole, 19.300*g/cm3);
  //cmaterial  = new G4Material("Pb",     82., 207.1900*g/mole, 11.350*g/cm3);
  //cmaterial  = new G4Material("Bi",     83., 208.9800*g/mole, 12.000*g/cm3);
  //cmaterial  = new G4Material("U",      92., 238.0300*g/mole, 18.950*g/cm3);

  //cmaterial = new G4Material("O2", 8., 16.00*g/mole, 1.1*g/cm3);
  //cmaterial->SetChemicalFormula("O_2");

  //cmaterial = new G4Material ("Water" , 1.*g/cm3, 2);
  //cmaterial->AddElement(H,2);
  //cmaterial->AddElement(O,1);
  //cmaterial->SetChemicalFormula("H_2O");

  //cmaterial = new G4Material ("Ethane" , 0.4241*g/cm3, 2);
  //cmaterial->AddElement(H,6);
  //cmaterial->AddElement(C,2);
  //cmaterial->SetChemicalFormula("C_2H_6");

  //cmaterial  = new G4Material("Graphite",6.,  12.0000*g/mole,  2.265*g/cm3);
  //cmaterial->SetChemicalFormula("Graphite");

  //cmaterial = new G4Material ("CsI" , 4.53*g/cm3, 2);
  //cmaterial->AddElement(Cs,1);
  //cmaterial->AddElement(I,1);
  //cmaterial->SetChemicalFormula("CsI");

  //cmaterial = new G4Material("Air"  , 1.290*mg/cm3, 2);
  //cmaterial->AddElement(N, 0.7);
  //cmaterial->AddElement(O, 0.3);

  // ----------- Hadronic Physics definition ---------------------
  //Test29Physics*   phys = new Test29Physics(); // ?
  // ---------- Define material for the simulation ------------------
  //G4double tgA     = 26.98; // @@ Important? Can it be just tgZ+tgN?
  G4double tgA     = tgZ+tgN; // @@ Temporary, not good
  G4double tgR     = 2.7;   // @@ Not important for the thin target example. Can be any
  G4String nameMat = "Tg";  // @@ Not important can be an arbitrary name
  // The material can be copied from the commented cMaterial Factory above
  G4Material* material = new G4Material(nameMat, tgZ*1., tgA*g/mole, tgR*g/cm3);
  // For the Material Factory case
  //G4String  nameMat  = "Si";
  //G4Material* material = G4Material::GetMaterial(nameMat); // Get definition of Material
  if(!material)
  {
    G4cout<<"Test29:Material "<<nameMat<<" is not defined in the Test29Material"<<G4endl;
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

  G4cout<<"Test29:###### Start new run #####" << G4endl; // only one run

  // Particles Definitions Factory for incident and secondary particles
  //G4ParticleDefinition* proton = G4Proton::Proton();
  //G4ParticleDefinition* neutron = G4Neutron::Neutron();
  //G4ParticleDefinition* lambda = G4Lambda::Lambda();
  //G4ParticleDefinition* sigmam = G4SigmaMinus::SigmaMinus();
  //G4ParticleDefinition* sigmaz = G4SigmaZero::SigmaZero();
  //G4ParticleDefinition* sigmap = G4SigmaPlus::SigmaPlus();
  //G4ParticleDefinition* gamma  = G4Gamma::Gamma();
  //G4ParticleDefinition* elect  = G4Electron::Electron();
  //G4ParticleDefinition* negmu  = G4MuonMinus::MuonMinus();
  //G4ParticleDefinition* pin    = G4PionMinus::PionMinus();
  //G4ParticleDefinition* pip    = G4PionPlus::PionPlus();
  //G4ParticleDefinition* pi0    = G4PionZero::PionZero();
  //G4ParticleDefinition* deu    = G4Deuteron::DeuteronDefinition();
  //G4ParticleDefinition* tri    = G4Triton::TritonDefinition();
  //G4ParticleDefinition* he3    = G4He3::He3Definition();
  //G4ParticleDefinition* alp    = G4Alpha::AlphaDefinition();

  //G4ParticleDefinition* part=G4MuonMinus::MuonMinus(); // Definition of the projectile
  G4ParticleDefinition* part=G4AntiProton::AntiProton(); // Definition of the projectile
  G4double pMass = part->GetPDGMass();                 // Mass of the projectile
  // @@ G4double totCN  = tgZ+part->GetPDGCharge();
  // @@ G4int    totBNN = tgZ+tgN+part->GetBaryonNumber();
#ifdef pdebug
  G4cout<<"Test29:tC="<<totC<<"?="<<totCN<<",tB="<<totBN<<"?="<<totBNN<<G4endl;
#endif
  G4double theStep   = 0.01*micrometer;

  // G4int maxz = (G4int)((*(material->GetElementVector()))[0]->GetZ()) + 1;
  G4cout<<"Test29: The particle: "<<part->GetParticleName()<<G4endl;
  G4cout<<"Test29: The material: "<<material->GetName()<<G4endl;
  G4cout<<"Test29: The step:     "<<theStep/mm<<" mm"<<G4endl;
  G4cout<<"Test29: The position: "<<aPosition/mm<<" mm"<<G4endl;
  G4cout<<"Test29: The direction:"<<aDirection<<G4endl;
  G4cout<<"Test29: The time:     "<<aTime/ns<<" ns"<<G4endl;

  // Different Process Managers are used for the atRest and onFlight processes
  G4ProcessManager* man = new G4ProcessManager(part);
  //G4VDiscreteProcess* proc = new G4QMuonNuclearOnFlight;
  G4VRestProcess* proc = new G4QCaptureAtRest;
  if(!proc)
  {
    G4cout<<"Tst29:For"<<part->GetParticleName()<<" no G4QMuonNuclearOnFlight"<<G4endl;
	exit(1);
  }
  //man->AddDiscreteProcess(proc);
  man->AddRestProcess(proc);

  // Create a DynamicParticle
  //G4double  energy   = 1000.*MeV;                            // 1 GeV particle energy
  G4double  energy   = 0.*MeV;                               // 0 GeV particle energy
  G4DynamicParticle dParticle(part,aDirection,energy);

  // General Track definition
  G4Track* gTrack = new G4Track(&dParticle,aTime,aPosition); // Track definition (NoTarg)

  // Step Definition
  G4Step* step = new G4Step();
  step->SetTrack(gTrack);          // Step is initialized by the Track (?)

  G4StepPoint *aPoint = new G4StepPoint(); // It cant be initialized right away (!?)
  aPoint->SetPosition(aPosition);
  aPoint->SetMaterial(material);
  G4double safety = 10000.*cm;
  aPoint->SetSafety(safety);
  step->SetPreStepPoint(aPoint);   // Begin of the step

  G4StepPoint* bPoint = aPoint;
  G4ThreeVector bPosition = aDirection*theStep+aPosition; // aDirection is defined byCard
  bPoint->SetPosition(bPosition);
  step->SetPostStepPoint(bPoint);  // End of the step
  step->SetStepLength(theStep);    // Step is set byCard above

  G4RotationMatrix* rot  = new G4RotationMatrix(); //Is it necessary in this application?
  G4double phi0 = aDirection.phi();
  G4double theta0 = aDirection.theta();
  rot->rotateZ(-phi0);
  rot->rotateY(-theta0);

  G4cout<<"Test29: rotation= "<<(*rot)*(aDirection)<<G4endl;

  G4Timer* timer = new G4Timer();
  timer->Start();
  const G4DynamicParticle* sec = 0;
  G4ParticleDefinition* pd;
  G4ThreeVector  mom;
  G4LorentzVector totSum, lorV;
  G4double e, p, m;
  // @@ G4double px, py, pz, pt;
  G4VParticleChange* aChange = 0;
  G4double e0 = energy+pMass;
  G4double pmax=sqrt(e0*e0-pMass*pMass);
  G4double et=e0+mt;
  //G4int nEvt=100;
  // Randomization loop: cycle random generator, using 2 lower digits in nEvt
  G4int    iRandCount = nEvt%100;
  G4double vRandCount = 0.;
  while (iRandCount>0)
  {
    vRandCount = G4UniformRand();
    iRandCount--;
  }
  for (G4int iter=0; iter<nEvt; iter++)
  {
#ifdef debug
    G4cout<<"Test29: ### "<<iter<< "-th event starts.###"<<G4endl;
#endif
    if(!(iter%1000) && iter) G4cout<<"TEST29: "<<iter<<" events are simulated"<<G4endl;
    dParticle.SetKineticEnergy(energy);// Fill the Kinetic Energy of the projectile

    gTrack->SetStep(step);            // Now step is included in the Track (see above)
    gTrack->SetKineticEnergy(energy); // Duplication of Kin. Energy for the Track (?!)

    totSum = G4LorentzVector(0., 0., pmax, et);

    G4double        totCharge = totC;
    G4int           totBaryN = totBN;

    aChange = proc->PostStepDoIt(*gTrack,*step);

    //G4double de = aChange->GetLocalEnergyDeposit(); // Init Total Energy by EnergyDeposit
    G4int nSec = aChange->GetNumberOfSecondaries();

    //G4int nbar = 0;

    // @@ G4int npt=0;
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
    G4cout<<"Test29:----DONE^^^^^^^************^^^^^^^^:ir="<<iter<<": #ofH="<<tNH<<G4endl;
    if(!(iter%100)) G4cerr<<"#"<<iter<<G4endl;
#endif
    G4bool alarm=false;
    // @@ G4bool rad=false;
    // @@ G4bool hyp=false;
    // @@ G4bool badPDG=false;
    // ------- LOOP over secondary particles -------
    for(G4int i=0; i<nSec; i++)
    {
      sec = aChange->GetSecondary(i)->GetDynamicParticle();
      pd  = sec->GetDefinition();
      c   = pd->GetPDGEncoding();
      m   = pd->GetPDGMass();
      mom = sec->GetMomentumDirection();
      e   = sec->GetKineticEnergy();
	  if (e < 0.0)
      {
	    G4cerr<<"**Test29:Event#"<<iter<<",Hadron#"<<i<<", E="<<e<<" <0 (Set 0)"<<G4endl;
        e = 0.0;
      }

	  // for exclusive reaction 2 particles in final state
	  p = sqrt(e*(e + m + m));
	  mom *= p;
      lorV = G4LorentzVector(mom, e + m);    // "e" is a Kinetic energy!
      totSum -= lorV;
      if(abs(m-lorV.m())>.005)
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
      if(c==223) nOmega++;
      if(c==22) nPhotons++;
      if(c==311||c==321||c==-311||c==-321) nKaons++; // kaons
      if(c==221) nEta++;                             // etas
      //if(c==90002002) nAlphas++;                     // Alphas
      if(c==2212) nProtons++;                        // Protons
      if(c==2112) nNeutrons++;                       // Neutrons
      if(c==2112 && abs(e-1005.)<3.) nSpNeut++;// Dibar-Neutrons
      //if(c==90002002 && e-m<7.) nSpAlph++;     // Special Alphas
      if(c==111) nP0++;                              // Neutral  pions
      if(c==-211) nPN++;                             // Negative pions
      if(c==211) nPP++;                              // Positive pions
      if(c==22) nGamma++;                            // Gammas
      if(c==22) EGamma+=e;                           // Energy of gammas
      totCharge-=pd->GetPDGCharge();
      totBaryN-=pd->GetBaryonNumber();
#ifdef pdebug
      G4cout<<"Test29:#"<<i<<"PDG="<<c<<",4M="<<lorV<<m<<",T="<<lorV.e()-m<<G4endl;
#endif
      delete aChange->GetSecondary(i);
	} // End of the LOOP over secondaries
	//	delete secondaries in the end of the event       	 
#ifdef pdebug
    G4cout<<"TEST29: 4M="<<totSum<<", Charge="<<totCharge<<", BaryN="<<totBaryN<<G4endl;
#endif
    G4double ss=abs(totSum.t())+abs(totSum.x())+abs(totSum.y())+abs(totSum.z());    
	if (totCharge ||totBaryN || !(ss<.01) || alarm || nGamma&&!EGamma)
    {
      G4cerr<<"***Test29:#"<<iter<<":n="<<nSec<<",4M="<<totSum<<",Charge="<<totCharge
            <<",BaryN="<<totBaryN<<G4endl;
      if(nGamma&&!EGamma)G4cerr<<"***Test29: Egamma=0"<<G4endl;
      totSum = G4LorentzVector(0., 0., pmax, et);
      for (G4int indx=0; indx<nSec; indx++)
      {
        sec = aChange->GetSecondary(indx)->GetDynamicParticle();
        pd  = sec->GetDefinition();
        c   = pd->GetPDGEncoding();
        m   = pd->GetPDGMass();
        mom = sec->GetMomentumDirection();
        e   = sec->GetKineticEnergy();
	    p = sqrt(e*(e + m + m));
	    mom *= p;
        lorV = G4LorentzVector(mom, e + m);    // "e" is a Kinetic energy!
        totSum -= lorV;
        G4cerr<<"Test29:#"<<indx<<",PDG="<<c<<",m="<<m<<",4M="<<lorV<<",T="<<e
              <<", d4M="<<totSum<<G4endl;
      }
      G4Exception("***Test29: ALARM or baryn/charge/energy/momentum is not conserved");
    }
    // =============== May be print it here if it is not zero...
    aChange->Clear();
  } // End of the LOOP over events

  // Stop the timer to estimate the speed of the generation
  timer->Stop();
  G4cout<<"Test29: Calculation time="<<*timer<<G4endl;
  delete timer;

  G4cerr<<"Test29: ########## End of run ##########"<<G4endl;


  delete material;
  //delete phys;
  //delete cmaterial;                   // Temporary material definition pointer (?)

  G4cout << "###### End of Test29 #####" << G4endl;
  exit(1);
}
