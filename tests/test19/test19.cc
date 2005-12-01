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

//#define nout
//#define pscan
//#define smear
//#define escan
//#define idebug
//#define tdebug
//#define debug
//#define pdebug
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
#include "G4QCollision.hh"
//#include "G4QuasmonString.hh"

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

#include "Test19Physics.hh"
#include "Test19PhysicsList.hh"

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

// Temporary AIDA can be not used in test19
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
  const G4int nTg=5;   // Length of the target list for the Performance test
  G4int tli[nTg]={90001000,90002002,90007007,90027032,90092146}; // PDG Codes of targets
  G4String tnm[nTg]={"Hydrogen","Helium","Nitrogen","Cobalt","Uranium"}; // Target names
  G4String tsy[nTg]={"1H","He","14N","59Co","238U"}; // Target symbols for the Target Loop
  G4Material* mat[nTg]={0,0,0,0,0}; // Material pointers for the Target Loop
  const G4int nPr=4;  // Length of the projectile list for the Performance test
  G4int pli[nPr] = {22, 11, 13, 15}; // PDG Codes of the projectile particles
  const G4int nEn=3;  // Length of the kin. energy list for the Performance test
  G4double eli[nEn] = {27., 227., 999.}; //
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
  if(mAZ>=nAZ) G4cout<<"***Test19: Too big mAZ="<<mAZ<<", nAZ="<<nAZ<<G4endl;

  // Run manager
  G4RunManager* runManager = new G4RunManager;
  runManager->SetUserInitialization(new Test19PhysicsList);
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
  G4QCollision::SetParameters(temperature,ssin2g,eteps,fN,fD,cP,rM,nop,sA);
  G4QCollision::SetManual();
  // ********** Now momb is a momentum of the incident particle, if =0 => LOOP ************
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
  // ---------- Define material for the simulation ------------------
  G4int tgA        = tgZ+tgN; // Mass number - fake
  G4double tgR     = 2.7;   // @@ Not important for the thin target example. Can be any
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
      isotope = new G4Isotope(tsy[tgi], tgZ, tgA, tgR*g/mole);
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
    isotope = new G4Isotope("Isotop", tgZ, tgA, tgR*g/mole);
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
  // ---------------------- PSEUDO-TRACK DEFINISIONS -------------------------
  G4double nx = 0.;
  G4double ny = 0.;
  G4double nz = 0.;
  G4int npart=1;                                      // By default only one particle
  if(!pPDG) npart=nPr;                                // Make a LOOP ove all particles
  // Different Process Managers are used for the atRest and onFlight processes
		//G4ProcessManager* man = new G4ProcessManager(part); //Does not help to go out
  //G4VDiscreteProcess* proc = new G4QCollision;
  G4QCollision* proc = new G4QCollision;
  if(!proc)
  {
    G4cout<<"Tst19: there is no G4QCollision process"<<G4endl;
	   exit(1);
  }
#ifdef debug
  G4cout<<"Test19:--***-- process is created --***--" << G4endl; // only one run
#endif
  proc->SetParameters(temperature, ssin2g, eteps, fN, fD, cP, rM, nop, sA);
  //man->AddDiscreteProcess(proc); //Does not help to go out
  G4int nTot=npart*tgm*cnE;
  G4int nCur=0;
  for(G4int pnb=0; pnb<npart; pnb++) // LOOP over particles
  {
   if (npart>1) pPDG=pli[pnb];
   G4QContent pQC=G4QPDGCode(pPDG).GetQuarkContent();
   G4int        cp = pQC.GetCharge();
   if(pPDG==22) cp = 0;                        // QuarkContent of photon is not supported
   if(pPDG==11 || pPDG==13 || pPDG==15) cp=-1; // QuarkContent of leptons is not supported
   if(pPDG==-11|| pPDG==-13|| pPDG==-15) cp=1; // QuarkContent of antilept isn't supported
   G4ParticleDefinition* part=0;               // DefaultProj is particle is empty
   if(pPDG==2212) part=G4Proton::Proton();           // Definition of Proton projectile
   else if(pPDG==22) part=G4Gamma::Gamma();          // Definition of gamma projectile
   else if(pPDG==11) part=G4Electron::Electron();    // Definition of electron projectile
   else if(pPDG==-11) part=G4Positron::Positron();   // Definition of positron projectile
   else if(pPDG==13) part=G4MuonMinus::MuonMinus();  // Definition of mu- projectile
   else if(pPDG==-13) part=G4MuonPlus::MuonPlus();   // Definition of the mu+ projectile
   else if(pPDG==15) part=G4TauMinus::TauMinus();    // Definition of tau- projectile
   else if(pPDG==-15) part=G4TauPlus::TauPlus();     // Definition of tau+ projectile
   else if(pPDG==2112) part=G4Neutron::Neutron();    // Definition of Neutron projectile
   else if(pPDG==-211) part=G4PionMinus::PionMinus();// Definition of Pi- projectile
   else if(pPDG==211)  part=G4PionMinus::PionMinus();// Definition of Pi+ projectile
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
     G4cerr<<"***Test19: "<<pPDG<<" is a PDG code of not supported particle"<<G4endl;
     G4Exception("***Test19: OnFlight Process is called for not supported particle");
   }
   G4double pMass = part->GetPDGMass();                 // Mass of the projectile
   //
   G4ThreeVector aPosition(nx*mm, ny*mm, nz*mm);
   G4DynamicParticle* dParticle = new G4DynamicParticle(part,aDirection,0.);// Dummy Energy

   // General Track definition
   G4Track* gTrack = new G4Track(dParticle,aTime,aPosition); // Track definition (NoTarg)

   G4int    bnp=pQC.GetBaryonNumber();
   G4LorentzVector proj4Mom(0.,0.,momb,ep);
   // ---------- Define material for the simulation ------------------
   //G4double tgA     = 26.98; // @@ Important? Can it be just tgZ+tgN?
   G4int tgA        = tgZ+tgN; // @@ Temporary, not good
   G4double tgR     = 1.;   // @@ Not important for the thin target example. Can be any
   G4String nameMat = "Thin Target";  // @@ Not important can be an arbitrary name
   // The material can be copied from the commented cMaterial Factory above
   G4Isotope* isotope = new G4Isotope("Isotop", tgZ, tgA, tgR*g/mole);
   G4Element* element = new G4Element("ZA_Isotop", "ZA", 1);
   element->AddIsotope(isotope, 100.*perCent);
   G4Material* material = new G4Material("ZA_Isomer", 1.*g/cm3, 1);
   material->AddElement(element, 1);
   //G4Material* material = new G4Material(nameMat, tgZ*1., tgA*g/mole, tgR*g/cm3);
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
   G4int    totBNN = tgZ+tgN+part->GetBaryonNumber();
#endif
   G4double theStep   = 0.01*micrometer;
   // G4int maxz = (G4int)((*(material->GetElementVector()))[0]->GetZ()) + 1;
#ifdef pdebug
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
   for(G4int nen=0; nen<cnE; nen++)                         // LOOP over projectile energy
		 {
    G4double  energy = (ep-mp)*MeV;                          // kinetic particle energy
    if(cnE>1) energy = eli[nen]*MeV;
#ifdef pdebug
    G4cout<<"Test19: M="<<mp<<", T="<<energy<<", MeV"<<G4endl;
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
      material=mat[tgi];
      G4Element* curEl=(*(material->GetElementVector()))[0];
      G4cout<<"Test19: Material="<<material->GetName()<<", Element[0]="<<curEl->GetName()
												<<",A[0]="<<(*(curEl->GetIsotopeVector()))[0]->GetN()<<" is selected."<<G4endl;
     }
  			G4cout<<"Test19:NewRun:Targ="<<tPDG<<",Proj="<<pPDG<<", "<<nCur<<" of "<<nTot<<G4endl;
     G4double mt=G4QPDGCode(tPDG).GetMass();             // @@ just for check
     G4QContent tQC=G4QPDGCode(tPDG).GetQuarkContent();
     G4int    ct=tQC.GetCharge();
     G4int    bnt=tQC.GetBaryonNumber();
#ifdef debug
     G4cout<<"Test19: pQC"<<pQC<<", pch="<<cp<<", tQC"<<tQC<<", tch="<<ct<<G4endl;
#endif

     G4int    totC=cp+ct;
     G4int    totBN=bnp+bnt;
#ifdef debug
     G4cout<<"Test19:tC="<<totC<<"?="<<totCN<<",tB="<<totBN<<"?="<<totBNN<<G4endl;
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
     G4ThreeVector bPosition = aDirection*theStep+aPosition; // aDirection is defined byCard
     bPoint->SetPosition(bPosition);
     step->SetPostStepPoint(bPoint);  // End of the step
     step->SetStepLength(theStep);    // Step is set byCard above
#ifdef pverb
     G4cout<<"Test19: The end point is defined and filled in the step "<<G4endl;
#endif
     G4Timer* timer = new G4Timer();
     timer->Start();
#ifdef idebug
     G4cout<<"Test19: Run is started, timer is started, kinEnergy="<<energy<<G4endl;
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
     G4double et=e0+mt;
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
     G4cout<<"Test19: Before the event loop, nEvents= "<<nEvt<<G4endl;
#endif
     G4double dTot=0.;
     G4double dEl=0.;
     for (G4int iter=0; iter<nEvt; iter++)
     {
#ifdef debug
      G4cout<<"Test19: ### "<<iter<< "-th event starts.### energy="<<energy<<G4endl;
#endif

      if(!(iter%100)&&iter)G4cout<<"***=>TEST19: "<<iter<<" events are simulated"<<G4endl;

      gTrack->SetStep(step);            // Now step is included in the Track (see above)
      gTrack->SetKineticEnergy(energy); // Duplication of Kin. Energy for the Track (?!)
#ifdef debug
      G4cout<<"Test19: Before PostStepDoIt"<<G4endl;
#endif
      aChange = proc->PostStepDoIt(*gTrack,*step); // For On Flight
#ifdef debug
      G4cout<<"Test19: After PostStepDoIt"<<G4endl;
#endif
      G4int nSec = aChange->GetNumberOfSecondaries();
      //G4cout<<"Test19: "<<nSec<<" secondary particles are generated"<<G4endl;
      G4double totCharge = totC;
      G4int    curN=proc->GetNumberOfNeutronsInTarget();
      G4int    dBN = curN-tgN;
      G4int    totBaryN = totBN+dBN;
      G4int    curPDG=tPDG+dBN;
      G4double curM=G4QPDGCode(curPDG).GetMass(); // Update mass of the TargetNucleus
      totSum = G4LorentzVector(0., 0., pmax, et+curM-mt);
      G4LorentzVector Residual=proc->GetEnegryMomentumConservation();
#ifdef debug
      G4double de = aChange->GetLocalEnergyDeposit();// Init TotalEnergy by EnergyDeposit
      G4cout<<"Test19: "<<nSec<<" secondary particles are generated, dE="<<de<<G4endl;
#endif
      // @@ ----------------------- Begin
      G4double weight = aChange->GetSecondary(0)->GetDynamicParticle()->GetKineticEnergy();
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
      G4cout<<"Test19:----DONE^^^^^^^*******^^^^^^^^:ir="<<iter<<": #ofH="<<nSec<<G4endl;
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
        if(!c)
        {
          G4int chrg=static_cast<G4int>(pd->GetPDGCharge());
          G4int bary=static_cast<G4int>(pd->GetBaryonNumber());
          c=90000000+chrg*999+bary;
        }
        m   = pd->GetPDGMass();
        mom = sec->GetMomentumDirection();
        e   = sec->GetKineticEnergy();
	       if (e < 0.0)
        {
	         G4cerr<<"**Test19:Event#"<<iter<<",Hadron#"<<i<<", E="<<e<<" <0 (Set 0)"<<G4endl;
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
		        G4cerr<<"***Test19: m="<<lorV.m()<<" # "<<m<<", d="<<lorV.m()-m<<G4endl;
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
        G4cout<<"Test19:#"<<i<<",PDG="<<c<<",C="<<cCG<<",B="<<cBN<<",4M="<<lorV<<m<<",T="
              <<lorV.e()-m<<G4endl;
#endif
        //delete aChange->GetSecondary(i);
	     } // End of the LOOP over secondaries
	     //	delete secondaries in the end of the event       	 
      G4double ss=std::fabs(totSum.t())+std::fabs(totSum.x())+std::fabs(totSum.y())+
                  std::fabs(totSum.z());
      G4double sr=std::fabs(Residual.t())+std::fabs(Residual.x())+std::fabs(Residual.y())+
                  std::fabs(Residual.z());    
#ifdef pdebug
      G4cout<<">TEST19:r4M="<<totSum<<ss<<",rCh="<<totCharge<<",rBaryN="<<totBaryN<<G4endl;
#endif
	     //if (1>2) // @@ The check is temporary closed
						if (totCharge ||totBaryN || ss>.27 || alarm || nGamma&&!EGamma)
      {
        totSum = G4LorentzVector(0., 0., pmax, et);
        G4cerr<<"**Test19:#"<<iter<<":n="<<nSec<<",4M="<<totSum<<",Charge="<<totCharge
              <<",BaryN="<<totBaryN<<", R="<<Residual<<",D2="<<ss<<",nN="<<curN<<G4endl;
        if(nGamma&&!EGamma)G4cerr<<"***Test19: Egamma=0"<<G4endl;
        for (G4int indx=0; indx<nSec; indx++)
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
	         p = std::sqrt(e*(e + m + m));
	         mom *= p;
          lorV = G4LorentzVector(mom, e + m);    // "e" is a Kinetic energy!
          totSum -= lorV;
          G4cerr<<"Test19:#"<<indx<<",PDG="<<c<<",m="<<m<<",4M="<<lorV<<",T="<<e
                <<", d4M="<<totSum<<G4endl;
        }
        if(sr>.27)
          G4Exception("***Test19: ALARM or baryn/charge/energy/momentum is not conserved");
      }
#ifndef nout
	     ntp->FillEvt(aChange); // Fill the simulated event in the ASCII "ntuple"
#endif
      // =============== May be print it here if it is not zero...
      for(G4int ides=0; ides<nSec; ides++) delete aChange->GetSecondary(ides);
      aChange->Clear();
#ifdef debug
      G4cout<<"Test19:--->>> After ntp.FillEvt"<<G4endl;
#endif
     } // End of the LOOP over events
     // Stop the timer to estimate the speed of the generation
     timer->Stop();
     G4cout<<"Test19:CalculationTimePerEvent= "<<timer->GetUserElapsed()/nEvt<<" s"<<G4endl;
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
     G4cerr<<"Test19: ########## End of run ##########"<<G4endl;
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
  delete runManager;
#ifdef pverb
  G4cout << "###### End of Test19 #####" << G4endl;
#endif
  //exit(1); // Never do this !
  //return no_of_errors;
  //return 0;
  //abort();
  //return EXIT_SUCCESS;
}
