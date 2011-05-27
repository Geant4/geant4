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
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//
//      File name:     G4EvaporationProbabilityTest.cc 
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it), 
// 
//      Creation date: 27 October 1998
//
//      Modifications: 
//
// -------------------------------------------------------------------

#include "globals.hh"

#include "G4ios.hh"
#include <fstream>
#include <iomanip>
#include <iostream>

#include "CLHEP/Hist/TupleManager.h"
#include "CLHEP/Hist/HBookFile.h"
#include "CLHEP/Hist/Histogram.h"
#include "CLHEP/Hist/Tuple.h"
#include "CLHEP/String/Strings.h"
#include "G4VEvaporationChannel.hh"
#include "G4CompetitiveFission.hh"

#include "G4VEvaporationChannel.hh"
#include "G4EvaporationChannel.hh"
#include "G4PhotonEvaporation.hh"
#include "G4DataVector.hh"
#include "G4LorentzVector.hh"
#include "G4ThreeVector.hh"
#include "G4NucleiProperties.hh"
#include "G4Fragment.hh"
#include "Randomize.hh"
#include "G4FragmentVector.hh"

#include "g4rw/tvvector.h"

int main()
{
  // MGP ---- HBOOK initialization
  HepTupleManager* hbookManager;
  hbookManager = new HBookFile("prob.hbook", 58);
  assert (hbookManager != 0);

  // MGP ---- Book histograms

  HepHistogram* hGammaE;
  hGammaE = hbookManager->histogram("Gamma energy", 100,0.,10.);
  assert (hGammaE != 0);  

  HepHistogram* hNProducts;
  hNProducts = hbookManager->histogram("Number of products", 20,0.,20.);
  assert (hNProducts != 0);  

  HepHistogram* hProb;
  hProb = hbookManager->histogram("Probability * 1.e25", 100,0.,10.);
  assert (hProb != 0);  

  // MGP ---- Book a ntuple
  HepTuple* ntuple;
  ntuple = hbookManager->ntuple("G4EvaporationProbability");
  assert (ntuple != 0);

  G4int Z;
  G4int A;

  G4cout << "Enter Z and A" << G4endl;
  G4cin >> Z >> A;

  assert (Z > 0);
  assert (A > 0);
  assert (A > Z);
  
  G4int iter;
  G4cout << "Enter number of iterations " << G4endl;
  G4cin >> iter;
  if (iter <1) iter = 1;

  G4int nProbs;
  G4cout << "Enter max number of channels (0 - 32) " << G4endl;
  G4cin >> nProbs;
  if (nProbs < 0 || nProbs > 32) nProbs = 32;

  G4double excMin;
  G4double excMax;    
  G4cout << "Enter initial min an max excitation energy" << G4endl;
  G4cin >> excMin >> excMax;
  assert (excMin >= 0.);
  assert (excMax > 0.);
  assert (excMax >= excMin);


  // G4EvaporationProbability for this (Z,A) material

  // Excitation energy levels for each channel

  G4int nExcitedStates = 10;

  G4double zero = 0.;

  G4RWTValVector<G4double> ExcitEnergyChann00(nExcitedStates,zero);
  G4RWTValVector<G4double> ExcitEnergyChann01(nExcitedStates,zero);
  G4RWTValVector<G4double> ExcitEnergyChann02(nExcitedStates,zero);
  G4RWTValVector<G4double> ExcitEnergyChann03(nExcitedStates,zero);
  G4RWTValVector<G4double> ExcitEnergyChann04(nExcitedStates,zero);
  G4RWTValVector<G4double> ExcitEnergyChann05(nExcitedStates,zero);
  G4RWTValVector<G4double> ExcitEnergyChann06(nExcitedStates,zero);
  G4RWTValVector<G4double> ExcitEnergyChann07(nExcitedStates,zero);
  G4RWTValVector<G4double> ExcitEnergyChann08(nExcitedStates,zero);
  G4RWTValVector<G4double> ExcitEnergyChann09(nExcitedStates,zero);
  G4RWTValVector<G4double> ExcitEnergyChann10(nExcitedStates,zero);
  G4RWTValVector<G4double> ExcitEnergyChann11(nExcitedStates,zero);
  G4RWTValVector<G4double> ExcitEnergyChann12(nExcitedStates,zero);
  G4RWTValVector<G4double> ExcitEnergyChann13(nExcitedStates,zero);
  G4RWTValVector<G4double> ExcitEnergyChann14(nExcitedStates,zero);
  G4RWTValVector<G4double> ExcitEnergyChann15(nExcitedStates,zero);
  G4RWTValVector<G4double> ExcitEnergyChann16(nExcitedStates,zero);
  G4RWTValVector<G4double> ExcitEnergyChann17(nExcitedStates,zero);
  G4RWTValVector<G4double> ExcitEnergyChann18(nExcitedStates,zero);
  G4RWTValVector<G4double> ExcitEnergyChann19(nExcitedStates,zero);
  G4RWTValVector<G4double> ExcitEnergyChann20(nExcitedStates,zero);
  G4RWTValVector<G4double> ExcitEnergyChann21(nExcitedStates,zero);
  G4RWTValVector<G4double> ExcitEnergyChann22(nExcitedStates,zero);
  G4RWTValVector<G4double> ExcitEnergyChann23(nExcitedStates,zero);
  G4RWTValVector<G4double> ExcitEnergyChann24(nExcitedStates,zero);
  G4RWTValVector<G4double> ExcitEnergyChann25(nExcitedStates,zero);
  G4RWTValVector<G4double> ExcitEnergyChann26(nExcitedStates,zero);
  G4RWTValVector<G4double> ExcitEnergyChann27(nExcitedStates,zero);
  G4RWTValVector<G4double> ExcitEnergyChann28(nExcitedStates,zero);
  G4RWTValVector<G4double> ExcitEnergyChann29(nExcitedStates,zero);
  G4RWTValVector<G4double> ExcitEnergyChann30(nExcitedStates,zero);
  G4RWTValVector<G4double> ExcitEnergyChann31(nExcitedStates,zero);


  // Spin of excitation energy levels for each channel

  G4int izero = 0;

  G4RWTValVector<G4int> ExcitSpinChann00(nExcitedStates,izero);
  G4RWTValVector<G4int> ExcitSpinChann01(nExcitedStates,izero);
  G4RWTValVector<G4int> ExcitSpinChann02(nExcitedStates,izero);
  G4RWTValVector<G4int> ExcitSpinChann03(nExcitedStates,izero);
  G4RWTValVector<G4int> ExcitSpinChann04(nExcitedStates,izero);
  G4RWTValVector<G4int> ExcitSpinChann05(nExcitedStates,izero);
  G4RWTValVector<G4int> ExcitSpinChann06(nExcitedStates,izero);
  G4RWTValVector<G4int> ExcitSpinChann07(nExcitedStates,izero);
  G4RWTValVector<G4int> ExcitSpinChann08(nExcitedStates,izero);
  G4RWTValVector<G4int> ExcitSpinChann09(nExcitedStates,izero);
  G4RWTValVector<G4int> ExcitSpinChann10(nExcitedStates,izero);
  G4RWTValVector<G4int> ExcitSpinChann11(nExcitedStates,izero);
  G4RWTValVector<G4int> ExcitSpinChann12(nExcitedStates,izero);
  G4RWTValVector<G4int> ExcitSpinChann13(nExcitedStates,izero);
  G4RWTValVector<G4int> ExcitSpinChann14(nExcitedStates,izero);
  G4RWTValVector<G4int> ExcitSpinChann15(nExcitedStates,izero);
  G4RWTValVector<G4int> ExcitSpinChann16(nExcitedStates,izero);
  G4RWTValVector<G4int> ExcitSpinChann17(nExcitedStates,izero);
  G4RWTValVector<G4int> ExcitSpinChann18(nExcitedStates,izero);
  G4RWTValVector<G4int> ExcitSpinChann19(nExcitedStates,izero);
  G4RWTValVector<G4int> ExcitSpinChann20(nExcitedStates,izero);
  G4RWTValVector<G4int> ExcitSpinChann21(nExcitedStates,izero);
  G4RWTValVector<G4int> ExcitSpinChann22(nExcitedStates,izero);
  G4RWTValVector<G4int> ExcitSpinChann23(nExcitedStates,izero);
  G4RWTValVector<G4int> ExcitSpinChann24(nExcitedStates,izero);
  G4RWTValVector<G4int> ExcitSpinChann25(nExcitedStates,izero);
  G4RWTValVector<G4int> ExcitSpinChann26(nExcitedStates,izero);
  G4RWTValVector<G4int> ExcitSpinChann27(nExcitedStates,izero);
  G4RWTValVector<G4int> ExcitSpinChann28(nExcitedStates,izero);
  G4RWTValVector<G4int> ExcitSpinChann29(nExcitedStates,izero);
  G4RWTValVector<G4int> ExcitSpinChann30(nExcitedStates,izero);
  G4RWTValVector<G4int> ExcitSpinChann31(nExcitedStates,izero);


  // (in MeV)

  ExcitEnergyChann09(0) = 3.56;

  ExcitEnergyChann10(0) = 0.48;

  ExcitEnergyChann11(0) = 0.98;

  ExcitEnergyChann12(0) = 0.43;
  
  ExcitEnergyChann15(0) = 3.37;
  ExcitEnergyChann15(1) = 5.96;
  ExcitEnergyChann15(2) = 6.18;
  ExcitEnergyChann15(3) = 6.26;

  ExcitEnergyChann17(0) = 0.72;
  ExcitEnergyChann17(1) = 1.74;
  ExcitEnergyChann17(2) = 2.15;
  ExcitEnergyChann17(3) = 3.59;

  ExcitEnergyChann18(0) = 2.13;
  ExcitEnergyChann18(1) = 4.44;
  ExcitEnergyChann18(2) = 5.02;
  ExcitEnergyChann18(3) = 6.76;
  ExcitEnergyChann18(4) = 7.29;
  ExcitEnergyChann18(5) = 7.98;
  ExcitEnergyChann18(6) = 8.56;

  ExcitEnergyChann19(0) = 0.95;
  ExcitEnergyChann19(1) = 1.67;
  ExcitEnergyChann19(2) = 6.25;

  ExcitEnergyChann20(0) = 2.00;
  ExcitEnergyChann20(1) = 4.32;
  ExcitEnergyChann20(2) = 4.80;
  ExcitEnergyChann20(3) = 6.34;
  ExcitEnergyChann20(4) = 6.48;
  ExcitEnergyChann20(5) = 6.90;
  ExcitEnergyChann20(6) = 7.50;
  ExcitEnergyChann20(7) = 8.10;
  ExcitEnergyChann20(8) = 8.42;
  ExcitEnergyChann20(9) = 8.66;

  ExcitEnergyChann21(0) = 4.44;

  ExcitEnergyChann22(0) = 3.09;
  ExcitEnergyChann22(1) = 3.68;
  ExcitEnergyChann22(2) = 3.85;

  ExcitEnergyChann23(0) = 6.09;
  ExcitEnergyChann23(1) = 6.69;
  ExcitEnergyChann23(2) = 6.96;
  ExcitEnergyChann23(3) = 7.34;

  ExcitEnergyChann25(0) = 2.31;
  ExcitEnergyChann25(1) = 3.95;
  ExcitEnergyChann25(2) = 4.92;
  ExcitEnergyChann25(3) = 5.11;
  ExcitEnergyChann25(4) = 5.69;
  ExcitEnergyChann25(5) = 5.83;
  ExcitEnergyChann25(6) = 6.20;
  ExcitEnergyChann25(7) = 6.44;
  ExcitEnergyChann25(8) = 7.03;

  ExcitEnergyChann26(0) = 5.28;
  ExcitEnergyChann26(1) = 6.32;
  ExcitEnergyChann26(2) = 7.22;
  ExcitEnergyChann26(3) = 7.57;
  ExcitEnergyChann26(4) = 8.31;
  ExcitEnergyChann26(5) = 8.57;
  ExcitEnergyChann26(6) = 9.15;
  ExcitEnergyChann26(7) = 9.79;
  ExcitEnergyChann26(8) = 10.0;

  ExcitEnergyChann27(0) = 0.12;
  ExcitEnergyChann27(1) = 0.30;
  ExcitEnergyChann27(2) = 0.40;

  ExcitEnergyChann28(0) = 5.22;
  ExcitEnergyChann28(1) = 6.18;
  ExcitEnergyChann28(2) = 6.83;
  ExcitEnergyChann28(3) = 7.28;

  ExcitEnergyChann29(0) = 6.10;
  ExcitEnergyChann29(1) = 6.92;
  ExcitEnergyChann29(2) = 7.12;

  ExcitEnergyChann30(0) = 0.87;
  ExcitEnergyChann30(1) = 3.06;
  ExcitEnergyChann30(2) = 3.84;

  ExcitEnergyChann31(0) = 1.98;
  ExcitEnergyChann31(1) = 3.57;
  ExcitEnergyChann31(2) = 3.92;
  ExcitEnergyChann31(3) = 4.46;
  ExcitEnergyChann31(4) = 5.10;
  ExcitEnergyChann31(5) = 5.33;
  ExcitEnergyChann31(6) = 5.53;
  ExcitEnergyChann31(7) = 6.20;
  ExcitEnergyChann31(8) = 6.38;
  ExcitEnergyChann31(9) = 6.88;


  ExcitSpinChann09(0) = 1;

  ExcitSpinChann10(0) = 2;

  ExcitSpinChann11(0) = 3;

  ExcitSpinChann12(0) = 2;

  ExcitSpinChann15(0) = 5;
  ExcitSpinChann15(1) = 8;
  ExcitSpinChann15(2) = 1;
  ExcitSpinChann15(3) = 5;

  ExcitSpinChann17(0) = 3;
  ExcitSpinChann17(1) = 1;
  ExcitSpinChann17(2) = 3;
  ExcitSpinChann17(3) = 5;

  ExcitSpinChann18(0) = 2;
  ExcitSpinChann18(1) = 6;
  ExcitSpinChann18(2) = 4;
  ExcitSpinChann18(3) = 10;
  ExcitSpinChann18(4) = 6;
  ExcitSpinChann18(5) = 4;
  ExcitSpinChann18(6) = 6;

  ExcitSpinChann19(0) = 5;
  ExcitSpinChann19(1) = 5;
  ExcitSpinChann19(2) = 4;

  ExcitSpinChann20(0) = 2;
  ExcitSpinChann20(1) = 6;
  ExcitSpinChann20(2) = 4;
  ExcitSpinChann20(3) = 2;
  ExcitSpinChann20(4) = 8;
  ExcitSpinChann20(5) = 6;
  ExcitSpinChann20(6) = 4;
  ExcitSpinChann20(7) = 4;
  ExcitSpinChann20(8) = 6;
  ExcitSpinChann20(9) = 8;

  ExcitSpinChann21(0) = 5;

  ExcitSpinChann22(0) = 2;
  ExcitSpinChann22(1) = 4;
  ExcitSpinChann22(2) = 6;

  ExcitSpinChann23(0) = 3;
  ExcitSpinChann23(1) = 8;
  ExcitSpinChann23(2) = 6;
  ExcitSpinChann23(3) = 5;

  ExcitSpinChann25(0) = 1;
  ExcitSpinChann25(1) = 3;
  ExcitSpinChann25(2) = 1;
  ExcitSpinChann25(3) = 5;
  ExcitSpinChann25(4) = 3;
  ExcitSpinChann25(5) = 7;
  ExcitSpinChann25(6) = 3;
  ExcitSpinChann25(7) = 7;
  ExcitSpinChann25(8) = 5;

  ExcitSpinChann26(0) = 8;
  ExcitSpinChann26(1) = 4;
  ExcitSpinChann26(2) = 10;
  ExcitSpinChann26(3) = 8;
  ExcitSpinChann26(4) = 2;
  ExcitSpinChann26(5) = 4;
  ExcitSpinChann26(6) = 14;
  ExcitSpinChann26(7) = 14;
  ExcitSpinChann26(8) = 8;

  ExcitSpinChann27(0) = 1;
  ExcitSpinChann27(1) = 7;
  ExcitSpinChann27(2) = 3;

  ExcitSpinChann28(0) = 8;
  ExcitSpinChann28(1) = 4;
  ExcitSpinChann28(2) = 10;
  ExcitSpinChann28(3) = 8;

  ExcitSpinChann29(0) = 8;
  ExcitSpinChann29(1) = 5;
  ExcitSpinChann29(2) = 3;

  ExcitSpinChann30(0) = 2;
  ExcitSpinChann30(1) = 2;
  ExcitSpinChann30(2) = 6;

  ExcitSpinChann31(0) = 5;
  ExcitSpinChann31(1) = 10;
  ExcitSpinChann31(2) = 5;
  ExcitSpinChann31(3) = 3;
  ExcitSpinChann31(4) = 7;
  ExcitSpinChann31(5) = 13;
  ExcitSpinChann31(6) = 5;
  ExcitSpinChann31(7) = 3;
  ExcitSpinChann31(8) = 12;
  ExcitSpinChann31(9) = 1;

  G4RWTPtrOrderedVector<G4VEvaporationChannel>* theChannels = new G4RWTPtrOrderedVector<G4VEvaporationChannel>;

  //                                        |Gamma|A| Z| 
  //                                        +-----+-+--+

  theChannels->insert(new G4EvaporationChannel(  2,  1, 0, &ExcitEnergyChann00, &ExcitSpinChann00)); // n
  theChannels->insert(new G4EvaporationChannel(  2,  1, 1, &ExcitEnergyChann01, &ExcitSpinChann01)); // p
  theChannels->insert(new G4EvaporationChannel(  6,  2, 1, &ExcitEnergyChann02, &ExcitSpinChann02)); // H2
  theChannels->insert(new G4EvaporationChannel(  6,  3, 1, &ExcitEnergyChann03, &ExcitSpinChann03)); // H3
  theChannels->insert(new G4EvaporationChannel(  6,  3, 2, &ExcitEnergyChann04, &ExcitSpinChann04)); // He3
  theChannels->insert(new G4EvaporationChannel(  4,  4, 2, &ExcitEnergyChann05, &ExcitSpinChann05)); // He4
  theChannels->insert(new G4EvaporationChannel( 20,  5, 2, &ExcitEnergyChann06, &ExcitSpinChann06)); // He5
  theChannels->insert(new G4EvaporationChannel( 30,  6, 2, &ExcitEnergyChann07, &ExcitSpinChann07)); // He6
  theChannels->insert(new G4EvaporationChannel( 20,  5, 3, &ExcitEnergyChann08, &ExcitSpinChann08)); // Li5
  theChannels->insert(new G4EvaporationChannel( 54,  6, 3, &ExcitEnergyChann09, &ExcitSpinChann09)); // Li6
  theChannels->insert(new G4EvaporationChannel( 73,  7, 3, &ExcitEnergyChann10, &ExcitSpinChann10)); // Li7
  theChannels->insert(new G4EvaporationChannel(101,  8, 3, &ExcitEnergyChann11, &ExcitSpinChann11)); // Li8
  theChannels->insert(new G4EvaporationChannel( 73,  7, 4, &ExcitEnergyChann12, &ExcitSpinChann12)); // Be7
  theChannels->insert(new G4EvaporationChannel(  8,  8, 4, &ExcitEnergyChann13, &ExcitSpinChann13)); // Be8
  theChannels->insert(new G4EvaporationChannel(146,  9, 4, &ExcitEnergyChann14, &ExcitSpinChann14)); // Be9
  theChannels->insert(new G4EvaporationChannel(100, 10, 4, &ExcitEnergyChann15, &ExcitSpinChann15)); // Be10 
  theChannels->insert(new G4EvaporationChannel(100,  9, 5, &ExcitEnergyChann16, &ExcitSpinChann16)); // B9
  theChannels->insert(new G4EvaporationChannel(343, 10, 5, &ExcitEnergyChann17, &ExcitSpinChann17)); // B10
  theChannels->insert(new G4EvaporationChannel(174, 11, 5, &ExcitEnergyChann18, &ExcitSpinChann18)); // B11
  theChannels->insert(new G4EvaporationChannel(393, 12, 5, &ExcitEnergyChann19, &ExcitSpinChann19)); // B12
  theChannels->insert(new G4EvaporationChannel(186, 11, 6, &ExcitEnergyChann20, &ExcitSpinChann20)); // C11
  theChannels->insert(new G4EvaporationChannel( 61, 12, 6, &ExcitEnergyChann21, &ExcitSpinChann21)); // C12
  theChannels->insert(new G4EvaporationChannel(202, 13, 6, &ExcitEnergyChann22, &ExcitSpinChann22)); // C13
  theChannels->insert(new G4EvaporationChannel(113, 14, 6, &ExcitEnergyChann23, &ExcitSpinChann23)); // C14
  theChannels->insert(new G4EvaporationChannel(213, 13, 7, &ExcitEnergyChann24, &ExcitSpinChann24)); // N13
  theChannels->insert(new G4EvaporationChannel(233, 14, 7, &ExcitEnergyChann25, &ExcitSpinChann25)); // N14
  theChannels->insert(new G4EvaporationChannel(180, 15, 7, &ExcitEnergyChann26, &ExcitSpinChann26)); // N15
  theChannels->insert(new G4EvaporationChannel(696, 16, 7, &ExcitEnergyChann27, &ExcitSpinChann27)); // N16
  theChannels->insert(new G4EvaporationChannel(194, 15, 8, &ExcitEnergyChann28, &ExcitSpinChann28)); // O15
  theChannels->insert(new G4EvaporationChannel(120, 16, 8, &ExcitEnergyChann29, &ExcitSpinChann29)); // O16
  theChannels->insert(new G4EvaporationChannel(458, 17, 8, &ExcitEnergyChann30, &ExcitSpinChann30)); // O17
  theChannels->insert(new G4EvaporationChannel(590, 18, 8, &ExcitEnergyChann31, &ExcitSpinChann31)); // O18
  theChannels->insert(new G4CompetitiveFission()); // Fission Channel
  theChannels->insert(new G4PhotonEvaporation()); // Photon Channel

  G4int nChannels = theChannels->entries();

  G4PhotonEvaporation* photonEvaporation = new G4PhotonEvaporation;

  G4int i;
  for (i=0; i<iter; i++)
    {
      G4double excitation = excMin + G4UniformRand() * (excMax - excMin);

      G4cout << G4endl << "TEST >>>>>>>>> Iteration " << i 
		 << " <<<<<<<<< Initial excitation " << excitation << G4endl;
	  
      G4LorentzVector p4(0.,0.,0.,G4NucleiProperties::GetAtomicMass(A,Z)+excitation);
      G4Fragment nucleus(A,Z,p4);

      ntuple->column("exc",excitation);

      // Photon evaporation probability
      photonEvaporation->Initialize(nucleus);
      G4FragmentVector* products = photonEvaporation->BreakUp(nucleus);
      G4double probPhoton = photonEvaporation->GetEmissionProbability();
      G4cout << "TEST: Photon probability = " << probPhoton << G4endl;
      ntuple->column("probg",probPhoton);

      // Probability for each channel

      G4int n;
      for (n=0; n<=nProbs; n++)
	{
	  G4Fragment initialNucleus(A,Z,p4);
          theChannels->at(n)->Initialize(initialNucleus);
	  G4double prob = theChannels->at(n)->GetEmissionProbability();
	  G4cout << "TEST: probability(" << n << ") = " << prob << G4endl;
	  ntuple->column("prob"+HepString(n),prob);
	}


      // Fill histograms 

      G4int nProducts = 0;
      if (products !=0) nProducts = products->entries();

      G4cout << "TEST: " << nProducts << " products generated: gammaE ";
      
      G4int ig = 0;
      for (ig=0; ig<nProducts; ig++)
	{
	  G4double productE = products->at(ig)->GetMomentum().e();
	  G4ThreeVector pProd(products->at(ig)->GetMomentum());
	  if (products->at(ig)->GetA() < 1)
	    {
	      G4cout << productE << " ";
	      hGammaE->accumulate(productE);
              ntuple->column("egamma",productE);
	    }
	  else
	    { G4cout << "(" << products->at(ig)->GetZ() << "," << products->at(ig)->GetA() << ") "; }
	}
      G4cout << G4endl;
      
      ntuple->dumpData();
      
      delete products;

    }

  hbookManager->write();

  theChannels->clearAndDestroy();
  delete theChannels;
  delete photonEvaporation;
  
  return EXIT_SUCCESS;
}



















