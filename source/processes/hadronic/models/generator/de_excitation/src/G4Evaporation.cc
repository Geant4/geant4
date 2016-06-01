// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)

#include "G4Evaporation.hh"


G4Evaporation::G4Evaporation()
{
  ExcitEnergyChann00.reshape(NumExcitedStates);
  ExcitEnergyChann01.reshape(NumExcitedStates);
  ExcitEnergyChann02.reshape(NumExcitedStates);
  ExcitEnergyChann03.reshape(NumExcitedStates);
  ExcitEnergyChann04.reshape(NumExcitedStates);
  ExcitEnergyChann05.reshape(NumExcitedStates);
  ExcitEnergyChann06.reshape(NumExcitedStates);
  ExcitEnergyChann07.reshape(NumExcitedStates);
  ExcitEnergyChann08.reshape(NumExcitedStates);
  ExcitEnergyChann09.reshape(NumExcitedStates);
  ExcitEnergyChann10.reshape(NumExcitedStates);
  ExcitEnergyChann11.reshape(NumExcitedStates);
  ExcitEnergyChann12.reshape(NumExcitedStates);
  ExcitEnergyChann13.reshape(NumExcitedStates);
  ExcitEnergyChann14.reshape(NumExcitedStates);
  ExcitEnergyChann15.reshape(NumExcitedStates);
  ExcitEnergyChann16.reshape(NumExcitedStates);
  ExcitEnergyChann17.reshape(NumExcitedStates);
  ExcitEnergyChann18.reshape(NumExcitedStates);
  ExcitEnergyChann19.reshape(NumExcitedStates);
  ExcitEnergyChann20.reshape(NumExcitedStates);
  ExcitEnergyChann21.reshape(NumExcitedStates);
  ExcitEnergyChann22.reshape(NumExcitedStates);
  ExcitEnergyChann23.reshape(NumExcitedStates);
  ExcitEnergyChann24.reshape(NumExcitedStates);
  ExcitEnergyChann25.reshape(NumExcitedStates);
  ExcitEnergyChann26.reshape(NumExcitedStates);
  ExcitEnergyChann27.reshape(NumExcitedStates);
  ExcitEnergyChann28.reshape(NumExcitedStates);
  ExcitEnergyChann29.reshape(NumExcitedStates); 
  ExcitEnergyChann30.reshape(NumExcitedStates);
  ExcitEnergyChann31.reshape(NumExcitedStates);

  ExcitSpinChann00.reshape(NumExcitedStates);
  ExcitSpinChann01.reshape(NumExcitedStates);
  ExcitSpinChann02.reshape(NumExcitedStates);
  ExcitSpinChann03.reshape(NumExcitedStates);
  ExcitSpinChann04.reshape(NumExcitedStates);
  ExcitSpinChann05.reshape(NumExcitedStates);
  ExcitSpinChann06.reshape(NumExcitedStates);
  ExcitSpinChann07.reshape(NumExcitedStates);
  ExcitSpinChann08.reshape(NumExcitedStates);
  ExcitSpinChann09.reshape(NumExcitedStates);
  ExcitSpinChann10.reshape(NumExcitedStates);
  ExcitSpinChann11.reshape(NumExcitedStates);
  ExcitSpinChann12.reshape(NumExcitedStates);
  ExcitSpinChann13.reshape(NumExcitedStates);
  ExcitSpinChann14.reshape(NumExcitedStates);
  ExcitSpinChann15.reshape(NumExcitedStates);
  ExcitSpinChann16.reshape(NumExcitedStates);
  ExcitSpinChann17.reshape(NumExcitedStates);
  ExcitSpinChann18.reshape(NumExcitedStates);
  ExcitSpinChann19.reshape(NumExcitedStates);
  ExcitSpinChann20.reshape(NumExcitedStates);
  ExcitSpinChann21.reshape(NumExcitedStates);
  ExcitSpinChann22.reshape(NumExcitedStates);
  ExcitSpinChann23.reshape(NumExcitedStates);
  ExcitSpinChann24.reshape(NumExcitedStates);
  ExcitSpinChann25.reshape(NumExcitedStates);
  ExcitSpinChann26.reshape(NumExcitedStates);
  ExcitSpinChann27.reshape(NumExcitedStates);
  ExcitSpinChann28.reshape(NumExcitedStates);
  ExcitSpinChann29.reshape(NumExcitedStates); 
  ExcitSpinChann30.reshape(NumExcitedStates);
  ExcitSpinChann31.reshape(NumExcitedStates);

  for (G4int i = 0; i < NumExcitedStates; i++) {
    ExcitEnergyChann00(i) = 0.0;
    ExcitEnergyChann01(i) = 0.0;
    ExcitEnergyChann02(i) = 0.0;
    ExcitEnergyChann03(i) = 0.0;
    ExcitEnergyChann04(i) = 0.0;
    ExcitEnergyChann05(i) = 0.0;
    ExcitEnergyChann06(i) = 0.0;
    ExcitEnergyChann07(i) = 0.0;
    ExcitEnergyChann08(i) = 0.0;
    ExcitEnergyChann09(i) = 0.0;
    ExcitEnergyChann10(i) = 0.0;
    ExcitEnergyChann11(i) = 0.0;
    ExcitEnergyChann12(i) = 0.0;
    ExcitEnergyChann13(i) = 0.0;
    ExcitEnergyChann14(i) = 0.0;
    ExcitEnergyChann15(i) = 0.0;
    ExcitEnergyChann16(i) = 0.0;
    ExcitEnergyChann17(i) = 0.0;
    ExcitEnergyChann18(i) = 0.0;
    ExcitEnergyChann19(i) = 0.0;
    ExcitEnergyChann20(i) = 0.0;
    ExcitEnergyChann21(i) = 0.0;
    ExcitEnergyChann22(i) = 0.0;
    ExcitEnergyChann23(i) = 0.0;
    ExcitEnergyChann24(i) = 0.0;
    ExcitEnergyChann25(i) = 0.0;
    ExcitEnergyChann26(i) = 0.0;
    ExcitEnergyChann27(i) = 0.0;
    ExcitEnergyChann28(i) = 0.0;
    ExcitEnergyChann29(i) = 0.0; 
    ExcitEnergyChann30(i) = 0.0;
    ExcitEnergyChann31(i) = 0.0;

    ExcitSpinChann00(i) = 0;
    ExcitSpinChann01(i) = 0;
    ExcitSpinChann02(i) = 0;
    ExcitSpinChann03(i) = 0;
    ExcitSpinChann04(i) = 0;
    ExcitSpinChann05(i) = 0;
    ExcitSpinChann06(i) = 0;
    ExcitSpinChann07(i) = 0;
    ExcitSpinChann08(i) = 0;
    ExcitSpinChann09(i) = 0;
    ExcitSpinChann10(i) = 0;
    ExcitSpinChann11(i) = 0;
    ExcitSpinChann12(i) = 0;
    ExcitSpinChann13(i) = 0;
    ExcitSpinChann14(i) = 0;
    ExcitSpinChann15(i) = 0;
    ExcitSpinChann16(i) = 0;
    ExcitSpinChann17(i) = 0;
    ExcitSpinChann18(i) = 0;
    ExcitSpinChann19(i) = 0;
    ExcitSpinChann20(i) = 0;
    ExcitSpinChann21(i) = 0;
    ExcitSpinChann22(i) = 0;
    ExcitSpinChann23(i) = 0;
    ExcitSpinChann24(i) = 0;
    ExcitSpinChann25(i) = 0;
    ExcitSpinChann26(i) = 0;
    ExcitSpinChann27(i) = 0;
    ExcitSpinChann28(i) = 0;
    ExcitSpinChann29(i) = 0; 
    ExcitSpinChann30(i) = 0;
    ExcitSpinChann31(i) = 0;  
  }

  // (in MeV)
  // neutrons
  ExcitEnergyChann00( 9) = 3.56;
  ExcitEnergyChann00(10) = 0.48;
  ExcitEnergyChann00(11) = 0.98;
  ExcitEnergyChann00(12) = 0.43;
  ExcitEnergyChann00(15) = 3.37;
  ExcitEnergyChann00(17) = 0.72;
  ExcitEnergyChann00(18) = 2.13;
  ExcitEnergyChann00(19) = 0.95;
  ExcitEnergyChann00(20) = 2.00;
  ExcitEnergyChann00(21) = 4.44;
  ExcitEnergyChann00(22) = 3.09;
  ExcitEnergyChann00(23) = 6.09;
  ExcitEnergyChann00(25) = 2.31;
  ExcitEnergyChann00(26) = 5.28;
  ExcitEnergyChann00(27) = 0.12;
  ExcitEnergyChann00(28) = 5.22;
  ExcitEnergyChann00(29) = 6.10;
  ExcitEnergyChann00(30) = 0.87;
  ExcitEnergyChann00(31) = 1.98;
  
  // protons
  ExcitEnergyChann01(15) = 5.96;
  ExcitEnergyChann01(17) = 1.74;
  ExcitEnergyChann01(18) = 4.44;
  ExcitEnergyChann01(19) = 1.67;
  ExcitEnergyChann01(20) = 4.32;
  ExcitEnergyChann01(22) = 3.68;
  ExcitEnergyChann01(23) = 6.69;
  ExcitEnergyChann01(25) = 3.95;
  ExcitEnergyChann01(26) = 6.32;
  ExcitEnergyChann01(27) = 0.30;
  ExcitEnergyChann01(28) = 6.18;
  ExcitEnergyChann01(29) = 6.92;
  ExcitEnergyChann01(30) = 3.06;
  ExcitEnergyChann01(31) = 3.57;

  // deuterons
  ExcitEnergyChann02(15) = 6.18;
  ExcitEnergyChann02(17) = 2.15;
  ExcitEnergyChann02(18) = 5.02;
  ExcitEnergyChann02(19) = 2.65;
  ExcitEnergyChann02(20) = 4.80;
  ExcitEnergyChann02(22) = 3.85;
  ExcitEnergyChann02(23) = 6.96;
  ExcitEnergyChann02(25) = 4.92;
  ExcitEnergyChann02(26) = 7.22;
  ExcitEnergyChann02(27) = 0.40;
  ExcitEnergyChann02(28) = 6.83;
  ExcitEnergyChann02(29) = 7.12;
  ExcitEnergyChann02(30) = 3.84;
  ExcitEnergyChann02(31) = 3.92;

  // tritons
  ExcitEnergyChann03(15) = 6.26;
  ExcitEnergyChann03(17) = 3.59;
  ExcitEnergyChann03(18) = 6.76;
  ExcitEnergyChann03(20) = 6.34;
  ExcitEnergyChann03(23) = 7.34;
  ExcitEnergyChann03(25) = 5.11;
  ExcitEnergyChann03(26) = 7.57;
  ExcitEnergyChann03(28) = 7.28;
  ExcitEnergyChann03(31) = 4.46;


  // He3
  ExcitEnergyChann04(18) = 7.29;
  ExcitEnergyChann04(20) = 6.48;
  ExcitEnergyChann04(25) = 5.69;
  ExcitEnergyChann04(26) = 8.31;
  ExcitEnergyChann04(31) = 5.10;

  // alphas
  ExcitEnergyChann05(18) = 7.98;
  ExcitEnergyChann05(20) = 6.90;
  ExcitEnergyChann05(25) = 5.83;
  ExcitEnergyChann05(26) = 8.57;
  ExcitEnergyChann05(31) = 5.33;
  
  // He5
  ExcitEnergyChann06(18) = 8.56;
  ExcitEnergyChann06(20) = 7.50;
  ExcitEnergyChann06(25) = 6.20;
  ExcitEnergyChann06(26) = 9.15;
  ExcitEnergyChann06(31) = 5.53;

  // He6
  ExcitEnergyChann07(20) = 8.10;
  ExcitEnergyChann07(25) = 6.44;
  ExcitEnergyChann07(26) = 9.79;
  ExcitEnergyChann07(31) = 6.20;


  // Li5 
  ExcitEnergyChann08(20) = 8.42;
  ExcitEnergyChann08(25) = 7.03;
  ExcitEnergyChann08(26) = 10.0;
  ExcitEnergyChann08(31) = 6.38;


  // Li6
  ExcitEnergyChann09(20) = 8.66;
  ExcitEnergyChann09(31) = 6.88;

  
  // Spin (2s+1)

  // neutrons
  ExcitSpinChann00( 9) = 1;
  ExcitSpinChann00(10) = 2;
  ExcitSpinChann00(11) = 3;
  ExcitSpinChann00(12) = 2;
  ExcitSpinChann00(15) = 5;
  ExcitSpinChann00(17) = 3;
  ExcitSpinChann00(18) = 2;
  ExcitSpinChann00(19) = 5;
  ExcitSpinChann00(20) = 2;
  ExcitSpinChann00(21) = 5;
  ExcitSpinChann00(22) = 2;
  ExcitSpinChann00(23) = 3;
  ExcitSpinChann00(25) = 1;
  ExcitSpinChann00(26) = 8;
  ExcitSpinChann00(27) = 1;
  ExcitSpinChann00(28) = 8;
  ExcitSpinChann00(29) = 8;
  ExcitSpinChann00(30) = 2;
  ExcitSpinChann00(31) = 5;

  // protons
  ExcitSpinChann01(15) = 8;
  ExcitSpinChann01(17) = 1;
  ExcitSpinChann01(18) = 6;
  ExcitSpinChann01(19) = 5;
  ExcitSpinChann01(20) = 6;
  ExcitSpinChann01(22) = 4;
  ExcitSpinChann01(23) = 8;
  ExcitSpinChann01(25) = 3;
  ExcitSpinChann01(26) = 4;
  ExcitSpinChann01(27) = 7;
  ExcitSpinChann01(28) = 4;
  ExcitSpinChann01(29) = 5;
  ExcitSpinChann01(30) = 2;
  ExcitSpinChann01(31) = 10;

  // deuterons
  ExcitSpinChann02(15) = 1;
  ExcitSpinChann02(17) = 3;
  ExcitSpinChann02(18) = 4;
  ExcitSpinChann02(19) = 4;
  ExcitSpinChann02(20) = 4;
  ExcitSpinChann02(22) = 6;
  ExcitSpinChann02(23) = 6;
  ExcitSpinChann02(25) = 1;
  ExcitSpinChann02(26) = 10;
  ExcitSpinChann02(27) = 3;
  ExcitSpinChann02(28) = 10;
  ExcitSpinChann02(29) = 3;
  ExcitSpinChann02(30) = 6;
  ExcitSpinChann02(31) = 5;

  // tritons
  ExcitSpinChann03(15) = 5;
  ExcitSpinChann03(17) = 5;
  ExcitSpinChann03(18) = 10;
  ExcitSpinChann03(20) = 2;
  ExcitSpinChann03(23) = 5;
  ExcitSpinChann03(25) = 5;
  ExcitSpinChann03(26) = 8;
  ExcitSpinChann03(28) = 8;
  ExcitSpinChann03(31) = 3;


  // He3
  ExcitSpinChann04(18) = 6;
  ExcitSpinChann04(20) = 8;
  ExcitSpinChann04(25) = 3;
  ExcitSpinChann04(26) = 2;
  ExcitSpinChann04(31) = 7;

  // alphas
  ExcitSpinChann05(18) = 4;
  ExcitSpinChann05(20) = 6;
  ExcitSpinChann05(25) = 7;
  ExcitSpinChann05(26) = 4;
  ExcitSpinChann05(31) = 13;
  
  // He5
  ExcitSpinChann06(18) = 6;
  ExcitSpinChann06(20) = 4;
  ExcitSpinChann06(25) = 3;
  ExcitSpinChann06(26) = 14;
  ExcitSpinChann06(31) = 5;

  // He6 
  ExcitSpinChann07(20) = 4;
  ExcitSpinChann07(25) = 7;
  ExcitSpinChann07(26) = 14;
  ExcitSpinChann07(31) = 3;

  // Li5 
  ExcitSpinChann08(20) = 6;
  ExcitSpinChann08(25) = 5;
  ExcitSpinChann08(26) = 8;
  ExcitSpinChann08(31) = 12;

  // Li6
  ExcitSpinChann09(20) = 8;
  ExcitSpinChann09(31) = 1;
  

  //                                        |Gamma|A| Z| 
  //                                        +-----+-+--+
  theChannels[ 0] = new G4EvaporationChannel(  2,  1, 0, &ExcitEnergyChann00, &ExcitSpinChann00); // n
  theChannels[ 1] = new G4EvaporationChannel(  2,  1, 1, &ExcitEnergyChann01, &ExcitSpinChann01); // p
  theChannels[ 2] = new G4EvaporationChannel(  6,  2, 1, &ExcitEnergyChann02, &ExcitSpinChann02); // H2
  theChannels[ 3] = new G4EvaporationChannel(  6,  3, 1, &ExcitEnergyChann03, &ExcitSpinChann03); // H3
  theChannels[ 4] = new G4EvaporationChannel(  6,  3, 2, &ExcitEnergyChann04, &ExcitSpinChann04); // He3
  theChannels[ 5] = new G4EvaporationChannel(  4,  4, 2, &ExcitEnergyChann05, &ExcitSpinChann05); // He4
  theChannels[ 6] = new G4EvaporationChannel( 20,  5, 2, &ExcitEnergyChann06, &ExcitSpinChann06); // He5
  theChannels[ 7] = new G4EvaporationChannel( 30,  6, 2, &ExcitEnergyChann07, &ExcitSpinChann07); // He6
  theChannels[ 8] = new G4EvaporationChannel( 20,  5, 3, &ExcitEnergyChann08, &ExcitSpinChann08); // Li5
  theChannels[ 9] = new G4EvaporationChannel( 54,  6, 3, &ExcitEnergyChann09, &ExcitSpinChann09); // Li6
  theChannels[10] = new G4EvaporationChannel( 73,  7, 3, &ExcitEnergyChann10, &ExcitSpinChann10); // Li7
  theChannels[11] = new G4EvaporationChannel(101,  8, 3, &ExcitEnergyChann11, &ExcitSpinChann11); // Li8
  theChannels[12] = new G4EvaporationChannel( 73,  7, 4, &ExcitEnergyChann12, &ExcitSpinChann12); // Be7
  theChannels[13] = new G4EvaporationChannel(  8,  8, 4, &ExcitEnergyChann13, &ExcitSpinChann13); // Be8
  theChannels[14] = new G4EvaporationChannel(146,  9, 4, &ExcitEnergyChann14, &ExcitSpinChann14); // Be9
  theChannels[15] = new G4EvaporationChannel(100, 10, 4, &ExcitEnergyChann15, &ExcitSpinChann15); // Be10 
  theChannels[16] = new G4EvaporationChannel(100,  9, 5, &ExcitEnergyChann16, &ExcitSpinChann16); // B9
  theChannels[17] = new G4EvaporationChannel(343, 10, 5, &ExcitEnergyChann17, &ExcitSpinChann17); // B10
  theChannels[18] = new G4EvaporationChannel(174, 11, 5, &ExcitEnergyChann18, &ExcitSpinChann18); // B11
  theChannels[19] = new G4EvaporationChannel(393, 12, 5, &ExcitEnergyChann19, &ExcitSpinChann19); // B12
  theChannels[20] = new G4EvaporationChannel(186, 11, 6, &ExcitEnergyChann20, &ExcitSpinChann20); // C11
  theChannels[21] = new G4EvaporationChannel( 61, 12, 6, &ExcitEnergyChann21, &ExcitSpinChann21); // C12
  theChannels[22] = new G4EvaporationChannel(202, 13, 6, &ExcitEnergyChann22, &ExcitSpinChann22); // C13
  theChannels[23] = new G4EvaporationChannel(113, 14, 6, &ExcitEnergyChann23, &ExcitSpinChann23); // C14
  theChannels[24] = new G4EvaporationChannel(213, 13, 7, &ExcitEnergyChann24, &ExcitSpinChann24); // N13
  theChannels[25] = new G4EvaporationChannel(233, 14, 7, &ExcitEnergyChann25, &ExcitSpinChann25); // N14
  theChannels[26] = new G4EvaporationChannel(180, 15, 7, &ExcitEnergyChann26, &ExcitSpinChann26); // N15
  theChannels[27] = new G4EvaporationChannel(696, 16, 7, &ExcitEnergyChann27, &ExcitSpinChann27); // N16
  theChannels[28] = new G4EvaporationChannel(194, 15, 8, &ExcitEnergyChann28, &ExcitSpinChann28); // O15
  theChannels[29] = new G4EvaporationChannel(120, 16, 8, &ExcitEnergyChann29, &ExcitSpinChann29); // O16
  theChannels[30] = new G4EvaporationChannel(458, 17, 8, &ExcitEnergyChann30, &ExcitSpinChann30); // O17
  theChannels[31] = new G4EvaporationChannel(590, 18, 8, &ExcitEnergyChann31, &ExcitSpinChann31); // O18 
  theChannels[32] = new G4CompetitiveFission(); // Fission Channel
  theChannels[33] = new G4PhotonEvaporation(); // Photon Channel

}

G4Evaporation::G4Evaporation(const G4Evaporation &right)
{
  G4Exception("G4Evaporation::copy_constructor meant to not be accessable.");
}


G4Evaporation::~G4Evaporation()
{
  for (G4int i = 0; i < TotNumberOfChannels; i++) 
    delete theChannels[i];

}

const G4Evaporation & G4Evaporation::operator=(const G4Evaporation &right)
{
  G4Exception("G4Evaporation::operator= meant to not be accessable.");
  return *this;
}


G4bool G4Evaporation::operator==(const G4Evaporation &right) const
{
  return false;
}

G4bool G4Evaporation::operator!=(const G4Evaporation &right) const
{
  return true;
}


G4FragmentVector * G4Evaporation::BreakItUp(const G4Fragment &theNucleus)
{
  G4FragmentVector * theResult = new G4FragmentVector;

  // CHECK that Excitation Energy != 0
  if (theNucleus.GetExcitationEnergy() == 0) {
    theResult->insert(new G4Fragment(theNucleus));
    return theResult;
  }

  // The residual nucleus (after evaporation of each fragment)
  G4Fragment theResidualNucleus = theNucleus;


  // Starts loop over evaporated particles
  for (;;) {
    // loop over evaporation channels
    G4int i;
    for (i=0; i < TotNumberOfChannels; i++)
      theChannels[i]->Initialize(theResidualNucleus);
    
    // Work out total decay probability by summing over channels 
    G4double TotalProbability = 0;
    for (i=0; i < TotNumberOfChannels; i++)
      TotalProbability += theChannels[i]->GetEmissionProbability();

//     G4cout << "---------------- " << theResidualNucleus.GetExcitationEnergy()/MeV << "-----------------------" << endl;
//     G4cout << "Prob of neutron: " << theChannels[0]->GetEmissionProbability()/TotalProbability << endl;
//     G4cout << "Prob of proton : " << theChannels[1]->GetEmissionProbability()/TotalProbability<< endl;
//     G4cout << "Prob of alpha  : " << theChannels[5]->GetEmissionProbability()/TotalProbability<< endl;
//     G4cout << "Prob of fission: " << theChannels[NumberOfFissionChannel]->GetEmissionProbability()/TotalProbability<< endl;


    if (TotalProbability <= 0.0) {
      // Will be no evaporation more
      // write information about residual nucleus
      theResult->insert(new G4Fragment(theResidualNucleus));
      break; 
    } else {
      // Selection of evaporation channel, fission or gamma
      G4double EmissionProbChannel[TotNumberOfChannels];
	
      
      EmissionProbChannel[0] = theChannels[0]->GetEmissionProbability();


      for (i=1; i < TotNumberOfChannels; i++) 
	EmissionProbChannel[i] = EmissionProbChannel[i-1] + theChannels[i]->GetEmissionProbability();


      G4double shoot = G4UniformRand() * TotalProbability;

      for (i=0; i < TotNumberOfChannels; i++)
	if (shoot < EmissionProbChannel[i]) 
	  break;

      if( i == TotNumberOfChannels )
	G4Exception( "Can't define emission probability of the channels (G4Evaporation::BreakItUp)" );
      else if (i == NumberOfFissionChannel) {
	// Fission has to be performed
	G4FragmentVector * theFissionResult = theChannels[i]->BreakUp(theResidualNucleus);
	while (theFissionResult->entries() > 0)
	  theResult->insert(theFissionResult->removeFirst());

	theFissionResult->clearAndDestroy();
	delete theFissionResult;
	break;

      } else if (i == NumberOfGammaChannel) {
	// Gamma evaporation has to be performed
	G4FragmentVector * theGammaResult = theChannels[i]->BreakUp(theResidualNucleus);
	while (theGammaResult->entries() > 0)
	  theResult->insert(theGammaResult->removeFirst());

	theGammaResult->clearAndDestroy();
	delete theGammaResult;
	break;

      } else {
	// Evaporation has to be performed
	G4FragmentVector * theEvaporationResult = theChannels[i]->BreakUp(theResidualNucleus);
	while (theEvaporationResult->entries() > 1)
	  theResult->insert(theEvaporationResult->removeFirst());

	theResidualNucleus = *(theEvaporationResult->at(0));
	theEvaporationResult->clearAndDestroy();
	delete theEvaporationResult;
      }
    }
  }

  return theResult;
}



