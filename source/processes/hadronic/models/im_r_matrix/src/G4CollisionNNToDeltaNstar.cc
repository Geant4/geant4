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
// $Id: G4CollisionNNToDeltaNstar.cc,v 1.2 2003-12-12 15:38:22 hpw Exp $ //

#include "globals.hh"
#include "G4CollisionNNToDeltaNstar.hh"
#include "G4KineticTrack.hh"
#include "G4VCrossSectionSource.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4LorentzRotation.hh"
#include "G4KineticTrackVector.hh"
#include "G4ParticleTable.hh"
#include "G4CollisionVector.hh"
#include "G4CollisionNNToDeltaNstar.hh"
#include "G4ConcreteNNToDeltaNstar.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4HadParticleCodes.hh"
#include "G4Pair.hh"

typedef G4ConcreteNNToDeltaNstar channelType;

typedef INT4(channelType, NeutronPC, NeutronPC, Delta0PC,  N1400nPC)  theC1;
typedef INT4(channelType, NeutronPC, NeutronPC, DeltamPC,  N1400pPC)  theC2;
typedef INT4(channelType, ProtonPC,  ProtonPC,  DeltapPC,  N1400pPC)  theC3;
typedef INT4(channelType, ProtonPC,  ProtonPC,  DeltappPC, N1400nPC)  theC4;
typedef INT4(channelType, NeutronPC, ProtonPC,  Delta0PC,  N1400pPC)  theC5;
typedef INT4(channelType, NeutronPC, ProtonPC,  DeltapPC,  N1400nPC)  theC6;

typedef INT4(channelType, NeutronPC, NeutronPC, Delta0PC,  N1520nPC)  theC7;
typedef INT4(channelType, NeutronPC, NeutronPC, DeltamPC,  N1520pPC)  theC8;
typedef INT4(channelType, ProtonPC,  ProtonPC,  DeltapPC,  N1520pPC)  theC9;
typedef INT4(channelType, ProtonPC,  ProtonPC,  DeltappPC, N1520nPC)  theC10;
typedef INT4(channelType, NeutronPC, ProtonPC,  Delta0PC,  N1520pPC)  theC11;
typedef INT4(channelType, NeutronPC, ProtonPC,  DeltapPC,  N1520nPC)  theC12;

typedef INT4(channelType, NeutronPC, NeutronPC, Delta0PC,  N1535nPC)  theC13;
typedef INT4(channelType, NeutronPC, NeutronPC, DeltamPC,  N1535pPC)  theC14;
typedef INT4(channelType, ProtonPC,  ProtonPC,  DeltapPC,  N1535pPC)  theC15;
typedef INT4(channelType, ProtonPC,  ProtonPC,  DeltappPC, N1535nPC)  theC16;
typedef INT4(channelType, NeutronPC, ProtonPC,  Delta0PC,  N1535pPC)  theC17;
typedef INT4(channelType, NeutronPC, ProtonPC,  DeltapPC,  N1535nPC)  theC18;

typedef INT4(channelType, NeutronPC, NeutronPC, Delta0PC,  N1650nPC)  theC19;
typedef INT4(channelType, NeutronPC, NeutronPC, DeltamPC,  N1650pPC)  theC20;
typedef INT4(channelType, ProtonPC,  ProtonPC,  DeltapPC,  N1650pPC)  theC21;
typedef INT4(channelType, ProtonPC,  ProtonPC,  DeltappPC, N1650nPC)  theC22;
typedef INT4(channelType, NeutronPC, ProtonPC,  Delta0PC,  N1650pPC)  theC23;
typedef INT4(channelType, NeutronPC, ProtonPC,  DeltapPC,  N1650nPC)  theC24;

typedef INT4(channelType, NeutronPC, NeutronPC, Delta0PC,  N1675nPC)  theC25;
typedef INT4(channelType, NeutronPC, NeutronPC, DeltamPC,  N1675pPC)  theC26;
typedef INT4(channelType, ProtonPC,  ProtonPC,  DeltapPC,  N1675pPC)  theC27;
typedef INT4(channelType, ProtonPC,  ProtonPC,  DeltappPC, N1675nPC)  theC28;
typedef INT4(channelType, NeutronPC, ProtonPC,  Delta0PC,  N1675pPC)  theC29;
typedef INT4(channelType, NeutronPC, ProtonPC,  DeltapPC,  N1675nPC)  theC30;

typedef INT4(channelType, NeutronPC, NeutronPC, Delta0PC,  N1680nPC)  theC31;
typedef INT4(channelType, NeutronPC, NeutronPC, DeltamPC,  N1680pPC)  theC32;
typedef INT4(channelType, ProtonPC,  ProtonPC,  DeltapPC,  N1680pPC)  theC33;
typedef INT4(channelType, ProtonPC,  ProtonPC,  DeltappPC, N1680nPC)  theC34;
typedef INT4(channelType, NeutronPC, ProtonPC,  Delta0PC,  N1680pPC)  theC35;
typedef INT4(channelType, NeutronPC, ProtonPC,  DeltapPC,  N1680nPC)  theC36;

typedef INT4(channelType, NeutronPC, NeutronPC, Delta0PC,  N1700nPC)  theC37;
typedef INT4(channelType, NeutronPC, NeutronPC, DeltamPC,  N1700pPC)  theC38;
typedef INT4(channelType, ProtonPC,  ProtonPC,  DeltapPC,  N1700pPC)  theC39;
typedef INT4(channelType, ProtonPC,  ProtonPC,  DeltappPC, N1700nPC)  theC40;
typedef INT4(channelType, NeutronPC, ProtonPC,  Delta0PC,  N1700pPC)  theC41;
typedef INT4(channelType, NeutronPC, ProtonPC,  DeltapPC,  N1700nPC)  theC42;

typedef INT4(channelType, NeutronPC, NeutronPC, Delta0PC,  N1710nPC)  theC43;
typedef INT4(channelType, NeutronPC, NeutronPC, DeltamPC,  N1710pPC)  theC44;
typedef INT4(channelType, ProtonPC,  ProtonPC,  DeltapPC,  N1710pPC)  theC45;
typedef INT4(channelType, ProtonPC,  ProtonPC,  DeltappPC, N1710nPC)  theC46;
typedef INT4(channelType, NeutronPC, ProtonPC,  Delta0PC,  N1710pPC)  theC47;
typedef INT4(channelType, NeutronPC, ProtonPC,  DeltapPC,  N1710nPC)  theC48;

typedef INT4(channelType, NeutronPC, NeutronPC, Delta0PC,  N1720nPC)  theC49;
typedef INT4(channelType, NeutronPC, NeutronPC, DeltamPC,  N1720pPC)  theC50;
typedef INT4(channelType, ProtonPC,  ProtonPC,  DeltapPC,  N1720pPC)  theC51;
typedef INT4(channelType, ProtonPC,  ProtonPC,  DeltappPC, N1720nPC)  theC52;
typedef INT4(channelType, NeutronPC, ProtonPC,  Delta0PC,  N1720pPC)  theC53;
typedef INT4(channelType, NeutronPC, ProtonPC,  DeltapPC,  N1720nPC)  theC54;

typedef INT4(channelType, NeutronPC, NeutronPC, Delta0PC,  N1900nPC)  theC55;
typedef INT4(channelType, NeutronPC, NeutronPC, DeltamPC,  N1900pPC)  theC56;
typedef INT4(channelType, ProtonPC,  ProtonPC,  DeltapPC,  N1900pPC)  theC57;
typedef INT4(channelType, ProtonPC,  ProtonPC,  DeltappPC, N1900nPC)  theC58;
typedef INT4(channelType, NeutronPC, ProtonPC,  Delta0PC,  N1900pPC)  theC59;
typedef INT4(channelType, NeutronPC, ProtonPC,  DeltapPC,  N1900nPC)  theC60;

typedef INT4(channelType, NeutronPC, NeutronPC, Delta0PC,  N1990nPC)  theC61;
typedef INT4(channelType, NeutronPC, NeutronPC, DeltamPC,  N1990pPC)  theC62;
typedef INT4(channelType, ProtonPC,  ProtonPC,  DeltapPC,  N1990pPC)  theC63;
typedef INT4(channelType, ProtonPC,  ProtonPC,  DeltappPC, N1990nPC)  theC64;
typedef INT4(channelType, NeutronPC, ProtonPC,  Delta0PC,  N1990pPC)  theC65;
typedef INT4(channelType, NeutronPC, ProtonPC,  DeltapPC,  N1990nPC)  theC66;

typedef INT4(channelType, NeutronPC, NeutronPC, Delta0PC,  N2090nPC)  theC67;
typedef INT4(channelType, NeutronPC, NeutronPC, DeltamPC,  N2090pPC)  theC68;
typedef INT4(channelType, ProtonPC,  ProtonPC,  DeltapPC,  N2090pPC)  theC69;
typedef INT4(channelType, ProtonPC,  ProtonPC,  DeltappPC, N2090nPC)  theC70;
typedef INT4(channelType, NeutronPC, ProtonPC,  Delta0PC,  N2090pPC)  theC71;
typedef INT4(channelType, NeutronPC, ProtonPC,  DeltapPC,  N2090nPC)  theC72;

typedef INT4(channelType, NeutronPC, NeutronPC, Delta0PC,  N2190nPC)  theC73;
typedef INT4(channelType, NeutronPC, NeutronPC, DeltamPC,  N2190pPC)  theC74;
typedef INT4(channelType, ProtonPC,  ProtonPC,  DeltapPC,  N2190pPC)  theC75;
typedef INT4(channelType, ProtonPC,  ProtonPC,  DeltappPC, N2190nPC)  theC76;
typedef INT4(channelType, NeutronPC, ProtonPC,  Delta0PC,  N2190pPC)  theC77;
typedef INT4(channelType, NeutronPC, ProtonPC,  DeltapPC,  N2190nPC)  theC78;

typedef INT4(channelType, NeutronPC, NeutronPC, Delta0PC,  N2220nPC)  theC79;
typedef INT4(channelType, NeutronPC, NeutronPC, DeltamPC,  N2220pPC)  theC80;
typedef INT4(channelType, ProtonPC,  ProtonPC,  DeltapPC,  N2220pPC)  theC81;
typedef INT4(channelType, ProtonPC,  ProtonPC,  DeltappPC, N2220nPC)  theC82;
typedef INT4(channelType, NeutronPC, ProtonPC,  Delta0PC,  N2220pPC)  theC83;
typedef INT4(channelType, NeutronPC, ProtonPC,  DeltapPC,  N2220nPC)  theC84;

typedef INT4(channelType, NeutronPC, NeutronPC, Delta0PC,  N2250nPC)  theC85;
typedef INT4(channelType, NeutronPC, NeutronPC, DeltamPC,  N2250pPC)  theC86;
typedef INT4(channelType, ProtonPC,  ProtonPC,  DeltapPC,  N2250pPC)  theC87;
typedef INT4(channelType, ProtonPC,  ProtonPC,  DeltappPC, N2250nPC)  theC88;
typedef INT4(channelType, NeutronPC, ProtonPC,  Delta0PC,  N2250pPC)  theC89;
typedef INT4(channelType, NeutronPC, ProtonPC,  DeltapPC,  N2250nPC)  theC90;

typedef GROUP90(theC1, theC2, theC3, theC4, theC5, theC6, theC7, theC8, theC9, theC10,
              theC11, theC12, theC13, theC14, theC15, theC16, theC17, theC18, theC19, theC20,
              theC21, theC22, theC23, theC24, theC25, theC26, theC27, theC28, theC29, theC30,
              theC31, theC32, theC33, theC34, theC35, theC36, theC37, theC38, theC39, theC40,
              theC41, theC42, theC43, theC44, theC45, theC46, theC47, theC48, theC49, theC50,
              theC51, theC52, theC53, theC54, theC55, theC56, theC57, theC58, theC59, theC60,
              theC61, theC62, theC63, theC64, theC65, theC66, theC67, theC68, theC69, theC70,
              theC71, theC72, theC73, theC74, theC75, theC76, theC77, theC78, theC79, theC80,
              theC81, theC82, theC83, theC84, theC85, theC86, theC87, theC88, theC89, theC90) theChannels;

G4CollisionNNToDeltaNstar::G4CollisionNNToDeltaNstar()
{ 
  G4ForEach<theChannels, Resolve>::Apply(this);
}

